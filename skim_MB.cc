#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <TSpline.h>

#include <vector>
#include <math.h>
#include <TRandom3.h>

#define n_event_vars 4

TTree * get_tree(TString file_name){

  std::cout << "Opening: " << file_name << std::endl;
  TFile *file = TFile::Open(file_name);

  if (file == NULL) {
    std::cout << " failed" << std::endl;
    exit(EXIT_FAILURE);
  }

  TTree * tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));
  if ( tree_event == NULL) {
    tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
    if (tree_event == NULL) {
      std::cout << " tree fail " << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  return tree_event;
}

TH1F ** get_trigger_histos(TTree * tree){
  //function loops through trigger ntuple
  //and fills event data histograms

  TSpline3 ** trig_splines = new TSpline3*[n_event_vars];
  TH1F ** trig_histos = new TH1F*[n_event_vars];

  Float_t centrality;
  Double_t primary_vertex[3];
  Float_t event_plane_angle[3];
  Float_t multiplicity_v0[64];

  tree->SetBranchAddress("centrality_v0m", &centrality);
  tree->SetBranchAddress("primary_vertex", primary_vertex);
  tree->SetBranchAddress("multiplicity_v0", &multiplicity_v0[0]);
  tree->SetBranchAddress("event_plane_psi_v0", &event_plane_angle[0]);

  TH1F* h_centrality = new TH1F("trigger_centrality", "centrality (Trigger)", 200, 0, 100);
  TH1F* h_z_vertices = new TH1F("trigger_Primary_Vertex", "Z-vertex (Trigger)", 240, -12, 12);
  TH1F* h_flow = new TH1F("trigger_Flow", "Flow (Trigger)", 500, -2, 2);
  TH1F* h_multiplicity = new TH1F("trigger_multiplicity", "multiplicity (Trigger)", 1000, 0, 1000);

  UInt_t nEvents = tree->GetEntries();
  /* nEvents = 10000; */
  for(Long64_t ievent = 0; ievent < nEvents ; ievent++)
  {
    fprintf(stderr, "\r%s:%d: Entry %llu / %llu", __FILE__, __LINE__, ievent, nEvents);
    tree->GetEntry(ievent);
    h_centrality->Fill(centrality);
    h_z_vertices->Fill(primary_vertex[0]);
    h_flow->Fill(event_plane_angle[1]); //elliptic flow
    h_multiplicity->Fill(multiplicity_v0[0]);
  }

  //normalize histograms to 1.0
  //this way, 1/spline can be used as probability
  //We want the frequency probabilyt, not the prob. density function. 
  /* h_centrality->Scale(1.0/h_centrality->Integral("width"),"width"); //This  gives PDF, but isn't working. Understand div by Sum() and without width. */
  /* h_centrality->Scale(1.0/h_centrality->Sum(),"width"); //This  gives PDF, but isn't working. Understand div by Sum() and without width. */
  h_centrality->Scale(1.0/h_centrality->GetMaximum());
  h_z_vertices->Scale(1.0/h_z_vertices->Integral());
  h_flow->Scale(1.0/h_flow->Integral());
  h_multiplicity->Scale(1.0/h_multiplicity->Integral());

  trig_histos[0] = h_centrality;
  trig_histos[1] = h_z_vertices;
  trig_histos[2] = h_flow;
  trig_histos[3] = h_multiplicity;

  return trig_histos;
}



TSpline3 ** get_trigger_splines(TH1F ** trig_histos){

  TSpline3 ** trig_splines = new TSpline3*[n_event_vars];
  //FIXME: use iterators and auto
  trig_splines[0] = new TSpline3(trig_histos[0]); trig_splines[0]->SetName("spline_centrality");
  trig_splines[1] = new TSpline3(trig_histos[1]); trig_splines[1]->SetName("spline_z_vertex");
  trig_splines[2] = new TSpline3(trig_histos[2]); trig_splines[2]->SetName("spline_flow");
  trig_splines[3] = new TSpline3(trig_histos[3]); trig_splines[3]->SetName("spline_multiplicity");

  return trig_splines;
}

/* std::array<TSpline3*> get_splines_arrays(TH1F ** trig_histos){ */

/*   TSpline3 ** trig_splines = new TSpline3*[n_event_vars]; */
/*   std::array<TSpline3*> */   
/*   //FIXME: use iterators and auto */
/*   trig_splines[0] = new TSpline3(trig_histos[0]); trig_splines[0]->SetName("spline_centrality"); */
/*   trig_splines[1] = new TSpline3(trig_histos[1]); trig_splines[1]->SetName("spline_z_vertex"); */
/*   trig_splines[2] = new TSpline3(trig_histos[2]); trig_splines[2]->SetName("spline_flow"); */
/*   trig_splines[3] = new TSpline3(trig_histos[3]); trig_splines[3]->SetName("spline_multiplicity"); */

/*   return trig_splines; */
/* } */

bool skim(TRandom3 * TRand, TSpline3 * cent_spline, Float_t centrality)
{
  float r = TRand->Rndm();
  float prob = 1.0-cent_spline->Eval(centrality);
  std::cout<<std::endl<<r<<" "<<prob<<"bool = "<<(r < prob)<<std::endl;
  return(r < prob);
}

TTree * skim_tree(TTree *tree, TSpline3 ** trig_splines){

  //function loops through trigger ntuple
  //and fills event data histograms

  Float_t centrality;
  Double_t primary_vertex[3];
  Float_t event_plane_angle[3];
  Float_t multiplicity_v0[64];

  tree->SetBranchAddress("centrality_v0m", &centrality);
  tree->SetBranchAddress("primary_vertex", primary_vertex);
  tree->SetBranchAddress("event_plane_psi_v0", &event_plane_angle[0]);
  tree->SetBranchAddress("multiplicity_v0", &multiplicity_v0[0]);

  TFile* fout = new TFile("Skim.root","RECREATE");
  fout->cd();
  TTree *skimmed_tree = tree->CloneTree(0);
  TRandom3 * rando = new TRandom3();

  UInt_t nEvents = tree->GetEntries();
  /* nEvents = 1000; */
  for(Long64_t ievent = 0; ievent < nEvents ; ievent++)
  {
    /* fprintf(stderr, "\r%s:%d: Entry %llu / %llu", __FILE__, __LINE__, ievent, nEvents); */
    tree->GetEntry(ievent);
    if (skim(rando,trig_splines[0],centrality))
      skimmed_tree->Fill();
  }

  /* skimmed_tree->AutoSave(); */
  //normalize histograms to 1.0
  //this way, 1/spline can be used as probability

  skimmed_tree->Write();
  fout->Close();
  return skimmed_tree;
}


TCanvas ** get_canvases(TH1F ** trigger_histos, TSpline3 ** trigger_splines)
{
  TCanvas ** histo_spline = new TCanvas*[n_event_vars];
  histo_spline[0] = new TCanvas("centrality_canv", "centrality_canv",800,800);
  trigger_histos[0]->Draw();
  trigger_splines[0]->SetLineColor(2);
  trigger_splines[0]->Draw("same");
  histo_spline[0]->Update();
  histo_spline[0]->Draw();

  histo_spline[1] = new TCanvas("vertex_canv", "vertex_canv",800,800);
  trigger_histos[1]->Draw();
  trigger_splines[1]->SetLineColor(2);
  trigger_splines[1]->Draw("same");
  histo_spline[1]->Update();
  histo_spline[1]->Draw();

  histo_spline[2] = new TCanvas("flow_canv", "flow_canv",800,800);
  trigger_histos[2]->Draw();
  trigger_splines[2]->SetLineColor(2);
  trigger_splines[2]->Draw("same");
  histo_spline[2]->Update();
  histo_spline[2]->Draw();
}

int main(int argc, char *argv[])
{

  if (argc < 3)
  {
    fprintf(stderr,"Batch Syntax is [Gamma-Triggered Paired Root]");
    exit(EXIT_FAILURE);
  }

  TTree *trigger_Tree = get_tree((TString)argv[1]);
  TTree *mb_tree= get_tree((TString)argv[2]);

  TH1F ** trigger_histos = get_trigger_histos(trigger_Tree);
  TSpline3 ** trigger_splines = get_trigger_splines(trigger_histos);

  /* TCanvas ** canvases = get_canvases(trigger_histos,trigger_splines); */
  TTree * skimmed_tree = skim_tree(mb_tree,trigger_splines);

  TFile* fout = new TFile("Skim.root","UPDATE");
  fout->cd();
  TCanvas ** histo_spline = new TCanvas*[n_event_vars];
  histo_spline[0] = new TCanvas("centrality_canv", "centrality_canv",800,800);
  trigger_histos[0]->Draw();
  trigger_splines[0]->SetLineColor(2);
  trigger_splines[0]->Draw("same");
  histo_spline[0]->Update();
  histo_spline[0]->Draw();

  histo_spline[1] = new TCanvas("vertex_canv", "vertex_canv",800,800);
  trigger_histos[1]->Draw();
  trigger_splines[1]->SetLineColor(2);
  trigger_splines[1]->Draw("same");
  histo_spline[1]->Update();
  histo_spline[1]->Draw();

  histo_spline[2] = new TCanvas("flow_canv", "flow_canv",800,800);
  trigger_histos[2]->Draw();
  trigger_splines[2]->SetLineColor(2);
  trigger_splines[2]->Draw("same");
  histo_spline[2]->Update();

  histo_spline[0]->Write();
  histo_spline[1]->Write();
  histo_spline[2]->Write();

  trigger_histos[0]->Write();
  trigger_histos[1]->Write();
  trigger_histos[2]->Write();
  trigger_histos[3]->Write();
  trigger_splines[0]->Write();
  trigger_splines[1]->Write();
  trigger_splines[2]->Write();
  trigger_splines[3]->Write();



  fout->Close();
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}
/* TString trigger_name = (TString)argv[1]; */
/* std::cout << "Opening: " << (TString)argv[1] << std::endl; */
/* TFile *file = TFile::Open(trigger_name); */

/* if (file == NULL) { */
/*   std::cout << " fail" << std::endl; */
/*   exit(EXIT_FAILURE); */
/* } */

/* TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event")); */
/* if (_tree_event == NULL) { */
/*   _tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event")); */
/*   if (_tree_event == NULL) { */
/*     std::cout << " tree fail " << std::endl; */
/*     exit(EXIT_FAILURE); */
/*   } */
/* } */

/* //MB File and Tree Check */
/* TString mb_name = (TString)argv[2]; */
/* std::cout << "Opening: " << (TString)argv[2] << std::endl; */
/* TFile *mb_file = TFile::Open(mb_name); */

/* if (file == NULL) { */
/*   std::cout << " fail" << std::endl; */
/*   exit(EXIT_FAILURE); */
/* } */

/* TTree * mb_tree = dynamic_cast<TTree *>(file->Get("_tree_event")); */
/* if (mb_tree == NULL) { */
/*   mb_tree = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event")); */
/*   if (mb_tree == NULL) { */
/*     std::cout << " tree fail " << std::endl; */
/*     exit(EXIT_FAILURE); */
/*   } */
/* } */
/* TSpline3 *spline_z = new TSpline3(h_z_vertices); */
/* TSpline3 *spline_centrality = new TSpline3(h_centrality); */
/* TSpline3 *spline_flow = new TSpline3(h_flow); */

/* TCanvas *cent_canv = new TCanvas("centrality_canv", "centrality_canv",800,800); */
/* h_centrality->Draw(); */
/* spline_centrality->SetLineColor(2); */
/* spline_centrality->Draw("same"); */
/* cent_canv->Update(); */
/* cent_canv->Draw(); */
/* cent_canv->Write(); */

/* TCanvas *flow_canv = new TCanvas("flow_canv", "flow_canv",800,800); */
/* h_flow->Draw(); */
/* spline_flow->SetLineColor(2); */
/* spline_flow->Draw("same"); */
/* flow_canv->Update(); */
/* flow_canv->Draw(); */
/* flow_canv->Write(); */

/* TCanvas *vertex_canv = new TCanvas("vertex_canv", "vertex_canv",800,800); */
/* h_flow->Draw(); */
/* spline_z->SetLineColor(2); */
/* spline_z->Draw("same"); */
/* vertex_canv->Update(); */
/* vertex_canv->Draw(); */
/* vertex_canv->Write(); */
