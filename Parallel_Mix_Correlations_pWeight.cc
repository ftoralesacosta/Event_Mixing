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
#include "H5Cpp.h"
#include <chrono>
#include "omp.h"


//#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sched.h>
#include <omp.h>

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_ISO_ITS_04_SUB, CLUSTER_ISO_TPC_04_SUB, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};
using namespace H5;

//int main(int argc, const char rootfilename, const char hdf5filename, const char *mixing_start, const char *mixing_end)

Float_t Get_Purity_ErrFunction(Float_t pT_GeV, std::string deviation,bool Is_pp=false) {

  Float_t purity_val = 0;                                                             

  //p-Pb                                                                                                                                                                                                    
  Float_t par[3] = {0.564400011365,
                    8.81589638183,
                    12.5209151159};

  if (Is_pp){
    par[0] = 0.507999134652;
    par[1] = 8.33474628579;
    par[2] = 13.3714887905;
  }

  if (strcmp(deviation.data(),"Plus")==0){
    par[0] = 0.60750016509;
    par[1] = 7.05184155403;
    par[2] = 13.6116163603;
  }

  if (strcmp(deviation.data(),"Minus")==0){
    par[0] = 0.479958593235;
    par[1] = 9.05392932723;
    par[2] = 10.2061359452;
  }


  purity_val = par[0]*TMath::Erf((pT_GeV-par[1])/par[2]);
  return purity_val;
}


int main(int argc ,char *argv[])
{
  if (argc < 7) {
    fprintf(stderr,"\n Batch Syntax is [Gamma-Triggered Paired Root], [Min-Bias HDF5] [Mix Start] [Mix End] [Track Skim GeV] [pp or p-Pb]\n");
    exit(EXIT_FAILURE);
  }


  bool Is_pp = false;

  
  std::string coll_system = argv[6];
  if (strcmp(coll_system.c_str(), "pp") == 0)
    Is_pp = true;
    
  if (Is_pp)
      fprintf(stderr,"\n PROTON PROTON SELECTED \n \n");

  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");
  
  TString root_file = (TString)argv[1];
  std::cout << "Opening: " << (TString)argv[1] << std::endl;

  std::string file_str = argv[2];
  const H5std_string hdf5_file_name(file_str.c_str()); //explicitly requires c_string
  TString hdf5_file = (TString)argv[2];
  fprintf(stderr,"%d: HDF5 File = ",__LINE__);
  fprintf(stderr,hdf5_file);
  std::cout<<std::endl;
  
  size_t mix_start = stoull(std::string(argv[3]));
  size_t mix_end = stoull(std::string(argv[4]));

  int GeV_Track_Skim = atoi(argv[5]);
  std::cout<<std::endl<<"mix start is "<<mix_start<<std::endl;
  std::cout<<"mix end is "<<mix_end<<std::endl;
  fprintf(stderr,"Using %iGeV Track Skimmed from batch Script \n",GeV_Track_Skim);

  size_t nmix = 300;
  fprintf(stderr,"Number of Mixed Events: %i \n",nmix);

  //Config File
  FILE* config = fopen("../Corr_config.yaml", "r");
  std::cout<<config<<std::endl;
  double DNN_min = 0;
  double DNN_max = 0;
  double Lambda0_cut = 0;
  double pT_min = 0;
  double pT_max = 0;
  double Eta_max = 0;
  double Cluster_min = 0;
  float Cluster_DtoBad = 0;
  UChar_t Cluster_NLocal_Max = 0;
  double EcrossoverE_min = 0;
  double cluster_time = 20;
  bool do_pile = false;
  
  float track_pT_min = 0.0;
  float track_pT_max = 0.0;
  int Track_Cut_Bit = 0;
  int track_chi_max = 1;
  double iso_max = 0;
  double noniso_min = 0;
  double noniso_max = 0;
  double deta_max = 0;
  isolationDet determiner = CLUSTER_ISO_ITS_04; //replaced by config file. Check on Print
  int n_eta_bins = 0;
  int n_phi_bins = 0;
  std::string shower_shape = "DNN";
  std::string purity_deviation = "None";
  // zT & pT bins
  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins+1];
  ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
   
  int nptbins = 3;
  float* ptbins;
  ptbins = new float[nptbins+1];
  ptbins[0] = 10.0; ptbins[1] = 11; ptbins[2] = 12.5; ptbins[3] = 16;

  // Read/Loop through config file
  char line[MAX_INPUT_LENGTH];

  while (fgets(line, MAX_INPUT_LENGTH, config) != NULL) {
      if (line[0] == '#') continue;

      char key[MAX_INPUT_LENGTH];
      char dummy[MAX_INPUT_LENGTH];
      char value[MAX_INPUT_LENGTH];

      // Cap off key[0] and value[0] with null characters and load the key, dummy-characters, and value of the line into their respective arrays
      key[0] = '\0';
      value[0] = '\0';
      sscanf(line, "%[^:]:%[ \t]%100[^\n]", key, dummy, value);
      
      //Read Config File: Detect Keys 
      if (strcmp(key, "DNN_min") == 0) {
          DNN_min = atof(value);
          std::cout << "DNN_min: " << DNN_min << std::endl; }

      else if (strcmp(key, "DNN_max") == 0) {
          DNN_max = atof(value);
          std::cout << "DNN_max: " << DNN_max << std::endl; }

      else if (strcmp(key,"Lambda0_cut") == 0){
	Lambda0_cut = atof(value);
	std::cout << "Lambda0 cut at: " <<Lambda0_cut <<std::endl;}

      else if (strcmp(key, "pT_min") == 0) {
          pT_min = atof(value);
          std::cout << "pT_min: " << pT_min << std::endl; }

      else if (strcmp(key, "pT_max") == 0) {
          pT_max = atof(value);
          std::cout << "pT_max: " << pT_max << std::endl; }

      else if (strcmp(key, "Eta_max") == 0) {
          Eta_max = atof(value);
          std::cout << "Eta_max: " << Eta_max << std::endl;
      }
      else if (strcmp(key, "Cluster_min") == 0) {
          Cluster_min = atof(value);
          std::cout << "Cluster_min: " << Cluster_min << std::endl; }

      else if (strcmp(key, "Cluster_dist_to_bad_channel") == 0){
        Cluster_DtoBad = atof(value);
	std::cout << "Cluster_DtoBad: "<< Cluster_DtoBad << std::endl;}

      else if (strcmp(key, "Cluster_N_Local_Maxima") == 0){
        Cluster_NLocal_Max = atof(value);
	std::cout << "Cluster_NLocal_Max: "<< Cluster_NLocal_Max << std::endl;}

      else if (strcmp(key, "EcrossoverE_min") == 0) {
          EcrossoverE_min = atof(value);
          std::cout << "EcrossoverE_min; " << EcrossoverE_min << std::endl; }

      else if (strcmp(key, "iso_max") == 0) {
          iso_max = atof(value);
          std::cout << "iso_max: " << iso_max << std::endl; }

      else if (strcmp(key, "noniso_min") == 0) {
          noniso_min = atof(value);
          std::cout << "noniso_min: " << noniso_min << std::endl; }

      else if (strcmp(key, "noniso_max") == 0) {
          noniso_max = atof(value);
          std::cout << "noniso_max: " << noniso_max << std::endl; }

      else if (strcmp(key, "deta_max") == 0) {
          deta_max = atof(value);
          std::cout << "deta_max: " << deta_max << std::endl; }

      else if (strcmp(key, "do_pileup_cut") == 0) {
        if (strcmp(value,"true") == 0)
          do_pile = true;
        std::cout << "do_pileup_cut: " << do_pile << std::endl; }
      
      else if (strcmp(key, "Cluster_Time") == 0){
	cluster_time = atof(value);
        std::cout << "Cluster_Time: "<< cluster_time << std::endl;}
      
      else if (strcmp(key, "N_Phi_Bins") == 0) {
	n_phi_bins = atoi(value);
	std::cout << "Number of Phi Bins: " << n_phi_bins << std::endl; }

      else if (strcmp(key, "N_Eta_Bins") == 0) {
        n_eta_bins = atoi(value);
	std::cout << "Number of Eta Bins: " << n_eta_bins << std::endl; }

      else if (strcmp(key, "Track_pT_Min") == 0) {
          track_pT_min = atof(value);
          std::cout << "Track pT Min: " << track_pT_min << std::endl; }

      else if (strcmp(key, "Track_pT_Max") == 0) {
          track_pT_max = atof(value);
          std::cout << "Track pT Max: " << track_pT_max << std::endl; }

      else if (strcmp(key, "Track_Cut_Bit") == 0) {
          Track_Cut_Bit = atoi(value);
	  if (Is_pp) {
	    fprintf(stderr,"%d: Overwrite Track Cut Bit for PROTON PROTON COLLISIONS \n",__LINE__);
	    Track_Cut_Bit = 3;
	  }
          std::cout << "Track Cut Bit: " << Track_Cut_Bit << std::endl;}

      else if (strcmp(key, "Track_Chi_Max") == 0) {
          track_chi_max = atoi(value);
          std::cout << "Track Chi Max: " << track_chi_max << std::endl; }

      
      else if (strcmp(key, "Zt_bins") == 0) {
          nztbins = -1;
          for (const char *v = value; *v != ']';) {
              while (*v != ']' && !isdigit(*v)) v++;
	      nztbins++;
              while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }

          ztbins = new float[nztbins + 1];
          int i = 0;
          for (const char *v = value; *v != ']' ;) {
              while (*v != ']' && !isdigit(*v)) v++;
              ztbins[i] = atof(v);
              i++;              
              while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }

          std::cout << "Number of Zt bins: " << nztbins << std::endl << "Zt bins: {";
          for (int i = 0; i <= nztbins; i++)
	    std::cout << ztbins[i] << ", ";
          std::cout << "}\n"; 
      }

      else if (strcmp(key, "Pt_bins") == 0) {
          nptbins = -1;
          for (const char *v = value; *v != ']';) {
            while (*v != ']' && !isdigit(*v)) v++;
	    nptbins++;
	    while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
	
          ptbins = new float[nptbins + 1];
	  int i = 0;
	  for (const char *v = value; *v != ']' ;) {
	     while (*v != ']' && !isdigit(*v))  v++;
	     ptbins[i] = atof(v);
	     i++;
	     while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
	
	  std::cout << "Number of Pt bins: " << nptbins << std::endl << "Pt bins: {";
          for (int i = 0; i <= nptbins; i++)
	    std::cout << ptbins[i] << ", ";
	  std::cout << "}\n";
      }

      else if (strcmp(key, "Cluster_isolation_determinant") == 0) {
          if (strcmp(value, "cluster_iso_tpc_04") == 0){
              determiner = CLUSTER_ISO_TPC_04;
              std::cout << "Isolation Variable: cluster_iso_tpc_04" << std::endl; }

          else if (strcmp(value, "cluster_iso_its_04") == 0){
              determiner = CLUSTER_ISO_ITS_04;
              std::cout << "Isolation Variable: cluster_iso_its_04" << std::endl; }

          else if (strcmp(value, "cluster_iso_its_04_sub") == 0){
              determiner = CLUSTER_ISO_ITS_04_SUB;
              std::cout << "Isolation Variable: cluster_iso_its_04_sub" << std::endl; }

          else if (strcmp(value, "cluster_iso_tpc_04_sub") == 0){
              determiner = CLUSTER_ISO_TPC_04_SUB;
              std::cout << "Isolation Variable: cluster_iso_tpc_04_sub" << std::endl; }

	  else if (strcmp(value, "cluster_frixione_tpc_04_02") == 0){
              determiner = CLUSTER_FRIXIONE_TPC_04_02;
              std::cout << "Isolation Variable: cluster_frixione_tpc_04_02" << std::endl; }

          else if (strcmp(value, "cluster_frixione_its_04_02") == 0){
              determiner = CLUSTER_FRIXIONE_ITS_04_02;
              std::cout << "Isolation Variable: cluster_frixione_its_04_02" << std::endl; }

          else {
              std::cout << "ERROR: Cluster_isolation_determinant in configuration file must be \"cluster_iso_tpc_04\", \"cluster_iso_its_04\", \"cluster_frixione_tpc_04_02\", or \"cluster_frixione_its_04_02\"" << std::endl << "Aborting the program" << std::endl;
              exit(EXIT_FAILURE); }
      }

      else if (strcmp(key, "Shower_Shape") == 0){
	shower_shape = value;
	std::cout<<"Shower Shape: "<<shower_shape.data()<<std::endl;
      }
      
      else if (strcmp(key, "Purity_Dev") == 0){
        purity_deviation = value;
	std::cout<<"Purity Deviation Change: "<<purity_deviation.data()<<std::endl;
      }

      else std::cout << "WARNING: Unrecognized keyvariable " << key << std::endl;
  
  }
  //end Config Loop

  fclose(config);
  
  for (int i = 0; i <= nztbins; i++)
    std::cout << "zt bound: " << ztbins[i] << std::endl;
  for (int i = 0; i <= nptbins; i++)
    std::cout << "pt bound: " << ptbins[i] << std::endl;


  //HISTOGRAMS
  TCanvas canvas("canvas", "");

  TH1D h_ntrig("h_ntrig", "", 2, -0.5,1.0);

  TH2D* Signal_pT_Dist = new TH2D("Signal_pT_Dist","Cluster Pt Spectrum For Isolation (its_04) bins 0.55 < DNN < 0.85",24,10,16,5,-0.5,2);
  TH2D* BKGD_pT_Dist = new TH2D("BKGD_pT_Dist","Cluster Pt Spectrum For Isolation (its_04) bins 0.0 < DNN < 0.3",24,10,16,5,-0.5,2);

  TH2D* Corr[nztbins*nptbins];
  TH2D* Signal_Corr[nztbins*nptbins];
  TH2D* BKGD_Corr[nztbins*nptbins];
  TH2D* Inclusive_Corr[nztbins*nptbins];  

  TH1D* H_Signal_Triggers[nptbins];
  TH1D* H_BKGD_Triggers[nptbins];
  float N_Signal_Triggers = 0;
  float N_BKGD_Triggers = 0;
  

  TH1D* z_Vertices_individual = new TH1D("Primary_Vertex_root", "Z-vertex (ROOT)",240, -12, 12);
  TH1D* z_Vertices_hdf5 = new TH1D("Primary_Vertex_hdf5", "Z-vertex (hdf5)", 240,-12, 12);
  TH1D* z_Vertices = new TH1D("Delta_Primary_Vertex", "#Delta V_z Distribution", 240,-12, 12);

  TH1D* Multiplicity_individual = new TH1D("Multiplicity_root", "Multiplicity (ROOT)", 1000, 0, 1000);
  TH1D* Multiplicity_hdf5 = new TH1D("Multplicity_hdf5", "Multiplicity (hdf5)", 500, 0, 1000);
  TH1D* Multiplicity = new TH1D("Delta_Multiplicity", "#Delta Multiplicit Distribution", 500, 0, 1000);

  TH2D* N_ME = new TH2D("N_ME", "Distribution No. Mixed Events Passed",300,0,300,500,0,1000);



    for (int ipt = 0; ipt <nptbins; ipt++) {
      H_Signal_Triggers[ipt] = new TH1D(
      Form("N_DNN%i_Triggers_pT%1.0f_%1.0f",1,ptbins[ipt],ptbins[ipt+1]),
      "Number of Isolated Photon Triggers", 2, -0.5,1.0);

      H_BKGD_Triggers[ipt] = new TH1D(
      Form("N_DNN%i_Triggers_pT%1.0f_%1.0f",2,ptbins[ipt],ptbins[ipt+1]),
      "Number of Isolated Low DNN Photon Triggers", 2, -0.5,1.0);

      for (int izt = 0; izt<nztbins; izt++){

      Corr[izt+ipt*nztbins] = new TH2D(Form("Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H [all] Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);

      Corr[izt+ipt*nztbins]->Sumw2();
      Corr[izt+ipt*nztbins]->SetMinimum(0.);

      Signal_Corr[izt+ipt*nztbins] = new TH2D(Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",1,ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H [Iso] Correlation", n_phi_bins,0,M_PI,n_eta_bins, -1.4, 1.4);

      Signal_Corr[izt+ipt*nztbins]->Sumw2();
      Signal_Corr[izt+ipt*nztbins]->SetMinimum(0.);

      BKGD_Corr[izt+ipt*nztbins] = new TH2D(Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H [AntiIso] Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);

      BKGD_Corr[izt+ipt*nztbins]->Sumw2();
      BKGD_Corr[izt+ipt*nztbins]->SetMinimum(0.);

      Inclusive_Corr[izt+ipt*nztbins] = new TH2D(Form("Inclusive_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",ptbins[ipt],ptbins[ipt+1],
      100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H Isclusive Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);

      Inclusive_Corr[izt+ipt*nztbins]->Sumw2();
      Inclusive_Corr[izt+ipt*nztbins]->SetMinimum(0.);

    }//zt bins
  }//pt bins                           

    TFile *file = TFile::Open(root_file);

    if (file == NULL) {
      std::cout << " file fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();
    
    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));
    if (_tree_event == NULL) {
      _tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
      if (_tree_event == NULL) {
	std::cout << " tree fail " << std::endl;
	exit(EXIT_FAILURE);
      }  
    }

    //variables
    Bool_t is_pileup_from_spd_5_08;
    Double_t primary_vertex[3];
    Float_t multiplicity_v0[64];
    Float_t ue_estimate_its_const;
    Float_t ue_estimate_tpc_const;
    
    UInt_t ntrack;
    Float_t track_e[NTRACK_MAX];
    Float_t track_pt[NTRACK_MAX];
    Float_t track_eta[NTRACK_MAX];
    Float_t track_phi[NTRACK_MAX];
    UChar_t track_quality[NTRACK_MAX];

    UInt_t ncluster;
    Float_t cluster_e[NTRACK_MAX];
    Float_t cluster_e_cross[NTRACK_MAX];
    Float_t cluster_pt[NTRACK_MAX];
    Float_t cluster_eta[NTRACK_MAX];
    Float_t cluster_phi[NTRACK_MAX]; 
    Float_t cluster_iso_tpc_04[NTRACK_MAX];
    Float_t cluster_iso_its_04[NTRACK_MAX];
    Float_t cluster_frixione_tpc_04_02[NTRACK_MAX];
    Float_t cluster_frixione_its_04_02[NTRACK_MAX];
    Float_t cluster_s_nphoton[NTRACK_MAX][4];
    unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
    Int_t cluster_ncell[NTRACK_MAX];
    UShort_t  cluster_cell_id_max[NTRACK_MAX];
    Float_t cluster_lambda_square[NTRACK_MAX][2];   
    Float_t cell_e[17664];
    Float_t cluster_distance_to_bad_channel[NTRACK_MAX];
    UChar_t cluster_nlocal_maxima[NTRACK_MAX];

    Float_t cluster_tof[NTRACK_MAX];
    Float_t cluster_iso_its_04_ue[NTRACK_MAX];
    Float_t cluster_iso_tpc_04_ue[NTRACK_MAX];
    
    fprintf(stderr,"Initializing Mixing Branch to %i ME",nmix);
    Long64_t mix_events[300];

    //MC
    unsigned int nmc_truth;
    Float_t mc_truth_pt[NTRACK_MAX];
    Float_t mc_truth_eta[NTRACK_MAX];  
    Float_t mc_truth_phi[NTRACK_MAX];
    short mc_truth_pdg_code[NTRACK_MAX];
    short mc_truth_first_parent_pdg_code[NTRACK_MAX];
    char mc_truth_charge[NTRACK_MAX];

    Float_t mc_truth_first_parent_e[NTRACK_MAX];
    Float_t mc_truth_first_parent_pt[NTRACK_MAX];
    Float_t mc_truth_first_parent_eta[NTRACK_MAX];
    Float_t mc_truth_first_parent_phi[NTRACK_MAX];
    UChar_t mc_truth_status[NTRACK_MAX];
     
    // Set the branch addresses of the branches in the TTrees
    _tree_event->SetBranchStatus("*mc*", 0);
  
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("multiplicity_v0", multiplicity_v0);
    _tree_event->SetBranchAddress("ue_estimate_its_const", &ue_estimate_its_const);
    _tree_event->SetBranchAddress("ue_estimate_tpc_const", &ue_estimate_tpc_const);
    _tree_event->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);
    
    _tree_event->SetBranchAddress("ntrack", &ntrack);
    _tree_event->SetBranchAddress("track_e", track_e);
    _tree_event->SetBranchAddress("track_pt", track_pt);
    _tree_event->SetBranchAddress("track_eta", track_eta);
    _tree_event->SetBranchAddress("track_phi", track_phi);
    _tree_event->SetBranchAddress("track_quality", track_quality);

    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_e", cluster_e);
    _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
    _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
    _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
    _tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
    _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
    _tree_event->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
    _tree_event->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);
    _tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);
    _tree_event->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);

    _tree_event->SetBranchAddress("cluster_tof", cluster_tof);
    _tree_event->SetBranchAddress("cluster_iso_its_04_ue",cluster_iso_its_04_ue);
    _tree_event->SetBranchAddress("cluster_iso_tpc_04_ue",cluster_iso_tpc_04_ue);
    
    _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    _tree_event->SetBranchAddress("cell_e", cell_e);

    _tree_event->SetBranchAddress("mixed_events", mix_events);
    //_tree_event->SetBranchAddress("Mix_Events", mix_events);
    //_tree_event->SetBranchAddress("LimitUse_Mixed_Events", Mix_Events);
    
    std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;

    //Using low level hdf5 API
    //open hdf5: Define size of data from file, explicitly allocate memory in hdf5 space and array size

  
    //H5File h5_file( hdf5_file_name, H5F_ACC_RDONLY );
    H5File h5_file( file_str.c_str(), H5F_ACC_RDONLY );

    const std::string track_ds_name("track");
    DataSet track_dataset = h5_file.openDataSet( track_ds_name.c_str());
    DataSpace track_dataspace = track_dataset.getSpace();
    
    const std::string cluster_ds_name( "cluster" );
    DataSet cluster_dataset = h5_file.openDataSet( cluster_ds_name.c_str() );
    DataSpace cluster_dataspace = cluster_dataset.getSpace();

    const std::string event_ds_name( "event" );
    DataSet event_dataset = h5_file.openDataSet( event_ds_name.c_str() );
    DataSpace event_dataspace = event_dataset.getSpace();

    //Initialize Track Dimensions
    const int track_ndims = track_dataspace.getSimpleExtentNdims();
    hsize_t track_maxdims[track_ndims];
    hsize_t trackdims[track_ndims];
    track_dataspace.getSimpleExtentDims(trackdims, track_maxdims);

    UInt_t ntrack_max = trackdims[1];
    UInt_t NTrack_Vars = trackdims[2];

    //Initialize Cluster Dimensions
    const int cluster_ndims = cluster_dataspace.getSimpleExtentNdims();
    hsize_t cluster_maxdims[cluster_ndims];
    hsize_t clusterdims[cluster_ndims];
    cluster_dataspace.getSimpleExtentDims(clusterdims, cluster_maxdims);

    UInt_t ncluster_max = clusterdims[1];
    UInt_t NCluster_Vars = clusterdims[2];

    //Initialize Event Dimensions
    const int event_ndims = event_dataspace.getSimpleExtentNdims();
    hsize_t event_maxdims[event_ndims];
    hsize_t eventdims[event_ndims];
    event_dataspace.getSimpleExtentDims(eventdims, event_maxdims);

    //Event dataspace rank only rank 2
    UInt_t NEvent_Vars = eventdims[1];

    fprintf(stderr, "\n%s:%d: n track variables:%i n cluster variables:%i n event variables:%i\n", __FILE__, __LINE__, NTrack_Vars,NCluster_Vars,NEvent_Vars);
    fprintf(stderr, "\n%s:%d: maximum tracks:%i maximum clusters:%i\n", __FILE__, __LINE__, ntrack_max,ncluster_max);

    //Define array hyperslab will be read into
    hsize_t block_size = eventdims[0];

    //block_size = eventdims[0]/2; //FIXME: RAM on cori recently reduced.
    if (eventdims[0] > 1000000){
      block_size = 1000000;
      fprintf(stderr,"%d: USING limited block size for large hdf5 file (block size = %d)\n",__LINE__,block_size);                                                                                                                                        
    }
    // if((strstr(root_file,"13f") !=NULL) && (strstr(root_file,"0") !=NULL)){
    //   block_size = eventdims[0]/2; //FIXME: RAM on cori recently reduced.
    //   fprintf(stderr,"%d: USING Half block size for 13f (block size = %d)\n",__LINE__,block_size);
    // }
    
    // if (strcmp("13f",basic_name)==0){

    // }
    //const Long64_t block_size = 2000;
    
    float track_data_out[block_size][ntrack_max][NTrack_Vars];
    float cluster_data_out[block_size][ncluster_max][NCluster_Vars];
    float event_data_out[block_size][NEvent_Vars];

    //Define hyperslab size and offset in  FILE;
    hsize_t track_offset[3] = {0, 0, 0};
    hsize_t track_count[3] = {block_size, ntrack_max, NTrack_Vars};
    hsize_t cluster_offset[3] = {0, 0, 0};
    hsize_t cluster_count[3] = {block_size, ncluster_max, NCluster_Vars};
    hsize_t event_offset[2] = {0,0};
    hsize_t event_count[2] = {block_size, NEvent_Vars};

    track_dataspace.selectHyperslab( H5S_SELECT_SET, track_count, track_offset );
    cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
    event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");

    //Define the memory dataspace to place hyperslab
    const int RANK_OUT = 3; //# of Dimensions
    const int Event_RANK_OUT = 2; //Different for Events
    DataSpace track_memspace( RANK_OUT, trackdims );
    DataSpace cluster_memspace( RANK_OUT, clusterdims );    
    DataSpace event_memspace( Event_RANK_OUT, eventdims );

    //Define memory offset for hypreslab starting at begining:
    hsize_t track_offset_out[3] = {0};
    hsize_t cluster_offset_out[3] = {0};
    hsize_t event_offset_out[2] = {0};

    //define Dimensions of array, for writing slab to array
    hsize_t track_count_out[3] = {block_size, ntrack_max, NTrack_Vars};
    hsize_t cluster_count_out[3] = {block_size, ncluster_max, NCluster_Vars};
    hsize_t event_count_out[2] = {block_size, NEvent_Vars};


    //READ Min Bias MIXED EVENTS
    //define space in memory for hyperslab, then write ENTIRE file to memory
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    track_memspace.selectHyperslab( H5S_SELECT_SET, track_count_out, track_offset_out );
    track_dataset.read( track_data_out, PredType::NATIVE_FLOAT, track_memspace, track_dataspace );
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "track dataset read into array: OK"); //Test only

    cluster_memspace.selectHyperslab( H5S_SELECT_SET, cluster_count_out, cluster_offset_out );
    cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "cluster dataset read into array: OK");

    event_memspace.selectHyperslab( H5S_SELECT_SET, event_count_out, event_offset_out );
    event_dataset.read( event_data_out, PredType::NATIVE_FLOAT, event_memspace, event_dataspace);

    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    std::cout << "HDF5 READ Time in micro seconds = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;



    Long64_t nentries = _tree_event->GetEntries();   
    UInt_t N_mix_Histos = 50; //Temporary Histograms for THREADING

    // fprintf(stderr,"Obtainng max no. clusters, used for threading\n\n");

    //UInt_t N_mix_Histos = 0; //Temporary Histograms for THREADING
    // for(Long64_t ievent = 0; ievent < nentries ; ievent++) 
    //   {
    // 	fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);
    // 	_tree_event->GetEntry(ievent);
    // 	N_mix_Histos = std::max(N_mix_Histos,ncluster);
    //   }
    
    TH2D* MixCorr[N_mix_Histos][nztbins*nptbins];
    TH2D* MixCorr_Signal[N_mix_Histos][nztbins*nptbins];
    TH2D* MixCorr_BKGD[N_mix_Histos][nztbins*nptbins];
    TH2D* MixInclusiveCorr[N_mix_Histos][nztbins*nptbins];

    for (Long64_t imix = 0; imix < N_mix_Histos; imix++){
      for (int ipt = 0; ipt <nptbins; ipt++) {
	for (int izt = 0; izt<nztbins; izt++){
	  MixCorr[imix][izt+ipt*nztbins] = new TH2D(Form("Mix_Correlation__ME_%i_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",imix,ptbins[ipt],ptbins[ipt+1],
	  100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H [all] Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);
	  MixCorr[imix][izt+ipt*nztbins]->Sumw2();

	  MixCorr_Signal[imix][izt+ipt*nztbins] = new TH2D(Form("DNN%i_Mix_Correlation__ME_%i_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",1,imix,ptbins[ipt],ptbins[ipt+1],
	  100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H [all] Signal Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);
	  MixCorr_Signal[imix][izt+ipt*nztbins]->Sumw2();

	  MixCorr_BKGD[imix][izt+ipt*nztbins] = new TH2D(Form("DNN%i_Mix_Correlation__ME_%i_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,imix,ptbins[ipt],ptbins[ipt+1],
	  100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H [all] Background Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);
	  MixCorr_BKGD[imix][izt+ipt*nztbins]->Sumw2();

	  MixInclusiveCorr[imix][izt+ipt*nztbins] = new TH2D(Form("Mix_Inclusive_Correlation__ME_%i_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",imix,ptbins[ipt],ptbins[ipt+1],
	  100*ztbins[izt],100*ztbins[izt+1]),"#gamma-H [all] Correlation", n_phi_bins,0,M_PI, n_eta_bins, -1.4, 1.4);
	  MixInclusiveCorr[imix][izt+ipt*nztbins]->Sumw2();

	}
      }
    }


    Float_t purity_weight_sum[nptbins] = {0};
    Float_t BR_purity_weight_sum[nptbins] = {0};
    bool TPC_Iso_Flag = false;
    
    //MAIN LOOP
    auto event_begin = std::chrono::steady_clock::now();
    for(Long64_t ievent = 0; ievent < nentries ; ievent++){     

      fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);
      _tree_event->GetEntry(ievent);

      float multiplicity_sum = 0;
      for (int k = 0; k < 64; k++)  multiplicity_sum += multiplicity_v0[k];

      int ME_pass_Counter = 0;
      bool first_cluster = true;

      if(do_pile && is_pileup_from_spd_5_08) continue;
            
#pragma omp parallel 
   {   
#pragma omp for schedule(dynamic,1)

      for (ULong64_t n = 0; n < ncluster; n++) {
	if( not(cluster_pt[n]>pT_min and cluster_pt[n]<pT_max)) continue;   //select pt of photons
	if( not(TMath::Abs(cluster_eta[n])<Eta_max)) continue;              //cut edges of detector                                                                          
        if( not(cluster_ncell[n]>Cluster_min)) continue;                    //removes clusters with 1 or 2 cells 
	if( not(cluster_e_cross[n]/cluster_e[n]>EcrossoverE_min)) continue; //removes "spiky" clusters
        if( not(cluster_distance_to_bad_channel[n]>=Cluster_DtoBad)) continue; //removes clusters near bad channels
        //if( not(cluster_nlocal_maxima[n] < Cluster_NLocal_Max)) continue; //require to have at most 2 local maxima.
        if( not(cluster_nlocal_maxima[n] < 3)) continue; //require to have at most 2 local maxima.
	if (not(abs(cluster_tof[n]) < cluster_time)) continue;
	//fprintf(stderr,"\n%s:%d: Cluster Number: %i  THREAD %i\n", __FILE__, __LINE__, n, omp_get_thread_num());

        float isolation;
        if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04[n];
        else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04[n];
        else if (determiner == CLUSTER_ISO_ITS_04_SUB) isolation = cluster_iso_its_04[n]
                                                         + cluster_iso_its_04_ue[n]
                                                         - ue_estimate_its_const*3.1416*0.4*0.4;

        else if (determiner == CLUSTER_ISO_TPC_04_SUB) {
          isolation = cluster_iso_tpc_04[n] + cluster_iso_tpc_04_ue[n]
            - ue_estimate_tpc_const*3.1416*0.4*0.4;
          TPC_Iso_Flag = true;
        }

	Bool_t Signal = false;
        Bool_t Background = false;
	Bool_t Isolated = isolation<iso_max;

        if (strcmp(shower_shape.data(),"Lambda")== 0) {
          if ((cluster_lambda_square[n][0] < Lambda0_cut))
            Signal = true;

          if ((cluster_lambda_square[n][0] > Lambda0_cut))
            Background = true;
        }

        else if (strcmp(shower_shape.data(),"DNN")==0){
          if ((cluster_s_nphoton[n][1] > DNN_min) && (cluster_s_nphoton[n][1]<DNN_max))
            Signal = true;
          if (cluster_s_nphoton[n][1] > 0.0 && cluster_s_nphoton[n][1] < 0.3) //0.3 = DNN_Bkgd
            Background = true;
        }

	Float_t BR_purity_weight;
	if (Background and Isolated){
	  BR_purity_weight = (1.0/Get_Purity_ErrFunction(cluster_pt[n],purity_deviation,Is_pp) - 1);
	  for (int ipt = 0; ipt < nptbins; ipt++){
            if (cluster_pt[n] >= ptbins[ipt] && cluster_pt[n] < ptbins[ipt+1]){
	      BR_purity_weight_sum[ipt] += BR_purity_weight;
	    }
	  }
	}

	Float_t purity_weight;
	if (Signal and Isolated){
          purity_weight = 1.0/Get_Purity_ErrFunction(cluster_pt[n],purity_deviation,Is_pp);
          for (int ipt = 0; ipt < nptbins; ipt++){
            if (cluster_pt[n] >= ptbins[ipt] && cluster_pt[n] < ptbins[ipt+1]){
	      purity_weight_sum[ipt] += purity_weight;
            }
          }
          //fprintf(stderr,"%d: Purity From Fit = %f\n",__LINE__,purity_weight);
        }

	if(first_cluster){
	  z_Vertices_individual->Fill(primary_vertex[2]);
	  Multiplicity_individual->Fill(multiplicity_sum);
	}

	begin = std::chrono::steady_clock::now();

	for (Long64_t mix_counter = mix_start; mix_counter < mix_end+1; mix_counter++){
	  
	  Long64_t mix_event = mix_events[mix_counter];
	  //fprintf(stderr,"%s: %d: Mixed Event Number %lu FROM BRANCH \n",__FILE__,__LINE__,mix_event);

	  if(mix_event >= 9999999) continue;  
	  Long64_t imix = mix_event;

	  //fprintf(stderr,"\n%s:%d: Mixed Event: %imix  THREAD %i\n", __FILE__, __LINE__, imix, omp_get_thread_num());

	  //Cut Paring Tails
	  if (std::isnan(event_data_out[imix][0])) continue;
	  if (std::isnan(event_data_out[imix][1])) continue;
	  if (std::abs(primary_vertex[2]-event_data_out[imix][0]) > 2) continue;
	  if (std::abs(multiplicity_sum - event_data_out[imix][1]) > 40) continue;
	  if (first_cluster) ME_pass_Counter ++;  //simple conuter for ME that pass
	  //MIX with Associated Tracks

	  for (ULong64_t itrack = 0; itrack < ntrack_max; itrack++) {
	    if (std::isnan(track_data_out[imix][itrack][1])) continue;
	    if ((int(track_data_out[imix][itrack][4]+ 0.5)&Track_Cut_Bit)==0) continue; //selection 16
	    if (track_data_out[imix][itrack][1] < track_pT_min) continue; //less than 1GeV or 0.5GeV (YAML)
	    if (track_data_out[imix][itrack][1] > track_pT_max) continue;
	    if (abs(track_data_out[imix][itrack][2]) > 0.8) continue;
	    if (track_data_out[imix][itrack][7] < 4) continue;
	    if ((track_data_out[imix][itrack][8]/track_data_out[imix][itrack][7]) > track_chi_max) continue;
	    //if( not(TMath::Abs(track_data_out[imix][itrack][9])<0.0231+0.0315/TMath::Power(track_data_out[imix][itrack][4],1.3 ))) continue;
	    //	    if (not(TMath::Abs(track_dca_xy[itrack]<2.4))) continue;
	    if (not(TMath::Abs(track_data_out[imix][itrack][9]<2.4))) continue; //dcaxy
	    if (not(TMath::Abs(track_data_out[imix][itrack][10]<3.2))) continue; //dcaz
	    
	    
	    double dRmin = 0.02;
	    //veto charged particles from mixed event tracks
	    bool MixTrack_HasMatch = false;
	    for (unsigned int l = 0; l < ncluster_max; l++){
	      if (std::isnan(cluster_data_out[imix][l][0])) break;
		float dphi = TMath::Abs(cluster_data_out[imix][l][2] - track_data_out[imix][itrack][5]);
		float deta = TMath::Abs(cluster_data_out[imix][l][3] - track_data_out[imix][itrack][6]);
		float dR = sqrt(dphi*dphi + deta*deta);
		if (dR < dRmin)	MixTrack_HasMatch = true;
		break; 
	    }
	    //if (MixTrack_HasMatch) continue;

	    //Observables
 	    Double_t zt = track_data_out[imix][itrack][1]/cluster_pt[n];
	    Float_t DeltaPhi = TMath::Abs(TVector2::Phi_mpi_pi(cluster_phi[n] - track_data_out[imix][itrack][3]));
	    Float_t DeltaEta = cluster_eta[n] - track_data_out[imix][itrack][2];
	    if ((TMath::Abs(DeltaPhi) < 0.005) && (TMath::Abs(DeltaEta) < 0.005)) continue;
	    for (int ipt = 0; ipt < nptbins; ipt++){
	      if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1]){
		for(int izt = 0; izt<nztbins ; izt++){
		  if(zt>ztbins[izt] and  zt<ztbins[izt+1]){

		    if(isolation< iso_max){
		      if (Signal){
			//Signal_Corr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
			MixCorr_Signal[n][izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta,purity_weight);
		      }
		    }
		    if(isolation<iso_max){
		      if (Background)
			//BKGD_Corr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);}
			MixCorr_BKGD[n][izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta,BR_purity_weight);
		    }
		    if(isolation<iso_max){
		      //Corr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
		      MixCorr[n][izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta); 
		    }
		    MixInclusiveCorr[n][izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
		    
		  }//if in zt bin
		} // end loop over zt bins
	      }//end if in pt bin
	    }//end pt loop bin
	  }//end loop over tracks
	}//end loop over mixed events
   } //Parallel

   // end= std::chrono::steady_clock::now();
   // std::cout << "PARALLEL LOOP Time difference in micro seconds = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;

   first_cluster = false;
      }//end loop on clusters.
	//N_ME->Fill(ME_pass_Counter,multiplicity_sum);       
    } //end loop over events   

    auto event_end = std::chrono::steady_clock::now();
    std::cout << "\nEVENT LOOP Time difference in micro seconds = " << std::chrono::duration_cast<std::chrono::microseconds>(event_end - event_begin).count() <<std::endl;

    
    // Write to fout    
    std::string basic_name = std::string(root_file);
    std::string::size_type pos = basic_name.find('_');
    if (pos != std::string::npos)
      basic_name = basic_name.substr(0, pos);

    TFile* fout;
    fout = new TFile(Form("%s_%luGeVTracks_Mixed_%lu_Correlation_%s.root",basic_name.data(),GeV_Track_Skim,(mix_end-mix_start),shower_shape.data()),"RECREATE");
    //File name something like: 13d_0GeVTracks_Mixed_Correlation_Lambda.root

    //Combine and Write Histograms to File
      for (int ipt = 0; ipt<nptbins; ipt++){    
	for (int izt = 0; izt<nztbins; izt++){
	  for (Long64_t imix =0; imix < N_mix_Histos; imix++ ){
	    Corr[izt+ipt*nztbins]->Add(MixCorr[imix][izt+ipt*nztbins],1);
	    Signal_Corr[izt+ipt*nztbins]->Add(MixCorr_Signal[imix][izt+ipt*nztbins],1);
	    BKGD_Corr[izt+ipt*nztbins]->Add(MixCorr_BKGD[imix][izt+ipt*nztbins],1);
	  }
	  Corr[izt+ipt*nztbins]->Write();	  
	}
	for (int izt = 0; izt<nztbins; izt++){
	  Signal_Corr[izt+ipt*nztbins]->Write();
	}
	for (int izt = 0; izt<nztbins; izt++){
	  BKGD_Corr[izt+ipt*nztbins]->Write();
	}
      }
      

    z_Vertices->Write();
    Multiplicity->Write();
    z_Vertices_individual->Write();
    z_Vertices_hdf5->Write();
    Multiplicity_individual->Write();
    Multiplicity_hdf5->Write();
    N_ME->Write();

    fout->Close();     
    
    //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() <<std::endl;

    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}
