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

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};

using namespace H5;

int main(int argc, char *argv[])
{
  if (argc < 5) {
    fprintf(stderr,"Batch Syntax is [Gamma-Triggered Paired Root], [Min-Bias HDF5] [Mix Start] [Mix End]");
    exit(EXIT_FAILURE);
  }

  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");

  TString root_file = (TString)argv[1];
  std::cout << "Opening: " << (TString)argv[1] << std::endl;
  TFile *file = TFile::Open(root_file);

  if (file == NULL) {
    std::cout << " fail" << std::endl;
    exit(EXIT_FAILURE);
  }

  const H5std_string hdf5_file_name(argv[2]);
  H5File h5_file( hdf5_file_name, H5F_ACC_RDONLY ); //hdf5_file_name from argv[2]
  TString hdf5_file = (TString)argv[2];
  std::cout << "Opening: " << hdf5_file << std::endl;

  size_t mix_start = atoi(argv[3]);
  size_t mix_end = atoi(argv[4]);
  fprintf(stderr,"Mixing Event range is %i to %i \n",mix_start,mix_end);

  TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));
  if (_tree_event == NULL) {
    _tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
    if (_tree_event == NULL) {
      std::cout << " tree fail " << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  //------------------------------------- ROOT ------------------------------------

  Double_t primary_vertex[3];
  Float_t multiplicity_v0[64];
  Float_t event_plane_angle[3];
  Float_t centrality;
  Long64_t mix_events[300];

  _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
  _tree_event->SetBranchAddress("multiplicity_v0", &multiplicity_v0[0]);
  _tree_event->SetBranchAddress("event_plane_psi_v0", &event_plane_angle[0]);
  _tree_event->SetBranchAddress("mixed_events", mix_events);
  _tree_event->SetBranchAddress("centrality_v0m", &centrality);


  //------------------------------------- HISTOGRAMS ------------------------------------
  TFile* fout = new TFile("event_distances.root","RECREATE");
  TH1D* z_vertices_root = new TH1D("Primary_Vertex_root", "Z-vertex (ROOT)", 240, -12, 12);
  TH1D* z_vertices_hdf5 = new TH1D("Primary_Vertex_hdf5", "Z-vertex (hdf5)", 240, -12, 12);
  TH1D* delta_z_vertices = new TH1D("Delta_Primary_Vertex", "#Delta V_z Distribution", 240, -12, 12);

  TH1D* multiplicity_root = new TH1D("multiplicity_root", "multiplicity (ROOT)", 1000, 0, 1000);
  TH1D* multiplicity_hdf5 = new TH1D("Multplicity_hdf5", "multiplicity (hdf5)", 500, 0, 1000);
  TH1D* delta_multiplicity = new TH1D("Delta_multiplicity", "#Delta Multiplicit Distribution", 500, 0, 1000);

  TH1D* flow_root = new TH1D("Flow_root", "Flow (ROOT)", 500, -2, 2);
  TH1D* flow_hdf5 = new TH1D("Flow_hdf5", "Flow (hdf5)", 500, -2, 2);
  TH1D* delta_flow= new TH1D("Delta_Flow", "#Delta Flow Distribution", 500, 0, 4);

  TH1D* centrality_root = new TH1D("centrality_root", "centrality (ROOT)", 200, 0, 100);
  TH1D* centrality_hdf5 = new TH1D("centrality_hdf5", "centrality (hdf5)", 200, 0, 100);
  TH1D* delta_centrality = new TH1D("Delta_centrality", "#Delta centrality Distribution", 200, 0, 100);

  TH2D* N_ME = new TH2D("N_ME", "Distribution No. Mixed Events Passed",300,0,300,500,0,1000);

  //-------------------------- Using low level hdf5 API ------------------------------------------------

  //open hdf5: Define size of data from file, explicitly allocate memory in hdf5 space and array size
  const H5std_string event_ds_name( "event" );
  DataSet event_dataset = h5_file.openDataSet( event_ds_name );
  DataSpace event_dataspace = event_dataset.getSpace();

  //Load the dimensions of dataset from file, to be used in array/hyperslab
  const int event_ndims = event_dataspace.getSimpleExtentNdims();
  hsize_t event_maxdims[event_ndims];
  hsize_t eventdims[event_ndims];
  event_dataspace.getSimpleExtentDims(eventdims, event_maxdims);

  UInt_t NEvent_Vars = eventdims[1];
  fprintf(stderr, "\n%s:%d: n event variables\n", __FILE__, __LINE__, NEvent_Vars);

  //Define array hyperslab will be fed into
  float event_data_out[1][NEvent_Vars];

  //Define hyperslab size and offset in  FILE;
  hsize_t event_offset[2] = {0, 0};
  hsize_t event_count[2] = {1, NEvent_Vars};

  event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");

  //Define the memory dataspace in which to place hyperslab
  const int RANK_OUT = 2; //# of Dimensions
  DataSpace event_memspace( RANK_OUT, eventdims );

  //Define memory offset for hypreslab starting at begining:
  hsize_t event_offset_out[2] = {0};

  //define Dimensions of array, for writing slab to array
  hsize_t event_count_out[2] = {1, NEvent_Vars};

  //define space in memory for hyperslab, then write from file to memory
  event_memspace.selectHyperslab( H5S_SELECT_SET, event_count_out, event_offset_out );
  event_dataset.read( event_data_out, PredType::NATIVE_FLOAT, event_memspace, event_dataspace );
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "event dataset read into array: OK");

  //MONEY MAKING LOOP

  //Loop structure logic:
  //hdf5 file is set to "chunk" data in blocks of 2000. Putting the 
  //mixing layer on the outside, makes the hdf5 read happen much faster
  //root files are optimized for sequential read, so the file is read
  //through 300x sequentially. hdf5-hdf5 would be the fastest for parallel.

  UInt_t nEvents = _tree_event->GetEntries();
  nEvents = 100000;
  for (Long64_t imix = mix_start; imix < mix_end; imix++){
    for(Long64_t ievent = 0; ievent < nEvents ; ievent++){
      _tree_event->GetEntry(ievent);

      fprintf(stderr, "\r%s:%d: mixed event %llu / %llu : Entry %llu / %llu", __FILE__, __LINE__, imix, mix_end,ievent, nEvents);

      size_t mix_event = mix_events[imix];

      if(mix_event < 0 ) continue; //Unpaired events have mix=-999
      event_offset[0]=mix_event;
      event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
      event_dataset.read( event_data_out, PredType::NATIVE_FLOAT, event_memspace, event_dataspace );

      z_vertices_hdf5->Fill(event_data_out[0][0]);
      z_vertices_root->Fill(primary_vertex[2]);
      delta_z_vertices->Fill(TMath::Abs(event_data_out[0][0]-primary_vertex[2]));

      flow_hdf5->Fill(event_data_out[0][2]);
      flow_root->Fill(event_plane_angle[1]);
      delta_flow->Fill(TMath::Abs(event_data_out[0][2] - event_plane_angle[1]));

      if (centrality){
        centrality_hdf5->Fill(event_data_out[0][3]);
        centrality_root->Fill(centrality);
        delta_centrality->Fill(TMath::Abs(event_data_out[0][3] - centrality));
      }
    }
  }//End loop over events

  file->Close();


  z_vertices_root->Write();
  z_vertices_hdf5->Write();
  delta_z_vertices->Write();

  centrality_root->Write();
  centrality_hdf5->Write();
  delta_centrality->Write();

  flow_root->Write();
  flow_hdf5->Write();
  delta_flow->Write();

  fout->Close();
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}
