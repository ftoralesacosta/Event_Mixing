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
  if (argc < 2) {
    fprintf(stderr,"Please Indicate HDF5 fil");
    exit(EXIT_FAILURE);
  }
    
    
  int dummyc = 1;
  char **dummyv = new char *[1];
    
  dummyv[0] = strdup("main");
    
  const H5std_string hdf5_file_name(argv[1]);
  TString hdf5_file = (TString)argv[1];
  fprintf(stderr,hdf5_file);
    
  //HISTOGRAMS
  TCanvas canvas("canvas", "");
    
  TH1D* z_Vertices_individual = new TH1D("Primary_Vertex_root", "Z-vertex (ROOT)", 240, -12, 12);
  TH1D* z_Vertices_hdf5 = new TH1D("Primary_Vertex_hdf5", "Z-vertex (hdf5)", 240, -12, 12);
  TH1D* z_Vertices = new TH1D("Delta_Primary_Vertex", "#Delta V_z Distribution", 240, -12, 12);
    
  TH1D* Multiplicity_individual = new TH1D("Multiplicity_root", "Multiplicity (ROOT)", 1000, 0, 1000);
  TH1D* Multiplicity_hdf5 = new TH1D("Multplicity_hdf5", "Multiplicity (hdf5)", 500, 0, 1000);
  TH1D* Multiplicity = new TH1D("Delta_Multiplicity", "#Delta Multiplicit Distribution", 500, 0, 1000);
    
  TH2D* N_ME = new TH2D("N_ME", "Distribution No. Mixed Events Passed",300,0,300,500,0,1000);
    
  TH1D* jet_pT_hdf5 = new TH1D("jet_pt_distribution", "Jet p_{T} distribution (HDF5)", 60, -15, 15);
  TH1D* jet_eta_hdf5 = new TH1D("jet_eta_distribution", "Jet #eta distribution (HDF5)", 60, -0.8, 0.8);
  TH1D* jet_phi_hdf5 = new TH1D("jet_phi_distribution", "Jet #phi distribution (HDF5)", 60, -3.1415926, 3.1415926);
  TH1D* jet_pTD_hdf5 = new TH1D("jet_ptD_distribution", "Jet p_{T}D distribution (HDF5)", 60, -1, 1);
  TH1D* jet_Multiplicity_hdf5 = new TH1D("jet_ptD_distribution", "Jet p_{T}D distribution (HDF5)", 60, -15, 15);


  //Using low level hdf5 API -------------------------------------------------------------------------------
    

  //open hdf5: Define size of data from file, explicitly allocate memory in hdf5 space and array size
  const H5std_string event_ds_name( "event" );
  H5File h5_file( hdf5_file_name, H5F_ACC_RDONLY ); //hdf5_file_name from argv[2]
  DataSet event_dataset = h5_file.openDataSet( event_ds_name );
  DataSpace event_dataspace = event_dataset.getSpace();
  //Load the dimensions of dataset from file, to be used in array/hyperslab
  const int event_ndims = event_dataspace.getSimpleExtentNdims();
  hsize_t event_maxdims[event_ndims];
  hsize_t eventdims[event_ndims];
  event_dataspace.getSimpleExtentDims(eventdims, event_maxdims);
  //UInt_t nevent_max = eventdims[1];
  UInt_t NEvent_Vars = eventdims[1];
  fprintf(stderr, "\n%s:%d: n track variables\n", __FILE__, __LINE__, NEvent_Vars);
    
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
  std::cout << "Made it to line " << __LINE__  << std::endl;
  event_dataset.read( event_data_out, PredType::NATIVE_FLOAT, event_memspace, event_dataspace );
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "event dataset read into array: OK");
    
  // Get jet
  const H5std_string jet_ds_name( "jet" );
  H5File h5_file_jet( hdf5_file_name, H5F_ACC_RDONLY ); //hdf5_file_name from argv[2]
  DataSet jet_dataset = h5_file_jet.openDataSet( jet_ds_name );
  DataSpace jet_dataspace = jet_dataset.getSpace();
    
  //test
  //Load the dimensions of dataset from file, to be used in array/hyperslab
  const int jet_ndims = jet_dataspace.getSimpleExtentNdims();
  hsize_t jet_maxdims[jet_ndims];
  hsize_t jetdims[jet_ndims];
  jet_dataspace.getSimpleExtentDims(jetdims, jet_maxdims);

  UInt_t nEvents = jetdims[0];
  UInt_t njet_max = jetdims[1];
  UInt_t Njet_Vars = jetdims[2];
  fprintf(stderr, "\n%s:%d: n jet variables\n", __FILE__, __LINE__, Njet_Vars);
    
  //Define array hyperslab will be fed into
  float jet_data_out[1][njet_max][Njet_Vars];
    
  //Define hyperslab size and offset in  FILE;
  hsize_t jet_offset[3] = {0, 0, 0};
  hsize_t jet_count[3] = {1, njet_max, Njet_Vars};
    
  /*
     The Offset is how we iterate over the entire hdf5 file.
     For example, To obtain data for jet 68, set the
     offset's to {68, njet_max, Njet_Vars}.
  */
    
    
  jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");
    
  //Define the memory dataspace in which to place hyperslab
  const int RANK_OUT_jet = 3; //# of Dimensions
  DataSpace jet_memspace( RANK_OUT_jet, jetdims );
    
  //Define memory offset for hypreslab starting at begining:
  hsize_t jet_offset_out[3] = {0};
    
  //define Dimensions of array, for writing slab to array
  hsize_t jet_count_out[3] = {1, njet_max, Njet_Vars};
  std::cout<< "Made it to line 398" << std::endl;

  //define space in memory for hyperslab, then write from file to memory
  jet_memspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
  jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "jet dataset read into array: OK");

  //MONEY MAKING LOOP
    
  for(Long64_t ievent = 0; ievent < nEvents ; ievent++){
    //for(Long64_t ievent = 0; ievent < 200; ievent++){                                                                                                                         
    
    fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nEvents);

    //Cuts/Variables from the ROOT file go here
            
      //if (mix_event == ievent) continue; //not needed for gamma-MB pairing: Different Triggers
            
      //adjust offset for next mixed event
      event_offset[0]=ievent;
      event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
      event_dataset.read( event_data_out, PredType::NATIVE_FLOAT, event_memspace, event_dataspace );

      jet_offset[0]=ievent;      //adjust offset for next mixed event                                                                                                    
      jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
      jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );

      z_Vertices_hdf5->Fill(event_data_out[0][0]);
      Multiplicity_hdf5->Fill(event_data_out[0][3]);
    
            
      // Loop over jets
      for (Long64_t ijet = 0; ijet < njet_max; ijet++) {
	double jet_pT = -9000;
	double jet_phi = -9000;
	double jet_eta = -9000;
	double jet_pTD = -9000;
	double jet_Multiplicity = -9000;
                
	jet_pT = jet_data_out[0][ijet][0];
	jet_phi = jet_data_out[0][ijet][2];
	jet_eta = jet_data_out[0][ijet][1];
	jet_pTD = jet_data_out[0][ijet][3];
	jet_Multiplicity = jet_data_out[0][ijet][4];
                
	// Cut on pT jet = 0 and pT jet NaN
	if(TMath::IsNaN(jet_pT)) continue;
	if(TMath::IsNaN(jet_eta)) continue;
	if(TMath::IsNaN(jet_phi)) continue;
	if(TMath::IsNaN(jet_pTD)) continue;
	if(TMath::IsNaN(jet_Multiplicity)) continue;
	//if(not(jet_data_out[0][ijet][0] == 0)) {continue;}
                
	// Fill the histogram
	jet_pT_hdf5->Fill(jet_pT);
	jet_eta_hdf5->Fill(jet_eta);
	jet_phi_hdf5->Fill(jet_phi);
	jet_pTD_hdf5->Fill(jet_pTD);
	jet_Multiplicity_hdf5->Fill(jet_Multiplicity);
      }
  }//End loop over events

  // if(ievent % 10000 == 0)
  //     std::cout << "Event " << ievent << std::endl;
  
  
  //very particular about file names to ease scripting
  // jet_Multiplicity_hdf5->Write();
  
  TCanvas* c = new TCanvas();
  jet_pT_hdf5->Draw();
  c->SaveAs("HDF5_jetpT_Distribution.pdf");
  c->Clear();
    
  jet_eta_hdf5->Draw();
  c->SaveAs("HDF5_jeteta_Distribution.pdf");
  c->Clear();
    
  jet_phi_hdf5->Draw();
  c->SaveAs("HDF5_jetphi_Distribution.pdf");
  c->Clear();
    
  jet_pTD_hdf5->Draw();
  c->SaveAs("HDF5_jetMultiplicity_Distribution.pdf");
  c->Clear();
    
  //fout->Close();
    
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}
