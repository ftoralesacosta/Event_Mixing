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
  if (argc < 6) {
    fprintf(stderr,"Batch Syntax is [Gamma-Triggered Paired Root], [Min-Bias HDF5] [Mix Start] [Mix End] [jet Skim GeV]");
    exit(EXIT_FAILURE);
  }


  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");

  TString root_file = (TString)argv[1];
  std::cout << "Opening: " << (TString)argv[1] << std::endl;

  const H5std_string hdf5_file_name(argv[2]);
  TString hdf5_file = (TString)argv[2];
  fprintf(stderr,hdf5_file);

  size_t mix_start = atoi(argv[3]);
  size_t mix_end = atoi(argv[4]);

  int GeV_jet_Skim = atoi(argv[5]);
  std::cout<<"mix start is "<<mix_start<<std::endl;
  std::cout<<"mix end is "<<mix_end<<std::endl;
  fprintf(stderr,"Using %iGeV jet Skimmed from batch Script \n",GeV_jet_Skim);

  size_t nmix = 300;
  fprintf(stderr,"Number of Mixed Events: %i \n",nmix);

  //Config File ---------------------------------------------------------------------------

  //Declaration and Initialize Variables TO BE SET BY [Corr_config.yaml]
  FILE* config = fopen("Corr_config.yaml", "r");
  double Cluster_pT_min = 0;
  double Cluster_pT_max = 0;
  double Eta_max = 0;
  double Cluster_min = 0;
  isolationDet determiner = CLUSTER_ISO_ITS_04;
  int n_eta_bins = 0;
  int n_phi_bins = 0;  

  // zT & pT bins, overwritten by config
  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins+1];
  ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;

  int nptbins = 3;
  float* ptbins;
  ptbins = new float[nptbins+1];
  ptbins[0] = 10.0; ptbins[1] = 11; ptbins[2] = 12.5; ptbins[3] = 16;


  //READ CONFIG
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
    if (strcmp(key, "Cluster_pT_min") == 0) {
      Cluster_pT_min = atof(value);
      std::cout << "Cluster_pT_min: " << Cluster_pT_min << std::endl; }

    else if (strcmp(key, "Cluster_pT_max") == 0) {
      Cluster_pT_max = atof(value);
      std::cout << "Cluster_pT_max: " << Cluster_pT_max << std::endl; }

    else if (strcmp(key, "Eta_max") == 0) {
      Eta_max = atof(value);
      std::cout << "Eta_max: " << Eta_max << std::endl;
    }
    else if (strcmp(key, "Cluster_min") == 0) {
      Cluster_min = atof(value);
      std::cout << "Cluster_min: " << Cluster_min << std::endl; }

    else if (strcmp(key, "N_Phi_Bins") == 0) {
      n_phi_bins = atoi(value);
      std::cout << "Number of Phi Bins: " << n_phi_bins << std::endl; }

    else if (strcmp(key, "N_Eta_Bins") == 0) {
      n_eta_bins = atoi(value);
      std::cout << "Number of Eta Bins: " << n_eta_bins << std::endl; }

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

      else {
        std::cout << "ERROR: Cluster_isolation_determinant in configuration file must be \"cluster_iso_tpc_04\", \"cluster_iso_its_04\", \"cluster_frixione_tpc_04_02\", or \"cluster_frixione_its_04_02\"" << std::endl << "Aborting the program" << std::endl;
        exit(EXIT_FAILURE); }
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

  TH1D* Jet_pT_Dist = new TH1D("Jet_pT_Dist","Jet Pt Spectrum",400,0,100 );
  //For this example, fill with pt and energy

  //ROOT --------------------------------------------------------------------------------------

  TFile *file = TFile::Open(root_file);

  if (file == NULL) {
    std::cout << " fail" << std::endl;
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
  /* UInt_t ntrack; */
  /* Float_t track_e[NTRACK_MAX]; */
  /* Float_t track_pt[NTRACK_MAX]; */
  /* Float_t track_eta[NTRACK_MAX]; */
  /* Float_t track_phi[NTRACK_MAX]; */
  /* UChar_t track_quality[NTRACK_MAX]; */

  UInt_t ncluster;
  Float_t cluster_e[NTRACK_MAX];

  Float_t cluster_pt[NTRACK_MAX];
  Float_t cluster_eta[NTRACK_MAX];
  Float_t cluster_phi[NTRACK_MAX]; 
  Float_t cluster_e_max[NTRACK_MAX]; 
  Int_t cluster_ncell[NTRACK_MAX];
  Float_t cluster_tof[NTRACK_MAX]; 
  Int_t cluster_distance_to_bad_channel[NTRACK_MAX];
  Float_t cluster_lambda_square[NTRACK_MAX][2];   
  Float_t cluster_iso_tpc_02[NTRACK_MAX];
  Float_t cluster_iso_tpc_02_ue[NTRACK_MAX];
  Float_t ue_estimate_tpc_const;

  Float_t jet_ak04tpc_pt_raw[NTRACK_MAX];
  Float_t jet_ak04tpc_phi[NTRACK_MAX];
  Float_t jet_ak04tpc_eta[NTRACK_MAX];

  Double_t primary_vertex[3];
  Float_t centrality_v0m[64];

  fprintf(stderr,"Initializing Mixing Branch to %i ME",nmix);
  Long64_t mix_events[300];

  /* _tree_event->SetBranchAddress("ntrack", &ntrack); */
  /* _tree_event->SetBranchAddress("track_e", track_e); */
  /* _tree_event->SetBranchAddress("track_pt", track_pt); */
  /* _tree_event->SetBranchAddress("track_eta", track_eta); */
  /* _tree_event->SetBranchAddress("track_phi", track_phi); */
  /* _tree_event->SetBranchAddress("ncluster", &ncluster); */
  /* _tree_event->SetBranchAddress("cluster_e", cluster_e); */
  /* _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton); */

  _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
  _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
  _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
  _tree_event->SetBranchAddress("cluster_e_max", cluster_e_max);
  _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
  _tree_event->SetBranchAddress("cluster_tof", cluster_tof);
  _tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);
  _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
  _tree_event->SetBranchAddress("cluster_iso_tpc_02",cluster_iso_tpc_02);
  _tree_event->SetBranchAddress("cluster_iso_tpc_02_ue",cluster_iso_tpc_02_ue);
  _tree_event->SetBranchAddress("ue_estimate_tpc_const",&ue_estimate_tpc_const);

  _tree_event->SetBranchAddress("jet_ak04tpc_pt_raw",jet_ak04tpc_pt_raw);
  _tree_event->SetBranchAddress("jet_ak04tpc_phi",jet_ak04tpc_phi);
  _tree_event->SetBranchAddress("jet_ak04tpc_eta",jet_ak04tpc_eta);

  _tree_event->SetBranchAddress("centrality_v0m",centrality_v0m);
  _tree_event->SetBranchAddress("mixed_events", mix_events);

  std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;


  //Using low level hdf5 API -------------------------------------------------------------------------------

  /* See to_hdf5.cc for jet_vars. Useful map below: */
  /* Jet Variariable 0 = jet_ak04tpc_pt_raw */
  /* Jet Variariable 1 = jet_ak04tpc_eta_raw */
  /* Jet Variariable 2 = jet_ak04tpc_phi */

  //open hdf5: Define size of data from file, explicitly allocate memory in hdf5 space and array size
  const H5std_string jet_ds_name( "jet" );
  H5File h5_file( hdf5_file_name, H5F_ACC_RDONLY ); //hdf5_file_name from argv[2]
  DataSet jet_dataset = h5_file.openDataSet( jet_ds_name );
  DataSpace jet_dataspace = jet_dataset.getSpace();

  //Load the dimensions of dataset from file, to be used in array/hyperslab
  const int jet_ndims = jet_dataspace.getSimpleExtentNdims();
  hsize_t jet_maxdims[jet_ndims];
  hsize_t jetdims[jet_ndims];
  jet_dataspace.getSimpleExtentDims(jetdims, jet_maxdims);
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
     For example, To obtain data for event 68, set the
     offset's to {68, njet_max, Njet_Vars}.
     */


  jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");

  //Define the memory dataspace in which to place hyperslab
  const int RANK_OUT = 3; //# of Dimensions
  DataSpace jet_memspace( RANK_OUT, jetdims );

  //Define memory offset for hypreslab starting at begining:
  hsize_t jet_offset_out[3] = {0};

  //define Dimensions of array, for writing slab to array
  hsize_t jet_count_out[3] = {1, njet_max, Njet_Vars};

  //define space in memory for hyperslab, then write from file to memory
  jet_memspace.selectHyperslab( H5S_SELECT_SET, jet_count_out, jet_offset_out );
  jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "jet dataset read into array: OK");


  //MONEY MAKING LOOP
  Long64_t nentries = _tree_event->GetEntries();    

  for(Long64_t ievent = 0; ievent < nentries ; ievent++){     
    _tree_event->GetEntry(ievent);

    //Cuts/Variables from the ROOT file go here

    for (Long64_t imix = mix_start; imix < mix_end+1; imix++){
      Long64_t mix_event = mix_events[imix];
      fprintf(stderr,"\n %s:%d: Mixed event = %lu",__FILE__,__LINE__,mix_event);

      //if (mix_event == ievent) continue; //not needed for gamma-MB pairing: Different Triggers
      if(mix_event >= 9999999) continue;  

      //adjust offset for next mixed event
      jet_offset[0]=mix_event;
      jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
      jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );

      /* IMPORTANT */
      /* jet_data_out[0][ijet][n_variables] */
      /* Jet Variariable 0 = jet_ak04tpc_pt_raw */
      /* Jet Variariable 1 = jet_ak04tpc_eta_raw */
      /* Jet Variariable 2 = jet_ak04tpc_phi */
      
      for (ULong64_t ijet = 0; ijet < njet_max; ijet++) {
        if (std::isnan(jet_data_out[0][ijet][0])) continue;
        if (jet_data_out[0][ijet][0] < 5.) continue; //greater than 5.0 GeV
        Jet_pT_Dist->Fill(jet_data_out[0][ijet][0]);

      }//end loop over jets
    }//end loop over mixed events
  } //end loop over events

  // Write to fout    
  size_t lastindex = std::string(root_file).find_last_of("."); 
  std::string rawname = std::string(root_file).substr(0, lastindex);
  //std::string rawname = std::string(argv[1]);
  TFile* fout = new TFile(Form("%s_%luGeVjets_Correlation_%1.1lu_to_%1.1lu.root",rawname.data(),GeV_jet_Skim,mix_start,mix_end),"RECREATE");    
  //Write histograms here    

  Jet_pT_Dist->Write();    

  fout->Close();     

  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}
