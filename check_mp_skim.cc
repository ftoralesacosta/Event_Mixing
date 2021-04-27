/*
  This code is a variation of Fernando's Skeleton_Mix_Correlations.cc code, this time for events
*/
// Author: Ivan Chernyshev; Creator of template code: Fernando Torales-Acosta

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
    // if (argc < 6) {
    //     fprintf(stderr,"Batch Syntax is [Gamma-Triggered Paired Root], [Min-Bias HDF5] [Mix Start] [Mix End] [Track Skim GeV]");
    //     exit(EXIT_FAILURE);
    // }
    
  std::cout<<"made it to main"<<std::endl;    
    int dummyc = 1;
    char **dummyv = new char *[1];
    
    dummyv[0] = strdup("main");
    
    TString root_file = (TString)argv[1];
    std::cout << "Opening: " << (TString)argv[1] << std::endl;
    
    // const H5std_string hdf5_file_name(argv[2]);
    // TString hdf5_file = (TString)argv[2];
    // fprintf(stderr,hdf5_file);
    
    size_t mix_start = atoi(argv[2]);
    size_t mix_end = atoi(argv[3]);
    
    int GeV_Track_Skim = atoi(argv[4]);
    std::cout<<"mix start is "<<mix_start<<std::endl;
    std::cout<<"mix end is "<<mix_end<<std::endl;
    fprintf(stderr,"Using %iGeV Track Skimmed from batch Script \n",GeV_Track_Skim);
    
    size_t nmix = 300;
    fprintf(stderr,"Number of Mixed Events: %i \n",nmix);
    
    //Config File ---------------------------------------------------------------------------
    
    //Declaration and Initialize Variables TO BE SET BY [Corr_config.yaml]
    FILE* config = fopen("Corr_config.yaml", "r");
    double DNN_min = 0;
    double DNN_max = 0;
    double pT_min = 0;
    double pT_max = 0;
    double Eta_max = 0;
    double Cluster_min = 0;
    double EcrossoverE_min = 0;
    int Track_Cut_Bit = 0;
    double iso_max = 0;
    double noniso_min = 0;
    double noniso_max = 0;
    double deta_max = 0;
    isolationDet determiner = CLUSTER_ISO_ITS_04;
    int n_eta_bins = 0;
    int n_phi_bins = 0;
    
    // zT & pT bins
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
        if (strcmp(key, "DNN_min") == 0) {
            DNN_min = atof(value);
            std::cout << "DNN_min: " << DNN_min << std::endl; }
        
        else if (strcmp(key, "DNN_max") == 0) {
            DNN_max = atof(value);
            std::cout << "DNN_max: " << DNN_max << std::endl; }
        
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
        
        else if (strcmp(key, "N_Phi_Bins") == 0) {
            n_phi_bins = atoi(value);
            std::cout << "Number of Phi Bins: " << n_phi_bins << std::endl; }
        
        else if (strcmp(key, "N_Eta_Bins") == 0) {
            n_eta_bins = atoi(value);
            std::cout << "Number of Eta Bins: " << n_eta_bins << std::endl; }
        
        else if (strcmp(key, "Track_Cut_Bit") == 0) {
            Track_Cut_Bit = atoi(value);
            std::cout << "Track Cut Bit: " << Track_Cut_Bit << std::endl; }
        
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
    
    TH1D* z_Vertices_individual = new TH1D("Primary_Vertex_root", "Z-vertex (ROOT)", 240, -12, 12);
    TH1D* z_Vertices_hdf5 = new TH1D("Primary_Vertex_hdf5", "Z-vertex (hdf5)", 240, -12, 12);
    TH1D* z_Vertices = new TH1D("Delta_Primary_Vertex", "#Delta V_z Distribution", 240, -12, 12);
    
    TH1D* Multiplicity_individual = new TH1D("Multiplicity_root", "Multiplicity (ROOT)", 1000, 0, 1000);
    TH1D* Multiplicity_hdf5 = new TH1D("Multplicity_hdf5", "Multiplicity (hdf5)", 500, 0, 1000);
    TH1D* Multiplicity = new TH1D("Delta_Multiplicity", "#Delta Multiplicit Distribution", 500, 0, 1000);
    
    TH2D* N_ME = new TH2D("N_ME", "Distribution No. Mixed Events Passed",300,0,300,500,0,1000);

    //TH2D* Signal_pT_Dist = new TH2D("Signal_pT_Dist","Cluster Pt Spectrum For Isolation (its_04) bins 0.55 < DNN < 0.85",59,0.5,30,59,0.5,30);
    //For this example, fill with pt and energy
    
    //ROOT --------------------------------------------------------------------------------------
    
    TFile *file = TFile::Open(root_file);
    
    if (file == NULL) {
        std::cout << " fail" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout<<"check"<<std::endl;
    //file->Print();
    
    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));
    if (_tree_event == NULL) {
      std::cout<<"first TTree attempt failed"<<std::endl;
      _tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
        if (_tree_event == NULL) {
            std::cout << " tree fail " << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    
    //variables
    Double_t primary_vertex[3];
    Float_t multiplicity_v0[64];
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
    
    _tree_event->SetBranchStatus("*mc*", 0);
    
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("multiplicity_v0", multiplicity_v0);
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
    
    _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    _tree_event->SetBranchAddress("cell_e", cell_e);
    
    _tree_event->SetBranchAddress("mixed_events", mix_events);
    
    std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;    
    
    
    //MONEY MAKING LOOP
    Long64_t nentries = _tree_event->GetEntries();
    //TF1 *land = new TF1("land","TMath::Landau(x,[0],[1],0)*17502.3",0,200);
    //land->SetParameters(133.019,48.6502);
    
    for(Long64_t ievent = 0; ievent < nentries ; ievent++){
        _tree_event->GetEntry(ievent);

    	float multiplicity_sum = 0;
    	for (int k = 0; k < 64; k++)  multiplicity_sum += multiplicity_v0[k];
    	//fprintf(stderr,"multp sum = %f\n",multiplicity_sum);
	
	z_Vertices_individual->Fill(primary_vertex[2]);
	Multiplicity_individual->Fill(multiplicity_sum);


	if(ievent % 10000 == 0)
	  std::cout << "Event " << ievent << std::endl;

    } //end loop over events
    //very particular about file names to ease scripting
    // Write to fout
    std::string filepath = argv[1];
    std::string opened_files = "_" + filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
    //std::string rawname = std::string(argv[1]);
    TFile* fout = new TFile(Form("%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu.root", opened_files.c_str(),GeV_Track_Skim,mix_start,mix_end),"RECREATE");
    std::cout<< "Created ROOT file " << Form("%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu.root",opened_files.c_str(),GeV_Track_Skim,mix_start,mix_end) << std::endl;
    
    
    //Write histograms here
    
    z_Vertices->Write();
    Multiplicity->Write();
    z_Vertices_individual->Write();
    z_Vertices_hdf5->Write();
    Multiplicity_individual->Write();
    Multiplicity_hdf5->Write();
    N_ME->Write();
    
    
    fout->Close();
    
    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}
