#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <emcal.h>
#include <H5Cpp.h>
#define NTRACK_MAX (1U << 14)

// This is chosen to be the CPU L2 cache size, which should exceed 512
// kB for many years now
#ifndef HDF5_DEFAULT_CACHE
#define HDF5_DEFAULT_CACHE (512 * 1024)
#endif // HDF5_DEFAULT_CACHE

#ifndef HDF5_USE_DEFLATE
#define HDF5_USE_DEFLATE
#endif // HDF5_USE_DEFLATE

// Rank of the tensor written in HDF5, rank 3 being (index_event,
// index_track, index_properties)
#define RANK 3
#define Event_RANK 2


void find_ntrack_ncluster_max(char *argv_first[], char *argv_last[], UInt_t &nevent_max, UInt_t &ntrack_max, UInt_t &ncluster_max, UInt_t &njet_ak04tpc_max, UInt_t &njet_ak02tpc_max)
{
  for (char **p = argv_first; p != argv_last; p++) {
    // Cautious opening of the TTree, capturing all modes of
    // failure, and keep the TDirectoryFile (to be deleted later)
    // to avoid memory leak

    std::string filename = *p;

    // hard-code the ones we know about to save time
    if (filename == "/global/project/projectdirs/alice/NTuples/PbPb/18q_pass3_cluster15_celltrack.root") {
      nevent_max = 585601;
      ntrack_max = 3786;
      ncluster_max = 2031;
      njet_ak04tpc_max = 52;
      njet_ak02tpc_max = 167;
    } else {
      TFile *file = TFile::Open(*p);
      std::string auxfilename = filename.replace(filename.find(".root"), 5, "_AUX.root");
      TFile *auxfile = TFile::Open((TString) auxfilename);

      if (file == NULL) {
        std::cout << " fail" << std::endl;
        exit(EXIT_FAILURE);
      }
      file->Print();

      TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));
      TTree *auxtree = dynamic_cast<TTree*> (auxfile->Get("ntupleaux"));

      if (_tree_event == NULL) {
        _tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
        if (_tree_event == NULL) {
          std::cout << " fail " << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      UInt_t ntrack;
      UInt_t ncluster;
      UInt_t njet_ak04tpc;
      UInt_t njet_ak02tpc;
      if (nevent_max == 0)
        nevent_max = UInt_t(_tree_event->GetEntries());

      _tree_event->SetBranchAddress("ntrack", &ntrack);
      _tree_event->SetBranchAddress("ncluster", &ncluster);
      _tree_event->SetBranchAddress("njet_ak04tpc", &njet_ak04tpc);
      auxtree->SetBranchAddress("njet_ak02tpc", &njet_ak02tpc);

      fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Obtaining ntrack, ncluster, and njet max for hdf5 file");

      for (Long64_t i = 0; i < nevent_max; i++) {
        /* for (Long64_t i = 0; i < 10000; i++) { */
        _tree_event->GetEntry(i);
        auxtree->GetEntry(i);

        ntrack_max = std::max(ntrack_max, ntrack);
        ncluster_max = std::max(ncluster_max, ncluster);
        njet_ak04tpc_max = std::max(njet_ak04tpc_max, njet_ak04tpc);
        njet_ak02tpc_max = std::max(njet_ak02tpc_max, njet_ak02tpc);
        fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, i, nevent_max);
      }

      fprintf(stderr, "\n");

      // Fully delete everything
      _tree_event->Delete();
      auxfile->Close();
      delete auxfile;
      file->Close();
      delete file;
    }
  }
}

void write_track_cluster(H5::DataSet &event_data_set, H5::DataSet &track_data_set, H5::DataSet &cluster_data_set,
                         H5::DataSet &jet_ak04tpc_data_set, H5::DataSet &jet_ak02tpc_data_set,
                         hsize_t *event_offset, hsize_t *offset,
                         const hsize_t *event_dim_extend, const hsize_t *track_dim_extend, const hsize_t *cluster_dim_extend,
                         const hsize_t *jet_ak04tpc_dim_extend, const hsize_t *jet_ak02tpc_dim_extend,
                         const UInt_t nevent_max, const UInt_t ntrack_max, const UInt_t ncluster_max,
                         const UInt_t njet_ak04tpc_max, const UInt_t njet_ak02tpc_max, const UInt_t block_size,
                         char *argv_first[], char *argv_last[])
{
  for (char **p = argv_first; p != argv_last; p++) {
    TFile *file = TFile::Open(*p);
    std::string filename = *p;
    std::string auxfilename = filename.replace(filename.find(".root"), 5, "_AUX.root");
    TFile *auxfile = TFile::Open((TString)auxfilename);

    if (file == NULL) {
      std::cout << " fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();

    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));
    TTree *auxtree = dynamic_cast<TTree*> (auxfile->Get("ntupleaux"));

    if (_tree_event == NULL) {
      _tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
      if (_tree_event == NULL) {
        std::cout << " fail " << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    //Events
    UInt_t nevent;
    Float_t centrality;
    Float_t ue_estimate_its_const;
    Float_t ue_estimate_tpc_const;
    Bool_t is_pileup_from_spd_5_08;
    std::vector<Double_t> primary_vertex(3, NAN);
    std::vector<Float_t> multiplicity_v0(64, NAN);//64 channels for v0 detector, to be summed
    std::vector<Float_t> event_plane_angle(3, NAN); //directed/eliptic/triangular


    //Tracks
    UInt_t ntrack;
    std::vector<Float_t> track_e(ntrack_max, NAN);
    std::vector<Float_t> track_pt(ntrack_max, NAN);
    std::vector<Float_t> track_eta(ntrack_max, NAN);
    std::vector<Float_t> track_phi(ntrack_max, NAN);
    std::vector<UChar_t> track_quality(ntrack_max, 64);
    std::vector<Float_t> track_eta_emcal(ntrack_max, NAN);
    std::vector<Float_t> track_phi_emcal(ntrack_max, NAN);
    std::vector<UChar_t> track_tpc_ncluster(ntrack_max, 0);
    //std::vector<Float_t> track_tpc_chi_square(ntrack_max, NAN);
    std::vector<Float_t> track_dca_xy(ntrack_max, NAN);
    std::vector<Float_t> track_dca_z(ntrack_max, NAN);

    //Clusters
    UInt_t ncluster;
    Float_t cluster_e[NTRACK_MAX];
    Float_t cluster_pt[NTRACK_MAX];
    Float_t cluster_eta[NTRACK_MAX];
    Float_t cluster_phi[NTRACK_MAX];
    Float_t cluster_lambda_square[NTRACK_MAX][2];
    Float_t cluster_e_max[NTRACK_MAX];
    Float_t cluster_e_cross[NTRACK_MAX];
    Float_t cluster_iso_tpc_02[NTRACK_MAX];
    Float_t cluster_iso_tpc_03[NTRACK_MAX];
    Float_t cluster_iso_tpc_04[NTRACK_MAX];
    Float_t cluster_iso_its_02[NTRACK_MAX];
    Float_t cluster_iso_its_03[NTRACK_MAX];
    Float_t cluster_iso_its_04[NTRACK_MAX];
    Float_t cluster_frixione_tpc_04_02[NTRACK_MAX];
    Float_t cluster_frixione_its_04_02[NTRACK_MAX];
    Float_t cluster_s_nphoton[NTRACK_MAX][4];
    unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
    Int_t cluster_ncell[NTRACK_MAX];
    UShort_t  cluster_cell_id_max[NTRACK_MAX];
    Float_t cell_e[17664];
    Float_t cluster_distance_to_bad_channel[NTRACK_MAX];
    UChar_t cluster_nlocal_maxima[NTRACK_MAX];
    Float_t cluster_tof[NTRACK_MAX];
    Float_t cluster_iso_its_02_ue[NTRACK_MAX];
    Float_t cluster_iso_its_03_ue[NTRACK_MAX];
    Float_t cluster_iso_its_04_ue[NTRACK_MAX];
    Float_t cluster_iso_tpc_02_ue[NTRACK_MAX];
    Float_t cluster_iso_tpc_03_ue[NTRACK_MAX];
    Float_t cluster_iso_tpc_04_ue[NTRACK_MAX];
    Float_t cluster_5x5all[NTRACK_MAX];

    //Jets (R=0.4)
    UInt_t njet_ak04tpc;
    std::vector<Float_t> jet_ak04tpc_pt_raw(njet_ak04tpc_max, NAN);
    std::vector<Float_t> jet_ak04tpc_eta_raw(njet_ak04tpc_max, NAN);
    std::vector<Float_t> jet_ak04tpc_phi(njet_ak04tpc_max, NAN);
    std::vector<Float_t> jet_ak04tpc_area(njet_ak04tpc_max, NAN);
    std::vector<UShort_t> jet_ak04tpc_multiplicity_raw(njet_ak04tpc_max, 0);

    // Jets (R=0.2)
    UInt_t njet_ak02tpc;
    std::vector<Float_t> jet_ak02tpc_pt_raw(njet_ak02tpc_max, NAN);
    std::vector<Float_t> jet_ak02tpc_eta_raw(njet_ak02tpc_max, NAN);
    std::vector<Float_t> jet_ak02tpc_phi(njet_ak02tpc_max, NAN);
    std::vector<Float_t> jet_ak02tpc_area(njet_ak02tpc_max, NAN);
    std::vector<UShort_t> jet_ak02tpc_multiplicity_raw(njet_ak02tpc_max, 0);

    _tree_event->SetBranchAddress("primary_vertex", &primary_vertex[0]);
    _tree_event->SetBranchAddress("multiplicity_v0", &multiplicity_v0[0]);
    _tree_event->SetBranchAddress("event_plane_psi_v0", &event_plane_angle[0]);
    _tree_event->SetBranchAddress("centrality_v0m", &centrality);
    _tree_event->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);
    _tree_event->SetBranchAddress("ue_estimate_its_const", &ue_estimate_its_const);
    _tree_event->SetBranchAddress("ue_estimate_tpc_const", &ue_estimate_tpc_const);
    /* _tree_event->SetBranchAddress("mixed_events",mix_events); */

    _tree_event->SetBranchAddress("ntrack", &ntrack);
    _tree_event->SetBranchAddress("track_e", &track_e[0]);
    _tree_event->SetBranchAddress("track_pt", &track_pt[0]);
    _tree_event->SetBranchAddress("track_eta", &track_eta[0]);
    _tree_event->SetBranchAddress("track_phi", &track_phi[0]);
    _tree_event->SetBranchAddress("track_quality", &track_quality[0]);
    _tree_event->SetBranchAddress("track_eta_emcal", &track_eta_emcal[0]);
    _tree_event->SetBranchAddress("track_phi_emcal", &track_phi_emcal[0]);
    _tree_event->SetBranchAddress("track_tpc_ncluster", &track_tpc_ncluster[0]);
    //_tree_event->SetBranchAddress("track_tpc_chi_square", &track_tpc_chi_square[0]);
    _tree_event->SetBranchAddress("track_dca_xy", &track_dca_xy[0]);
    _tree_event->SetBranchAddress("track_dca_z", &track_dca_z[0]);

    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_e", cluster_e);
    _tree_event->SetBranchAddress("cluster_e_max", cluster_e_max);
    _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
    _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
    _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
    _tree_event->SetBranchAddress("cluster_iso_tpc_02", cluster_iso_tpc_02);
    _tree_event->SetBranchAddress("cluster_iso_tpc_03", cluster_iso_tpc_03);
    _tree_event->SetBranchAddress("cluster_iso_tpc_04", cluster_iso_tpc_04);
    _tree_event->SetBranchAddress("cluster_iso_its_02", cluster_iso_its_02);
    _tree_event->SetBranchAddress("cluster_iso_its_03", cluster_iso_its_03);
    _tree_event->SetBranchAddress("cluster_iso_its_04", cluster_iso_its_04);
    _tree_event->SetBranchAddress("cluster_frixione_tpc_04_02", cluster_frixione_tpc_04_02);
    _tree_event->SetBranchAddress("cluster_frixione_its_04_02", cluster_frixione_its_04_02);
    _tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);
    _tree_event->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);
    _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    _tree_event->SetBranchAddress("cell_e", cell_e);
    _tree_event->SetBranchAddress("cluster_tof", cluster_tof);
    _tree_event->SetBranchAddress("cluster_iso_its_02_ue", cluster_iso_its_02_ue);
    _tree_event->SetBranchAddress("cluster_iso_its_03_ue", cluster_iso_its_03_ue);
    _tree_event->SetBranchAddress("cluster_iso_its_04_ue", cluster_iso_its_04_ue);
    _tree_event->SetBranchAddress("cluster_iso_tpc_02_ue", cluster_iso_tpc_02_ue);
    _tree_event->SetBranchAddress("cluster_iso_tpc_03_ue", cluster_iso_tpc_03_ue);
    _tree_event->SetBranchAddress("cluster_iso_tpc_04_ue", cluster_iso_tpc_04_ue);
    auxtree->SetBranchAddress("cluster_5x5all", cluster_5x5all);


    _tree_event->SetBranchAddress("njet_ak04tpc", &njet_ak04tpc);
    _tree_event->SetBranchAddress("jet_ak04tpc_pt_raw", &jet_ak04tpc_pt_raw[0]);
    _tree_event->SetBranchAddress("jet_ak04tpc_eta_raw", &jet_ak04tpc_eta_raw[0]);
    _tree_event->SetBranchAddress("jet_ak04tpc_phi", &jet_ak04tpc_phi[0]);//for some reason, phi_raw is nan, but phi is not. Must be a mix up.
    _tree_event->SetBranchAddress("jet_ak04tpc_area", &jet_ak04tpc_area[0]);
    _tree_event->SetBranchAddress("jet_ak04tpc_multiplicity_raw", &jet_ak04tpc_multiplicity_raw[0]);

    auxtree->SetBranchAddress("njet_ak02tpc", &njet_ak02tpc);
    auxtree->SetBranchAddress("jet_ak02tpc_pt_raw", &jet_ak02tpc_pt_raw[0]);
    auxtree->SetBranchAddress("jet_ak02tpc_eta", &jet_ak02tpc_eta_raw[0]);
    auxtree->SetBranchAddress("jet_ak02tpc_phi", &jet_ak02tpc_phi[0]);//for some reason, phi_raw is nan, but phi is not. Must be a mix up.
    auxtree->SetBranchAddress("jet_ak02tpc_area", &jet_ak02tpc_area[0]);
    auxtree->SetBranchAddress("jet_ak02tpc_multiplicity_raw", &jet_ak02tpc_multiplicity_raw[0]);

    //Changes here should be also be done on line 596
    static const size_t event_row_size = 7;
    static const size_t track_row_size = 10;
    static const size_t cluster_row_size = 32;
    static const size_t jet_row_size = 5;

    fprintf(stderr, "\n %d: HDF5 Chunk Cache Size = %u\n", __LINE__, block_size);

    std::vector<float> event_data (block_size * event_row_size, NAN); //4 variables, multp, vertx, 2 event angles
    std::vector<float> track_data(block_size * ntrack_max * track_row_size, NAN);
    std::vector<float> cluster_data(block_size * ncluster_max * cluster_row_size, NAN);
    std::vector<float> jet_ak04tpc_data(block_size * njet_ak04tpc_max * jet_row_size, NAN);
    std::vector<float> jet_ak02tpc_data(block_size * njet_ak02tpc_max * jet_row_size, NAN);

    for (Long64_t i = 0; i < nevent_max; i++) {
      /* for (Long64_t i = 0; i <10000; i++) { */
      /* for (Long64_t i = 0; i <4001; i++) { */
      _tree_event->GetEntry(i);
      auxtree->GetEntry(i);

      int iblock = i % block_size;
      //writing to file is done every [block_size] number of events
      //this variable keeps track of the current increment within a block
      //, as opposed to [i] which is looping through all events

      //fprintf(stderr,"\n %d: iblock = %i \n",__LINE__,iblock);

      float multiplicity_sum = 0;
      for (int k = 0; k < 64; k++) multiplicity_sum += multiplicity_v0[k];
      event_data[iblock * event_row_size + 0] = primary_vertex[2]; //xyz, choose 3rd element, z
      event_data[iblock * event_row_size + 1] = multiplicity_sum;
      event_data[iblock * event_row_size + 2] = event_plane_angle[1]; //elliptic flow
      event_data[iblock * event_row_size + 3] = centrality;
      event_data[iblock * event_row_size + 4] = is_pileup_from_spd_5_08;
      event_data[iblock * event_row_size + 5] = ue_estimate_its_const;
      event_data[iblock * event_row_size + 6] = ue_estimate_tpc_const;

      //event_data[iblock*nEventVariables + N] = event_plane_angle[0]; //directed flow

      for (Long64_t j = 0; j < ntrack; j++) {
        track_data[(iblock * ntrack_max + j)*track_row_size + 0] = track_e[j];
        track_data[(iblock * ntrack_max + j)*track_row_size + 1] = track_pt[j];
        track_data[(iblock * ntrack_max + j)*track_row_size + 2] = track_eta[j];
        track_data[(iblock * ntrack_max + j)*track_row_size + 3] = track_phi[j];
        track_data[(iblock * ntrack_max + j)*track_row_size + 4] = track_quality[j];
        track_data[(iblock * ntrack_max + j)*track_row_size + 5] = track_eta_emcal[j];
        track_data[(iblock * ntrack_max + j)*track_row_size + 6] = track_phi_emcal[j];
        track_data[(iblock * ntrack_max + j)*track_row_size + 7] = track_tpc_ncluster[j];
        track_data[(iblock * ntrack_max + j)*track_row_size + 8] = track_dca_xy[j];
        track_data[(iblock * ntrack_max + j)*track_row_size + 9] = track_dca_z[j];
      }

      for (Long64_t n = 0; n < ncluster; n++) {
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 0] = cluster_e[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 1] = cluster_pt[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 2] = cluster_eta[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 3] = cluster_phi[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 4] = cluster_lambda_square[n][0];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 5] = cluster_e_max[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 6] = cluster_e_cross[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 7] = cluster_iso_tpc_02[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 8] = cluster_iso_tpc_03[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 9] = cluster_iso_tpc_04[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 10] = cluster_iso_its_02[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 11] = cluster_iso_its_03[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 12] = cluster_iso_its_04[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 13] = cluster_frixione_tpc_04_02[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 14] = cluster_frixione_its_04_02[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 15] = cluster_s_nphoton[n][0];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 16] = cluster_mc_truth_index[n][0]; //32. 0 is placholder.FIXME
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 17] = cluster_ncell[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 18] = cluster_cell_id_max[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 19] = cell_e[0]; //length 17664. Placholder
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 20] = cluster_distance_to_bad_channel[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 21] = cluster_nlocal_maxima[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 22] = cluster_tof[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 23] = cluster_iso_its_02_ue[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 24] = cluster_iso_its_03_ue[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 25] = cluster_iso_its_04_ue[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 26] = cluster_iso_tpc_02_ue[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 27] = cluster_iso_tpc_03_ue[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 28] = cluster_iso_tpc_04_ue[n];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 29] = cluster_lambda_square[n][1];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 30] = cluster_s_nphoton[n][1];
        cluster_data[(iblock * ncluster_max + n)*cluster_row_size + 31] = cluster_5x5all[n];
      }

      for (Long64_t j = 0; j < njet_ak04tpc; j++) {
        jet_ak04tpc_data[(iblock * njet_ak04tpc_max + j)*jet_row_size + 0] = jet_ak04tpc_pt_raw[j];
        jet_ak04tpc_data[(iblock * njet_ak04tpc_max + j)*jet_row_size + 1] = jet_ak04tpc_eta_raw[j];
        jet_ak04tpc_data[(iblock * njet_ak04tpc_max + j)*jet_row_size + 2] = jet_ak04tpc_phi[j];
        jet_ak04tpc_data[(iblock * njet_ak04tpc_max + j)*jet_row_size + 3] = jet_ak04tpc_area[j];
        jet_ak04tpc_data[(iblock * njet_ak04tpc_max + j)*jet_row_size + 4] = jet_ak04tpc_multiplicity_raw[j];
      }

      for (Long64_t j = 0; j < njet_ak02tpc; j++) {
        jet_ak02tpc_data[(iblock * njet_ak02tpc_max + j)*jet_row_size + 0] = jet_ak02tpc_pt_raw[j];
        jet_ak02tpc_data[(iblock * njet_ak02tpc_max + j)*jet_row_size + 1] = jet_ak02tpc_eta_raw[j];
        jet_ak02tpc_data[(iblock * njet_ak02tpc_max + j)*jet_row_size + 2] = jet_ak02tpc_phi[j];
        jet_ak02tpc_data[(iblock * njet_ak02tpc_max + j)*jet_row_size + 3] = jet_ak02tpc_area[j];
        jet_ak02tpc_data[(iblock * njet_ak02tpc_max + j)*jet_row_size + 4] = jet_ak02tpc_multiplicity_raw[j];
      }

      bool print_stuff = true;
      if (print_stuff) {

        if (i % 100000 == 0) {

          std::cout << "Event Info:\n";
          std::cout << "z-vertex = " << primary_vertex[2] << "\n";
          std::cout << "multilplicity = " << multiplicity_sum << "\n";
          std::cout << "v2 flow = " << event_plane_angle[1] << "\n";
          std::cout << "centrality = " << centrality << "\n";
          std::cout << "is_pileup_from_spd_5_08 = " << is_pileup_from_spd_5_08 << "\n";
          std::cout << "ue_estimate_its_const = " << ue_estimate_its_const << "\n";
          std::cout << "\n";

          std::cout << "Track Info:\n";
          std::cout << "track E = " << track_e[0] << "\n";
          std::cout << "track pT = " << track_pt[0] << "\n";
          std::cout << "track Eta = " << track_eta[0] << "\n";
          std::cout << "track Phi = " << track_phi[0] << "\n";
          std::cout << "track quality = " << track_quality[0] << "\n";
          std::cout << "track Eta_EMCal = " << track_eta_emcal[0] << "\n";
          std::cout << "track Phi_EMCal = " << track_phi_emcal[0] << "\n";
          std::cout << "track TPC NCluster = " << track_tpc_ncluster[0] << "\n";
          std::cout << "track DCA XY = " << track_dca_xy[0] << "\n";
          std::cout << "track DCA Z = " << track_dca_z[0] << "\n";
          std::cout << "\n";

          std::cout << "Cluster Info: \n";
          std::cout << "cluster E = " << cluster_e[0] << std::endl;
          std::cout << "cluster Eta = " << cluster_eta[0] << std::endl;
          std::cout << "cluster Phi = " << cluster_phi[0] << std::endl;
          std::cout << "cluster Pt = " << cluster_pt[0] << std::endl;
          std::cout << "cluster Lamda_0^2 = " << cluster_lambda_square[0][0] << std::endl;
          //std::cout<<cluster_e_cross[0]<<std::endl;
          std::cout << "\n";

          std::cout << "Jet R=0.4 Info: \n";
          std::cout << "jet pT raw = " << jet_ak04tpc_pt_raw[0] << std::endl;
          std::cout << "jet Eta raw = " << jet_ak04tpc_eta_raw[0] << std::endl;
          std::cout << "jet Phi = " << jet_ak04tpc_phi[0] << std::endl;
          std::cout << "jet area = " << jet_ak04tpc_area[0] << std::endl;
          std::cout << "jet Multp raw = " << jet_ak04tpc_multiplicity_raw[0] << std::endl;
          std::cout << "\n";

          std::cout << "Jet R=0.2 Info: \n";
          std::cout << "jet pT raw = " << jet_ak02tpc_pt_raw[0] << std::endl;
          std::cout << "jet Eta raw = " << jet_ak02tpc_eta_raw[0] << std::endl;
          std::cout << "jet Phi = " << jet_ak02tpc_phi[0] << std::endl;
          std::cout << "jet area = " << jet_ak02tpc_area[0] << std::endl;
          std::cout << "jet Multp raw = " << jet_ak02tpc_multiplicity_raw[0] << std::endl;
          std::cout << "\n";
        }
      }//print

      if (iblock == (block_size - 1)) { //writes 1 block (2000 events) at a time. Faster/less memory

        //write first event
        if (event_offset[0] == 0) {
          event_data_set.write(&event_data[0], H5::PredType::NATIVE_FLOAT);
          std::fill(event_data.begin(), event_data.end(), NAN);
        }

        //extend dataset otherwise
        else {
          const hsize_t event_dim_extended[Event_RANK] = {
            event_offset[0] + event_dim_extend[0], event_dim_extend[1]
          };

          //Extend to new Dimension
          event_data_set.extend(event_dim_extended);
          H5::DataSpace event_file_space = event_data_set.getSpace();

          event_file_space.selectHyperslab(
            H5S_SELECT_SET, event_dim_extend, event_offset);

          H5::DataSpace event_memory_space(Event_RANK, event_dim_extend, NULL);
          event_data_set.write(&event_data[0], H5::PredType::NATIVE_FLOAT,
                               event_memory_space, event_file_space);

          std::fill(event_data.begin(), event_data.end(), NAN);
        }

        //put in block logic
        event_offset[0] += block_size;


        //---------Event Stuff Done. Moving onto track,clustr, and jets------------//

        if (offset[0] == 0) {
          // Writing the first event. The track_data space is already
          // created with space for one event (see when
          // file.createDataSet() was called)
          track_data_set.write(&track_data[0], H5::PredType::NATIVE_FLOAT);
          cluster_data_set.write(&cluster_data[0], H5::PredType::NATIVE_FLOAT);
          jet_ak04tpc_data_set.write(&jet_ak04tpc_data[0], H5::PredType::NATIVE_FLOAT);
          jet_ak02tpc_data_set.write(&jet_ak02tpc_data[0], H5::PredType::NATIVE_FLOAT);

          //Make sure to clear previous arrays
          std::fill(track_data.begin(), track_data.end(), NAN);
          std::fill(cluster_data.begin(), cluster_data.end(), NAN);
          std::fill(jet_ak04tpc_data.begin(), jet_ak04tpc_data.end(), NAN);
          std::fill(jet_ak02tpc_data.begin(), jet_ak02tpc_data.end(), NAN);
        }


        else {
          // The new, extended-by-1 dimension
          const hsize_t track_dim_extended[RANK] = {
            offset[0] + track_dim_extend[0], track_dim_extend[1],
            track_dim_extend[2]
          };

          const hsize_t cluster_dim_extended[RANK] = {
            offset[0] + cluster_dim_extend[0], cluster_dim_extend[1],
            cluster_dim_extend[2]
          };

          const hsize_t jet_ak04tpc_dim_extended[RANK] = {
            offset[0] + jet_ak04tpc_dim_extend[0], jet_ak04tpc_dim_extend[1],
            jet_ak04tpc_dim_extend[2]
          };

          const hsize_t jet_ak02tpc_dim_extended[RANK] = {
            offset[0] + jet_ak02tpc_dim_extend[0], jet_ak02tpc_dim_extend[1],
            jet_ak02tpc_dim_extend[2]
          };

          // Extend to the new dimension
          track_data_set.extend(track_dim_extended);
          cluster_data_set.extend(cluster_dim_extended);
          jet_ak04tpc_data_set.extend(jet_ak04tpc_dim_extended);
          jet_ak02tpc_data_set.extend(jet_ak02tpc_dim_extended);

          // Select the hyperslab that only encompass the
          // difference from extending the data space (i.e. the
          // new event, but offset at the existing event)
          H5::DataSpace track_file_space = track_data_set.getSpace();
          H5::DataSpace cluster_file_space = cluster_data_set.getSpace();
          H5::DataSpace jet_ak04tpc_file_space = jet_ak04tpc_data_set.getSpace();
          H5::DataSpace jet_ak02tpc_file_space = jet_ak02tpc_data_set.getSpace();

          track_file_space.selectHyperslab(
            H5S_SELECT_SET, track_dim_extend, offset);

          cluster_file_space.selectHyperslab(
            H5S_SELECT_SET, cluster_dim_extend, offset);

          jet_ak04tpc_file_space.selectHyperslab(
            H5S_SELECT_SET, jet_ak04tpc_dim_extend, offset);

          jet_ak02tpc_file_space.selectHyperslab(
            H5S_SELECT_SET, jet_ak02tpc_dim_extend, offset);

          // The memory space is the difference only (i.e. also
          // the new event, but at offset 0)
          H5::DataSpace track_memory_space(RANK, track_dim_extend, NULL);
          H5::DataSpace cluster_memory_space(RANK, cluster_dim_extend, NULL);
          H5::DataSpace jet_ak04tpc_memory_space(RANK, jet_ak04tpc_dim_extend, NULL);
          H5::DataSpace jet_ak02tpc_memory_space(RANK, jet_ak02tpc_dim_extend, NULL);

          // Write the data from memory space to file space
          track_data_set.write(&track_data[0], H5::PredType::NATIVE_FLOAT,
                               track_memory_space, track_file_space);

          cluster_data_set.write(&cluster_data[0], H5::PredType::NATIVE_FLOAT,
                                 cluster_memory_space, cluster_file_space);

          jet_ak04tpc_data_set.write(&jet_ak04tpc_data[0], H5::PredType::NATIVE_FLOAT,
                             jet_ak04tpc_memory_space, jet_ak04tpc_file_space);

          jet_ak02tpc_data_set.write(&jet_ak02tpc_data[0], H5::PredType::NATIVE_FLOAT,
                             jet_ak02tpc_memory_space, jet_ak02tpc_file_space);

          std::fill(track_data.begin(), track_data.end(), NAN);
          std::fill(cluster_data.begin(), cluster_data.end(), NAN);
          std::fill(jet_ak04tpc_data.begin(), jet_ak04tpc_data.end(), NAN);
          std::fill(jet_ak02tpc_data.begin(), jet_ak02tpc_data.end(), NAN);
        }
        //put in block log

        offset[0] += block_size;

        fprintf(stderr, "\r%s:%d: %llu / %lld\n", __FILE__, __LINE__, offset[0],
                _tree_event->GetEntries());

      }//BLOCK CHECK

    }//Events

    _tree_event->Delete();
    auxfile->Close();
    delete auxfile;
    file->Close();
    delete file;
  }
}

int main(int argc, char *argv[])
{
  if (argc < 3) {
    fprintf(stderr, "%s", "Syntax is [command] [root_file] [new hdf5 file name]");
    exit(EXIT_FAILURE);
  }

  UInt_t nevent_max = 0;
  /* UInt_t nevent_max = 20000; */
  UInt_t ntrack_max = 0;
  UInt_t ncluster_max = 0;
  UInt_t njet_ak04tpc_max = 0;
  UInt_t njet_ak02tpc_max = 0;
  UInt_t block_size = 2000; //affects chunk size, used from pairing

  find_ntrack_ncluster_max(argv + 1, argv + argc - 1, nevent_max, ntrack_max, ncluster_max, njet_ak04tpc_max, njet_ak02tpc_max);
  /* nevent_max = 529683; ntrack_max = 3786; ncluster_max = 2022; njet_max = 52; //18q_pass3_cluster15 */
  /* nevent_max = 908225; ntrack_max = 3176; ncluster_max = 381; njet_max = 50; //18q_mb */
  /* nevent_max = 908000; ntrack_max = 3176; ncluster_max = 381; njet_max = 50; //18q_mb */
  /* nevent_max = 400000; ntrack_max = 3786; ncluster_max = 2022; njet_max = 52; //18q_pass3_cluster15 */

  fprintf(stderr, "%s:%d: nevents = %u, ntrack_max = %u, ncluster_max = %u, njet_ak04tpc_max = %u, njet_ak02tpc_max = %u\n", __FILE__, __LINE__, nevent_max, ntrack_max, ncluster_max, njet_ak04tpc_max, njet_ak02tpc_max);

  // Access mode H5F_ACC_TRUNC truncates any existing file, while
  // not throwing any exception (unlike H5F_ACC_RDWR)
  std::string file_str = argv[argc - 1];
  H5::H5File file( file_str.c_str(), H5F_ACC_TRUNC );
  //H5::H5File file(argv[argc - 1], H5F_ACC_TRUNC);

  // How many properties per event is written
  //Shoudl be Same as Line 225:
  static const size_t event_row_size = 7;
  static const size_t track_row_size = 10;
  static const size_t cluster_row_size = 32;
  static const size_t jet_row_size = 5;

  // The tensor dimension increment for each new chunk of events
  //The chunking of data can be edited for performance. 2000 is used in mixing.
  hsize_t event_dim_extend[Event_RANK] = {2000, event_row_size};
  hsize_t track_dim_extend[RANK] = { 2000, ntrack_max, track_row_size };
  hsize_t cluster_dim_extend[RANK] = { 2000, ncluster_max, cluster_row_size };
  hsize_t jet_ak04tpc_dim_extend[RANK] = { 2000, njet_ak04tpc_max, jet_row_size };
  hsize_t jet_ak02tpc_dim_extend[RANK] = { 2000, njet_ak02tpc_max, jet_row_size };

  // The maximum tensor dimension, for unlimited number of events
  hsize_t event_dim_max[Event_RANK] = {H5S_UNLIMITED, event_row_size};
  hsize_t track_dim_max[RANK] = { H5S_UNLIMITED, ntrack_max, track_row_size };
  hsize_t cluster_dim_max[RANK] = { H5S_UNLIMITED, ncluster_max, cluster_row_size };
  hsize_t jet_ak04tpc_dim_max[RANK] = { H5S_UNLIMITED, njet_ak04tpc_max, jet_row_size };
  hsize_t jet_ak02tpc_dim_max[RANK] = { H5S_UNLIMITED, njet_ak02tpc_max, jet_row_size };

  // The extensible HDF5 data space
  H5::DataSpace event_data_space(Event_RANK, event_dim_extend, event_dim_max);
  H5::DataSpace track_data_space(RANK, track_dim_extend, track_dim_max);
  H5::DataSpace cluster_data_space(RANK, cluster_dim_extend, cluster_dim_max);
  H5::DataSpace jet_ak04tpc_data_space(RANK, jet_ak04tpc_dim_extend, jet_ak04tpc_dim_max);
  H5::DataSpace jet_ak02tpc_data_space(RANK, jet_ak02tpc_dim_extend, jet_ak02tpc_dim_max);

  // To enable zlib compression (there will be many NANs) and
  // efficient chunking (splitting of the tensor into contingous
  // hyperslabs), a HDF5 property list is needed
  H5::DSetCreatPropList event_property = H5::DSetCreatPropList();
  H5::DSetCreatPropList track_property = H5::DSetCreatPropList();
  H5::DSetCreatPropList cluster_property = H5::DSetCreatPropList();
  H5::DSetCreatPropList jet_ak04tpc_property = H5::DSetCreatPropList();
  H5::DSetCreatPropList jet_ak02tpc_property = H5::DSetCreatPropList();

#ifdef HDF5_USE_DEFLATE
  // Check for zlib (deflate) availability and enable only if
  // present
  if (!H5Zfilter_avail(H5Z_FILTER_DEFLATE)) {
    fprintf(stderr, "%s:%d: warning: deflate filter not "
            "available\n", __FILE__, __LINE__);
  }
  else {
    unsigned int filter_info;

    H5Zget_filter_info(H5Z_FILTER_DEFLATE, &filter_info);
    if (!(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED)) {
      fprintf(stderr, "%s:%d: warning: deflate filter not "
              "available for encoding\n", __FILE__, __LINE__);
    }
    else {
      event_property.setDeflate(1);
      track_property.setDeflate(1);
      cluster_property.setDeflate(1);
      jet_ak04tpc_property.setDeflate(1);
      jet_ak02tpc_property.setDeflate(1);
    }
  }
#endif // HDF5_USE_DEFLATE

  // Activate chunking, while observing the HDF5_DEFAULT_CACHE being
  // the CPU L2 cache size
  hsize_t event_dim_chunk[Event_RANK] = {
    event_dim_extend[0],
    event_dim_extend[1]
  };


  hsize_t track_dim_chunk[RANK] = {
    track_dim_extend[0],
    track_dim_extend[1],
    track_dim_extend[2]
  };

  hsize_t cluster_dim_chunk[RANK] = {
    cluster_dim_extend[0],
    cluster_dim_extend[1],
    cluster_dim_extend[2]
  };

  hsize_t jet_ak04tpc_dim_chunk[RANK] = {
    jet_ak04tpc_dim_extend[0],
    jet_ak04tpc_dim_extend[1],
    jet_ak04tpc_dim_extend[2]
  };

  hsize_t jet_ak02tpc_dim_chunk[RANK] = {
    jet_ak02tpc_dim_extend[0],
    jet_ak02tpc_dim_extend[1],
    jet_ak02tpc_dim_extend[2]
  };

  event_property.setChunk(Event_RANK, event_dim_chunk);
  track_property.setChunk(RANK, track_dim_chunk);
  cluster_property.setChunk(RANK, cluster_dim_chunk);
  jet_ak04tpc_property.setChunk(RANK, jet_ak04tpc_dim_chunk);
  jet_ak02tpc_property.setChunk(RANK, jet_ak02tpc_dim_chunk);

  // Create the data set, which will have space for the first event (chunk)
  H5::DataSet event_data_set =
    file.createDataSet("event", H5::PredType::NATIVE_FLOAT,
                       event_data_space, event_property);
  fprintf(stderr, "%s:%d: CREATED EVENT DATASET", __FILE__, __LINE__);
  H5::DataSet track_data_set =
    file.createDataSet("track", H5::PredType::NATIVE_FLOAT,
                       track_data_space, track_property);

  H5::DataSet cluster_data_set =
    file.createDataSet("cluster", H5::PredType::NATIVE_FLOAT,
                       cluster_data_space, cluster_property);

  H5::DataSet jet_ak04tpc_data_set =
    file.createDataSet("jet_ak04tpc", H5::PredType::NATIVE_FLOAT,
                       jet_ak04tpc_data_space, jet_ak04tpc_property);

  H5::DataSet jet_ak02tpc_data_set =
    file.createDataSet("jet_ak02tpc", H5::PredType::NATIVE_FLOAT,
                       jet_ak02tpc_data_space, jet_ak02tpc_property);

  hsize_t event_offset [Event_RANK] = {0, 0};
  hsize_t offset[RANK] = {0, 0, 0};

  write_track_cluster(event_data_set, track_data_set, cluster_data_set,
                      jet_ak04tpc_data_set, jet_ak02tpc_data_set,
                      event_offset, offset, event_dim_extend, track_dim_extend, cluster_dim_extend,
                      jet_ak04tpc_dim_extend, jet_ak02tpc_dim_extend,
                      nevent_max, ntrack_max, ncluster_max, njet_ak04tpc_max, njet_ak02tpc_max, block_size, argv + 1, argv + argc - 1);
  fprintf(stderr, "%s:%d: After Write Function\n", __func__, __LINE__);

  file.close();

  return EXIT_SUCCESS;
}
