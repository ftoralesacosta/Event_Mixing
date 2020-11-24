#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <H5Cpp.h>

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

void find_ntrack_ncluster_max(char *argv_first[], char *argv_last[], UInt_t &nevent_max, UInt_t &ntrack_max, UInt_t &ncluster_max, UInt_t &njet_max)
{ 
    for (char **p = argv_first; p != argv_last; p++) {
        // Cautious opening of the TTree, capturing all modes of
        // failure, and keep the TDirectoryFile (to be deleted later)
        // to avoid memory leak
        TFile *file = TFile::Open(*p);

        if (file == NULL) {
            continue;
        }

        TDirectoryFile *df = dynamic_cast<TDirectoryFile *>
            (file->Get("AliAnalysisTaskNTGJ"));

        if (df == NULL) {
            fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Cannot open TFile");
            continue;
        }

        TTree *hi_tree = dynamic_cast<TTree *>
            (df->Get("_tree_event"));

        if (hi_tree == NULL) {
            fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Cannot open _tree_event");
            continue;
        }

        // Get the maximum of the "ntrack"

        UInt_t ntrack;
        UInt_t ncluster;
        UInt_t njet_ak04its;
        nevent_max = UInt_t(hi_tree->GetEntries());
	if (nevent_max > 1000000)
	  nevent_max = 1000000;

        hi_tree->SetBranchAddress("ntrack", &ntrack);
        hi_tree->SetBranchAddress("ncluster", &ncluster);
        hi_tree->SetBranchAddress("njet_ak04its", &njet_ak04its);

        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Obtaining ntrack, ncluster, and njet max for hdf5 file");

        //for (Long64_t i = 0; i < hi_tree->GetEntries(); i++) {
	for (Long64_t i = 0; i < nevent_max; i++) {
            hi_tree->GetEntry(i);
            ntrack_max = std::max(ntrack_max, ntrack);
            ncluster_max = std::max(ncluster_max, ncluster);
            njet_max = std::max(njet_max, njet_ak04its);
            fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, i, hi_tree->GetEntries());
        }

        fprintf(stderr, "\n");

        // Fully delete everything
        hi_tree->Delete();
        delete df;
        file->Close();
        delete file;
    }
}

void write_track_cluster(H5::DataSet &event_data_set, H5::DataSet &track_data_set, H5::DataSet &cluster_data_set, H5::DataSet &jet_data_set,
             hsize_t *event_offset, hsize_t *offset,
             const hsize_t *event_dim_extend, const hsize_t *track_dim_extend, const hsize_t *cluster_dim_extend, const hsize_t *jet_dim_extend,
	     const UInt_t nevent_max, const UInt_t ntrack_max, const UInt_t ncluster_max, const UInt_t njet_max, const UInt_t block_size,
             char *argv_first[], char *argv_last[])
{
    for (char **p = argv_first; p != argv_last; p++) {
        TFile *file = TFile::Open(*p);

        if (file == NULL) {
            continue;
        }

        TDirectoryFile *df = dynamic_cast<TDirectoryFile *>
            (file->Get("AliAnalysisTaskNTGJ"));

        if (df == NULL) {
            continue;
        }

        TTree *hi_tree = dynamic_cast<TTree *>
            (df->Get("_tree_event"));

        if (hi_tree == NULL) {
            continue;
        }

        UInt_t nevent;
        std::vector<Double_t> primary_vertex(3, NAN);
        std::vector<Float_t> multiplicity_v0(64, NAN);//64 channels for v0 detector, to be summed

        UInt_t ntrack;
        std::vector<Float_t> track_e(ntrack_max, NAN);
        std::vector<Float_t> track_pt(ntrack_max, NAN);
        std::vector<Float_t> track_eta(ntrack_max, NAN);
        std::vector<Float_t> track_phi(ntrack_max, NAN);
        std::vector<UChar_t> track_quality(ntrack_max, NAN);
        std::vector<Float_t> track_eta_emcal(ntrack_max, NAN);
        std::vector<Float_t> track_phi_emcal(ntrack_max, NAN);
        std::vector<UChar_t> track_its_ncluster(ntrack_max, NAN);
        std::vector<Float_t> track_its_chi_square(ntrack_max, NAN);
        std::vector<Float_t> track_dca_xy(ntrack_max, NAN);
	std::vector<Float_t> track_dca_z(ntrack_max, NAN);

        UInt_t ncluster;
        std::vector<Float_t> cluster_e(ncluster_max, NAN);
        std::vector<Float_t> cluster_pt(ncluster_max, NAN);
        std::vector<Float_t> cluster_eta(ncluster_max, NAN);
        std::vector<Float_t> cluster_phi(ncluster_max, NAN);
        std::vector<Float_t> cluster_e_cross(ncluster_max, NAN);
        //std::vector<std::vector<Float_t> > cluster_s_nphoton(ncluster_max,std::vector <Float_t> (4, NAN) );

        UInt_t njet_ak04its;
        std::vector<Float_t> jet_ak04its_pt_raw(njet_max, NAN);
        std::vector<Float_t> jet_ak04its_eta_raw(njet_max, NAN);
        std::vector<Float_t> jet_ak04its_phi(njet_max, NAN);
        std::vector<Float_t> jet_ak04its_ptd_raw(njet_max, NAN);
        std::vector<UShort_t> jet_ak04its_multiplicity(njet_max, NAN);
        // std::vector<Float_t> jet_ak04its_width_sigma(njet_max, NAN);

        hi_tree->SetBranchAddress("primary_vertex", &primary_vertex[0]);
        hi_tree->SetBranchAddress("multiplicity_v0", &multiplicity_v0[0]);

        hi_tree->SetBranchAddress("ntrack", &ntrack);
        hi_tree->SetBranchAddress("track_e", &track_e[0]);
        hi_tree->SetBranchAddress("track_pt", &track_pt[0]);
        hi_tree->SetBranchAddress("track_eta", &track_eta[0]);
        hi_tree->SetBranchAddress("track_phi", &track_phi[0]);
        hi_tree->SetBranchAddress("track_quality", &track_quality[0]);
        hi_tree->SetBranchAddress("track_eta_emcal", &track_eta_emcal[0]);
        hi_tree->SetBranchAddress("track_phi_emcal", &track_phi_emcal[0]);
        hi_tree->SetBranchAddress("track_its_ncluster", &track_its_ncluster[0]);
        hi_tree->SetBranchAddress("track_its_chi_square", &track_its_chi_square[0]);
        hi_tree->SetBranchAddress("track_dca_xy", &track_dca_xy[0]);
	hi_tree->SetBranchAddress("track_dca_z", &track_dca_z[0]);

        hi_tree->SetBranchAddress("ncluster", &ncluster);
        hi_tree->SetBranchAddress("cluster_e", &cluster_e[0]);
        hi_tree->SetBranchAddress("cluster_pt", &cluster_pt[0]);
        hi_tree->SetBranchAddress("cluster_eta", &cluster_eta[0]);
        hi_tree->SetBranchAddress("cluster_phi", &cluster_phi[0]);
        hi_tree->SetBranchAddress("cluster_e_cross", &cluster_e_cross[0]);
        //hi_tree->SetBranchAddress("cluster_s_nphoton", &cluster_s_nphoton[0]);

        hi_tree->SetBranchAddress("njet_ak04its", &njet_ak04its);
        hi_tree->SetBranchAddress("jet_ak04its_pt_raw", &jet_ak04its_pt_raw[0]);
        hi_tree->SetBranchAddress("jet_ak04its_eta_raw", &jet_ak04its_eta_raw[0]);
        hi_tree->SetBranchAddress("jet_ak04its_phi", &jet_ak04its_phi[0]);
        hi_tree->SetBranchAddress("jet_ak04its_ptd_raw", &jet_ak04its_ptd_raw[0]);
        hi_tree->SetBranchAddress("jet_ak04its_multiplicity_raw", &jet_ak04its_multiplicity[0]);
        // hi_tree->SetBranchAddress("jet_ak04its_width_sigma", &jet_ak04its_width_sigma[0]);

	int nEventVariables = 2;
	int nClusterVariables = 5;
	int nTrackVariables = 10;
        int nJetVariables = 5;

	fprintf(stderr,"\n %d: BLOCK SIZE = %u\n",__LINE__,block_size);

	std::vector<float> event_data (block_size * nEventVariables, NAN); //2 variables, multp, vertx
	std::vector<float> track_data(block_size * ntrack_max * 11, NAN);
	std::vector<float> cluster_data(block_size * ncluster_max * 5, NAN);
	std::vector<float> jet_data(block_size * njet_max * nJetVariables, NAN);

        //for (Long64_t i = 0; i < hi_tree->GetEntries(); i++) {
	for (Long64_t i = 0; i < 1000000; i++) {
            hi_tree->GetEntry(i);
	    
	    int iblock = i % block_size;
	    //fprintf(stderr,"\n %d: iblock = %i \n",__LINE__,iblock);
	   
            float multiplicity_sum = 0;
            for (int k = 0; k < 64; k++) multiplicity_sum += multiplicity_v0[k];        
            event_data[iblock*nEventVariables + 0] = primary_vertex[2]; //xyz, choose 3rd element, z
            event_data[iblock*nEventVariables + 1] = multiplicity_sum;
        
            for (Long64_t j = 0; j < ntrack; j++) {
                // Note HDF5 is always row-major (C-like)
                track_data[iblock*ntrack_max*11 + j*11 + 0] = track_e[j];
                track_data[iblock*ntrack_max*11 + j*11 + 1] = track_pt[j];
                track_data[iblock*ntrack_max*11 + j*11 + 2] = track_eta[j];
                track_data[iblock*ntrack_max*11 + j*11 + 3] = track_phi[j];
                track_data[iblock*ntrack_max*11 + j*11 + 4] = track_quality[j];
                track_data[iblock*ntrack_max*11 + j*11 + 5] = track_eta_emcal[j];
                track_data[iblock*ntrack_max*11 + j*11 + 6] = track_phi_emcal[j];
                track_data[iblock*ntrack_max*11 + j*11 + 7] = track_its_ncluster[j];
                track_data[iblock*ntrack_max*11 + j*11 + 8] = track_its_chi_square[j];
                track_data[iblock*ntrack_max*11 + j*11 + 9] = track_dca_xy[j];
		track_data[iblock*ntrack_max*11 + j*11 + 10] = track_dca_z[j];
            }

            for(Long64_t n = 0; n < ncluster; n++){
                cluster_data[iblock*ncluster_max*5 + n*5 + 0] = cluster_e[n];
                cluster_data[iblock*ncluster_max*5 + n*5 + 1] = cluster_pt[n];
                cluster_data[iblock*ncluster_max*5 + n*5 + 2] = cluster_eta[n];
                cluster_data[iblock*ncluster_max*5 + n*5 + 3] = cluster_phi[n];
                cluster_data[iblock*ncluster_max*5 + n*5 + 4] = cluster_e_cross[n];
                //cluster_data[j * 7 + 5] = cluster_s_nphoton[j][0];
                //cluster_data[j * 7 + 6] = cluster_s_nphoton[j][1];
            }

            for (Long64_t j = 0; j < njet_ak04its; j++) {
                jet_data[iblock*njet_max*nJetVariables + j*nJetVariables + 0] = jet_ak04its_pt_raw[j];
                jet_data[iblock*njet_max*nJetVariables + j*nJetVariables + 1] = jet_ak04its_eta_raw[j];
                jet_data[iblock*njet_max*nJetVariables + j*nJetVariables + 2] = jet_ak04its_phi[j];
                jet_data[iblock*njet_max*nJetVariables + j*nJetVariables + 3] = jet_ak04its_ptd_raw[j];
                jet_data[iblock*njet_max*nJetVariables + j*nJetVariables + 4] = jet_ak04its_multiplicity[j];
            }
	
	    if (iblock == (block_size-1)) {

	      // for (int i = 0; i < block_size; i++){
	      // 	fprintf(stderr,"\n %d: event %i: track pT = %f\n",__LINE__,i,track_data[i*ntrack_max*10 + 0*ntrack_max + 1]); //look at track pT of first track per event
	      // } //check array

	      if (event_offset[0] == 0) {
		//write first event
		event_data_set.write(&event_data[0], H5::PredType::NATIVE_FLOAT);
	      }
	      
	      else{ 
		//set up extension
		const hsize_t event_dim_extended[2] = {
		  event_offset[0] + event_dim_extend[0], event_dim_extend[1] };
        
		//Extend to new Dimension
		event_data_set.extend(event_dim_extended);
		H5::DataSpace event_file_space = event_data_set.getSpace();

		event_file_space.selectHyperslab(
		   H5S_SELECT_SET, event_dim_extend, event_offset);

		H5::DataSpace event_memory_space(2, event_dim_extend, NULL);
		event_data_set.write(&event_data[0], H5::PredType::NATIVE_FLOAT,
				     event_memory_space, event_file_space);
	      }

	      //put in block logic
	      event_offset[0] += block_size;

	      if (offset[0] == 0) {
		// Writing the first event. The track_data space is already
		// created with space for one event (see when
		// file.createDataSet() was called)
		track_data_set.write(&track_data[0], H5::PredType::NATIVE_FLOAT);
		cluster_data_set.write(&cluster_data[0], H5::PredType::NATIVE_FLOAT);
		jet_data_set.write(&jet_data[0], H5::PredType::NATIVE_FLOAT);
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
		
                const hsize_t jet_dim_extended[RANK] = {
		  offset[0] + jet_dim_extend[0], jet_dim_extend[1],
		  jet_dim_extend[2]
                };

                // Extend to the new dimension
                track_data_set.extend(track_dim_extended);
                cluster_data_set.extend(cluster_dim_extended);
                jet_data_set.extend(jet_dim_extended);

                // Select the hyperslab that only encompass the
                // difference from extending the data space (i.e. the
                // new event, but offsetted at the existing event)
                H5::DataSpace track_file_space = track_data_set.getSpace();
                H5::DataSpace cluster_file_space = cluster_data_set.getSpace();
                H5::DataSpace jet_file_space = jet_data_set.getSpace();

                track_file_space.selectHyperslab(
		  H5S_SELECT_SET, track_dim_extend, offset);

                cluster_file_space.selectHyperslab(
		  H5S_SELECT_SET, cluster_dim_extend, offset);

                jet_file_space.selectHyperslab(
		  H5S_SELECT_SET, jet_dim_extend, offset);

                // The memory space is the difference only (i.e. also
                // the new event, but at offset 0)
                H5::DataSpace track_memory_space(RANK, track_dim_extend, NULL);
                H5::DataSpace cluster_memory_space(RANK, cluster_dim_extend, NULL);
                H5::DataSpace jet_memory_space(RANK, jet_dim_extend, NULL);

                // Write the data from memory space to file space
                track_data_set.write(&track_data[0], H5::PredType::NATIVE_FLOAT,
				     track_memory_space, track_file_space);

                cluster_data_set.write(&cluster_data[0], H5::PredType::NATIVE_FLOAT,
				       cluster_memory_space, cluster_file_space);

                jet_data_set.write(&jet_data[0], H5::PredType::NATIVE_FLOAT,
				   jet_memory_space, jet_file_space);

	      }
	      //put in block log

	      offset[0]+=block_size;	 
	      
	      fprintf(stderr, "%s:%d: %llu / %lld\n", __FILE__,__LINE__, offset[0],
		      hi_tree->GetEntries());

	    }//BLOCK CHECK
	    
        }//Events
	
        hi_tree->Delete();
        delete df;
        file->Close();
        delete file;
    }
}

int main(int argc, char *argv[])
{
    if (argc < 3) {
      fprintf(stderr, "%s", "Syntax is [root_file] [new hdf5 file name]");
        exit(EXIT_FAILURE);
    }

    //UInt_t ntrack_max = 1927 ;
    
    // hard-coded values for 13d_minbias_0GeVTracks.root
    // UInt_t nevent_max = 1946854;
    // UInt_t ntrack_max = 432;
    // UInt_t ncluster_max = 29;
    // UInt_t njet_max = 42;

    UInt_t nevent_max = 0;
    UInt_t ntrack_max = 0;
    UInt_t ncluster_max = 0;
    UInt_t njet_max = 0;
    UInt_t block_size = 2000; //affects chunk size, used from pairing

    find_ntrack_ncluster_max(argv + 1, argv + argc - 1, nevent_max, ntrack_max, ncluster_max, njet_max);
    fprintf(stderr, "%sf:%d: nevents = %u, ntrack_max = %u, ncluster_max = %u, njet_max = %u\n", __FILE__, __LINE__, nevent_max, ntrack_max, ncluster_max, njet_max);

    // Access mode H5F_ACC_TRUNC truncates any existing file, while
    // not throwing any exception (unlike H5F_ACC_RDWR)
    std::string file_str = argv[argc-1];
    H5::H5File file( file_str.c_str(), H5F_ACC_TRUNC );
    //H5::H5File file(argv[argc - 1], H5F_ACC_TRUNC);

    // How many properties per track is written
    static const size_t event_row_size = 2;
    static const size_t track_row_size = 11;
    static const size_t cluster_row_size = 5;
    static const size_t jet_row_size = 5;
    //easier to just make cluster use same # properties, reuse dim_extend


    // The tensor dimension increment for each new chunk of events
    hsize_t event_dim_extend[Event_RANK] = {2000, event_row_size};
    hsize_t track_dim_extend[RANK] = { 2000, ntrack_max, track_row_size };
    hsize_t cluster_dim_extend[RANK] = { 2000, ncluster_max, cluster_row_size };
    hsize_t jet_dim_extend[RANK] = { 2000, njet_max, jet_row_size };

    // The maximum tensor dimension, for unlimited number of events
    hsize_t event_dim_max[Event_RANK] = {H5S_UNLIMITED, event_row_size};
    hsize_t track_dim_max[RANK] = { H5S_UNLIMITED, ntrack_max, track_row_size };
    hsize_t cluster_dim_max[RANK] = { H5S_UNLIMITED, ncluster_max, cluster_row_size };
    hsize_t jet_dim_max[RANK] = { H5S_UNLIMITED, njet_max, jet_row_size };

    // The extensible HDF5 data space
    H5::DataSpace event_data_space(Event_RANK, event_dim_extend, event_dim_max);
    H5::DataSpace track_data_space(RANK, track_dim_extend, track_dim_max);
    H5::DataSpace cluster_data_space(RANK, cluster_dim_extend, cluster_dim_max);
    H5::DataSpace jet_data_space(RANK, jet_dim_extend, jet_dim_max);
    //might need two data_spaces: track_data_space & cluster_data_space


    // To enable zlib compression (there will be many NANs) and
    // efficient chunking (splitting of the tensor into contingous
    // hyperslabs), a HDF5 property list is needed
    H5::DSetCreatPropList event_property = H5::DSetCreatPropList();
    H5::DSetCreatPropList track_property = H5::DSetCreatPropList();
    H5::DSetCreatPropList cluster_property = H5::DSetCreatPropList();
    H5::DSetCreatPropList jet_property = H5::DSetCreatPropList();

// #ifdef HDF5_USE_DEFLATE
//     // Check for zlib (deflate) availability and enable only if
//     // present
//     if (!H5Zfilter_avail(H5Z_FILTER_DEFLATE)) {
//         fprintf(stderr, "%s:%d: warning: deflate filter not "
//                 "available\n", __FILE__, __LINE__);
//     }
//     else {
//         unsigned int filter_info;

//         H5Zget_filter_info(H5Z_FILTER_DEFLATE, &filter_info);
//         if (!(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED)) {
//             fprintf(stderr, "%s:%d: warning: deflate filter not "
//                     "available for encoding\n", __FILE__, __LINE__);
//         }
//         else {
//             event_property.setDeflate(1);
//             track_property.setDeflate(1);
//             cluster_property.setDeflate(1);
//             jet_property.setDeflate(1);
//         }
//     }
// #endif // HDF5_USE_DEFLATE

    // Activate chunking, while observing the HDF5_DEFAULT_CACHE being
    // the CPU L2 cache size
    hsize_t event_dim_chunk[Event_RANK] = {
        // std::max(static_cast<unsigned long long>(1),
        //          HDF5_DEFAULT_CACHE /
        //          std::max(static_cast<unsigned long long>(1),
        //                   event_dim_extend[0]*
        //                   sizeof(float))),
	event_dim_extend[0],
        event_dim_extend[1]
    };
    

    hsize_t track_dim_chunk[RANK] = {
        // std::max(static_cast<unsigned long long>(1),
        //          HDF5_DEFAULT_CACHE /
        //          std::max(static_cast<unsigned long long>(1),
        //                   track_dim_extend[1] * track_dim_extend[2] *
        //                   sizeof(float))),
	track_dim_extend[0],
	track_dim_extend[1],
        track_dim_extend[2]
    };

    hsize_t cluster_dim_chunk[RANK] = {
        // std::max(static_cast<unsigned long long>(1),
        //          HDF5_DEFAULT_CACHE /
        //          std::max(static_cast<unsigned long long>(1),
        //                   cluster_dim_extend[1] * cluster_dim_extend[2] *
        //                   sizeof(float))),
	cluster_dim_extend[0],
        cluster_dim_extend[1],
        cluster_dim_extend[2]
    };

    hsize_t jet_dim_chunk[RANK] = {
        // std::max(static_cast<unsigned long long>(1),
        //          HDF5_DEFAULT_CACHE /
        //          std::max(static_cast<unsigned long long>(1),
        //                   jet_dim_extend[1] * jet_dim_extend[2] *
        //                   sizeof(float))),
        jet_dim_extend[0],
        jet_dim_extend[1],
        jet_dim_extend[2]
    };

    event_property.setChunk(Event_RANK, event_dim_chunk);
    track_property.setChunk(RANK, track_dim_chunk);
    cluster_property.setChunk(RANK, cluster_dim_chunk);
    jet_property.setChunk(RANK, jet_dim_chunk);

    // Create the data set, which will have space for the first event (chunk)
    H5::DataSet event_data_set =
      file.createDataSet("event", H5::PredType::NATIVE_FLOAT,
             event_data_space, event_property);
    fprintf(stderr,"%s:%d: CREATED EVENT DATASET",__FILE__,__LINE__);
    H5::DataSet track_data_set =
      file.createDataSet("track", H5::PredType::NATIVE_FLOAT,
             track_data_space, track_property);

    H5::DataSet cluster_data_set =
      file.createDataSet("cluster", H5::PredType::NATIVE_FLOAT,
             cluster_data_space, cluster_property);

    H5::DataSet jet_data_set =
      file.createDataSet("jet", H5::PredType::NATIVE_FLOAT,
             jet_data_space, jet_property);

    hsize_t event_offset [2] = {0,0};
    hsize_t offset[RANK] = {0, 0, 0};

    write_track_cluster(event_data_set, track_data_set, cluster_data_set, jet_data_set,
            event_offset, offset, event_dim_extend, track_dim_extend, cluster_dim_extend, jet_dim_extend,
	    nevent_max, ntrack_max, ncluster_max, njet_max, block_size,argv + 1, argv + argc - 1);

    file.close();

    return EXIT_SUCCESS;
}
