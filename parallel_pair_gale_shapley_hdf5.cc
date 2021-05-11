#include <cstdio>
#include <cmath>
#include <sys/time.h>

#include <vector>
#include <list>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>

#include <omp.h>
#include <H5Cpp.h>

//#define HI_TREE "hiEvtAnalyzer/HiTree"
#define HI_TREE "_tree_event"
#define HI_TREE_2 "AliAnalysisTaskNTGJ/_tree_event"
//using namespace std;
namespace {

	typedef unsigned short index_t;

	size_t nevent(const char *filename)
	{

		TFile *root_file = TFile::Open(filename);

		if (root_file == NULL) {
		  fprintf(stderr, "%s:%d: ROOT FILE FAIL\n",__FILE__, __LINE__);
			return 0;
		}

		TTree *hi_tree = dynamic_cast<TTree *>(root_file->Get(HI_TREE));

		if (hi_tree == NULL) {
		  hi_tree = dynamic_cast<TTree *>(root_file->Get(HI_TREE_2));		
		  if(hi_tree == NULL){
		    fprintf(stderr, "%s:%d: TREE FAIL\n",__FILE__, __LINE__);
		    return 0;
		  }
		}

		const size_t ret = hi_tree->GetEntries();

		root_file->Close();

		return ret;
	}

	std::vector<float> feature_extract(const char *filename,
					   const size_t event_start,
					   const size_t event_end,
					   const size_t nfeature)
	{
		TFile *root_file = TFile::Open(filename);

		if (root_file == NULL) {
		  fprintf(stderr, "%s:%d: ROOT FILE FAIL\n",__FILE__, __LINE__);
			return std::vector<float>();
		}

		TTree *hi_tree = dynamic_cast<TTree *>(root_file->Get(HI_TREE));

		if (hi_tree == NULL) {
		  hi_tree = dynamic_cast<TTree *>(root_file->Get(HI_TREE_2));
		  if(hi_tree == NULL){
		    fprintf(stderr, "%s:%d: TREE FAIL\n",__FILE__, __LINE__);
		    return std::vector<float>();
		  }
		}

		double vertex[3];

		hi_tree->SetBranchAddress("primary_vertex", vertex);

		float multiplicity_v0[64];
		float event_plane_angle[3];//v1,v2,v3 flow coefficients

		switch (nfeature) {
		case 2:
			hi_tree->SetBranchAddress("multiplicity_v0", &multiplicity_v0);
			break;
		case 3:
		        hi_tree->SetBranchAddress("multiplicity_v0", &multiplicity_v0);
		        hi_tree->SetBranchAddress("event_plane_psi_v0", &event_plane_angle[0]);
			break;

		default:
			fprintf(stderr, "%s:%d: illegal nfeature = %lu\n",
					__FILE__, __LINE__, nfeature);
			return std::vector<float>();
		}

		std::vector<float> ret;

		for (size_t i = event_start; i < event_end; i++) {
			hi_tree->GetEntry(i);

			ret.push_back(vertex[2]);

			if (nfeature >= 2) {
 			        float multp_sum = 0;
			        for (int k = 0; k < 64; k++) {
				  multp_sum += multiplicity_v0[k];
				}
				ret.push_back(multp_sum);

				if (nfeature >= 3)
				  ret.push_back(event_plane_angle[1]);
        if (i > event_start - 3)
        fprintf(stderr,"%s: %d: z-vertex = %1.2f mp = %f, v2 = %1.2f\n",__func__,__LINE__,vertex[2],multp_sum,event_plane_angle[1]);
			}

		}
//    for (int i = 0; i < ret.size();i+=3){
 //       fprintf(stderr,"%d: size =  %i, z-vertex = %1.2f mp = %f, v2 = %1.2f\n",
  //          __LINE__,ret.size(),ret[i+0],ret[i+1],ret[i+2]);
    //}
		root_file->Close();
		return ret;
	}


  size_t nevent_hdf5(const char *filename)
  {

    std::string file_str = filename;

    const H5std_string hdf5_file_name(file_str.c_str());
    H5::H5File h5_file( file_str.c_str(), H5F_ACC_RDONLY );

    const std::string event_ds_name( "event" );
    H5::DataSet event_dataset = h5_file.openDataSet( event_ds_name.c_str() );
    H5::DataSpace event_dataspace = event_dataset.getSpace();

    //Initialize Event Dimensions
    const int event_ndims = event_dataspace.getSimpleExtentNdims();
    hsize_t event_maxdims[event_ndims];
    hsize_t eventdims[event_ndims];
    event_dataspace.getSimpleExtentDims(eventdims, event_maxdims);
    fprintf(stderr,"Number of events in hdf5 = %i\n",eventdims[0]);

    size_t ret = eventdims[0];
    return ret;
  }

  std::vector<float> feature_extract_hdf5(const char *filename,
					   const size_t event_start,
					   const size_t event_end,
					   const size_t nfeature)
	{

	  //HAVE THE BLOCK SIZE BE EVENT_END-EVENT_START, AND PUT OFFSET AS EVEST_START
	  /* fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Befor HDF5 READ"); */
	  std::string file_str = filename;

    const int EVENT_RANK = 2; // Rank 2 because [No. Events, Event Variables]
	  
		const H5std_string hdf5_file_name(file_str.c_str());
		H5::H5File h5_file( file_str.c_str(), H5F_ACC_RDONLY );

		const std::string event_ds_name( "event" );
		H5::DataSet event_dataset = h5_file.openDataSet( event_ds_name.c_str() );
		H5::DataSpace event_dataspace = event_dataset.getSpace();

		//Initialize Event Dimensions
		const int event_ndims = event_dataspace.getSimpleExtentNdims();
		hsize_t event_maxdims[event_ndims];
		hsize_t eventdims[event_ndims];
		event_dataspace.getSimpleExtentDims(eventdims, event_maxdims);

		UInt_t NEvent_Vars = eventdims[1]; //# event properties
		//Define array hyperslab will be read into 

		const hsize_t block_size = event_end-event_start;
		float event_data_out[block_size][NEvent_Vars] = {0};
		/* fprintf(stderr,"%s:%d: HDF5 Block Size = %i, Event Start = %i, NVars = %i\n", */
		/* 	__FILE__,__LINE__,block_size,event_start,NEvent_Vars); */

		hsize_t event_offset[EVENT_RANK] = {event_start,0};
		hsize_t event_count[EVENT_RANK] = {block_size, NEvent_Vars};
		//const int Event_RANK_OUT = 2; //Event #, Event properties
		H5::DataSpace event_memspace( EVENT_RANK, eventdims );
		hsize_t event_offset_out[EVENT_RANK] = {0,0};
		hsize_t event_count_out[EVENT_RANK] = {block_size, NEvent_Vars};

		event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
		event_memspace.selectHyperslab( H5S_SELECT_SET, event_count_out, event_offset_out );
		event_dataset.read( event_data_out, H5::PredType::NATIVE_FLOAT, event_memspace, event_dataspace);
		/* fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "event dataset read into array: OK"); */

		std::vector<float> ret;

		for (size_t i = 0; i < block_size; i++) {
		  float z_vtx = event_data_out[i][0];
		  ret.push_back(z_vtx);
		  if (nfeature >= 2) {
		    float multp = event_data_out[i][1];
		    ret.push_back(multp);
		    if (nfeature >=3){
          float v2 = event_data_out[i][2];
          ret.push_back(v2);//index 2 is v1, 3 is v2, 4 is v3
          if (std::isnan(v2)||v2==0)
            fprintf(stderr,"%s: %d: NAN or 0 v2; i = %zu, z-vertex = %1.2f mp = %f, v2 = %1.2f\n",__func__,__LINE__,i,z_vtx,multp,v2);
          if (i==0)
            fprintf(stderr,"%s: %d: z-vertex = %1.2f mp = %f, v2 = %1.2f\n",__func__,__LINE__,z_vtx,multp,v2);
          }
		  }
		}

    for (int i = 0; i < ret.size();i+=3){
      if (i < ret.size() - 3) continue;
      fprintf(stderr,"%d: size =  %i, z-vertex = %1.2f mp = %f, v2 = %1.2f\n",
          __LINE__,ret.size(),ret[i+0],ret[i+1],ret[i+2]);
    }
		return ret;
	}
  
	void feature_normalize(std::vector<float> &u,
			       std::vector<float> &v, const size_t n)
  {
    //normalize vectors such that event properties with very different magnitudes can be equally weighted in pairing
    //ex: multiplicity can range from 0-1,000. But z-vertex only ranges from 0-10. Normalization is needed to pair
    //events with both criteria simoultaneously.
    /* "n" is the number of features. Likely to be 3: z-vert, multp, and v2 */

    //sum each event property for each vector (datasets)
    std::vector<float> s(n, 0);
    for (size_t j = 0; j < n; j++) {
      float s_j = 0;

      for (size_t i = 0; i < u.size(); i += n) {
        s_j += fabsf(u[i + j]);
        /* if (j==1 && i%500 == 0) */
        /*   fprintf(stderr,"%s: %d: i = %zu, mp_sum = %f\n",__func__,__LINE__,i/3,s_j); */
        /* if (i/3 > 1997) */
        /*   fprintf(stderr,"%s: %d: i = %zu z = %f, mp = %f, v2 = %f \n",__func__, __LINE__,i/3,u[i+0],u[i+1],u[i+2]); */
      }
      s[j] = s_j;
    }

    for (size_t j = 0; j < n; j++) {
      float s_j = 0;

      for (size_t i = 0; i < v.size(); i += n) {
        s_j += fabsf(v[i + j]);
        /* if (j==0 && i/3>1997){ */
        /*   fprintf(stderr,"%s: %d: i = %i, z = %f, mp = %f, v2 = %f \n",__func__,__LINE__,i/3,v[i+0],v[i+1],v[i+2]); */
          /* fprintf(stderr,"%s: %d: v[i + j] = %f\n",__func__,__LINE__,v[i + j]); */
        /* } */
      }
      s[j] += s_j;
    }

    //get normalization factor: N_Events/property_sum
    for (size_t j = 0; j < n; j++) {
      s[j] = (u.size() + v.size()) / s[j];
    }

    //apply normalization to each vector. you want the same normalization applied to each dataset.
    for (size_t i = 0; i < u.size(); i += n) {
      for (size_t j = 0; j < n; j++) {
        u[i + j] *= s[j];
      }
    }

    for (size_t i = 0; i < v.size(); i += n) {
      for (size_t j = 0; j < n; j++) {
        v[i + j] *= s[j];
      }
    }

  }
}

bool preference_compare(const std::pair<float, index_t> u,
    const std::pair<float, index_t> v)
{
  return u.first > v.first;
}

void order_preference(std::vector<std::list<index_t> > &up,
    std::vector<std::list<index_t> > &vp,
    std::vector<float> u, std::vector<float> v,
    const size_t n, const size_t nduplicate_v)
{
  const size_t u_size_n = u.size() / n; //feature vec size / n_features = num events
  const size_t v_size_n = v.size() / n;

  //fprintf(stderr, "%s:%d: %lux%lu Thread %i\n", __FILE__, __LINE__,
  //u_size_n, v_size_n, omp_get_thread_num());

  const size_t size_max = std::max(u_size_n, v_size_n);

  up.resize(u_size_n, std::list<index_t>());

  //Loopt through each event in dataset_0 block
  for (size_t i = 0; i < u_size_n; i++) {
    std::vector<std::pair<float, index_t> > l;
    /* if (i==0) */
      /* fprintf(stderr,"%s: %d:[i] z = %f, mp = %f, v2 = %f\n",__func__,__LINE__,u[i*n+0],u[i*n+1],u[i*n+2]); */
    for (size_t j = 0; j < v_size_n; j++) {
      float d = 0;

      /* if (j==0 && i==0) */
        /* fprintf(stderr,"%s: %d:[j] z = %f, mp = %f, v2 = %f\n",__func__,__LINE__,v[i*n+0],v[i*n+1],v[i*n+2]); */

      for (size_t k = 0; k < n; k++)//loop through featurs for distance metric d
        d += std::pow(u[i * n + k] - v[j * n + k], 2);
      if (d==0) d = 999999; //Avoid pairing identical events
      l.push_back(std::pair<float, size_t>(d, j));
      /* if (i%1000 == 0 || j%1000 == 0); */
      /* fprintf(stderr,"%d: Distance = %f\n",__LINE__,d); */
    }
    std::sort(l.begin(), l.end(), preference_compare);
    // up.push_back(std::list<index_t>());
    for (size_t j = 0; j < l.size(); j++) {
      for (size_t k = 0; k < nduplicate_v; k++) {
        up[i].push_front(l[j].second + k * v_size_n);
        //fprintf(stderr,"%d: v-event = %i,  u-event = %i, d = %f \n",__LINE__,i,l[j].second,l[j].first);
      }
    }
    up[i].resize(size_max, v_size_n);
    //if (i % 100 == 0) {
    //  fprintf(stderr, "%s:%d: %lu/%lu\n", __FILE__, __LINE__,i, u_size_n);
    //}
  }

  vp.resize(v_size_n * nduplicate_v, std::list<index_t>());

  //Loop through each event in dataset_1 block
  for (size_t j = 0; j < v_size_n; j++) {
    std::vector<std::pair<float, index_t> > l;

    for (size_t i = 0; i < u_size_n; i++) {
      float d = 0;

      for (size_t k = 0; k < n; k++) {
        d += std::pow(u[i * n + k] - v[j * n + k], 2);
      }
      if (d==0) d = 999999;//avoid event self-pairing
      l.push_back(std::pair<float, index_t>(d, i));
      //fprintf(stderr,"%d: Distance = %f\n",__LINE__,d);
    }
    std::sort(l.begin(), l.end(), preference_compare);

    std::list<index_t> b;

    for (size_t i = 0; i < l.size(); i++) {
      b.push_front(l[i].second);
    }
    b.resize(size_max, u_size_n);
    for (size_t k = 0; k < nduplicate_v; k++) {
      vp[j * nduplicate_v + k] = b;
    }

    //if (j % 100 == 0) {
    //fprintf(stderr, "%s:%d: %lu/%lu\n", __FILE__, __LINE__,j, v_size_n);}
  }

  /* fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Order Preference done for this block"); */
}

std::vector<index_t> gale_shapley(std::vector<std::list<index_t> > &mp,
    std::vector<std::list<index_t> > &fp)
{
  /* pass in male and female preference lists. Create vectors to track */ 
  /* which males are engaged to which females and vice-versa. */
  std::vector<index_t> m_to_f_engaged(mp.size(), fp.size());
  std::vector<index_t> f_to_m_engaged(fp.size(), mp.size());

  std::vector<std::vector<std::pair<
    std::vector<std::list<index_t> >::iterator,
    std::list<index_t>::iterator> > > mp_index;

  mp_index.resize(fp.size());

  /* iterate through all males (outer), and then that male's preference list (inner) */
  for (std::vector<std::list<index_t> >::iterator
      iterator_outer = mp.begin();
      iterator_outer != mp.end(); iterator_outer++) {

    for (std::list<index_t>::iterator
        iterator_inner = iterator_outer->begin();
        iterator_inner != iterator_outer->end();
        iterator_inner++) {

      mp_index[*iterator_inner].push_back(
          std::pair<std::vector<std::list<index_t> >::iterator,
          std::list<index_t>::iterator>(
            iterator_outer, iterator_inner));

    }

    //if ((iterator_outer - mp.begin()) % 100 == 0) {
    //  fprintf(stderr, "%s:%d: %lu/%lu\n", __FILE__, __LINE__, 
    //      (iterator_outer - mp.begin()), mp.size()); }


  }

  for (;;) {
    std::vector<index_t>::const_iterator m_iterator =
      std::find(m_to_f_engaged.begin(),
          m_to_f_engaged.end(), fp.size());

    if (m_iterator == m_to_f_engaged.end()) {
      break;
    }

    const index_t m = m_iterator - m_to_f_engaged.begin();
    const index_t w = mp[m].front();

    /* if (m % 500 == 0 || w % 500 == 0) */
    /* fprintf(stderr, "%s:%d: %hu<>%hu\n", __FILE__, __LINE__,m, w); */


    // Some man p is engaged to w
    index_t p = f_to_m_engaged[w];

    if (p != mp.size()) {
      // Assign p to be free
      m_to_f_engaged[p] = fp.size(); //clever way of assigning non-existent fp element index (fp index are 0 through size-1)
    }
    f_to_m_engaged[w] = m;
    m_to_f_engaged[m] = w;

    std::list<index_t>::iterator s =
      std::find(fp[w].begin(), fp[w].end(), m);

    s++;

    //once paired, delete partners lower on preference list -> faster subsequent pairings

    /* for (std::vector<std::pair<std::vector<std::list<index_t> >:: */
    /*     iterator, std::list<index_t>::iterator> >::iterator */
    /*     iterator = mp_index[w].begin(); */
    for (auto iterator = mp_index[w].begin();
        iterator != mp_index[w].end(); iterator++) {
      iterator->first->erase(iterator->second);
    }
    fp[w].erase(s, fp[w].end());
  }
  return m_to_f_engaged;
}

std::map<size_t,std::vector<Long64_t> > mix_gale_shapley(const char *filename_0, const char *filename_1, 
    const char *mixing_start, const char *mixing_end,
    const char *GeV_Track_Skim, const int nfeature, 
    const int nduplicate)
{

  std::map<size_t,std::vector<Long64_t> > Matches;

  double time_spent = 0.0;
  clock_t begin = clock();

  size_t mix_start = atoi(mixing_start);
  size_t mix_end = atoi(mixing_end);
  int Track_Skim = atoi(GeV_Track_Skim);

  /* const size_t nevent_0 = nevent(filename_0); */
  size_t nevent_0 = nevent_hdf5(filename_0);
  size_t nevent_1 = nevent_hdf5(filename_1);

  fprintf(stderr,"\n N EVENTS IN %s = %u \n",filename_0,nevent_0);

  size_t block_size = 2000;
  /* size_t block_size = 1999; */
  // block size is the size of the "pool" events are mixed in. 2000 is a good start as it reaches a stable
  // pairing solution with events that are very similar. Block size can be reduced to 1000 if limited N events.

  int max_remainder_events = block_size/4;
  //max number of events not used if nevents/block_size don't divide evenly. "N/4" is somewhat arbitrary 

  while (nevent_0 % block_size > max_remainder_events){
    block_size --;
    //ensure n_events/blocksize isn't too large. This will correspond to triggered events
    //at the end of the root file without pairings.
  }

  /* const size_t nblocks_0 = 2; */
  const size_t nblocks_0 = nevent_0 / block_size;
  std::vector<std::vector<float> >feature_0_vec;

  //events from dataset_0 (usually the triggered dataset) make up the "block" structure. The number of
  // blocks is roughly Nevents/block_size. Each block of dataset_0 is eventually mixed with "n_mix" blocks of
  // dataset_1 (MB/other data), such that each triggered event has n_mix pairings, and each event pairing 
  // comes from a unique block-pairing.

  for(size_t h = 0; h < nblocks_0; h++){

    size_t event_start_0 = h * block_size;
    size_t event_end_0 = event_start_0 + block_size;

    feature_0_vec.push_back(feature_extract_hdf5(filename_0, event_start_0, event_end_0,nfeature));
    /* feature_0_vec.push_back(feature_extract(filename_0, event_start_0, event_end_0,nfeature)); */

    fprintf(stderr,"%d: GAMMA EVENT START=%u || EVENT END=%u\n",
        __LINE__,event_start_0,event_end_0);
  }

  //n_mix dictates the number of blocks in MB data, since events in the second dataset can be reused for different
  //events in the first dataset.. So, each mixing iterator, i, indicates an
  //iteration through n_mix blocks of MB data. Each and every triggered block of events is mixed
  //with n_mix blocks of MB data.

  size_t nblocks_1 = mix_end-mix_start;
  fprintf(stderr,"\n%d: mix end = %zu, mix start = %zu\n",__LINE__,mix_end,mix_start);
  /* const size_t nblocks_1 = ((nevent_1 * nduplicate)/block_size); */

  /* if (nblocks_1 < ((nevent_1 * nduplicate)/block_size)){ */
  /*   fprintf(stderr,"Not enough events in dataset 1 (hdf5 file) to mix %zu blocks at a blocksize of %zu \n",nblocks_1,block_size); */
  /*   exit(EXIT_FAILURE); */
  /*   /1* FIXME: add in a while loop to iterate the blocksize, as we did above *1/ */
  /* } */

  //For even distribution of MB events
  //assumes dataset_1 (min-bias) is larger than triggered dataset, and spreads remainder throughout the dataset
  //Note: remainder is only used here. Afterwards, the block of data is taken from the mixing index in the main loop
  //where consecutive blocks of dataset_1 are alreay separated by this remainder.
  int remainder_1 = (nevent_1-block_size*nblocks_1)/nblocks_1;
  remainder_1=0;
  fprintf(stderr,"%d: nblocks_1 = %zu, remainder_1 = %i\n",__LINE__,nblocks_1,remainder_1);

  std::vector<std::vector<float> >feature_1_vec;

  for (size_t i = mix_start; i < mix_end; i++) {

    size_t event_start_1 = i * (block_size + remainder_1);
    /* size_t event_start_1 = i * nevent_1 / (nblocks_1);  //old */
    size_t event_end_1 = event_start_1 + block_size;
    /* size_t event_start_1 = i * nevent_1 / (nblocks_1+1); */
    /* size_t event_end_1 = (i + nduplicate) * nevent_1 / */
    /*       (nblocks_1+ nduplicate); */

    fprintf(stderr,"%d: MB EVENT START=%u || EVENT END=%u \n",
        __LINE__,event_start_1,event_end_1);   

    feature_1_vec.push_back(feature_extract_hdf5(filename_1, event_start_1, event_end_1,nfeature));
  }

  /* fprintf(stderr,"%d: Number of Threads = %i",__LINE__,omp_get_num_threads()); */
#pragma omp parallel for ordered schedule(dynamic)
  //each thread will independently pair data_0 a single block with all [n_mix] blocks from data_1.
  //all relevant data are at this point in thread-safe vectors, root is no longer used.

  for(size_t h = 0; h < nblocks_0; h++){

    fprintf(stderr,"%s:%d: %s %lu %s %lu\n",__FILE__,__LINE__,"Block",h,"of",nblocks_0);	   

    size_t event_start_0 = h * block_size;
    std::vector<std::vector<Long64_t> > k;		  

    std::vector<float> feature_0_scaled = feature_0_vec[h];

    for (size_t i = mix_start; i < mix_end; i++) {

      //normalize datasets
      std::vector<float> feature_1_scaled = feature_1_vec[i];
      feature_normalize(feature_0_scaled, feature_1_scaled,nfeature);

      //create preference list
      std::vector<std::list<index_t> > preference_0;
      std::vector<std::list<index_t> > preference_1;

      order_preference(preference_0, preference_1, feature_0_scaled,
          feature_1_scaled, nfeature, nduplicate);

      std::vector <index_t> m;
      m = gale_shapley(preference_0, preference_1);


      const size_t feature_1_size_nfeature = feature_1_vec[i].size() / nfeature; //basically block_size_1
      /* fprintf(stderr,"\n\n %s:%d feature_1_size_nfeature = %zu\n",__FILE__,__LINE__,feature_1_size_nfeature); */

      size_t event_start_1 = i * (block_size + remainder_1);
      size_t event_end_1 = event_start_1+block_size; 
      /* size_t event_start_1 = i * nevent_1 / (nblocks_1+1); */
      /* size_t event_end_1 = (i + nduplicate) * nevent_1 / */
      /*     (nblocks_1+ nduplicate); */
      /* size_t event_start_1 = i * nevent_1 / (nblocks_1);  //old */

      //output of m = gale_shapley are consecutive elements of dataset_1. Need to convert back to raw event#
      /* fprintf(stderr," %s:%d event_start_1 for mixing map output = %zu; Next = %zu\n\n",__FILE__,__LINE__,event_start_1,event_end_1); */
      std::vector <Long64_t>tempflat;

      //This basically restructures the output,
      //such that all the mixed events for the same triggered event 
      //are together in a vector of dimension: [n_events][n_mix]
      //
      for (size_t j = 0; j < m.size(); j++) {                                        
        //following line is crucial: converts element index in m = gale_shapeley
        //to raw event number from dataset_1 ("MB" dataset). Since m only is filled
        //with elements, the correct event start for that block is added).
        Long64_t q = event_start_1 + (m[j] % feature_1_size_nfeature);
        /* fprintf(stderr,"%d: Index = %zu, Event Start = %zu, m[j] = %zu, Mixed Event q  = %zu\n",__LINE__,i,event_start_1,m[j] % feature_1_size_nfeature, q); */
        tempflat.push_back(q);
      }		
      k.push_back(tempflat);

    }//mixed events

    for (size_t j = 0; j < k[0].size(); j++)
    {
      std::vector <Long64_t> P;

      size_t event_num =  event_start_0+j;

      for (size_t l = 0; l < k.size(); l++){
        P.push_back(k[l][j]);
      }	     
#pragma omp critical
      Matches[event_num] = P;
    }	    
  }

  clock_t end = clock();
  time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
  fprintf(stderr,"\n\n\n%d: Time elpased is %f seconds\n",__LINE__,time_spent);
  return Matches;
}

void write_txt(std::map<size_t,std::vector<Long64_t> > Matches,
    const char *filename_0, const char *mixing_start, 
    const char *mixing_end, const char *GeV_Track_Skim)
{

  fprintf(stderr,"\n%d: WRITING TO TEXT FILE, GALE DONE!\n",__LINE__);

  size_t mix_start = atoi(mixing_start);
  size_t mix_end = atoi(mixing_end);
  int Track_Skim = atoi(GeV_Track_Skim);

  size_t lastindex = std::string(filename_0).find_last_of("."); 
  std::string rawname = std::string(filename_0).substr(0, lastindex);
  FILE * txtfile = fopen (Form("./%s_%iGeVTrack_Pairs_%lu_to_%lu.txt",rawname.data(),Track_Skim,mix_start,mix_end),"w");
  for (size_t t=0; t<Matches.size();t++){ 
    for (size_t s=0; s<Matches[t].size();s++){
      fprintf(txtfile, "%lld\t", Matches[t][s]);
    }
    fprintf(txtfile, "%s\n\n","");
  }
  fclose (txtfile);
}

void write_root(std::map<size_t,std::vector<Long64_t> > Matches, 
    const char *filename_0, const char *GeV_Track_Skim, unsigned int n_mix_events)
{

  fprintf(stderr,"\n%d: Paring Done, writing to ROOT file \n",__LINE__);
  int Track_Skim = atoi(GeV_Track_Skim);

  TFile *root_file = new TFile(filename_0,"update");
  TTree *hi_tree = dynamic_cast<TTree *>(root_file->Get(HI_TREE));
  if (hi_tree == NULL) {
    hi_tree = dynamic_cast<TTree *>(root_file->Get(HI_TREE_2));		
    if(hi_tree == NULL){
      fprintf(stderr, "%s:%d: TREE FAIL\n",__FILE__, __LINE__);
      return;
    }
  }

  size_t lastindex = std::string(filename_0).find_last_of("."); 
  std::string rawname = std::string(filename_0).substr(0, lastindex);

  TFile *newfile = new TFile(Form("%s_%luGeVTrack_paired_hdf5.root",rawname.data(),Track_Skim),"recreate");
  /* fprintf(stderr,"%d: Right Before Clone\n",__LINE__); */
  TTree *newtree = hi_tree->CloneTree(0);

  ULong64_t nentries = hi_tree->GetEntries();    
  unsigned int n_mix = Matches[0].size();
  Long64_t Mix_Events[n_mix];

  /* fprintf(stderr,"%d: Before ME branch is made\n",__LINE__); */
  TBranch *MixE = newtree->Branch("mixed_events", Mix_Events, Form("&mixed_events[%ui]/L",n_mix));

  for (ULong64_t t = 0; t<nentries;t++){ //Event # is key used in map <Matches>

    /* fprintf(stderr,"%d: Writing for event number %i  \n",__LINE__,t); */
    hi_tree->GetEntry(t);

    if (t < Matches.size()) {
      for (size_t s=0; t<Matches[0].size();t++)
      { 
        fprintf(stderr,"%d: %i>%i Filling With Mixed Event %i\n",__LINE__,t,s,Matches[t][s]);
        Mix_Events[s]=Matches[t][s]; 
      }
    }

    else
      for(size_t s = 0; s<n_mix; s++)
      {
        /* fprintf(stderr,"%d: %i>%i Filling With NAN \n",__LINE__,t,s); */
        /* Mix_Events[s] = NAN; */
        Mix_Events[s] = -999;
      }
    /* if(t < Matches.size()){ */

    /*   for (size_t s=0; s<(Matches[t]).size();s++){ */
    /*     fprintf(stderr,"%d: %i Mixed Event Number = %lu \n",__LINE__,s,Matches[t][s]); */
    /*     Mix_Events[s]=Matches[t][s]; */ 
    /*   } */
    /*   fprintf(stderr,"\n"); */
    /* } */	    
    /* //small remainder of unpared events due to block structure */
    /* else if (t >= Matches.size()){ */
    /*   for(size_t u = 0; u<n_mix_events; u++){ */
    /*     fprintf(stderr,"%d: Right Before Fill \n",__LINE__); */
    /*     Mix_Events[u] = NAN; */
    /*   } */
    /* } */
    /* fprintf(stderr,"%d: Right Before Fill \n",__LINE__); */
    newtree->Fill();  

  }//entries

  fprintf(stderr,"%d: Right Before Write\n",__LINE__);
  newtree->Write();

  delete root_file;
  delete newfile;

  gSystem->Exit(0);
}

int main(int argc, char *argv[])
{
  if (argc < 6) {
    fprintf(stderr,"%s\n","Argument Syntax is [Command] [File] [File] [mix start] [mix end] [GeV Track Skim]");
    return EXIT_FAILURE;
  }

  const char *root_file = argv[1];
  const char *hdf5_file = argv[2];
  const char *mix_start = argv[3];
  const char *mix_end   = argv[4];
  const char *track_skim= argv[5];
  int n_event_properties = 3; //z-vertex, multiplicity, v2 (flow)
  int n_duplicate = 1;

  std::map<size_t,std::vector<Long64_t> > Matches = mix_gale_shapley(root_file, hdf5_file,
      mix_start, mix_end, track_skim, n_event_properties, n_duplicate);

  write_txt(Matches,root_file,mix_start,mix_end,track_skim);
  /* write_root(Matches,root_file,track_skim,mix_end-mix_start); */
}
