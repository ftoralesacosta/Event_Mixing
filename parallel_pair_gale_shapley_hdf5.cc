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
#include "H5Cpp.h"

//#define HI_TREE "hiEvtAnalyzer/HiTree"
#define HI_TREE "_tree_event"
#define HI_TREE_2 "AliAnalysisTaskNTGJ/_tree_event"
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

		double v[3];

		hi_tree->SetBranchAddress("primary_vertex", v);

		float multiplicity_v0[64];
		  //float centrality_v0m;

		switch (nfeature) {
		case 2:
			hi_tree->SetBranchAddress("multiplicity_v0", &multiplicity_v0);
			break;
		default:
			fprintf(stderr, "%s:%d: illegal nfeature = %lu\n",
					__FILE__, __LINE__, nfeature);
			return std::vector<float>();
		}

		std::vector<float> ret;

		for (size_t i = event_start; i < event_end; i++) {
			hi_tree->GetEntry(i);

			ret.push_back(v[2]);

			if (nfeature >= 2) {
 			        float multp_sum = 0;
			        for (int k = 0; k < 64; k++) {
				  multp_sum += multiplicity_v0[k];
				}
				ret.push_back(multp_sum);
				//fprintf(stderr,"%d: z-vertex = %1.2f mp = %f\n",__LINE__,v[2],multp_sum);
			}

		}

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

	  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Befor HDF5 READ");
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

		//Event dataspace rank only rank 2
		UInt_t NEvent_Vars = eventdims[1];
		//Define array hyperslab will be read into 

		//const hsize_t block_size = eventdims[0];
		const hsize_t block_size = event_end-event_start;
		float event_data_out[block_size][NEvent_Vars] = {0};
		fprintf(stderr,"%s:%d: Block Size = %i, Event Start = %i, NVars = %i\n",__FILE__,__LINE__,block_size,event_start,NEvent_Vars);

		//hsize_t event_offset[2] = {0,0};

		hsize_t event_offset[2] = {event_start,0};
		hsize_t event_count[2] = {block_size, NEvent_Vars};
		const int Event_RANK_OUT = 2; //Different for Events       
		H5::DataSpace event_memspace( Event_RANK_OUT, eventdims );
		hsize_t event_offset_out[2] = {0,0};
		hsize_t event_count_out[2] = {block_size, NEvent_Vars};

		event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
		event_memspace.selectHyperslab( H5S_SELECT_SET, event_count_out, event_offset_out );
		//event_dataset.read( event_data_out, H5::PredType::NATIVE_FLOAT, event_memspace, event_dataspace);
		fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "event dataset read into array: OK");

		std::vector<float> ret;

		//for (size_t i = event_start; i < event_end; i++) {
		  for (size_t i = 0; i < block_size; i++) {
		    float v2 = event_data_out[i][0];
		    ret.push_back(v2);
			if (nfeature >= 2) {
			  float multp = event_data_out[i][1];
			  ret.push_back(multp);
			  //fprintf(stderr,"%d: z-vertex = %1.2f mp = %f\n",__LINE__,v2,multp);
			}
		}
		return ret;
	}

	void feature_normalize(std::vector<float> &u,
			       std::vector<float> &v, const size_t n)
	{
		std::vector<float> s(n, 0);

		for (size_t j = 0; j < n; j++) {
			float s_j = 0;

			for (size_t i = 0; i < u.size(); i += n) {
				s_j += fabsf(u[i + j]);
			}
			s[j] = s_j;
		}
		for (size_t j = 0; j < n; j++) {
			float s_j = 0;

			for (size_t i = 0; i < v.size(); i += n) {
				s_j += fabsf(v[i + j]);
			}
			s[j] += s_j;
		}
		for (size_t j = 0; j < n; j++) {
			s[j] = (u.size() + v.size()) / s[j];

			//fprintf(stderr, "%s:%d: %lu %f\n", __FILE__, __LINE__, j, s[j]);
		}

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
	const size_t u_size_n = u.size() / n;
	const size_t v_size_n = v.size() / n;

	//fprintf(stderr, "%s:%d: %lux%lu Thread %i\n", __FILE__, __LINE__,
	//u_size_n, v_size_n, omp_get_thread_num());

	const size_t size_max = std::max(u_size_n, v_size_n);

	up.resize(u_size_n, std::list<index_t>());

	for (size_t i = 0; i < u_size_n; i++) {
		std::vector<std::pair<float, index_t> > l;
		
		for (size_t j = 0; j < v_size_n; j++) {
			float d = 0;

			for (size_t k = 0; k < n; k++) {
				d += std::pow(u[i * n + k] - v[j * n + k], 2);
			}
			if (d==0) d = 999999; //Avoid pairing identical events
			l.push_back(std::pair<float, size_t>(d, j));
		}
		std::sort(l.begin(), l.end(), preference_compare);
		// up.push_back(std::list<index_t>());
		for (size_t j = 0; j < l.size(); j++) {
			for (size_t k = 0; k < nduplicate_v; k++) {
				up[i].push_front(l[j].second + k * v_size_n);
			}
		}
		up[i].resize(size_max, v_size_n);
		// if (i % 100 == 0) {
		//   fprintf(stderr, "%s:%d: %lu/%lu\n", __FILE__, __LINE__,i, u_size_n);
		// }
	}

	vp.resize(v_size_n * nduplicate_v, std::list<index_t>());

	for (size_t j = 0; j < v_size_n; j++) {
		std::vector<std::pair<float, index_t> > l;

		for (size_t i = 0; i < u_size_n; i++) {
			float d = 0;

			for (size_t k = 0; k < n; k++) {
				d += std::pow(u[i * n + k] - v[j * n + k], 2);
			}
			if (d==0) d = 999999;
			l.push_back(std::pair<float, index_t>(d, i));
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

		// if (j % 100 == 0) {
		//   fprintf(stderr, "%s:%d: %lu/%lu\n", __FILE__, __LINE__,j, v_size_n);
		// }
	}

	fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Order Preference done for this block");
	//fprintf(stderr,"\n %d : Thread %i : orderfunction: HERE \n",__LINE__,omp_get_thread_num());
}

std::vector<index_t> gale_shapley(std::vector<std::list<index_t> > &mp,
				  std::vector<std::list<index_t> > &fp)
{
	std::vector<index_t> m_to_f_engaged(mp.size(), fp.size());
	std::vector<index_t> f_to_m_engaged(fp.size(), mp.size());

	std::vector<std::vector<std::pair<
	            std::vector<std::list<index_t> >::iterator,
		    std::list<index_t>::iterator> > > mp_index;

	mp_index.resize(fp.size());

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

		// if ((iterator_outer - mp.begin()) % 100 == 0) {
		//   fprintf(stderr, "%s:%d: %lu/%lu\n", __FILE__, __LINE__,
		// 	  iterator_outer - mp.begin(), mp.size()); }
		
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

		// if (m % 100 == 0 || w % 100 == 0) {
		//   fprintf(stderr, "%s:%d: %hu<>%hu\n", __FILE__, __LINE__,m, w);
		// }

		// Some man p is engaged to w
		index_t p = f_to_m_engaged[w];

		if (p != mp.size()) {
			// Assign p to be free
			m_to_f_engaged[p] = fp.size();
		}
		f_to_m_engaged[w] = m;
		m_to_f_engaged[m] = w;

		std::list<index_t>::iterator s =
			std::find(fp[w].begin(), fp[w].end(), m);

		s++;
		for (std::vector<std::pair<std::vector<std::list<index_t> >::
				 iterator, std::list<index_t>::iterator> >::iterator
				 iterator = mp_index[w].begin();
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

	const size_t nevent_0 = nevent(filename_0);
	size_t nevent_1 = nevent_hdf5(filename_1);

	//TEMPORARY MAX
	// if (nevent_1 > 582695)
	//   nevent_1 = 582695;
	
	fprintf(stderr,"\n N EVENTS IN %s = %u \n",filename_0,nevent_0);

	size_t block_size = 2000;
	// block size is the size of the "pool" events are mixed in. 2000 is a good start as it reaches a stable
        // pairing solution with events that are very similar. Block size can be reduced to 1000 if limited N events.

	int max_remainder_events = block_size/4;
	//max number of events not used if nevents/block_size don't divide evenly.  

	size_t nblocks_temp = (std::min(nevent_0, nevent_1 * nduplicate) / block_size);
	while (nevent_0%block_size > max_remainder_events);
	    block_size --;

	const size_t nblocks_0 = nevent_0 / block_size; 
	const size_t nblocks_1 = ((nevent_1 * nduplicate)/block_size);

	//For even distribution of MB events
	int remainder_1 = (nevent_1%block_size)/(mixing_end-mixing_start);

	const size_t n_mix = 300;
	
	std::vector<std::vector<float> >feature_0_vec;
	std::vector<std::vector<float> >feature_1_vec;

	for(size_t h = 0; h < nblocks_0; h++){
	  //for(size_t h = 0; h < 2; h++){

	  size_t event_start_0 = h * block_size;
	  size_t event_end_0 = event_start_0 + block_size;

	    feature_0_vec.push_back(feature_extract(filename_0, event_start_0, event_end_0,nfeature));
	    
	    fprintf(stderr,"\n %d: GAMMA EVENT START=%u || EVENT END=%u",
	    	    __LINE__,event_start_0,event_end_0);
	}

	for (size_t i = mix_start; i < mix_end; i++) {
	     	    
	  size_t event_start_1 = i * (block_size + remainder_1);
	  size_t event_end_1 = event_start_1 + block_size;
	  
	  fprintf(stderr,"\n %d: MB EVENT START=%u || EVENT END=%u \n",
		  __LINE__,event_start_1,event_end_1);   
	  feature_1_vec.push_back(feature_extract_hdf5(filename_1, event_start_1, event_end_1,nfeature));
	}

#pragma omp parallel for ordered schedule(dynamic)
	//each thread will independently pair events in each block 
	for(size_t h = 0; h < nblocks_0; h++){

	  fprintf(stderr,"%s:%d: %s %lu %s %lu\n",__FILE__,__LINE__,"Block",h,"of",nblocks_0);	   

	  size_t event_start_0 = h * block_size;
	  std::vector<std::vector<Long64_t> > k;		  

	  for (size_t i = mix_start; i < mix_end; i++) {

	    std::vector<float> feature_0_scaled = feature_0_vec[h];
	    std::vector<float> feature_1_scaled = feature_1_vec[i];
	    
	    feature_normalize(feature_0_scaled, feature_1_scaled,nfeature);
	    std::vector<std::list<index_t> > preference_0;
	    std::vector<std::list<index_t> > preference_1;

	    order_preference(preference_0, preference_1, feature_0_scaled,
			     feature_1_scaled, nfeature, nduplicate);

	    std::vector <index_t> m;

	    m = gale_shapley(preference_0, preference_1);

	    const size_t feature_1_size_nfeature = feature_1_vec[i].size() / nfeature;

	    size_t event_start_1 = i * nevent_1 / (nblocks_1+1);
	    std::vector <Long64_t>tempflat;

	    for (size_t j = 0; j < m.size(); j++) {                                        
	      Long64_t q = event_start_1 + (m[j] % feature_1_size_nfeature);
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
	FILE * txtfile = fopen (Form("../InputData/%s_%iGeVTrack_Pairs_%lu_to_%lu.txt",rawname.data(),Track_Skim,mix_start,mix_end),"w");
	for (size_t t=0; t<Matches.size();t++){ 
	  for (size_t s=0; s<Matches[t].size();s++){
	      fprintf(txtfile, "%lld\t", Matches[t][s]);
	    }
	    fprintf(txtfile, "%s\n","");
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
	  TTree *newtree = hi_tree->CloneTree(0);

	  ULong64_t nentries = hi_tree->GetEntries();    
	  Long64_t Mix_Events[n_mix_events];

	  TBranch *MixE = newtree->Branch("mixed_events", Mix_Events, "&mixed_events[300]/L");
	  
	  for (ULong64_t t = 0; t<nentries;t++){ //Event # is key used in map <Matches>

	    hi_tree->GetEntry(t);
	    
	    if(t < Matches.size()){

	      for (size_t s=0; s<(Matches[t]).size();s++){
	  	Mix_Events[s]=Matches[t][s]; 
	      }
	    }	    
	    //small remainder of unpared events due to block structure
	    else if (t >= Matches.size()){
	      for(size_t u = 0; u<n_mix_events; u++)
	  	Mix_Events[u] = 999999999;
	    }
	    
	    newtree->Fill();  
	    
	  }//entries
	  
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

	std::map<size_t,std::vector<Long64_t> >	Matches = mix_gale_shapley(argv[1], argv[2], argv[3], argv[4],argv[5], 2, 1);
	write_txt(Matches,argv[1],"0","300",argv[5]);
	write_root(Matches,argv[1],argv[5],300);
}
