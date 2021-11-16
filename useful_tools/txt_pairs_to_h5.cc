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

#define RANK 2

#include <vector>
#include <math.h>

using namespace H5;

// Code injects mixed event from txt file into triggered hdf5 file.
// This lets one edit hdf5 information and re-use the pairing, so long as
// the information does not affect the pairing process (i.e. does not edit
// event information).
int main(int argc, char *argv[])
{
  if (argc < 5) {
    fprintf(stderr,"Batch Syntax is [Triggered HDF5] [pairing.txt] [Mix Start] [Mix End]");
    exit(EXIT_FAILURE);
  }

  //Triggered HDF5 File for augmenting
  const H5std_string hdf5_file_name(argv[1]);
  H5File h5_file( hdf5_file_name,  H5F_ACC_RDWR);//read+write
  TString hdf5_file = (TString)argv[1];
  std::cout << "Opening: " << hdf5_file << std::endl;

  //Text File with Pairings 
  std::ifstream pairing_textfile;
  std::ostringstream filename;
  filename << argv[2];
  pairing_textfile.open(argv[2]);
  std::cout<<"Opened Text File: "<<argv[2]<<std::endl;

  //Mixing Range
  size_t mix_start = atoi(argv[3]);
  size_t mix_end = atoi(argv[4]);
  fprintf(stderr,"Mixing Event range is %i to %i \n",mix_start,mix_end);

  //-------------------------- Get Number of Events ------------------------------------------------
  
  fprintf(stderr,"Opening event dataset\n");
  const std::string event_ds_name( "event" );
  H5::DataSet event_dataset = h5_file.openDataSet( event_ds_name.c_str() );
  H5::DataSpace event_dataspace = event_dataset.getSpace();
  H5::DSetCreatPropList mixing_property = H5::DSetCreatPropList();

  //Initialize Event Dimensions
  const int event_ndims = event_dataspace.getSimpleExtentNdims();
  hsize_t event_maxdims[event_ndims];
  hsize_t eventdims[event_ndims];
  event_dataspace.getSimpleExtentDims(eventdims, event_maxdims);
  size_t nEvents = eventdims[0];
  fprintf(stderr,"Number of events in hdf5 = %i\n",eventdims[0]);

  //-------------------------- Create Mixed Event DataSet ------------------------------------------------
  size_t mixing_row_size = 300; //N Mixed Events
  UInt_t block_size = 2000;
  hsize_t mixing_dim_extend[RANK] = {block_size, mixing_row_size };
  hsize_t mixing_dim_max[RANK] = { H5S_UNLIMITED, mixing_row_size };
  H5::DataSpace mixing_data_space(RANK, mixing_dim_extend, mixing_dim_max);

  //-------------------------- Compress DataSet if Possible --------------------------
#ifdef HDF5_USE_DEFLATE
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
      mixing_property.setDeflate(1);
    }   
  }   
#endif // HDF5_USE_DEFLATE

  //----------------------- Finish Mixed Event DataSet --------------------------

  hsize_t mixing_dim_chunk[RANK] = {
    mixing_dim_extend[0],
    mixing_dim_extend[1],
  };

  mixing_property.setChunk(RANK, mixing_dim_chunk);

  H5::DataSet mixing_data_set = h5_file.createDataSet("mixing", 
      H5::PredType::NATIVE_FLOAT, mixing_data_space, mixing_property);
  fprintf(stderr,"%s:%d: CREATED EVENT DATASET",__FILE__,__LINE__);

  hsize_t offset[RANK] = {0, 0};

  //vector to temporarily hold mixed events
  std::vector<float> mixing_data(block_size * mixing_row_size, NAN);

  Long64_t mix_range = (mix_end - mix_start);


  //-------------------------- LOOP --------------------------
  for(Long64_t ievent = 0; ievent < nEvents ; ievent++){

    int iblock = ievent % block_size;
    fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nEvents);

    std::string eventline;
    if (ievent > 0){ 
      //skips \n that separates each triggered event's pairings list
      getline(pairing_textfile, eventline); 
    }   
    //Should now grab the list of MB events from the text file if not eof
    getline(pairing_textfile, eventline);

    if (eventline.size() == 0 ) 
    {   
      for (int m = 0; m < 300; m++)
        mixing_data[iblock * mixing_row_size + m] = -999;
    }

    else
    {
      std::string mixednum_string;
      long mixednum;
      std::istringstream parser[1];
      parser[0].str(eventline);
      int currentindex 0;
      // Loop over mixed events, and fill the mixed_events histogram too
      //getline auto-increments to next value, specified by separator "\t"
      for(int m = 0; m < mix_range; m++) {
        getline(parser[currentindex], mixednum_string, '\t');
        mixing_data[iblock * mixing_row_size + m] = stoul(mixednum_string);
        /* std::cout<<currentindex<<" "<<mixednum_string<<std::endl; */
      }   
    }   

    if (iblock == (block_size-1)) {
      if (offset[0] == 0)  
        mixing_data_set.write(&mixing_data[0], H5::PredType::NATIVE_FLOAT);

      //extend dataset otherwise
      else { 
        const hsize_t mixing_dim_extended[RANK] = { 
          offset[0] + mixing_dim_extend[0], mixing_dim_extend[1] };

        //Extend to new Dimension
        mixing_data_set.extend(mixing_dim_extended);
        H5::DataSpace mixing_file_space = mixing_data_set.getSpace();

        mixing_file_space.selectHyperslab(
            H5S_SELECT_SET, mixing_dim_extend, offset);

        H5::DataSpace mixing_memory_space(RANK, mixing_dim_extend, NULL);
        mixing_data_set.write(&mixing_data[0], H5::PredType::NATIVE_FLOAT,
            mixing_memory_space, mixing_file_space);
      }   

      //put in block logic
      offset[0] += block_size;
    }
  }//events

  return EXIT_SUCCESS;
}

