# Event_Mixing

This is a framework for event mixing using a Gale Shapely Stable Matching Algorithm to pair collider events that are similar. This framework is intended to be used with NTuples created using [this repo](https://github.com/alwina/ntuple-gj). Within the UCB group, this is known as the NTupilizer. There are two major advantages of using this over standard pooled or bin mixing:

1. Once set up, it is much, much faster (3+days -> 2hrs)
2. Events are paired to be more similar on average

The trade-off, however, is larger overhead for this performance. One of the NTuples (usually minimum bias dataset) must be converted to HDF5 file format. The nature of mixing us to randomly sample events to some extend. ROOT files are in general optimized for impressive serial performance (the documentation makes reference to the use of "buckets" of event). This optimization, however, greatly inhibits random access speeds. HDF5 is designed for very fast random access, and plays nicely with parallelization. The steps to run the event mixing, in order are:

1. Convert one dataset to HDF5
2. Pair events from the remaining ROOT file with events in the HDF5 file
3. Perform your mixed event analysis, using the event parings

##Pairing:
The pairing is the most complex aspect. The stable pairing algorithm, in parallel_pair_gale_shapley_hdf5.cc works as follows:
1. Creates a ranked list of all events. In the case of triggered and Min-Bias pairings, the triggered dataset has a ranked list of all MB events according to how similar those events are (multiplicity, z-vertex, event-plane angle). The number of criteria for this event pairing is _nfeature_, which is set in line 817
2.Once the list is created, it will iterate through all the triggered events. It will find the highest ranked MB event on the current triggered events list
3. If that MB is not already paired, the pairing is recorded
4. If that MB is already paired, it checks the MB pairing list to see which trigger event should go to that MB event
5. Continue iterating.

This algorithm is done in blocks of 2k or 4k events, because the algorithm is of complexity n^2, and would take too long to run over the entire dataset. If the pairing performance is not sufficient, increasing the block size should be the first step. The block structure makes the loops a little harded to follow, but basically the hdf5 data is read 1 block at a time, and as you iterate through events and reach the end of a block, the code will then pull the next chunk of data out of the hdf5 file.

##Converting root to hdf5:
to_hdf5: [command] [root_file] [new hdf5 file name]

This code converts the NTuples to hdf5 files. The dimensions of the hdf5 datasets are quite rigid, and can't be changed once initialized. For this reason, the dimensions are set to the max number of clusters, jets, or tracks depending on the dataset, and fills the data structure with NaN if it does not reach the max for that event. The conversion code has jet, cluster, and track data, as well as some event data. Please edit this code to include whichever physics object you need, but please keep in mind the rank of the dataset is 3 (event #, track #, track variable/value). This means that NTuples branches of rank 2 or higher won't fit in this dataset easily without being converted to flat arrays. An easier solution is to just save each element in the higher order branches as it's own variable in the hdf5 file. This is done for the cluster_lambda_square variable (aka sigma^2_long) variable in lines 363 and line 388. The latter was added much later, and so the two varibales are not adjacent.

When adding a varibale in general, make sure to change the **row_size** constants on lines 309 and 606.

##Running the pairing:
[Command] [Triggerd H5 File] [Min-Bias H5 File] [mix start] [mix end] [GeV Track Skim]
The paring should not be run on the login node on CORI. Please at least use the interactive node for best parallel performance.

###set up interactive node on CORI###
**salloc -N 1 -C haswell -q interactive -t 04:00:00**
this requests a single interactive node (64 threads). It may take a bit. The group allocation is 64 nodes, and any user can request up to two interactive sessions.

**export HDF5_USE_FILE_LOCKING=FALSE**
This makes the HDF5 files able to read by multiple threads at the same time (read, NOT write).

parallel_pair_gale_shapley_hdf5: [Command] [Triggerd H5 File] [Min-Bias H5 File] [mix start] [mix end] [GeV Track Skim]
  - the last argument can be set to zero, this is meant to skim the dataset for high pT tracks (or jets)


##Improvements##
A lot of annoying work arounds were done to ensure the pairing finishes within the 4 hour timeslot for an interactive node. But this is mostly due to the use of only a single node of 64 threads because OpenMP parallelization was used. Changing this to MPI parallelization, which does not requere the memory to be shared between nodes unlike OpenMP, would mean being able to run on 64nodes at 64 threads each. That means you could mix up to 4096 blocks of events simoultaneously!!! If you are attempting this, or have any other questions, please email me at fernando_tta@berkeley.edu . Yue Shi Lai will also be able to help if you ask him some precises questions :] .
