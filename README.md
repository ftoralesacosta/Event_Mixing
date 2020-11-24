# Event_Mixing

This is a framework for event mixing using a Gale Shapely Stable Matching Algorithm to pair events that are similar. This framework is intended to be used with NTuples created with __LINK__. There are two major advantages of using this over standard pooled or bin mixing:

1. Once set up, it is much, much faster
2. Events are paired to be more similar on average

The trade-off, however, is larger overhead for this performance. One of the NTuples (usually minimum bias dataset) must be converted to HDF5 file format. The nature of mixing us to randomly sample events to some extend. ROOT files are in general optimized for impressive serial performance (the documentation makes reference to the use of "buckets" of event). This optimization, however, greatly inhibits random access speeds. HDF5 is designed for very fast random access, and plays nicely with parallelization. The steps to run the event mixing, in order are:

1. Convert one dataset to HDF5
2. Pair events from the remaining ROOT file with events in the HDF5 file
3. Perform your mixed event analysis, using the event parings