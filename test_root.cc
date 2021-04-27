#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

int main(int argc, char *argv[])
{ 
  TString file_name = argv[1];

    UInt_t nevent_max = 0;
    UInt_t ntrack_max= 0;
    UInt_t ncluster_max= 0;
    UInt_t njet_max= 0;


        // Cautious opening of the TTree, capturing all modes of
        // failure, and keep the TDirectoryFile (to be deleted later)
        // to avoid memory leak

        TFile *file = TFile::Open(file_name);
        if (file == NULL) {
            return 1;
        }

        TDirectoryFile *df = dynamic_cast<TDirectoryFile *>
            (file->Get("AliAnalysisTaskNTGJ"));

        if (df == NULL) {
            fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Cannot open TFile");
            return 1;
        }

        TTree *hi_tree = dynamic_cast<TTree *>
            (df->Get("_tree_event"));

        if (hi_tree == NULL) {
            fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Cannot open _tree_event");
            return 1;
        }

        UInt_t ntrack;
        UInt_t ncluster;
        UInt_t njet_ak04tpc;
        if (nevent_max == 0)
          nevent_max = UInt_t(hi_tree->GetEntries());

        hi_tree->SetBranchAddress("ntrack", &ntrack);
        hi_tree->SetBranchAddress("ncluster", &ncluster);
        hi_tree->SetBranchAddress("njet_ak04tpc", &njet_ak04tpc);

        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Obtaining ntrack, ncluster, and njet max for hdf5 file");

        //for (Long64_t i = 0; i < hi_tree->GetEntries(); i++) {
        for (Long64_t i = 0; i < nevent_max; i++) {
          //if ( i>=7893 ) continue;
          hi_tree->GetEntry(i);
          std::cout<<"\nntrack_max = "<<ntrack_max<<std::endl;
          std::cout<<"ntrack = "<<ntrack<<std::endl;
          std::cout<<"ncluster = "<<ncluster<<std::endl;
          std::cout<<"njet = "<<njet_ak04tpc<<std::endl;

          ntrack_max = std::max(ntrack_max, ntrack);
          ncluster_max = std::max(ncluster_max, ncluster);
          njet_max = std::max(njet_max, njet_ak04tpc);
          fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, i,nevent_max);
        }

        fprintf(stderr, "\n");

        // Fully delete everything
        hi_tree->Delete();
        delete df;
        file->Close();
        delete file;
    return 0;
}
