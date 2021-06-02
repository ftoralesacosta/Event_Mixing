/**
  This program clones an NTuple, then uses data contained in text files to addmixed events to the clone
  */
// Author: Fernando & Ivan Chernyshev; Date: 6/18/2018

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
#include <sstream>

#define NTRACK_MAX (1U << 15)

#include <vector>
#include <math.h>

int fileArg = 1;
int runArg = 2;

int main(int argc, char *argv[])
{
  if (argc < 5) {
    fprintf(stderr,"Syntax is [Command] [root file] [text file] [Mix Start] [Mix End]\n");
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");

  size_t mix_start = stoull(std::string(argv[3]));
  size_t mix_end = stoull(std::string(argv[4]));

  std::cout << "Opening: " << (TString)argv[fileArg] << std::endl;
  TFile *file = TFile::Open((TString)argv[fileArg]);

  if (file == NULL) {
    std::cout << " fail" << std::endl;
    exit(EXIT_FAILURE);
  }
  file->Print();

  TTree *_tree_event = NULL;
  _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
  if (_tree_event == NULL) {
    std::cout << "Failed to grab tree, perhaps AliAnalysisTaskNTGJ does not exist, trying again" << std::endl;
    _tree_event = dynamic_cast<TTree *> (file->Get("_tree_event"));
    if (_tree_event == NULL) {
      std::cout << " fail " << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  //_tree_event->Print();
  std::cout<<"TTree successfully acquired" << std::endl;
  std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;

  // New file
  size_t lastindex = std::string(argv[fileArg]).find_last_of("."); //removes .root
  std::string rawname = std::string(argv[fileArg]).substr(0, lastindex);
  std::cout<<rawname<<std::endl;
  TFile *newfile = new TFile(Form("%s_paired.root", rawname.data()), "RECREATE");
  TTree *newtree = _tree_event->CloneTree(0);
  newtree->SetMaxTreeSize(1000000000000LL); //1TB max
  //new branch: mixed_events
  Long64_t mixed_events[300];
  newtree->Branch("mixed_events", mixed_events, "mixed_events[300]/L"); // One more entry needed for this to work

  std::cout<< "New branch successfully created " <<std::endl;

  // Get the mixed event textfiles
  std::ifstream mixed_textfile;

  std::ostringstream filename;
  filename << argv[2];
  mixed_textfile.open(argv[2]);
  std::cout<<"Opened Text File: "<<argv[2]<<std::endl;

  const Long64_t nevents = _tree_event->GetEntries();
  /* const Long64_t nevents = 10; */ 
  // Loop over events
  for(Long64_t ievent = 0; ievent < nevents ; ievent++){

    fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, _tree_event->GetEntries());
    _tree_event->GetEntry(ievent);
    
    std::string eventline;
    if (ievent > 0){
      //skips \n that separates each triggered event's pairings list
      getline(mixed_textfile, eventline); 
    }
    //Should now grab the list of MB events from the text file if not eof
    getline(mixed_textfile, eventline);

    if (eventline.size() == 0 )
    {
      for (int m = 0; m < 300; m++)
      mixed_events[m] = -999;
    }

    else
    { /*if eventline != "" */
      std::string mixednum_string;
      long mixednum;
      std::istringstream parser[1];
      parser[0].str(eventline);
      int currentindex;
      Long64_t mix_range = (mix_end - mix_start);
      // Loop over mixed events, and fill the mixed_events histogram too
      for(int m = 0; m < mix_range; m++) {
        currentindex = m/mix_range;
        getline(parser[currentindex], mixednum_string, '\t');
        mixed_events[m] = stoul(mixednum_string);
        /* std::cout<<currentindex<<" "<<mixednum_string<<std::endl; */
      }
    }
    newtree->Fill();
  }
  std::cout << "Successfully exited the eventloop" << std::endl;
  /* newtree->AutoSave(); */
  newtree->Write();
  std::cout << "Successful autosave" <<std::endl;
  delete newfile;
  /* std::cout << "Deleted newfile" << std::endl; */

  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}
