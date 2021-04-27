/**
  This program clones an NTuple, then uses data contained in text files to addmixed events to the clone
  */
// Author: Fernando + Ivan Chernyshev; Date: 6/18/2018

// Syntax: ./mixed_injector <ROOT file for mixed events to be injected into goes here> <Run number goes here (13d, 13e, 13f, etc.)> <Track pair energy, in GeV (must be an integer or the program will fail>

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

int num_of_files = 15;



int fileArg = 1;
int runArg = 2;
int trackpairenergyArg = 2;

int main(int argc, char *argv[])
{
  if (argc < 3) {
    fprintf(stderr,"\nSyntax is [Command] [root file] [TrackSkim GeV] \nThis Assumes 15 .txt files \n\n");
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");


  std::cout << "Opening: " << (TString)argv[fileArg] << std::endl;
  TFile *file = TFile::Open((TString)argv[fileArg]);

  if (file == NULL) {
    std::cout << " fail" << std::endl;
    exit(EXIT_FAILURE);
  }
  file->Print();

  TTree *_tree_event = dynamic_cast<TTree *> (file->Get("_tree_event"));
  if (_tree_event == NULL) {
    std::cout << "Failed to grab tree, perhaps AliAnalysisTaskNTGJ does not exist, trying again" << std::endl;
    _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
    if (_tree_event == NULL) {
      std::cout << " fail " << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  //_tree_event->Print();
  std::cout<<"TTree successfully acquired" << std::endl;
  std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;

  // New file
  size_t lastindex = std::string(argv[fileArg]).find_last_of("."); 
  std::string rawname = std::string(argv[fileArg]).substr(0, lastindex);
  std::cout<<rawname<<std::endl;
  TFile *newfile = new TFile(Form("%s_%iGeVTrack_paired.root", rawname.data(), std::stoi((std::string)argv[trackpairenergyArg])), "RECREATE");
  TTree *newtree = _tree_event->CloneTree(0);

  // TFile *newfile = new TFile(Form("%s_mixedadded_output.root", ((std::string)argv[runArg]).c_str()), "RECREATE");
  // TTree *newtree = _tree_event->CloneTree(0);

  //new branch: mixed_events
  Long64_t mixed_events[NTRACK_MAX];
  newtree->Branch("mixed_events", mixed_events, "mixed_events[300]/L"); // One more entry needed for this to work

  std::cout<< "New branch successfully created " <<std::endl;

  // Get the mixed event textfiles
  std::ifstream mixed_textfiles[num_of_files];
  for(int i = 0; i < num_of_files; i++) {
    std::ostringstream filename;
    filename << Form("%s_%iGeVTrack_Pairs_%i_to_%i.txt", rawname.data(), std::stoi((std::string)argv[trackpairenergyArg]), i*20, ((i+1)*20)-1);
    //filename << Form("%s_%iGeVTrack_Pairs_%i_to_%i.txt", rawname.data(), std::stoi((std::string)argv[trackpairenergyArg]), 0, 1);
    mixed_textfiles[i].open(filename.str());
    std::cout<<"Opening Text File: "<<Form("%s_%iGeVTrack_Pairs_%i_to_%i.txt", rawname.data(), std::stoi((std::string)argv[trackpairenergyArg]), i*20, ((i+1)*20)-1)<<std::endl;
  }


  const Long64_t nevents = _tree_event->GetEntries();
  // Loop over events
  for(Long64_t ievent = 0; ievent < nevents ; ievent++){
    //for(Long64_t ievent = 0; ievent < 2000; ievent++){
    fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, _tree_event->GetEntries());
    _tree_event->GetEntry(ievent);
    // Get the appropriate line from each file, break out of the loop if you hit an empty file
    std::string eventlines[num_of_files];
    bool event_end = false;
    for(int i = 0; i < num_of_files; i++) {
      getline(mixed_textfiles[i], eventlines[i]);
      if (eventlines[i] == "") {
        event_end = true;
        std::cout<<std::endl<<"reached end of file: "<<i<<std::endl;
        break;
      }
    }
    if(event_end)
      break;

    //try {
    std::string mixed_num_string;
    long mixed_num;
    std::istringstream parsers[num_of_files];
    for(int i = 0; i < num_of_files; i++) {
      parsers[i].str(eventlines[i]);
    }
    int currentindex;
    // Loop over mixed events and fill the mixed_events histogram 
    //for (int m = 0; m<2; m++){
    for(int m = 0; m <300; m++) {
      currentindex = m/20;
      getline(parsers[currentindex], mixed_num_string, '\t');
      mixed_events[m] = stoul(mixed_num_string);
      //fprintf(stderr,"%lu\n",mixed_events[m]);
    }
    newtree->Fill();
  }
  std::cout << "Successfully exited the eventloop" << std::endl;
  newtree->AutoSave();
  //newtree->Write();
  std::cout << "Successful autosave" <<std::endl;
  delete newfile;
  delete file;
  std::cout << "Deleted newfile" << std::endl;

  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
  }
