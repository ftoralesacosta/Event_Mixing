/**
   This program clones an NTuple, then uses data contained in text files to addmixed events to the clone
*/
// Author: Fernando & Ivan Chernyshev; Date: 6/18/2018

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
    if (argc < 5) {
      fprintf(stderr,"Syntax is [Command] [root file] [13d, 13e or 13f] [TrackSkim GeV] [Mix Start] [Mix End]\n");
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
	size_t lastindex = std::string(argv[fileArg]).find_last_of("."); 
	std::string rawname = std::string(argv[fileArg]).substr(0, lastindex);
	std::cout<<rawname<<std::endl;
	TFile *newfile = new TFile(Form("%s_%iGeVTrack_paired_test.root", rawname.data(), std::stoi((std::string)argv[trackpairenergyArg])), "RECREATE");
        TTree *newtree = _tree_event->CloneTree(0);
        
	// TFile *newfile = new TFile(Form("%s_mixedadded_output.root", ((std::string)argv[runArg]).c_str()), "RECREATE");
        // TTree *newtree = _tree_event->CloneTree(0);

        //new branch: mixed_events
        Long64_t mixed_events[NTRACK_MAX];
        newtree->Branch("mixed_events", mixed_events, "mixed_events[300]/L"); // One more entry needed for this to work
        
        std::cout<< "New branch successfully created " <<std::endl;
        
        // Get the mixed event textfiles
        std::ifstream mixed_textfile;

	std::ostringstream filename;
	filename << Form("%s_%iGeVTrack_Pairs_%lu_to_%lu.txt", rawname.data(), std::stoi((std::string)argv[trackpairenergyArg]), mix_start, mix_end);
	    mixed_textfile.open(filename.str());
	    std::cout<<"Opening Text File: "<<Form("%s_%iGeVTrack_Pairs_%lu_to_%lu.txt", rawname.data(), std::stoi((std::string)argv[trackpairenergyArg]), mix_start, mix_end)<<std::endl;

        
        
        const Long64_t nevents = _tree_event->GetEntries();
        // Loop over events
        for(Long64_t ievent = 0; ievent < nevents ; ievent++){
        //for(Long64_t ievent = 0; ievent < 2000; ievent++){
	  fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, _tree_event->GetEntries());
            _tree_event->GetEntry(ievent);
            // Get the appropriate line from each file, break out of the loop if you hit an empty file
            std::string eventline;
	    getline(mixed_textfile, eventline);
	    if (eventline == "") {
	      std::cout<<std::endl<<"reached end of file: "<<std::endl;
	      break;
	    }
	
            //try {
            std::string mixednum_string;
            long mixednum;
            std::istringstream parsers[1];
	    parsers[0].str(eventline);
            int currentindex;
	    Long64_t mix_range = (mix_end - mix_start);
            // Loop over mixed events, fill the mixed_events histogram while at it
            for(int m = 0; m < mix_range; m++) {
                currentindex = m/mix_range;
                getline(parsers[currentindex], mixednum_string, '\t');
                mixed_events[m] = stoul(mixednum_string);
		//fprintf(stderr,"%lu\n",mixed_events[m]);
            }
            //}
            //catch(std::invalid_argument) {
                //std::cout << "std_invalid_argument thrown at Event " << ievent << std::endl;
            //}
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
