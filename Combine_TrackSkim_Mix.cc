#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TLine.h>
#include <TLegend.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>

const int MAX_INPUT_LENGTH = 200;

int main(int argc, char *argv[])
//int main()
{
  if (argc < 2){
    std::cout<<"Syntax is [Command] [Same_Event_root_file] (Assumes Mixed File)"<<std::endl;
    exit(EXIT_FAILURE);
  }
  
  int nTrackSkims = 2;
  float* trackPtSkims;
  trackPtSkims = new float[nTrackSkims+1];

  trackPtSkims[0] = 0.0; //trackPtSkims[1] = 40.0;
  trackPtSkims[1] = 4.0; //trackPtSkims[2] = 6.0;
  trackPtSkims[2] = 40.0;

    std::string basic_name = argv[1];
    std::string::size_type pos = basic_name.find('_');
    if (pos != std::string::npos)
      basic_name = basic_name.substr(0, pos);
    fprintf(stderr,"%d: basic_name = %s \n",__LINE__,basic_name.c_str());

  TFile* MixFile[nTrackSkims];
  for (int iSkim = 0; iSkim < nTrackSkims; iSkim ++){
    // MixFile[iSkim] = TFile::Open(Form("%s_%1.0fGeV_Skim_Correlation.root",basic_name.c_str(),trackPtSkims[iSkim]));
    //MixFile[iSkim] = TFile::Open(Form("%s_%luGeVTracks_Mixed_%lu_Correlation_%s.root"),basic_name.c_str(),trackPtSkims[iSkim],300,shower);
    //MixFile[iSkim] = TFile::Open(Form("%s_%1.0fGeVTrack_paired_%1.0fGeVTracks_Correlation_Lambda_0_to_300.root",basic_name.c_str(),trackPtSkims[iSkim],trackPtSkims[iSkim]));
    MixFile[iSkim] = TFile::Open(Form("%s_MB_%1.0fGeV_Skim_Correlation_L0.root",basic_name.c_str(),trackPtSkims[iSkim]));
    //MixFile[iSkim] = TFile::Open(Form("%s_%1.0fGeVTrack_paired_%1.0fGeVTracks_Correlation_Lambda_0_to_300.root",basic_name.c_str(),trackPtSkims[iSkim],trackPtSkims[iSkim]));
    //std::cout<<Form("%s_%1.0fGeVTracks.root",basic_name.c_str(),trackPtSkims[iSkim])<<std::endl;

    //MixFile[iSkim] = TFile::Open((TString)argv[2]);

    if (MixFile[iSkim] == NULL) {
      std::cout<<MixFile[iSkim]<<std::endl;
      std::cout << " Mixed file "<<iSkim<<" failed" << std::endl;
      exit(EXIT_FAILURE);}
  }


  FILE* config = fopen("../Corr_config.yaml", "r");

  int n_eta_bins = 0;
  int n_phi_bins = 0;

  //FIXME add Corr_config.yaml support
  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins+1];
  ztbins[0] = 0.05; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
  //FIXME: Why was it written like this? why not zt[] = {a,b,c,...}

  int nptbins = 3;
  float* ptbins;
  ptbins = new float[nptbins+1];
  ptbins[0] = 10.0; ptbins[1] = 11; ptbins[2] = 12.5; ptbins[3] = 16;


  // Loop through config file
  char line[MAX_INPUT_LENGTH];
  while (fgets(line, MAX_INPUT_LENGTH, config) != NULL) {
    if (line[0] == '#') continue;

    char key[MAX_INPUT_LENGTH];
    char dummy[MAX_INPUT_LENGTH];
    char value[MAX_INPUT_LENGTH];

    // Cap off key[0] and value[0] with null characters and load the key, dummy-characters, and value of the line into their respective arrays       
    key[0] = '\0';
    value[0] = '\0';
    sscanf(line, "%[^:]:%[ \t]%100[^\n]", key, dummy, value);
    
    if (strcmp(key, "N_Phi_Bins") == 0) {
      n_phi_bins = atoi(value);
      std::cout << "Number of Phi Bins: " << n_phi_bins << std::endl; }

    else if (strcmp(key, "N_Eta_Bins") == 0) {
      n_eta_bins = atoi(value);
      std::cout << "Number of Eta Bins: " << n_eta_bins << std::endl; }

    else if (strcmp(key, "Zt_bins") == 0) {
      nztbins = -1;
      for (const char *v = value; *v != ']';) {
	while (*v != ']' && !isdigit(*v)) v++;
	nztbins++;
	while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
      
      ztbins = new float[nztbins + 1];
      int i = 0;
      for (const char *v = value; *v != ']' ;) {
	while (*v != ']' && !isdigit(*v)) v++;
	ztbins[i] = atof(v);
	i++;
	while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
      
      std::cout << "Number of Zt bins: " << nztbins << std::endl << "Zt bins: {";
      for (int i = 0; i <= nztbins; i++)
	std::cout << ztbins[i] << ", ";
      std::cout << "}\n";
    }
    
    else if (strcmp(key, "Pt_bins") == 0) {
      nptbins = -1;
      for (const char *v = value; *v != ']';) {
	while (*v != ']' && !isdigit(*v)) v++;
	nptbins++;
	while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
      
      ptbins = new float[nptbins + 1];
      int i = 0;
      for (const char *v = value; *v != ']' ;) {
	while (*v != ']' && !isdigit(*v))  v++;
	ptbins[i] = atof(v);
	i++;
	while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
      
      std::cout << "Number of Pt bins: " << nptbins << std::endl << "Pt bins: {";
      for (int i = 0; i <= nptbins; i++)
	std::cout << ptbins[i] << ", ";
      std::cout << "}\n";
    }
  }

  fclose(config);

  for (int i = 0; i <= nztbins; i++)
    std::cout << "zt bound: " << ztbins[i] << std::endl;
  for (int i = 0; i <= nptbins; i++)
    std::cout << "pt bound: " << ptbins[i] << std::endl;
  std::cout<<"Number of eta bins"<<n_eta_bins<<std::endl;
  std::cout<<"Number of phi bins"<<n_phi_bins<<std::endl;
  //Config File Done -------------//

  TH2F* Mix_Inclusive_Corr[nztbins*nptbins];
  TH2F* Mix_Isolated_Corr[nztbins*nptbins];
  TH2F* Mix_DNN1_Corr[nztbins*nptbins];
  TH2F* Mix_DNN2_Corr[nztbins*nptbins];

  TH1D* N_Incl_Triggers[nptbins];
  TH1D* N_Isolated_Triggers[nptbins];
  TH1D* N_Sig_Triggers[nptbins];
  TH1D* N_BKGD_Triggers[nptbins];
  TH1D* N_Signal_Overlap_Triggers[nptbins];
  TH1D* N_BKGD_Overlap_Triggers[nptbins];

  //TFile *MyFile = new TFile(Form("%s_GMB_Ratio_COMBINED_OLD.root",basic_name.c_str()),"RECREATE");
  TFile *MyFile = new TFile(Form("%s_GMB_RATIO_COMBINED_OLD.root",basic_name.c_str()),"RECREATE");
 
  Long64_t DNN1_entries;
  for(int ipt = 0; ipt < nptbins; ipt++){
    for (int izt = 0; izt<nztbins; izt++){

      //Implement due to Track pT Skimming min bias
      for (int iSkim = 0; iSkim < nTrackSkims; iSkim ++){
	fprintf(stderr,"%s%d: iSkim = %i \n",__FILE__,__LINE__,iSkim);


	//SELECT WHICH MB CORRELATION TO USE

	if ((ptbins[ipt]*ztbins[izt] >= trackPtSkims[iSkim]) && ((ptbins[ipt]*ztbins[izt]) < trackPtSkims[iSkim+1]+0.5)){
	  //a bit hacky, since ptbin*ztbin do not always correspond to integers...

	  fprintf(stderr,"pt min max: %1.2f %1.2f, zt min max %1.2f %1.2f\n",ptbins[ipt],ptbins[ipt+1],ztbins[izt],ztbins[izt+1]);
	  fprintf(stderr,"Track pT between %1.2f and %1.2f\n",ptbins[ipt]*ztbins[izt],ptbins[ipt+1]*ztbins[izt+1]);
	  

	  // Mix_Inclusive_Corr[izt+ipt*nztbins] = (TH2F*)MixFile[iSkim]->Get(
	  //     Form("Inclusive_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
	  //     ptbins[ipt],ptbins[ipt+1],100*ztbins[izt],100*ztbins[izt+1]));

	  // if (Mix_Inclusive_Corr[izt+ipt*nztbins] == NULL) {
	  //   std::cout << "Mix Incl TH2D fail" << std::endl;
	  //   exit(EXIT_FAILURE);}


	  Mix_Isolated_Corr[izt+ipt*nztbins] = (TH2F*)MixFile[iSkim]->Get(
	      Form("Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
	      ptbins[ipt],ptbins[ipt+1],100*ztbins[izt],100*ztbins[izt+1]));

	  fprintf(stderr,"Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
		  ptbins[ipt],ptbins[ipt+1],100*ztbins[izt],100*ztbins[izt+1]);

	  DNN1_entries = Mix_Isolated_Corr[izt+ipt*nztbins]->GetEntries();
	  fprintf(stderr, "%s:%d: Number of Entries in Mix Isolated %lld \n",__FILE__,__LINE__,DNN1_entries);
      
	  if (Mix_Isolated_Corr[izt+ipt*nztbins] == NULL) {
	    std::cout << "Mix Iso TH2D fail" << std::endl;
	    exit(EXIT_FAILURE);}

	  Mix_DNN1_Corr[izt+ipt*nztbins] = (TH2F*)MixFile[iSkim]->Get(
          Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
          1,ptbins[ipt],ptbins[ipt+1],100*ztbins[izt],100*ztbins[izt+1]));   

	  if (Mix_DNN1_Corr[izt+ipt*nztbins] == NULL) {
	    std::cout << " mix TH2D DNN1 read fail "<<"Mix file "<<iSkim<< std::endl;
	    exit(EXIT_FAILURE);}

	  DNN1_entries = Mix_DNN1_Corr[izt+ipt*nztbins]->GetEntries();
	  fprintf(stderr, "%s:%d: Number of Entries in Mix DNN1 %lu \n",__FILE__,__LINE__,DNN1_entries);

	  Mix_DNN2_Corr[izt+ipt*nztbins] = (TH2F*)MixFile[iSkim]->Get(
          Form("DNN%i_Correlation__pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",
	  2,ptbins[ipt],ptbins[ipt+1],100*ztbins[izt],100*ztbins[izt+1]));

	  if (Mix_DNN2_Corr[izt+ipt*nztbins] == NULL) {
	    std::cout << " mix TH2D DNN2 read fail "<<"Mix file "<<iSkim<<std::endl;
	    exit(EXIT_FAILURE);}     
	}

      } //iSkim

      Double_t mix_DNN1_intgrl = Mix_DNN1_Corr[izt+ipt*nztbins]->GetBinContent(Mix_DNN1_Corr[izt+ipt*nztbins]->GetMaximumBin());
      Double_t mix_DNN2_intgrl = Mix_DNN2_Corr[izt+ipt*nztbins]->GetBinContent(Mix_DNN2_Corr[izt+ipt*nztbins]->GetMaximumBin());
      Double_t mix_Isolated_intgrl = Mix_Isolated_Corr[izt+ipt*nztbins]->GetBinContent(Mix_Isolated_Corr[izt+ipt*nztbins]->GetMaximumBin());
      //    Double_t mix_Inc_intgrl = Mix_Inclusive_Corr[izt+ipt*nztbins]->GetBinContent(Mix_Inclusive_Corr[izt+ipt*nztbins]->GetMaximumBin());                                                                                               

      //Mix_Inclusive_Corr[izt+ipt*nztbins]->Scale(1.0/mix_Inc_intgrl);                                                                                                                                                                       
      Mix_Isolated_Corr[izt+ipt*nztbins]->Scale(1.0/mix_Isolated_intgrl);
      Mix_DNN1_Corr[izt+ipt*nztbins]->Scale(1.0/mix_DNN1_intgrl);
      Mix_DNN2_Corr[izt+ipt*nztbins]->Scale(1.0/mix_DNN2_intgrl);

    fprintf(stderr, "%s: %d: Division OK\n",__FILE__,__LINE__);

    Mix_DNN1_Corr[izt+ipt*nztbins]->Write();
    fprintf(stderr, "%s: %d: Write DNN1 OK\n",__FILE__,__LINE__);
    }//zt bin
    
    for (int izt = 0; izt < nztbins; izt++){
      Mix_DNN2_Corr[izt+ipt*nztbins]->Write();
    fprintf(stderr, "%s: %d: Write DNN2 OK\n",__FILE__,__LINE__);
    }

    for (int izt = 0; izt < nztbins; izt++){
      Mix_Isolated_Corr[izt+ipt*nztbins]->Write();
      fprintf(stderr, "%s: %d: Write Isolated OK\n",__FILE__,__LINE__);
    }

    // for (int izt = 0; izt < nztbins; izt++){
    //   Same_Inclusive_Corr[izt+ipt*nztbins]->Write();
    //   fprintf(stderr, "%s: %d: Write Inclusive OK\n",__FILE__,__LINE__);
    // }

    //    N_Incl_Triggers[ipt] = (TH1D*)corr->Get(Form("N_Inclusive_Triggers_pT%1.0f_%1.0f",1,ptbins[ipt],ptbins[ipt+1]));
    fprintf(stderr, "%s: %d: Write Trigger Get OK",__FILE__,__LINE__);

    //N_Incl_Triggers[ipt]->Write();
  }//ipt

  for (int iSkim = 0; iSkim < nTrackSkims; iSkim ++) MixFile[iSkim]->Close();
  MyFile->Close();
  //  theApp.Run();
  std::cout<<"Success"<<std::endl;

  return 0;
}
    
