
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "TRandom3.h"
#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TTree.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TNamed.h"
#include "TFile.h"
#include "TChain.h"
#include "TColor.h"
#include "TLine.h"

using namespace std;

// !!!! Don't insert any value for "first_iter" and "last_iter" to use the whole result file !!!!

int my_crew(const TString results_alias, const Int_t first_iter=100, const Int_t last_iter=100) {
  if(first_iter!=last_iter) {
    if(last_iter<first_iter) {
      cout<<"\n\n ERROR: 'last_iretation_value' should be bigger that 'first_iteration_value' \n\n";
      exit(-1);
    }
    else if(first_iter<100) {
      cout<<"\n\n ERROR: file include 100 seconds schedule !! Insert a bigger value for 'first_iretation_value' \n\n";
      exit(-1);
    }
    else if(last_iter<100) {
      cout<<"\n\n ERROR: file include 100 seconds schedule !! Insert a bigger value for 'last_iteration_value' \n\n";
      exit(-1);
    }
    else if((last_iter-first_iter)<100) {
      cout<<"\n\n ERROR: file include 100 seconds schedule !! Insert a bigger value for 'last_iteration_value' \n\n";
      exit(-1);
    }
  }
      
  Bool_t found_item=true;
  vector<TH2D *> h_crew;
  Int_t counter=0,tmp_fI;
  TString input_results_path = "../results/";
  input_results_path+=results_alias;
  input_results_path+="_maps_result.root";
  TFile *input_results_file = new TFile(input_results_path.Data());
  if(input_results_file->IsZombie()) {
    cout<<"\n\nError opening input file. Macro finished\n\n";
    exit(-1);
    }
  TCanvas *c1 = new TCanvas("c1","Crew Canvas - POI Time Development");
  while(found_item) {
    TString wanna_file="POI_Healpix_";
    if(first_iter==100 && last_iter==100) {
      tmp_fI=first_iter+counter*100;
      wanna_file+=tmp_fI;
      if(input_results_file->Get(wanna_file)) {
	h_crew.push_back((TH2D*)(input_results_file->Get(wanna_file)));
	h_crew.at(counter)->SetMarkerColor(counter);
	h_crew.at(counter)->Draw("same");
	counter++;
      }
      else
	found_item=false;
    }
    else {
	tmp_fI=first_iter+counter*100;
	if(tmp_fI>last_iter)
          break;
	else {
	  wanna_file+=tmp_fI;
	  if(input_results_file->Get(wanna_file)) {
	    h_crew.push_back((TH2D*)(input_results_file->Get(wanna_file)));
	    h_crew.at(counter)->SetMarkerColor(counter+1);
	    h_crew.at(counter)->Draw("same");
	    counter++;
	  }
	  else
	    found_item=false;
	}
    }
  }
  /*
    TH2 *h_crew[counter];
    for(Int_t idx_h=0; idx_h<counter; idx_h++) {
    TString wanna_file="POI_Healpix_";
    tmp_fI=100+idx_h*100;
    wanna_file+=tmp_fI;
    TH2* h_crew[idx_h] = (TH2D*)(input_results_file->Get(wanna_file));
    h_crew[idx_h]->Draw("same");
    }
  */
  return 0;
}
