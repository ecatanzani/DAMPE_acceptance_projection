///////////////////////////// Software description
//
//
//
//
/////////////////////////////////////////////////

#include "MyHead.h"

int main(int argc,char * argv[]) {
  
  ///////////// Variables for the TTrees

  Float_t g_lat = 0,g_lon = 0,geo_lat = 0,geo_lon = 0,sat_ra[3],sat_dec[3];
  UInt_t sec = 0;
  UShort_t n_ev = 0;
  bool good = true;
  Int_t chain_entries;                  // -> Variable that stores the whole number of tchain entries

  for(Int_t idx=0; idx<3; idx++) {
    sat_ra[idx] = 0;
    sat_dec[idx] = 0;
  }

    Float_t new_dec,new_ra;
    new_dec = 50;
    new_ra = 50;
    
  //////////////////////// Variable description

  /*
    
    sat_ra[3]     -> Is the array that stores the right ascension values of the satellite
    sat_dec[3]    -> Is the array that stores the celestial declination values of the satellite

    geo_lat       -> Is the geographic latitude
    geo_lon       -> Is the geographic longitude

    g_lat         -> Is the galactic latitude
    g_lon         -> Is the galactic longitude
    
                  These two are the "absolute coordinates" used to plot the final maps !!!!

   
    n_ev          -> Is the nuber of triggered events in a second
    sec           -> Is DAMPE's acquisition second number

    good          -> Is the status of the SBI

   */

  //////////////////////////////////////

  ///////////////////// Costheta flat binning variables
  
  Int_t n_bin_lon=360;                  // -> Number of bins along longitude axis

  Double_t lon_bin_min=-180.0;          // -> Set max and min for longitude binning
  Double_t lon_bin_max=180.0;

  Int_t n_bin_lat=180;                  // -> Number of bins along latitude axis
  
  Double_t lat_bin_min=-90.0;           // -> Set max and min for latitude binning
  Double_t lat_bin_max=90.0;

  Double_t* binning;                    // -> Array used to store the custom binning intervals !!!
  
  ///////////////////////////////////////////////////////////// 

  ///////////////////// Acceptance border POI vectors
  
  Bool_t first_call=true;
    
  vector<Double_t> peacks_costheta,valley_costheta;
  vector<Double_t> peacks_phi,valley_phi;

  TGraph *HEALPIX_peacks[2];
  TGraph *HEALPIX_valley[2];

    TGraph *s_HEALPIX_peacks[2];
    TGraph *s_HEALPIX_valley[2];
    
  //////////////////// Points vectors
    
  vector<Double_t> p_dNS,p_dEW,v_dNS,v_dEW;
    vector<Double_t> sp_dNS,sp_dEW,sv_dNS,sv_dEW;
    vector<Double_t> iterations;
  
   /////////////////// Slope variables
    
    Double_t slope_pNS,slope_pEW;
    Double_t slope_vNS,slope_vEW;
    
  /////////////////////////////////////////////////////////////
    
    
  Double_t perc=0;                         // -> Just used to store the percentage valu
  string log_path = output_path_creator(0),root_out_path = output_path_creator(1);
  
  ofstream output_log_file(log_path);     //log file creation !
  if(!output_log_file.is_open()) {
    cout<<"\n\nCannot create output file! Program finished !"<<endl;
    exit(-1);
  }

  log_file_init(output_log_file);

  TChain *tree= new TChain("SBItree");      //Defining a TTree to read SBI data fil
  read_SBI_data(g_lat,g_lon,geo_lat,geo_lon,sat_ra,sat_dec,sec,n_ev,good,tree,sbi_path,output_log_file);    //That function fills the tree reading from the SBI data files

  gRandom->SetSeed(random_seed);

  create_binning(n_bin_lat,lat_bin_min,lat_bin_max,binning,true);

  TFile *results_file = new TFile(root_out_path.c_str(),"RECREATE");

  /////////////////////////////// Opening DAMPE acceptance 2D-histo
  
  static TFile* acc_file= new TFile(acceptance_final_plot.Data());
  if(acc_file->IsZombie()) {
    cout << "\n\nError opening DAMPE acceptance TFile. Prorgram finished \n\n";
    output_log_file << "\n\nError opening DAMPE acceptance TFile. Prorgram finished \n\n";
    exit(-1);
  }

  static TH2D* acc = (TH2D*)acc_file->Get("Acceptance");
  static TH1D* acc_costheta = (TH1D*)acc_file->Get("Acceptance_X_Proj");
  static TH1D* acc_phi = (TH1D*)acc_file->Get("Acceptance_Y_Proj");
  
  /////////////////////////////////////////////////////////////////
    
  chain_entries = tree->GetEntries();
  results_file->cd();
    
  for(Int_t tree_idx=0; tree_idx<chain_entries; tree_idx++) {
    
    tree->GetEntry(tree_idx);

    if((sec%100000)==0)
      continue;                         //there was a bug in the SBI production
    if(!good)
      continue;                         //good second
    
    if (((Double_t)tree_idx/chain_entries)>(perc*0.01)) {
      cout<<"-> Tree events: [ "<<perc<<" % ]"<<endl;
      output_log_file<<"-> Tree events: [ "<<perc<<" % ]"<<endl;
      perc++;
    }

    if (g_lon>180)
      g_lon-=360;
    if (geo_lon>180)
      geo_lon-=360;

    for (int i=0; i<3; i++) {
      sat_ra[i]*=TMath::DegToRad();
      sat_dec[i]*=TMath::DegToRad();
    }

      /*
      if(wanna_just_project_acc) {
          if(MC_acc_proj==false)
              project_acceptance(results_file,tree_idx,sat_ra,sat_dec,acc,acc_costheta,acc_phi,output_log_file,n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
          else
              project_acceptance_edgeMC(results_file,tree_idx,sat_ra,sat_dec,acc,acc_costheta,acc_phi,output_log_file,n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
      }
      else {
          acceptance_projection_study(results_file,tree_idx,sat_ra,sat_dec,acc,acc_costheta,output_log_file,n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
          break; //This kind of study is done just for one second
      }
       */
    
      
      ////////////////////////////////// This function is used to plot the POI points each 100 seconds, to study shape deviations eventually present at different DAMPE's acquisition moments
      
      //shift_satellite_position(new_dec,new_ra,sat_ra,sat_dec,acc,n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);

      
    if(plot_POI_interval)
      if((tree_idx%100)==0)
	POI_discrete_plotting(sat_ra,sat_dec,acc,peacks_costheta,peacks_phi,valley_costheta,valley_phi,results_file,tree_idx,output_log_file,n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning,first_call,p_dNS,p_dEW,v_dNS,v_dEW,sp_dNS,sp_dEW,sv_dNS,sv_dEW,iterations,slope_pNS,slope_pEW,slope_vNS,slope_vEW);
 
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
  } //end loop on seconds

  ////////////////////////////////// Construst distance graphs

    if(plot_POI_interval) {
        
        HEALPIX_peacks[0] = new TGraph(p_dNS.size(),&iterations[0],&p_dNS[0]);
        HEALPIX_peacks[0]->SetTitle("NS peacks distance;#iterations;distance");
        HEALPIX_peacks[0]->SetName("peacks_distance_NS");
        HEALPIX_peacks[0]->SetMarkerStyle(31);
  
        HEALPIX_peacks[1] = new TGraph(p_dEW.size(),&iterations[0],&p_dEW[0]);
        HEALPIX_peacks[1]->SetTitle("EW peacks distance;#iterations;distance");
        HEALPIX_peacks[1]->SetName("peacks_distance_EW");
        HEALPIX_peacks[1]->SetMarkerStyle(31);
  
        HEALPIX_valley[0] = new TGraph(v_dNS.size(),&iterations[0],&v_dNS[0]);
        HEALPIX_valley[0]->SetTitle("NS valleys distance;#iterations;distance");
        HEALPIX_valley[0]->SetName("valleys_distance_NS");
        HEALPIX_valley[0]->SetMarkerStyle(31);
  
        HEALPIX_valley[1] = new TGraph(v_dEW.size(),&iterations[0],&v_dEW[0]);
        HEALPIX_valley[1]->SetTitle("EW valleys distance;#iterations;distance");
        HEALPIX_valley[1]->SetName("valleys_distance_EW");
        HEALPIX_valley[1]->SetMarkerStyle(31);
  
        ///////////////////////////////////////////////////////////////////////////

        HEALPIX_peacks[0]->Write();
        HEALPIX_peacks[1]->Write();
        HEALPIX_valley[0]->Write();
        HEALPIX_valley[1]->Write();
        
        ///////////////////////       Slope Calcolous       ///////////////////////
        
        s_HEALPIX_peacks[0] = new TGraph(sp_dNS.size(),&iterations[0],&sp_dNS[0]);
        s_HEALPIX_peacks[0]->SetTitle("NS peacks distance (slope);#iterations;distance");
        s_HEALPIX_peacks[0]->SetName("slope_peacks_distance_NS");
        s_HEALPIX_peacks[0]->SetMarkerStyle(31);
        
        s_HEALPIX_peacks[1] = new TGraph(sp_dEW.size(),&iterations[0],&sp_dEW[0]);
        s_HEALPIX_peacks[1]->SetTitle("EW peacks distance (slope);#iterations;distance");
        s_HEALPIX_peacks[1]->SetName("slope_peacks_distance_EW");
        s_HEALPIX_peacks[1]->SetMarkerStyle(31);
        
        s_HEALPIX_valley[0] = new TGraph(sv_dNS.size(),&iterations[0],&sv_dNS[0]);
        s_HEALPIX_valley[0]->SetTitle("NS valleys distance (slope);#iterations;distance");
        s_HEALPIX_valley[0]->SetName("slope_valleys_distance_NS");
        s_HEALPIX_valley[0]->SetMarkerStyle(31);
        
        s_HEALPIX_valley[1] = new TGraph(sv_dEW.size(),&iterations[0],&sv_dEW[0]);
        s_HEALPIX_valley[1]->SetTitle("EW valleys distance(slope);#iterations;distance");
        s_HEALPIX_valley[1]->SetName("slope_valleys_distance_EW");
        s_HEALPIX_valley[1]->SetMarkerStyle(31);
        
        ///////////////////////////////////////////////////////////////////////////
        
        s_HEALPIX_peacks[0]->Write();
        s_HEALPIX_peacks[1]->Write();
        s_HEALPIX_valley[0]->Write();
        s_HEALPIX_valley[1]->Write();
        
    }
        
  ///////////////////////////////////////// Manually shift satellite position ////////////////////////////
    
    
    
    if(shift_pointing)
        shift_satellite_position(new_dec,new_ra,sat_ra,sat_dec,acc,n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(save_output_file)
    results_file->Write();

  cout<<"\n\nSimulation completed !!\n\n";
  output_log_file<<"\n\nSimulation completed !!\n\n";
  
}
