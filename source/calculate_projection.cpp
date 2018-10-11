
#include "MyHead.h"

void project_acceptance(TFile *results,Int_t h_idx,Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,TH1D *acc_costheta,TH1D* acc_phi,ofstream &out_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning) {

  Double_t b=0,l=0;
  Double_t costheta=-999, phi=-999;
  TString h_name;
  TLine *h_line[4];

  Double_t acc_glat_min=-999,acc_glat_max=-999,acc_glon_min=-999,acc_glon_max=-999;
  
  results->cd();
  
  h_name = "Acceptance_Sky_Projection_";
  h_name+=h_idx;
  TH2D *h_projection = new TH2D(h_name,"Acceptance Projection; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
  
  from_local_to_galactic(0.62,2.6,acc_glon_min,acc_glat_min,sat_ra,sat_dec);
  from_local_to_galactic(0.68,6.,acc_glon_max,acc_glat_max,sat_ra,sat_dec);
  
  if(acc_glon_min>180)
    acc_glon_min-=360.;
  if(acc_glon_max>180)
    acc_glon_max-=360.;
  
  for(Int_t idx_ev=0; idx_ev<sky_events; idx_ev++) {
    acc->GetRandom2(costheta,phi);
    from_local_to_galactic(costheta,phi,l,b,sat_ra,sat_dec);
    
    if(isnan(l)) {
      idx_ev--;
      continue;
    }

    if (l>180.0)
      l-=360.0;

    
    if( (l>acc_glon_max) || (l<acc_glon_min) || (b>acc_glat_max) || (b<acc_glat_min) ) {
      cout<<"\nLongitude or latitude bandwith exceeded ! -> Longitude:"<<l<<" (min = "<<acc_glon_min<<" max = "<<acc_glon_max<<" ) Latitude: "<<b<<" (min = "<<acc_glat_min<<" max = "<<acc_glat_max<<")";
      //out_file<<"\nLongitude or latitude bandwith exceeded ! -> Longitude:"<<l<<" (min = "<<acc_glon_min<<" max = "<<acc_glon_max<<" ) Latitude: "<<b<<" (min = "<<acc_glat_min<<" max = "<<acc_glat_max<<")";
    }
      
    

    h_projection->Fill(l,b,1);     //Is important to note that all the events has the same weigh. Here we just want to compute the acceptance projection on the sky. Anisotropy study is not performed here
  } //end loop on sky events

  h_line[0] = new TLine(acc_glon_min,acc_glat_min,acc_glon_max,acc_glat_min);
  h_line[1] = new TLine(acc_glon_max,acc_glat_min,acc_glon_max,acc_glat_max);
  h_line[2] = new TLine(acc_glon_max,acc_glat_max,acc_glon_min,acc_glat_max);
  h_line[3] = new TLine(acc_glon_min,acc_glat_max,acc_glon_min,acc_glat_min);

  for(Int_t idx_l=0; idx_l<4; idx_l++) {
    h_line[idx_l]->SetLineColor(kRed);
    h_projection->GetListOfFunctions()->Add(h_line[idx_l]);
  }

  //for(Int_t idx_l=0; idx_l<4; idx_l++)
  //h_line[idx_l]->Delete();

  h_projection->Write();
  h_projection->Delete();
  
}


extern void project_acceptance_edgeMC(TFile *results,Int_t h_idx,Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,TH1D *acc_costheta,TH1D* acc_phi,ofstream &out_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning) {

  Double_t b=0,l=0;
  Double_t costheta=-999, phi=-999;
  TString h_name;
  TLine *h_line[4];

  Double_t acc_glat_min=-999,acc_glat_max=-999,acc_glon_min=-999,acc_glon_max=-999;
  
  results->cd();

  h_name = "Acceptance_Sky_Projection_";
  h_name+=h_idx;
  TH2D *h_projection = new TH2D(h_name,"Acceptance Projection; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);

  obtain_edge(sat_ra,sat_dec,acc,acc_glat_min,acc_glat_max,acc_glon_min,acc_glon_max);
  
  for(Int_t idx_ev=0; idx_ev<sky_events; idx_ev++) {
    acc->GetRandom2(costheta,phi);
    from_local_to_galactic(costheta,phi,l,b,sat_ra,sat_dec);
    
    if(isnan(l)) {
      idx_ev--;
      continue;
    }

    if (l>180.0)
      l-=360.0;
    
    if( (l>acc_glon_max) || (l<acc_glon_min) || (b>acc_glat_max) || (b<acc_glat_min) ) {
      cout<<"\nLongitude or latitude bandwith exceeded ! -> Longitude:"<<l<<" (min = "<<acc_glon_min<<" max = "<<acc_glon_max<<" ) Latitude: "<<b<<" (min = "<<acc_glat_min<<" max = "<<acc_glat_max<<")";
      out_file<<"\nLongitude or latitude bandwith exceeded ! -> Longitude:"<<l<<" (min = "<<acc_glon_min<<" max = "<<acc_glon_max<<" ) Latitude: "<<b<<" (min = "<<acc_glat_min<<" max = "<<acc_glat_max<<")";
    }
    
    h_projection->Fill(l,b,1);     //Is important to note that all the events has the same weigh. Here we just want to compute the acceptance projection on the sky. Anisotropy study is not performed here
  } //end loop on sky events

  h_line[0] = new TLine(acc_glon_min,acc_glat_min,acc_glon_max,acc_glat_min);
  h_line[1] = new TLine(acc_glon_max,acc_glat_min,acc_glon_max,acc_glat_max);
  h_line[2] = new TLine(acc_glon_max,acc_glat_max,acc_glon_min,acc_glat_max);
  h_line[3] = new TLine(acc_glon_min,acc_glat_max,acc_glon_min,acc_glat_min);

  for(Int_t idx_l=0; idx_l<4; idx_l++) {
    h_line[idx_l]->SetLineColor(kRed);
    h_projection->GetListOfFunctions()->Add(h_line[idx_l]);
  }

  h_projection->Write();
  h_projection->Delete();

}

void obtain_edge(Float_t sat_ra[],Float_t sat_dec[],TH2D *h_acc,Double_t &acc_glat_min,Double_t &acc_glat_max,Double_t &acc_glon_min,Double_t &acc_glon_max) {

  Double_t b=0,l=0;
  Double_t costheta=-999, phi=-999;
  
  for(Int_t idx_ev=0; idx_ev<n_gevents; idx_ev++) {
    h_acc->GetRandom2(costheta,phi);
    from_local_to_galactic(costheta,phi,l,b,sat_ra,sat_dec);

    if(isnan(l)) {
      idx_ev--;
      continue;
    }

    if(l>180.0)
      l-=360.0;

    if(idx_ev==0) {
      acc_glat_max=acc_glat_min=b;
      acc_glon_max=acc_glon_min=l;
    }
    else {
      if(l>acc_glon_max)
	acc_glon_max=l;
      if(l<acc_glon_min)
	acc_glon_min=l;
      if(b>acc_glat_max)
	acc_glat_max=b;
      if(b<acc_glat_min)
	acc_glat_min=b;
    }
  } //end for on generated events  
}

void acceptance_projection_study(TFile *results_file,Int_t tree_idx,Float_t sat_ra[],Float_t sat_dec[],TH2D *acc,TH1D *acc_costheta,ofstream &output_log_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning) {
    
    Double_t b=0,l=0;
    Double_t costheta=-999, phi=-999;
    Double_t perc=0;
    vector<Double_t> peacks_costheta,valley_costheta;
    vector<Double_t> peacks_phi,valley_phi;
    
    UInt_t seed=22;
    TRandom3 *rnd_gen = new TRandom3(seed);
    
    results_file->cd();
    
    /////// Extract acceptance border
    
    static TH2D *acc_border = (TH2D*)acc->Clone("acc_border");
    acc_border->Reset();
    get_acceptance_border(acc,acc_border);
    acc_border->Write();
    
    /////////////////////////////////
    
    ////// Mark rilevant acceptance points
    
    TString p_acc,p_healpix,p_lb_uniform,p_lcosb_uniform,p_lcos2b_uniform;
    TString v_acc,v_healpix,v_lb_uniform,v_lcosb_uniform,v_lcos2b_uniform;
    
    p_acc="Peacks_Acceptance_Sky_Projection_";
    p_healpix="Healpix_Peacks_Acceptance_Sky_Projection_";
    p_lb_uniform="lb_uniform_Peacks_Acceptance_Sky_Projection_";
    p_lcosb_uniform="lcosb_uniform_Peacks_Acceptance_Sky_Projection_";
    p_lcos2b_uniform="lcos2b_uniform_Peacks_Acceptance_Sky_Projection_";
    
    p_acc+=tree_idx;
    p_healpix+=tree_idx;
    p_lb_uniform+=tree_idx;
    p_lcosb_uniform+=tree_idx;
    p_lcos2b_uniform+=tree_idx;
    
    v_acc="Valley_Acceptance_Sky_Projection_";
    v_healpix="Healpix_Valley_Acceptance_Sky_Projection_";
    v_lb_uniform="lb_uniform_Valley_Acceptance_Sky_Projection_";
    v_lcosb_uniform="lcosb_uniform_Valley_Acceptance_Sky_Projection_";
    v_lcos2b_uniform="lcos2b_uniform_Valley_Acceptance_Sky_Projection_";
    
    v_acc+=tree_idx;
    v_healpix+=tree_idx;
    v_lb_uniform+=tree_idx;
    v_lcosb_uniform+=tree_idx;
    v_lcos2b_uniform+=tree_idx;
    
    get_relevant_accceptance_points(acc_border,peacks_costheta,valley_costheta,peacks_phi,valley_phi);
    
    TH2D *h_p_acc = new TH2D(p_acc,"Acceptance Peacks Points; cos(#theta); #phi; Entries",1000,0,1,1000,0,2*TMath::Pi());
    TH2D *h_p_healpix = new TH2D(p_healpix,"Acceptance Peacks Sky Projection; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
    TH2D *h_p_lb_uniform = new TH2D(p_lb_uniform,"Acceptance Peacks Sky Projection; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",1000,-180,180,1000,-90,90);
    TH2D *h_p_lcosb_uniform = new TH2D(p_lcosb_uniform,"Acceptance Peacks Sky Projection; l (#circ);  cos(b); Entries",1000,-180,180,1000,-1,1);
    TH2D *h_p_lcos2b_uniform = new TH2D(p_lcos2b_uniform,"Acceptance Peacks Sky Projection; l (#circ);  cos^{2}(b); Entries",1000,-180,180,1000,0,1);
    
    TH2D *h_v_acc = new TH2D(v_acc,"Acceptance Valley Points; cos(#theta); #phi; Entries",1000,0,1,1000,0,2*TMath::Pi());
    TH2D *h_v_healpix = new TH2D(v_healpix,"Acceptance Valley Sky Projection; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
    TH2D *h_v_lb_uniform = new TH2D(v_lb_uniform,"Acceptance Valley Sky Projection; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",1000,-180,180,1000,-90,90);
    TH2D *h_v_lcosb_uniform = new TH2D(v_lcosb_uniform,"Acceptance Valley Sky Projection; l (#circ);  cos(b); Entries",1000,-180,180,1000,-1,1);
    TH2D *h_v_lcos2b_uniform = new TH2D(v_lcos2b_uniform,"Acceptance Valley Sky Projection; l (#circ);  cos^{2}(b); Entries",1000,-180,180,1000,0,1);
    
    h_p_acc->SetMarkerColor(2);
    h_p_acc->SetMarkerStyle(41);
    
    h_p_healpix->SetMarkerColor(2);
    h_p_healpix->SetMarkerStyle(41);
    
    h_p_lb_uniform->SetMarkerColor(2);
    h_p_lb_uniform->SetMarkerStyle(41);
    
    h_p_lcosb_uniform->SetMarkerColor(2);
    h_p_lcosb_uniform->SetMarkerStyle(41);
    
    h_p_lcos2b_uniform->SetMarkerColor(2);
    h_p_lcos2b_uniform->SetMarkerStyle(41);
    
    
    for(Int_t idx_v=0; idx_v<peacks_costheta.size(); idx_v++) {
        from_local_to_galactic(peacks_costheta.at(idx_v),peacks_phi.at(idx_v),l,b,sat_ra,sat_dec);
        if(l>180.0)
            l-=360.0;
        
        h_p_acc->Fill(peacks_costheta.at(idx_v),peacks_phi.at(idx_v),1);
        h_p_healpix->Fill(l,b,1);
        h_p_lb_uniform->Fill(l,b,1);
        h_p_lcosb_uniform->Fill(l,cos(b*TMath::DegToRad()),1);
        h_p_lcos2b_uniform->Fill(l,TMath::Power(cos(b*TMath::DegToRad()),2),1);
    }
    
    h_v_acc->SetMarkerColor(4);
    h_v_acc->SetMarkerStyle(41);
    
    h_v_healpix->SetMarkerColor(4);
    h_v_healpix->SetMarkerStyle(41);
    
    h_v_lb_uniform->SetMarkerColor(4);
    h_v_lb_uniform->SetMarkerStyle(41);
    
    h_v_lcosb_uniform->SetMarkerColor(4);
    h_v_lcosb_uniform->SetMarkerStyle(41);
    
    h_v_lcos2b_uniform->SetMarkerColor(4);
    h_v_lcos2b_uniform->SetMarkerStyle(41);
    
    
    for(Int_t idx_v=0; idx_v<valley_costheta.size(); idx_v++) {
        from_local_to_galactic(valley_costheta.at(idx_v),valley_phi.at(idx_v),l,b,sat_ra,sat_dec);
        if(l>180.0)
            l-=360.0;
        
        h_v_acc->Fill(valley_costheta.at(idx_v),valley_phi.at(idx_v),1);
        h_v_healpix->Fill(l,b,1);
        h_v_lb_uniform->Fill(l,b,1);
        h_v_lcosb_uniform->Fill(l,cos(b*TMath::DegToRad()),1);
        h_v_lcos2b_uniform->Fill(l,TMath::Power(cos(b*TMath::DegToRad()),2),1);
    }
    
    h_p_acc->Write();
    h_p_healpix->Write();
    h_p_lb_uniform->Write();
    h_p_lcosb_uniform->Write();
    h_p_lcos2b_uniform->Write();
    
    h_v_acc->Write();
    h_v_healpix->Write();
    h_v_lb_uniform->Write();
    h_v_lcosb_uniform->Write();
    h_v_lcos2b_uniform->Write();
    
    h_p_acc->Delete();
    h_p_healpix->Delete();
    h_p_lb_uniform->Delete();
    h_p_lcosb_uniform->Delete();
    h_p_lcos2b_uniform->Delete();
    
    h_v_acc->Delete();
    h_v_healpix->Delete();
    h_v_lb_uniform->Delete();
    h_v_lcosb_uniform->Delete();
    h_v_lcos2b_uniform->Delete();
    
    
    ////////////////////////////////
    
    TString healpix_acc,lb_uniform_acc,lcosb_uniform_acc,lcos2b_uniform_acc;
    TString healpix_Eacc,lb_uniform_Eacc,lcosb_uniform_Eacc,lcos2b_uniform_Eacc;
    TString healpix_SQacc,lb_uniform_SQacc,lcosb_uniform_SQacc,lcos2b_uniform_SQacc;
    TString healpix_ESQacc,lb_uniform_ESQacc,lcosb_uniform_ESQacc,lcos2b_uniform_ESQacc;
    Double_t costheta_min=obtain_min_histo(acc_costheta);
    
    
    ///////////////////////////////////////////// Considering real DAMPE acceptance in costheta and phi and evaluate the acceptance projection and its border
    
    
    healpix_acc = "Healpix_Acceptance_Sky_Projection_";
    lb_uniform_acc = "lb_uniform_Acceptance_Sky_Projection_";
    lcosb_uniform_acc = "lcosb_uniform_Acceptance_Sky_Projection_";
    lcos2b_uniform_acc = "lcos2b_uniform_Acceptance_Sky_Projection_";
    
    healpix_acc+=tree_idx;
    lb_uniform_acc+=tree_idx;
    lcosb_uniform_acc+=tree_idx;
    lcos2b_uniform_acc+=tree_idx;
    
    TH2D *h_healpix_acc = new TH2D(healpix_acc,"Acceptance Sky Projection; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
    TH2D *h_lb_uniform_acc = new TH2D(lb_uniform_acc,"Acceptance Sky Projection; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",1000,-180,180,1000,-90,90);
    TH2D *h_lcosb_uniform_acc = new TH2D(lcosb_uniform_acc,"Acceptance Sky Projection; l (#circ);  cos(b); Entries",1000,-180,180,1000,-1,1);
    TH2D *h_lcos2b_uniform_acc = new TH2D(lcos2b_uniform_acc,"Acceptance Sky Projection; l (#circ);  cos^{2}(b); Entries",1000,-180,180,1000,0,1);
    
    for(Int_t idx_ev=0; idx_ev<n_points; idx_ev++) {
        
        if (((Double_t)idx_ev/n_points)>(perc*0.01)) {
            cout<<"\t-> AccStudy events: [ "<<perc<<" % ]"<<endl;
            output_log_file<<"\t-> AccStudy events: [ "<<perc<<" % ]"<<endl;
            perc++;
        }
        
        acc->GetRandom2(costheta,phi);
        from_local_to_galactic(costheta,phi,l,b,sat_ra,sat_dec);
        
        if(isnan(l)) {
            idx_ev--;
            continue;
        }
        
        if (l>180.0)
            l-=360.0;
     
        h_healpix_acc->Fill(l,b,1);
        h_lb_uniform_acc->Fill(l,b,1);
        h_lcosb_uniform_acc->Fill(l,cos(b*TMath::DegToRad()),1);
        h_lcos2b_uniform_acc->Fill(l,TMath::Power(cos(b*TMath::DegToRad()),2),1);
        
    }
    
    h_healpix_acc->Write();
    h_lb_uniform_acc->Write();
    h_lcosb_uniform_acc->Write();
    h_lcos2b_uniform_acc->Write();
    
    h_healpix_acc->Delete();
    h_lb_uniform_acc->Delete();
    h_lcosb_uniform_acc->Delete();
    h_lcos2b_uniform_acc->Delete();
    
    
    ////////////////// And now just the border
    
    healpix_Eacc = "Healpix_EDGE_Acceptance_Sky_Projection_";
    lb_uniform_Eacc = "lb_uniform_EDGE_Acceptance_Sky_Projection_";
    lcosb_uniform_Eacc = "lcosb_uniform_EDGE_Acceptance_Sky_Projection_";
    lcos2b_uniform_Eacc = "lcos2b_uniform_EDGE_Acceptance_Sky_Projection_";
    
    healpix_Eacc+=tree_idx;
    lb_uniform_Eacc+=tree_idx;
    lcosb_uniform_Eacc+=tree_idx;
    lcos2b_uniform_Eacc+=tree_idx;
    
    TH2D *h_healpix_Eacc = new TH2D(healpix_Eacc,"Edge Acceptance Sky Projection; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
    TH2D *h_lb_uniform_Eacc = new TH2D(lb_uniform_Eacc,"Edge Acceptance Sky Projection; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",1000,-180,180,1000,-90,90);
    TH2D *h_lcosb_uniform_Eacc = new TH2D(lcosb_uniform_Eacc,"Edge Acceptance Sky Projection; l (#circ);  cos(b); Entries",1000,-180,180,1000,-1,1);
    TH2D *h_lcos2b_uniform_Eacc = new TH2D(lcos2b_uniform_Eacc,"Edge Acceptance Sky Projection; l (#circ);  cos^{2}(b); Entries",1000,-180,180,1000,0,1);
    
    cout<<endl;
    perc=0;
    
    for(Int_t idx_ev=0; idx_ev<n_points; idx_ev++) {
        
        if (((Double_t)idx_ev/n_points)>(perc*0.01)) {
            cout<<"\t-> EdgeAccStudy events: [ "<<perc<<" % ]"<<endl;
            output_log_file<<"\t-> EdgeAccStudy events: [ "<<perc<<" % ]"<<endl;
            perc++;
        }
        
        acc_border->GetRandom2(costheta,phi);
        from_local_to_galactic(costheta,phi,l,b,sat_ra,sat_dec);
        
        if(isnan(l)) {
            idx_ev--;
            continue;
        }
        
        if (l>180.0)
            l-=360.0;
        
        h_healpix_Eacc->Fill(l,b,1);
        h_lb_uniform_Eacc->Fill(l,b,1);
        h_lcosb_uniform_Eacc->Fill(l,cos(b*TMath::DegToRad()),1);
        h_lcos2b_uniform_Eacc->Fill(l,TMath::Power(cos(b*TMath::DegToRad()),2),1);
        
    }
    
    h_healpix_Eacc->Write();
    h_lb_uniform_Eacc->Write();
    h_lcosb_uniform_Eacc->Write();
    h_lcos2b_uniform_Eacc->Write();
    
    h_healpix_Eacc->Delete();
    h_lb_uniform_Eacc->Delete();
    h_lcosb_uniform_Eacc->Delete();
    h_lcos2b_uniform_Eacc->Delete();
    
    
    ///////////////////////////////////////////// Considering acceptance flat in costheta and phi and evaluate the acceptance projection and its border
    
    healpix_SQacc = "Healpix_SQAcceptance_Sky_Projection_";
    lb_uniform_SQacc = "lb_uniform_SQAcceptance_Sky_Projection_";
    lcosb_uniform_SQacc = "lcosb_uniform_SQAcceptance_Sky_Projection_";
    lcos2b_uniform_SQacc = "lcos2b_uniform_SQAcceptance_Sky_Projection_";
    
    healpix_SQacc+=tree_idx;
    lb_uniform_SQacc+=tree_idx;
    lcosb_uniform_SQacc+=tree_idx;
    lcos2b_uniform_SQacc+=tree_idx;
    
    TH2D *h_healpix_SQacc = new TH2D(healpix_SQacc,"Uniform Acceptance Sky Projection; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
    TH2D *h_lb_uniform_SQacc = new TH2D(lb_uniform_SQacc,"Uniform Acceptance Sky Projection; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",1000,-180,180,1000,-90,90);
    TH2D *h_lcosb_uniform_SQacc = new TH2D(lcosb_uniform_SQacc,"Uniform Acceptance Sky Projection; l (#circ);  cos(b); Entries",1000,-180,180,1000,-1,1);
    TH2D *h_lcos2b_uniform_SQacc = new TH2D(lcos2b_uniform_SQacc,"Uniform Acceptance Sky Projection; l (#circ);  cos^{2}(b); Entries",1000,-180,180,1000,0,1);
    
    perc=0;
    
    cout<<endl;
    
    for(Int_t idx_ev=0; idx_ev<n_points; idx_ev++) {
        
        if (((Double_t)idx_ev/n_points)>(perc*0.01)) {
            cout<<"\t-> SQAccStudy events: [ "<<perc<<" % ]"<<endl;
            output_log_file<<"\t-> SQAccStudy events: [ "<<perc<<" % ]"<<endl;
            perc++;
        }
        
        costheta=rnd_gen->Uniform(costheta_min,1);
        phi=rnd_gen->Uniform(0,2*TMath::Pi());
        from_local_to_galactic(costheta,phi,l,b,sat_ra,sat_dec);
        
        if(isnan(l)) {
            idx_ev--;
            continue;
        }
        
        if (l>180.0)
            l-=360.0;
        
        h_healpix_SQacc->Fill(l,b,1);
        h_lb_uniform_SQacc->Fill(l,b,1);
        h_lcosb_uniform_SQacc->Fill(l,cos(b*TMath::DegToRad()),1);
        h_lcos2b_uniform_SQacc->Fill(l,TMath::Power(cos(b*TMath::DegToRad()),2),1);
        
    }

    h_healpix_SQacc->Write();
    h_lb_uniform_SQacc->Write();
    h_lcosb_uniform_SQacc->Write();
    h_lcos2b_uniform_SQacc->Write();
    
    h_healpix_SQacc->Delete();
    h_lb_uniform_SQacc->Delete();
    h_lcosb_uniform_SQacc->Delete();
    h_lcos2b_uniform_SQacc->Delete();

    ////////////////// And now just the border
    
    healpix_ESQacc = "Healpix_EDGE_SQAcceptance_Sky_Projection_";
    lb_uniform_ESQacc = "lb_uniform_EDGE_SQAcceptance_Sky_Projection_";
    lcosb_uniform_ESQacc = "lcosb_uniform_EDGE_SQAcceptance_Sky_Projection_";
    lcos2b_uniform_ESQacc = "lcos2b_uniform_EDGE_SQAcceptance_Sky_Projection_";
    
    healpix_ESQacc+=tree_idx;
    lb_uniform_ESQacc+=tree_idx;
    lcosb_uniform_ESQacc+=tree_idx;
    lcos2b_uniform_ESQacc+=tree_idx;
    
    TH2D *h_healpix_ESQacc = new TH2D(healpix_ESQacc,"Uniform Acceptance Sky Projection Edge; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
    TH2D *h_lb_uniform_ESQacc = new TH2D(lb_uniform_ESQacc,"Uniform Acceptance Sky Projection Edge; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",1000,-180,180,1000,-90,90);
    TH2D *h_lcosb_uniform_ESQacc = new TH2D(lcosb_uniform_ESQacc,"Uniform Acceptance Sky Projection Edge; l (#circ);  cos(b); Entries",1000,-180,180,1000,-1,1);
    TH2D *h_lcos2b_uniform_ESQacc = new TH2D(lcos2b_uniform_ESQacc,"Uniform Acceptance Sky Projection Edge; l (#circ);  cos^{2}(b); Entries",1000,-180,180,1000,0,1);
    
    perc=0;
    costheta=costheta_min; //Set min value of costheta for all the events
    cout<<endl;
    
    for(Int_t idx_ev=0; idx_ev<n_points; idx_ev++) {
        
        if (((Double_t)idx_ev/n_points)>(perc*0.01)) {
            cout<<"\t-> EdgeSQAccStudy events: [ "<<perc<<" % ]"<<endl;
            output_log_file<<"\t-> EdgeSQAccStudy events: [ "<<perc<<" % ]"<<endl;
            perc++;
        }
        
        
        phi=rnd_gen->Uniform(0,2*TMath::Pi());
        from_local_to_galactic(costheta,phi,l,b,sat_ra,sat_dec);
        
        if(isnan(l)) {
            idx_ev--;
            continue;
        }
        
        if (l>180.0)
            l-=360.0;
        
        h_healpix_ESQacc->Fill(l,b,1);
        h_lb_uniform_ESQacc->Fill(l,b,1);
        h_lcosb_uniform_ESQacc->Fill(l,cos(b*TMath::DegToRad()),1);
        h_lcos2b_uniform_ESQacc->Fill(l,TMath::Power(cos(b*TMath::DegToRad()),2),1);
        
    }
    
    h_healpix_ESQacc->Write();
    h_lb_uniform_ESQacc->Write();
    h_lcosb_uniform_ESQacc->Write();
    h_lcos2b_uniform_ESQacc->Write();
    
    h_healpix_ESQacc->Delete();
    h_lb_uniform_ESQacc->Delete();
    h_lcosb_uniform_ESQacc->Delete();
    h_lcos2b_uniform_ESQacc->Delete();
    
    //////////////////////// Find projection behaviour at the edges of the galactic map
    
    
    

}

    Double_t obtain_min_histo(TH1D *histo) {
        Double_t min_histo,bin_lenght;
        
        bin_lenght=(Double_t)1/histo->GetNbinsX();
        
        for(Int_t idx_b=1; idx_b<=histo->GetNbinsX(); idx_b++)
            if(histo->GetBinContent(idx_b)!=0) {
                min_histo=0.5*(idx_b*bin_lenght+(idx_b-1)*bin_lenght);
                break; //the first not-empty bin has been found. No need to go ahead
            }
        
        return min_histo;
    }

void get_acceptance_border(TH2D *acc,TH2D* acc_border) {
    
    for(Int_t y_bin=1; y_bin<=acc->GetNbinsY(); y_bin++)
        for(Int_t x_bin=1; x_bin<=acc->GetNbinsX(); x_bin++)
            if(acc->GetBinContent(x_bin,y_bin)!=0) {
                acc_border->SetBinContent(x_bin,y_bin,acc->GetBinContent(x_bin,y_bin));
                break; //I found the first not-empty bin regarding a such phi (or y) value. Now I have to choose a new phi bin and search for the first not-empty costheta bin
            }
}

void get_relevant_accceptance_points(TH2D* acc_border,vector<Double_t> &peacks_costheta,vector<Double_t> &valley_costheta,vector<Double_t> &peacks_phi,vector<Double_t> &valley_phi) {
    
    Double_t l_binX,l_binY,costheta,phi;
    Double_t new_costheta,new_phi;
    Int_t bunch=get_bunch_dimension(acc_border),n_bunch=acc_border->GetNbinsY()/bunch;
    Bool_t first_b_bin,first_b=true;
    
    Double_t tmp_min_costheta1,tmp_max_costheta1,tmp_min_phi1,tmp_max_phi1;
    Double_t tmp_min_costheta2,tmp_max_costheta2,tmp_min_phi2,tmp_max_phi2;
    Double_t tmp_min_costheta3,tmp_max_costheta3,tmp_min_phi3,tmp_max_phi3;
    
    l_binX=(Double_t)1/acc_border->GetNbinsX();
    l_binY=(Double_t)2*TMath::Pi()/acc_border->GetNbinsY();
    
    for(Int_t idx_b=0; idx_b<n_bunch; idx_b++) {
        first_b_bin=true;
        for(Int_t idx_bY=((idx_b*bunch)+1); idx_bY<=((idx_b+1)*bunch); idx_bY++) {
            for(Int_t idx_bX=1; idx_bX<=acc_border->GetNbinsX(); idx_bX++) {
                if(acc_border->GetBinContent(idx_bX,idx_bY)!=0) {
                    costheta=.5*(idx_bX*l_binX+(idx_bX-1)*l_binX);
                    phi=.5*(idx_bY*l_binY+(idx_bY-1)*l_binY);
                    if(first_b_bin) {
                        tmp_min_costheta1=tmp_max_costheta1=costheta;
                        tmp_min_phi1=tmp_max_phi1=phi;
                        first_b_bin=false;
                    }
                    else {
                        if(costheta>tmp_max_costheta1) {
                            tmp_max_costheta1=costheta;
                            tmp_max_phi1=phi;
                        }
                        if(costheta<tmp_min_costheta1) {
                            tmp_min_costheta1=costheta;
                            tmp_min_phi1=phi;
                        }
                    }
                }
            }
        }
        first_b_bin=true;
        for(Int_t idx_bY=(((idx_b+1)*bunch)+1); idx_bY<=((idx_b+2)*bunch); idx_bY++) {
        //for(Int_t idx_bY=(((idx_b+3)*bunch)+1); idx_bY<=((idx_b+4)*bunch); idx_bY++) {
            for(Int_t idx_bX=1; idx_bX<=acc_border->GetNbinsX(); idx_bX++) {
                if(acc_border->GetBinContent(idx_bX,idx_bY)!=0) {
                    new_costheta=.5*(idx_bX*l_binX+(idx_bX-1)*l_binX);
                    new_phi=.5*(idx_bY*l_binY+(idx_bY-1)*l_binY);
                    if(first_b_bin) {
                        tmp_min_costheta3=tmp_max_costheta3=new_costheta;
                        tmp_min_phi3=tmp_max_phi3=new_phi;
                        first_b_bin=false;
                    }
                    else {
                        if(costheta>tmp_max_costheta3) {
                            tmp_max_costheta3=new_costheta;
                            tmp_max_phi3=new_phi;
                        }
                        if(costheta<tmp_min_costheta3) {
                            tmp_min_costheta3=new_costheta;
                            tmp_min_phi3=new_phi;
                        }
                    }
                }
            }
        }
        if(first_b)
            first_b=false;
        else {
            if(tmp_min_costheta1<tmp_min_costheta2 && tmp_min_costheta1<tmp_min_costheta3) {
                peacks_costheta.push_back(tmp_min_costheta1);
                peacks_phi.push_back(tmp_min_phi1);
            }
            if(tmp_max_costheta1>tmp_max_costheta2 && tmp_max_costheta1>tmp_max_costheta3) {
                valley_costheta.push_back(tmp_max_costheta1);
                valley_phi.push_back(tmp_max_phi1);
            }
        }
        tmp_min_costheta2=tmp_min_costheta1;
        tmp_max_costheta2=tmp_max_costheta1;
        tmp_min_phi2=tmp_min_phi1;
        tmp_max_phi2=tmp_max_phi1;
    }
    
}

Int_t get_bunch_dimension(TH2D* acc_border) {
    Int_t bin_num=10,Ybins=acc_border->GetNbinsY();
    Bool_t found=false;
    
    while(found==false)
        if((Ybins%bin_num)==0)
            found=true;
        else
            bin_num++;
    
    return bin_num;
}
    
Bool_t check_Us(Float_t sat_ra[],Float_t sat_dec[],const Double_t &old_ang_xy,const Double_t &old_ang_xz,const Double_t &old_ang_yz) {
    
    Bool_t good_angle;
    Double_t ang_xy,ang_xz,ang_yz;
    Double_t eps=1e-6;
    
    Float_t ux1[3];
    Float_t uy1[3];
    Float_t uz1[3];
    Float_t rax = sat_ra[0];
    Float_t ray = sat_ra[1];
    Float_t raz = sat_ra[2];
    Float_t decx = sat_dec[0];
    Float_t decy = sat_dec[1];
    Float_t decz = sat_dec[2];
    
    ux1[0] = cos(decx)*cos(rax);
    ux1[1] = cos(decx)*sin(rax);
    ux1[2] = sin(decx);
    
    uy1[0] = cos(decy)*cos(ray);
    uy1[1] = cos(decy)*sin(ray);
    uy1[2] = sin(decy);
    
    uz1[0] = cos(decz)*cos(raz);
    uz1[1] = cos(decz)*sin(raz);
    uz1[2] = sin(decz);
    
    ang_xy=ux1[0]*uy1[0]+ux1[1]*uy1[1]+ux1[2]*uy1[2];
    ang_xz=ux1[0]*uz1[0]+ux1[1]*uz1[1]+ux1[2]*uz1[2];
    ang_yz=uy1[0]*uz1[0]+uy1[1]*uz1[1]+uy1[2]*uz1[2];
    
    if( abs(ang_xy-old_ang_xy)<=eps && abs(ang_xz-old_ang_xz)<=eps && abs(ang_yz-old_ang_yz)<=eps )
        good_angle=true;
    else
        good_angle=false;
    
    return good_angle;
    
}

void evaluate_Us(Double_t &old_ang_xy,Double_t &old_ang_xz,Double_t &old_ang_yz,Float_t sat_ra[],Float_t sat_dec[]) {
    
    Float_t ux1[3];
    Float_t uy1[3];
    Float_t uz1[3];
    Float_t rax = sat_ra[0];
    Float_t ray = sat_ra[1];
    Float_t raz = sat_ra[2];
    Float_t decx = sat_dec[0];
    Float_t decy = sat_dec[1];
    Float_t decz = sat_dec[2];
    
    ux1[0] = cos(decx)*cos(rax);
    ux1[1] = cos(decx)*sin(rax);
    ux1[2] = sin(decx);
    
    uy1[0] = cos(decy)*cos(ray);
    uy1[1] = cos(decy)*sin(ray);
    uy1[2] = sin(decy);
    
    uz1[0] = cos(decz)*cos(raz);
    uz1[1] = cos(decz)*sin(raz);
    uz1[2] = sin(decz);
    
    old_ang_xy=ux1[0]*uy1[0]+ux1[1]*uy1[1]+ux1[2]*uy1[2];
    old_ang_xz=ux1[0]*uz1[0]+ux1[1]*uz1[1]+ux1[2]*uz1[2];
    old_ang_yz=uy1[0]*uz1[0]+uy1[1]*uz1[1]+uy1[2]*uz1[2];
    
}


void shift_satellite_position(Float_t new_dec,Float_t new_ra,Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning) {
 
    Double_t b=0,l=0;
    Double_t costheta=-999, phi=-999;
    Double_t perc=0;
    
    Double_t old_ang_xy,old_ang_xz,old_ang_yz;
    evaluate_Us(old_ang_xy,old_ang_xz,old_ang_yz,sat_ra,sat_dec);
    
    /*
    cout<<endl;
    cout<<"XY: -->"<<old_ang_xy<<endl;
    cout<<"XZ: -->"<<old_ang_xz<<endl;
    cout<<"YZ: -->"<<old_ang_yz<<endl;
    cout<<endl;
    */
     
    TString file_name=satellite_shift_path;
    file_name+=new_ra;
    file_name+="_";
    file_name+=new_ra;
    file_name+="_satellite_shift.root";
    
    //////////////////// changing sat_ra[] and sat_dec[]
    
    new_dec*=TMath::DegToRad();
    new_ra*=TMath::DegToRad();
    
    
    TRandom3 *r_gen = new TRandom3();
    
    
    sat_ra[2]=new_ra;
    sat_dec[2]=new_dec;
    /*
    for(Int_t p_idx=0; p_idx<2; ++p_idx) {
        sat_ra[p_idx]=r_gen->Uniform(-TMath::Pi(),TMath::Pi());
        sat_dec[p_idx]=r_gen->Uniform(-TMath::Pi()/2.,TMath::Pi()/2.);
    }
    */
    
    do {
        for(Int_t p_idx=0; p_idx<2; ++p_idx) {
            sat_ra[p_idx]=r_gen->Uniform(-TMath::Pi(),TMath::Pi());
            sat_dec[p_idx]=r_gen->Uniform(-TMath::Pi()/2.,TMath::Pi()/2.);
        }
    } while(check_Us(sat_ra,sat_dec,old_ang_xy,old_ang_xz,old_ang_yz)==false);
    
   
    /*
    Double_t norm_ra = TMath::Sqrt( TMath::Power(sat_ra[0],2) + TMath::Power(sat_ra[1],2) + TMath::Power(sat_ra[2],2) );
    Double_t norm_dec = TMath::Sqrt( TMath::Power(sat_dec[0],2) + TMath::Power(sat_dec[1],2) + TMath::Power(sat_dec[2],2) );
    
    cout<<"\n\nsat_ra normalization: "<<norm_ra;
    cout<<"\nsat_dec normalization: "<<norm_dec;
    
    cout<<endl;
    for(Int_t p_idx=0; p_idx<3; p_idx++)
        cout<<"sat_ra["<<p_idx<<"]: "<<sat_ra[p_idx]<<endl;

    cout<<endl;
    for(Int_t p_idx=0; p_idx<3; p_idx++)
        cout<<"sat_dec["<<p_idx<<"]: "<<sat_dec[p_idx]<<endl;
    cout<<endl;
    */
    /*
    sat_ra[2]=new_ra;
    sat_dec[2]=new_dec;
    
    sat_ra[0]=0.435177;
    sat_dec[0]=2.49168;
    
    sat_ra[1]=2.17559;
    sat_dec[1]=2.92305;
    */
    //obtain_full_ra_dec_arrays(sat_ra,sat_dec);
    /*
    evaluate_Us(old_ang_xy,old_ang_xz,old_ang_yz,sat_ra,sat_dec);
    cout<<endl;
    cout<<"XY: -->"<<old_ang_xy<<endl;
    cout<<"XZ: -->"<<old_ang_xz<<endl;
    cout<<"YZ: -->"<<old_ang_yz<<endl;
    cout<<endl;
    */
    
    /////////////////////////////////////////////////

    TFile *out_file = new TFile(file_name.Data(),"RECREATE");
    out_file->cd();
    if(out_file->IsZombie()) {
        cout<<"\n\nError writing try output file\n\n";
        exit(-1);
    }
    
    TH2D *h_shift = new TH2D("h_shift","Acceptance Sky Projection Shift; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
    
    for(Int_t idx_ev=0; idx_ev<n_points; idx_ev++) {
        
        if (((Double_t)idx_ev/n_points)>(perc*0.01)) {
            cout<<"\t-> Shift Study events: [ "<<perc<<" % ]"<<endl;
            perc++;
        }
        
        acc->GetRandom2(costheta,phi);
        from_local_to_galactic(costheta,phi,l,b,sat_ra,sat_dec);
        
        if(isnan(l)) {
            idx_ev--;
            continue;
        }
        
        if (l>180.0)
            l-=360.0;
    
        
        h_shift->Fill(l,b,1);
        
    }
    h_shift->Write();
    h_shift->Delete();
    out_file->Close();
}

void obtain_full_ra_dec_arrays(Float_t sat_ra[],Float_t sat_dec[]) {
    
    Float_t ux1[3];
    Float_t uy1[3];
    Float_t uz1[3];
    
    Float_t raz = sat_ra[2];
    Float_t decz = sat_dec[2];
    
    uz1[0] = cos(decz)*cos(raz);
    uz1[1] = cos(decz)*sin(raz);
    uz1[2] = sin(decz);
    
    ///// Rotate respect X to obtain uy1
    
    uy1[0]=uz1[0];
    uy1[1]=-uz1[2];
    uy1[2]=uz1[1];
    
    ///// Rotate respect Z to obtain ux1
    
    ux1[0]=uz1[2];
    ux1[1]=uz1[0];
    ux1[2]=uz1[1];

    //////////////////////////////////////////////////////////
    
    Double_t norm_x = TMath::Sqrt( TMath::Power(ux1[0],2) + TMath::Power(ux1[1],2) + TMath::Power(ux1[2],2) );
    Double_t norm_y = TMath::Sqrt( TMath::Power(uy1[0],2) + TMath::Power(uy1[1],2) + TMath::Power(uy1[2],2) );
    Double_t norm_z = TMath::Sqrt( TMath::Power(uz1[0],2) + TMath::Power(uz1[1],2) + TMath::Power(uz1[2],2) );
    
    
    cout<<"\n\nX normalization: "<<norm_x;
    cout<<"\nY normalization: "<<norm_y;
    cout<<"\nY normalization: "<<norm_z<<endl;
    
    cout<<"\n";
    cout<<"cos(XY) --> "<<ux1[0]*uy1[0]+ux1[1]*uy1[1]+ux1[2]*uy1[2];
    cout<<"cos(XZ) --> "<<ux1[0]*uz1[0]+ux1[1]*uz1[1]+ux1[2]*uz1[2];
    cout<<"cos(YZ) --> "<<uy1[0]*uz1[0]+uy1[1]*uz1[1]+uy1[2]*uz1[2];
    
    
    //////////////////////////////////////////////////////////
    
    /////// Obtaining sat_ra and sat_dec

    sat_dec[1]=asin(uy1[2]);
    sat_ra[1]=acos(uy1[0]/cos(sat_dec[1]));
    
    sat_dec[0]=asin(ux1[2]);
    sat_ra[0]=acos(ux1[0]/cos(sat_dec[0]));
    
}

void plot_POI(vector<Double_t> &peacks_costheta,vector<Double_t> &peacks_phi,vector<Double_t> &valley_costheta,vector<Double_t> &valley_phi,TFile* results_file,Int_t tree_idx,Float_t sat_ra[],Float_t sat_dec[],ofstream &output_log_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning,vector<Double_t> &p_dNS,vector<Double_t> &p_dEW,vector<Double_t> &v_dNS,vector<Double_t> &v_dEW,vector<Double_t> &sp_dNS,vector<Double_t> &sp_dEW,vector<Double_t> &sv_dNS,vector<Double_t> &sv_dEW,Double_t &slope_pNS,Double_t &slope_pEW,Double_t &slope_vNS,Double_t &slope_vEW) {

  ////////////////////////////////////////////////////////// Creating vectors for max and min of galactic latitude and longitude (obviously depending on DAMPE position)
  
    Double_t l=0,b=0;
    Bool_t all_POI=true;
    Double_t distance_pNS,distance_pEW;
    Double_t distance_vNS,distance_vEW;
    
    Double_t s_distance_pNS,s_distance_pEW;
    Double_t s_distance_vNS,s_distance_vEW;
    
  vector<Float_t> peacks_glat,peacks_glon;
  vector<Float_t> valley_glat,valley_glon;
  Int_t n_peacks=peacks_costheta.size(),n_valley=valley_costheta.size();

  peacks_glat.resize(n_peacks);
  peacks_glon.resize(n_peacks);
  valley_glat.resize(n_valley);
  valley_glon.resize(n_valley);
  
  for(Int_t idx_p=0; idx_p<n_peacks; idx_p++) {
    from_local_to_galactic(peacks_costheta.at(idx_p),peacks_phi.at(idx_p),l,b,sat_ra,sat_dec);
    if(isnan(l)) {
      all_POI=false;
      break;
    }
    if (l>180.0)
      l-=360.0;
    peacks_glat[idx_p]=b;
    peacks_glon[idx_p]=l;
  }

  if(all_POI) {
    for(Int_t idx_v=0; idx_v<n_valley; idx_v++) {
      from_local_to_galactic(valley_costheta.at(idx_v),valley_phi.at(idx_v),l,b,sat_ra,sat_dec);
      if(isnan(l)) {
	all_POI=false;
	break;
      }
      if (l>180.0)
	l-=360.0;
      valley_glat[idx_v]=b;
      valley_glon[idx_v]=l;
    }  
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /*
    cout<<"\n\n Peacks: \n";
    for(Int_t idx=0; idx<5; idx++) {
        cout<<idx<<"\t"<<peacks_glon[idx]<<endl;
    }
    
    cout<<"\n\nValley: \n";
    for(Int_t idx=0; idx<3; idx++) {
        cout<<idx<<"\t"<<valley_glon[idx]<<endl;
    }
    
    exit(-1);
    */
     
  if(all_POI) {
    results_file->cd();
      
      if(points_are_inside(peacks_glon,valley_glon)) {
          
          //cout<<"\npoints inside: "<<tree_idx;
          
          ///////////////////////////////// Filling vectors for graphs
          
          distance_pNS=TMath::Sqrt(TMath::Power(peacks_glat[4]-peacks_glat[1],2)+TMath::Power(peacks_glon[4]-peacks_glon[1],2));
          distance_pEW=TMath::Sqrt(TMath::Power(peacks_glat[0]-peacks_glat[3],2)+TMath::Power(peacks_glon[0]-peacks_glon[3],2));
          distance_vNS=TMath::Sqrt(TMath::Power(valley_glat[1]-valley_glat[0],2)+TMath::Power(valley_glon[1]-valley_glon[0],2));
          distance_vEW=TMath::Sqrt(TMath::Power(valley_glat[2]-peacks_glat[2],2)+TMath::Power(valley_glon[2]-peacks_glon[2],2));
          
          s_distance_pNS=TMath::Sqrt(TMath::Power(peacks_glat[4]-peacks_glat[1],2)+TMath::Power(peacks_glon[4]-peacks_glon[1],2));
          s_distance_pEW=TMath::Sqrt(TMath::Power(peacks_glat[0]-peacks_glat[3],2)+TMath::Power(peacks_glon[0]-peacks_glon[3],2));
          s_distance_vNS=TMath::Sqrt(TMath::Power(valley_glat[1]-valley_glat[0],2)+TMath::Power(valley_glon[1]-valley_glon[0],2));
          s_distance_vEW=TMath::Sqrt(TMath::Power(valley_glat[2]-peacks_glat[2],2)+TMath::Power(valley_glon[2]-peacks_glon[2],2));
          
          p_dNS.push_back(distance_pNS);
          p_dEW.push_back(distance_pEW);
          v_dNS.push_back(distance_vNS);
          v_dEW.push_back(distance_vEW);
          
          sp_dNS.push_back(distance_pNS);
          sp_dEW.push_back(distance_pEW);
          sv_dNS.push_back(distance_vNS);
          sv_dEW.push_back(distance_vEW);
          
          
          //////////////////////////////// Calculate slopes
      
          slope_pNS=acos(abs(peacks_glat[4]-peacks_glat[1])/distance_pNS);
          slope_pEW=acos(abs(peacks_glat[0]-peacks_glat[3])/distance_pEW);
          slope_vNS=acos(abs(valley_glat[1]-valley_glat[0])/distance_vNS);
          slope_vEW=acos(abs(valley_glat[2]-peacks_glat[2])/distance_vEW);
    
      }
      else{
          
          //cout<<"\npoints outside: "<<tree_idx;
          
          ///////////////////////////////// Filling vectors for graphs
          
          s_distance_pNS=peacks_glon[1]/sin(slope_pNS);
          s_distance_pNS+=(180-peacks_glon[4])/sin(slope_pNS);
          
          s_distance_pEW=peacks_glon[3]/sin(slope_pEW);
          s_distance_pEW+=(180-peacks_glon[0])/sin(slope_pEW);
          
          s_distance_vNS=valley_glon[0]/sin(slope_vNS);
          s_distance_vNS+=(180-valley_glon[1])/sin(slope_vNS);
          
          s_distance_vEW=peacks_glon[2]/sin(slope_vEW);
          s_distance_vEW+=(180-valley_glon[2])/sin(slope_vEW);
          
          sp_dNS.push_back(distance_pNS);
          sp_dEW.push_back(distance_pEW);
          sv_dNS.push_back(distance_vNS);
          sv_dEW.push_back(distance_vEW);
          
          //////////////////////////////// Calculate slopes
          
          slope_pNS=acos(abs(peacks_glat[4]-peacks_glat[1])/distance_pNS);
          slope_pEW=acos(abs(peacks_glat[0]-peacks_glat[3])/distance_pEW);
          slope_vNS=acos(abs(valley_glat[1]-valley_glat[0])/distance_vNS);
          slope_vEW=acos(abs(valley_glat[2]-peacks_glat[2])/distance_vEW);
          
          ///////////////////////////////// Shifting all points by 180 degrees. So they will fit inside the whole graph's canvas
          
          for(Int_t idx_v=0; idx_v<5; idx_v++) {
              if(idx_v<3) {
                  if(peacks_glon[idx_v]<0)
                      peacks_glon[idx_v]+=180;
                  else
                      peacks_glon[idx_v]-=180;
                  if(valley_glon[idx_v]<0)
                      valley_glon[idx_v]+=180;
                  else
                      valley_glon[idx_v]-=180;
              }
              else {
                  if(peacks_glon[idx_v]<0)
                      peacks_glon[idx_v]+=180;
                  else
                      peacks_glon[idx_v]-=180;
              }
          }
          
          
          ///////////////////////////////// Filling vectors for graphs
          
          distance_pNS=TMath::Sqrt(TMath::Power(peacks_glat[4]-peacks_glat[1],2)+TMath::Power(peacks_glon[4]-peacks_glon[1],2));
          distance_pEW=TMath::Sqrt(TMath::Power(peacks_glat[0]-peacks_glat[3],2)+TMath::Power(peacks_glon[0]-peacks_glon[3],2));
          distance_vNS=TMath::Sqrt(TMath::Power(valley_glat[1]-valley_glat[0],2)+TMath::Power(valley_glon[1]-valley_glon[0],2));
          distance_vEW=TMath::Sqrt(TMath::Power(valley_glat[2]-peacks_glat[2],2)+TMath::Power(valley_glon[2]-peacks_glon[2],2));
          
          p_dNS.push_back(distance_pNS);
          p_dEW.push_back(distance_pEW);
          v_dNS.push_back(distance_vNS);
          v_dEW.push_back(distance_vEW);
          
      }
    
      
      
      
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    TString POI_healpix = "POI_Healpix_";
    TString POI_lb_uniform = "POI_lb_uniform_";
    TString POI_lcosb_uniform = "POI_lcosb_uniform_";
    TString POI_lcos2b_uniform = "POI_lcos2b_uniform_";
    
    POI_healpix+=tree_idx;
    POI_lb_uniform+=tree_idx;
    POI_lcosb_uniform+=tree_idx;
    POI_lcos2b_uniform+=tree_idx;
    
    TH2D *h_POI_healpix = new TH2D(POI_healpix,"POI Healpix; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
    TH2D *h_POI_lb_uniform = new TH2D(POI_lb_uniform,"Uniform POI; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",1000,-180,180,1000,-90,90);
    TH2D *h_POI_lcosb_uniform = new TH2D(POI_lcosb_uniform,"Uniform POI; l (#circ);  cos(b); Entries",1000,-180,180,1000,-1,1);
    TH2D *h_POI_lcos2b_uniform = new TH2D(POI_lcos2b_uniform,"Uniform POI; l (#circ);  cos^{2}(b); Entries",1000,-180,180,1000,0,1);
    
    h_POI_healpix->SetMarkerStyle(20);
    h_POI_lb_uniform->SetMarkerStyle(20);
    h_POI_lcosb_uniform->SetMarkerStyle(20);
    h_POI_lcos2b_uniform->SetMarkerStyle(20);
    
    for(Int_t idx_p=0; idx_p<n_peacks; idx_p++) {
      h_POI_healpix->Fill(peacks_glon[idx_p],peacks_glat[idx_p],1);
      h_POI_lb_uniform->Fill(peacks_glon[idx_p],peacks_glat[idx_p],1);
      h_POI_lcosb_uniform->Fill(peacks_glon[idx_p],cos(peacks_glat[idx_p]),1);
      h_POI_lcos2b_uniform->Fill(peacks_glon[idx_p],TMath::Power(cos(peacks_glat[idx_p]),2),1); 
    }

    for(Int_t idx_p=0; idx_p<n_valley; idx_p++) {
      h_POI_healpix->Fill(valley_glon[idx_p],valley_glat[idx_p],1);
      h_POI_lb_uniform->Fill(valley_glon[idx_p],valley_glat[idx_p],1);
      h_POI_lcosb_uniform->Fill(valley_glon[idx_p],cos(valley_glat[idx_p]),1);
      h_POI_lcos2b_uniform->Fill(valley_glon[idx_p],TMath::Power(cos(valley_glat[idx_p]),2),1);
    }
    
    h_POI_healpix->Write();
    h_POI_lb_uniform->Write();
    h_POI_lcosb_uniform->Write();
    h_POI_lcos2b_uniform->Write();
    
    h_POI_healpix->Delete();
    h_POI_lb_uniform->Delete();
    h_POI_lcosb_uniform->Delete();
    h_POI_lcos2b_uniform->Delete();
    
  }
}
  
void calculate_POI_stuff(Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,vector<Double_t> &peacks_costheta,vector<Double_t> &peacks_phi,vector<Double_t> &valley_costheta,vector<Double_t> &valley_phi) {
  
  // That's the first call ! Acceptance EDGE and POI should be calculated !
  
  Int_t n_peacks,n_valley;
  
  /////// Extract acceptance border
  
  static TH2D *acc_border = (TH2D*)acc->Clone("acc_border");
  acc_border->Reset();
  get_acceptance_border(acc,acc_border);
  acc_border->Write();
  
  /////////////////////////////////
  
  /////// Extract acceptance POI
  
  get_relevant_accceptance_points(acc_border,peacks_costheta,valley_costheta,peacks_phi,valley_phi);

}

void POI_discrete_plotting(Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,vector<Double_t> &peacks_costheta,vector<Double_t> &peacks_phi,vector<Double_t> &valley_costheta,vector<Double_t> &valley_phi,TFile *results_file,Int_t tree_idx,ofstream &output_log_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* binning,Bool_t &first_call,vector<Double_t> &p_dNS,vector<Double_t> &p_dEW,vector<Double_t> &v_dNS,vector<Double_t> &v_dEW,vector<Double_t> &sp_dNS,vector<Double_t> &sp_dEW,vector<Double_t> &sv_dNS,vector<Double_t> &sv_dEW,vector<Double_t> &iterations,Double_t &slope_pNS,Double_t &slope_pEW,Double_t &slope_vNS,Double_t &slope_vEW) {
  
  // This kind of study, at the countrary, should be done on all DAMPE acquisition period
  
  if(first_call) {
    calculate_POI_stuff(sat_ra,sat_dec,acc,peacks_costheta,peacks_phi,valley_costheta,valley_phi);
    plot_POI(peacks_costheta,peacks_phi,valley_costheta,valley_phi,results_file,tree_idx,sat_ra,sat_dec,output_log_file,n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning,p_dNS,p_dEW,v_dNS,v_dEW,sp_dNS,sp_dEW,sv_dNS,sv_dEW,slope_pNS,slope_pEW,slope_vNS,slope_vEW);
        first_call=false;
  }
  plot_POI(peacks_costheta,peacks_phi,valley_costheta,valley_phi,results_file,tree_idx,sat_ra,sat_dec,output_log_file,n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning,p_dNS,p_dEW,v_dNS,v_dEW,sp_dNS,sp_dEW,sv_dNS,sv_dEW,slope_pNS,slope_pEW,slope_vNS,slope_vEW);

  iterations.push_back(tree_idx);
  
}

Bool_t points_are_inside(vector<Float_t> &peacks_glon,vector<Float_t> &valley_glon) {
    Bool_t all_points_inside=true;
    
    /*
    for(Int_t idx_v=0; idx_v<5; idx_v++) {
        if(idx_v==0) {
            if(valley_glon[idx_v]>peacks_glon[0]) {
                all_points_inside=false;
                break;
            }
        }
        if(idx_v>0 && idx_v<3) {
            if(peacks_glon[idx_v]>peacks_glon[0] || valley_glon[idx_v]>peacks_glon[0]) {
                all_points_inside=false;
                break;
            }
        }
        else {
            if(peacks_glon[idx_v]>peacks_glon[0]) {
                all_points_inside=false;
                break;
            }
        }
    }
     */
    
    Double_t max=peacks_glon[0],min=max;
    
    for(Int_t idx_v=0; idx_v<3; idx_v++) {
        if(valley_glon[idx_v]>max)
            max=valley_glon[idx_v];
        if(valley_glon[idx_v]<min)
            min=valley_glon[idx_v];
    }
    
    for(Int_t idx_v=1; idx_v<5; idx_v++) {
        if(peacks_glon[idx_v]>max)
            max=peacks_glon[idx_v];
        if(peacks_glon[idx_v]<min)
            min=peacks_glon[idx_v];
    }
    
    if((max-min)>180.)
        all_points_inside=false;
        
    return all_points_inside;
}
