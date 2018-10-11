
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <ctime>
#include <vector>
#include <string>
#include <cmath>

///// ROOT libraries 

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

#include "orbitStruct.h"

using namespace std;

#define EPS 1.e-12

///////////////////////////////////// Simulation variables

//const static Int_t sky_events = 1e+9;
const static Int_t sky_events = 1e+3;
const static UInt_t random_seed = 22;
const static Int_t n_gevents=1e+9;              //Number of generated events for the MC calculation of galactic latitude and longitude edges
const static Int_t n_points=1e+6;               //Number of generated point for the high statistic acceptance sky projection
const static Int_t POI_interval=100;            //Interval, in seconds, between POI plots
const static Bool_t plot_POI_interval=true;    //Boolean variable regarding POI plot: if "TRUE" software plots POI each "POI_interval" seconds
const static Bool_t shift_pointing=false;        //Boolean variable regarding DAMPE pointing shift

const static time_t time_stamp = time(0);       //Setting timestamp for the out files

const static TString sbi_path = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/Stuff/acceptance-sky-projection-svn/trunk/SBI_data/";
const static TString acceptance_final_plot = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/dampe-gacceptance-svn/trunk/results/1516819290_acceptance_result.root";
const static TString sbi_subsample = "010";
const static string string_sbi_subsample = "010";
const static Int_t number_SBI_files = 3;       // To be precise, at the moment of sotware writing, they are 0102800000_SBI.root 0102900000_SBI.root 0103000000_SBI.root

const static string output_log = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/Stuff/acceptance-sky-projection-svn/trunk/logs/";
const static string output_root = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/Stuff/acceptance-sky-projection-svn/trunk/results/";

const static TString satellite_shift_path = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/Stuff/acceptance-sky-projection-svn/trunk/results/";

const static Bool_t save_output_file = true;   //Variable to save result ROOT file 
const static Bool_t MC_acc_proj = true;        //Switch to compute MC simulation for the acceptance projection on galactic coordinates
const static Bool_t wanna_just_project_acc = false; //Switch to project or not acceptance into galactic coordinates. If false just other study are done

////////////////////////////////////////////////////////////////////////

extern string output_path_creator(const Int_t out_choose);
extern void log_file_init(ofstream &out_file);
extern void read_SBI_data(Float_t &galactic_lat,Float_t &galactic_lon,Float_t &geographic_lat,Float_t &geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_events,bool &good,TChain* tree,TString sbi_data_path,ofstream &out_file);
extern Bool_t check_sbi_loading(Float_t galactic_lat,Float_t galactic_lon,Float_t geographic_lat,Float_t geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t sec,UShort_t n_events);
extern Bool_t chech_if_null_variable(Float_t in_variable);
extern void reinitialize_all_variables(Float_t &galactic_lat,Float_t &galactic_lon,Float_t &geographic_lat,Float_t &geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_events);
extern void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale);
extern void project_acceptance(TFile *results,Int_t h_idx,Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,TH1D *acc_costheta,TH1D* acc_phi,ofstream &out_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning);
extern void from_satellite_to_celestial(Float_t ra[],Float_t dec[],double vectorin[],AtPolarVect &vector_out);
extern void AtVect_To_AtPolarVect(double in_vector[],AtPolarVect &vector_out);
extern void invert_AtPolarVect_direction(AtPolarVect vector_out,AtPolarVect &vector_out_inv);
extern void AtPolarVect_to_vector(AtPolarVect &input_polar,double out_array[]);
extern void from_celestial_to_galactic(Double_t ra,Double_t dec,Double_t &l,Double_t &b);
extern void from_local_to_galactic(Double_t costheta,Double_t phi,Double_t &l,Double_t &b,Float_t sat_ra[],Float_t sat_dec[]);
extern void obtain_edge_coordinate(Int_t tree_idx,Float_t sat_ra[],Float_t sat_dec[],TH2D *h_acc);
extern void project_acceptance_edgeMC(TFile *results,Int_t h_idx,Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,TH1D *acc_costheta,TH1D* acc_phi,ofstream &out_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning);
extern void obtain_edge(Float_t sat_ra[],Float_t sat_dec[],TH2D *h_acc,Double_t &acc_glat_min,Double_t &acc_glat_max,Double_t &acc_glon_min,Double_t &acc_glon_max);
extern void acceptance_projection_study(TFile *results_file,Int_t tree_idx,Float_t sat_ra[],Float_t sat_dec[],TH2D *acc,TH1D *acc_costheta,ofstream &output_log_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning);
extern Double_t obtain_min_histo(TH1D *histo);
extern void get_acceptance_border(TH2D *acc,TH2D* acc_border);
extern void get_relevant_accceptance_points(TH2D* acc_border,vector<Double_t> &peacks_costheta,vector<Double_t> &valley_costheta,vector<Double_t> &peacks_phi,vector<Double_t> &valley_phi);
extern Int_t get_bunch_dimension(TH2D* acc_border);
extern void plot_POI(vector<Double_t> &peacks_costheta,vector<Double_t> &peacks_phi,vector<Double_t> &valley_costheta,vector<Double_t> &valley_phi,TFile* results_file,Int_t tree_idx,Float_t sat_ra[],Float_t sat_dec[],ofstream &output_log_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning,vector<Double_t> &p_dNS,vector<Double_t> &p_dEW,vector<Double_t> &v_dNS,vector<Double_t> &v_dEW,vector<Double_t> &sp_dNS,vector<Double_t> &sp_dEW,vector<Double_t> &sv_dNS,vector<Double_t> &sv_dEW,Double_t &slope_pNS,Double_t &slope_pEW,Double_t &slope_vNS,Double_t &slope_vEW);
extern void calculate_POI_stuff(Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,vector<Double_t> &peacks_costheta,vector<Double_t> &peacks_phi,vector<Double_t> &valley_costheta,vector<Double_t> &valley_phi);
extern void shift_satellite_position(Float_t new_dec,Float_t new_ra,Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* &binning);
extern void POI_discrete_plotting(Float_t sat_ra[],Float_t sat_dec[],TH2D* acc,vector<Double_t> &peacks_costheta,vector<Double_t> &peacks_phi,vector<Double_t> &valley_costheta,vector<Double_t> &valley_phi,TFile *results_file,Int_t tree_idx,ofstream &output_log_file,Int_t n_bin_lon,Double_t lon_bin_min,Double_t lon_bin_max,Int_t n_bin_lat,Double_t* binning,Bool_t &first_call,vector<Double_t> &p_dNS,vector<Double_t> &p_dEW,vector<Double_t> &v_dNS,vector<Double_t> &v_dEW,vector<Double_t> &sp_dNS,vector<Double_t> &sp_dEW,vector<Double_t> &sv_dNS,vector<Double_t> &sv_dEW,vector<Double_t> &iterations,Double_t &slope_pNS,Double_t &slope_pEW,Double_t &slope_vNS,Double_t &slope_vEW);
extern Bool_t check_Us(Float_t sat_ra[],Float_t sat_dec[],const Double_t &old_ang_xy,const Double_t &old_ang_xz,const Double_t &old_ang_yz);
extern void evaluate_Us(Double_t &old_ang_xy,Double_t &old_ang_xz,Double_t &old_ang_yz,Float_t sat_ra[],Float_t sat_dec[]);
extern void obtain_full_ra_dec_arrays(Float_t sat_ra[],Float_t sat_dec[]);
extern Bool_t points_are_inside(vector<Float_t> &peacks_glon,vector<Float_t> &valley_glon);
