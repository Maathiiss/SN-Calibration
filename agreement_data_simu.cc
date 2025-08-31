#include <math.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstring>
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TApplication.h"
#include "TMultiGraph.h"
#include "TFeldmanCousins.h"
#include "TGaxis.h"
#include "TLeaf.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include "TRandom.h"
#include <TLegend.h>
#include <TParameter.h>
#include <TSpectrum.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TSystem.h>
#include "/home/granjon/Documents/stage_radon/Stage/Mathis/sndisplay/sndisplay.cc"
using namespace std;

int suppr_bad_om(int om_number_data){
  std::ifstream infile("../../list_om_tab.txt");
  std::string line;
  std::getline(infile, line);
  vector<double> long_term_jump_bi_vec, long_term_jump_li_vec, Se_num_vec, b_num_vec;
  int long_term_jump_bi_num, long_term_jump_li_num, Se_num, b_num;
  while (std::getline(infile, line)) {
    std::stringstream ss(line);
    std::string val;
    std::getline(ss, val, ';');
    long_term_jump_bi_num = val.empty() ? -1 : std::stoi(val);
    std::getline(ss, val, ';');
    long_term_jump_li_num = val.empty() ? -1 : std::stoi(val);
    std::getline(ss, val, ';');
    Se_num = val.empty() ? -1 : std::stoi(val);
    std::getline(ss, val, ';');
    b_num = val.empty() ? -1 : std::stoi(val);
    if (long_term_jump_bi_num != -1) long_term_jump_bi_vec.push_back(long_term_jump_bi_num);
    if (long_term_jump_li_num != -1) long_term_jump_li_vec.push_back(long_term_jump_li_num);
    if (Se_num != -1) Se_num_vec.push_back(Se_num);
    if (b_num != -1) b_num_vec.push_back(b_num);
  }
  if (!(std::find(long_term_jump_bi_vec.begin(), long_term_jump_bi_vec.end(), om_number_data) != long_term_jump_bi_vec.end() ||
	std::find(long_term_jump_li_vec.begin(), long_term_jump_li_vec.end(), om_number_data) != long_term_jump_li_vec.end())) {
    return 0;
  }
  else{
    return 1;
  }    
}


int compute_num_om(int type, int module_useless, int side, int column, int row){
  if(type==1302){
    return side*260+column*13+row;
  }
  return -1;
}


double getGeomFWHM(int om_number) {
  std::ifstream infile("calorimeter_regime_database_v1.db");
  std::string geomid, fwhm_str;
  int d0, d1, d2, d3, d4;
  double fwhm;

  while (infile >> geomid >> fwhm_str) {
    if (geomid.empty() || geomid[0] != '[') continue;
    std::string all = geomid.substr(1, geomid.size()-2); // retire [ ]
    std::replace(all.begin(), all.end(), ':', ' ');
    std::replace(all.begin(), all.end(), '.', ' ');
    std::stringstream ss(all);
    if (!(ss >> d0 >> d1 >> d2 >> d3 >> d4)) continue;
    fwhm_str = fwhm_str.substr(0, fwhm_str.find('%'));
    fwhm = std::stod(fwhm_str);
    if(om_number == compute_num_om(d0,d1,d2,d3,d4)){
      return fwhm;
    }
  }
  return 0;
}



void draw_pb_om() {
  //sncalo = new sndisplay::calorimeter("sncalo_display",true);

    TFile *f = TFile::Open("../flat_pb_Final_merged_1000.root");
    TTree *tree = (TTree*)f->Get("MergedTree");
    int om_numberr, source_number;
    double FWHM_1MeV, FWHM_simu, FWHM_simu_error, FWHM_1MeV_error;
    tree->SetBranchAddress("om_number", &om_numberr);
    tree->SetBranchAddress("source_number", &source_number);
    tree->SetBranchAddress("FWHM_1MeV", &FWHM_1MeV);
    tree->SetBranchAddress("FWHM_1MeV_error", &FWHM_1MeV_error);
    tree->SetBranchAddress("FWHM_simu", &FWHM_simu);
    tree->SetBranchAddress("FWHM_simu_error", &FWHM_simu_error);
    std::map<int, std::pair<double, double>> bestRatio;
    
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if (FWHM_1MeV != 0 && FWHM_simu != 0) {
	  double ratio = fabs(FWHM_1MeV / FWHM_simu);
            double ratio_error = ratio * sqrt(pow(FWHM_1MeV_error/FWHM_1MeV,2) + pow(FWHM_simu_error/FWHM_simu,2));
	    //if(om_numberr == 1){
	    cout<<"om number = "<<om_numberr <<" ratio = "<<ratio<<" fwhm = "<< FWHM_1MeV <<endl;
	      //}
            if (bestRatio.find(om_numberr) == bestRatio.end() || fabs(ratio - 1.0) < fabs(bestRatio[om_numberr].first - 1.0)) {
                bestRatio[om_numberr] = std::make_pair(ratio, ratio_error);
		//if(om_numberr == 1){
		  cout<<"best ratio = "<<ratio<<endl;
		  //}
			    
            }
        }
    }

    TFile *outFile = new TFile("FWHM_ratios.root", "RECREATE");
    TTree *outTree = new TTree("BestRatios", "Ratios compatibles avec 1 sigma");

    int out_source, out_om;
    double out_FWHM_1MeV, out_FWHM_simu, out_FWHM_1MeV_error, out_FWHM_simu_error, out_ratio, FWHM_om_ageing, FWHM_om_real, fwhm_bordeaux;

    outTree->Branch("source_number", &out_source);
    outTree->Branch("om_number", &out_om);
    outTree->Branch("FWHM_1MeV", &out_FWHM_1MeV);
    outTree->Branch("FWHM_simu", &out_FWHM_simu);
    outTree->Branch("FWHM_1MeV_error", &out_FWHM_1MeV_error);
    outTree->Branch("FWHM_simu_error", &out_FWHM_simu_error);
    outTree->Branch("FWHM_om_ageing", &FWHM_om_ageing);
    outTree->Branch("FWHM_om_real", &FWHM_om_real);
    outTree->Branch("fwhm_bordeaux", &fwhm_bordeaux);
    outTree->Branch("ratio", &out_ratio);

    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if (FWHM_1MeV != 0 && FWHM_simu != 0 && suppr_bad_om(om_numberr)==0) {
	  double ratio = fabs(FWHM_1MeV / FWHM_simu);
            double ratio_error = ratio * sqrt(pow(FWHM_1MeV_error/FWHM_1MeV,2) + pow(FWHM_simu_error/FWHM_simu,2));
            double best = bestRatio[om_numberr].first;
            double best_err = bestRatio[om_numberr].second;
            //if (ratio >= best - best_err && ratio <= best + best_err) {
	    if(ratio==best){
                out_source = source_number;
                out_om = om_numberr;
                out_FWHM_1MeV = FWHM_1MeV;
                out_FWHM_simu = FWHM_simu;
                out_FWHM_1MeV_error = FWHM_1MeV_error;
                out_FWHM_simu_error = FWHM_simu_error;
                out_ratio = ratio;
		FWHM_om_ageing = sqrt(pow(out_FWHM_1MeV,2)-pow(out_FWHM_simu,2));
		fwhm_bordeaux = getGeomFWHM(om_numberr);
		FWHM_om_real = sqrt(pow(fwhm_bordeaux,2)+pow(FWHM_om_ageing,2));
		//sncalo->setcontent(om_numberr,ratio);
                outTree->Fill();
            }
        }
    }

    outTree->Write();
    outFile->Close();
    f->Close();
    //   sncalo->draw1();
    // sncalo->draw1();
    // sncalo->canvas->SaveAs("test.png");

}





void compute_best_fwhm_spectra(){
  TFile *file = new TFile("../simu_flat_70M.root", "READ");
  TTree* Result_tree = (TTree*)file->Get("Event");
  vector<int>* num_om = nullptr;
  vector<int>* number_of_kinks_per_track = nullptr;
  vector<int>* num_source = nullptr;
  vector<double>* ellipse_source = nullptr;
  vector<double>* evis_bcu = nullptr;
  int number_of_electrons, om_num;
  double fwhm_compute, E_new;
  
  gROOT->cd();
  Result_tree->SetBranchStatus("*",0);
  Result_tree->SetBranchStatus("*",0);
  Result_tree->SetBranchStatus("num_om",1);
  Result_tree->SetBranchAddress("num_om", &num_om);
  Result_tree->SetBranchStatus("number_of_kinks_per_track",1);
  Result_tree->SetBranchAddress("number_of_kinks_per_track", &number_of_kinks_per_track);
  Result_tree->SetBranchStatus("number_of_electrons",1);
  Result_tree->SetBranchAddress("number_of_electrons", &number_of_electrons);
  Result_tree->SetBranchStatus("num_source",1);
  Result_tree->SetBranchAddress("num_source", &num_source);
  Result_tree->SetBranchStatus("ellipse_source",1);
  Result_tree->SetBranchAddress("ellipse_source", &ellipse_source);
  Result_tree->SetBranchStatus("evis_bcu",1);
  Result_tree->SetBranchAddress("evis_bcu", &evis_bcu);

  TFile *file_create = new TFile("several_fwhm.root", "RECREATE");
  TTree final_tree("final_tree","");
  final_tree.Branch("fwhm",&fwhm_compute);
  final_tree.Branch("om_num",&om_num);
  final_tree.Branch("E_new",&E_new);

  TRandom rando;  
  for(int entry=0; entry<Result_tree->GetEntries();entry++){
    Result_tree->GetEntry(entry);
    if(entry%100000==0){
      cout<< 100*entry/Result_tree->GetEntries()<<" percent"<<endl;
    }
    for(int k=0; k<number_of_electrons; k++){
      if(num_om->at(k)>519){continue;}
      om_num = num_om->at(k);
      if(number_of_kinks_per_track->at(k)==0){
        if(ellipse_source->at(k)<1){
	  double fwhm = getGeomFWHM(num_om->at(k));
	  for(int i=0; i<100; i++){
	    double sigma_add = sqrt(pow((10+i/10)/235.5*evis_bcu->at(k), 2) - pow(fwhm/235.5*evis_bcu->at(k), 2));
	    E_new = rando.Gaus(evis_bcu->at(k), sigma_add); // conversion 8% → 12% FWHM
	    fwhm_compute = 10+i/235.5;
	    final_tree.Fill();
	  }
	}
      }
    }
  }

      

  final_tree.Write();
  file_create->Close();
}




int main(int argc, char const *argv[]){
  draw_pb_om();
  //compute_best_fwhm_spectra();
  return 0;
}
