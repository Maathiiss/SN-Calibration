#define _USE_MATH_DEFINES
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
#include <TLegend.h>
#include <TParameter.h>
#include <TSpectrum.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TSystem.h>
#include "/home/granjon/Documents/stage_radon/Stage/Mathis/sndisplay/sndisplay.cc"



using namespace std;

void compute_simulated_parameters(){
  double chi2,ndf, bcu_pic_1MeV, bcu_pic_1MeV_error, bcu_pic_0_5MeV, bcu_pic_0_5MeV_error, sigma_bcu_pic_1MeV, sigma_vis_pic_1MeV, vis_pic_1MeV,vis_pic_1MeV_error, vis_pic_0_5MeV, vis_pic_0_5MeV_error, bc_pic_1MeV, bc_pic_1MeV_error, bc_pic_0_5MeV, bc_pic_0_5MeV_error, u_pic_1MeV, u_pic_1MeV_error, u_pic_0_5MeV, u_pic_0_5MeV_error, time_bi;
  int om_number,nb_entries;
  TFile *file_create = new TFile("simu/simu_mean_extracted.root", "RECREATE");
  TTree final_tree("final_tree","");
  final_tree.Branch("Chi2",&chi2);
  final_tree.Branch("ndf",&ndf);
  final_tree.Branch("om_number",&om_number);
  final_tree.Branch("bcu_pic_1MeV",&bcu_pic_1MeV);
  final_tree.Branch("sigma_bcu_pic_1MeV",&sigma_bcu_pic_1MeV);
  final_tree.Branch("bcu_pic_1MeV_error",&bcu_pic_1MeV_error);
  final_tree.Branch("bcu_pic_0_5MeV",&bcu_pic_0_5MeV);
  final_tree.Branch("bcu_pic_0_5MeV_error",&bcu_pic_0_5MeV_error);
  
  final_tree.Branch("vis_pic_1MeV",&vis_pic_1MeV);
  final_tree.Branch("sigma_vis_pic_1MeV",&sigma_vis_pic_1MeV);
  final_tree.Branch("vis_pic_1MeV_error",&vis_pic_1MeV_error);
  final_tree.Branch("vis_pic_0_5MeV",&vis_pic_0_5MeV);
  final_tree.Branch("vis_pic_0_5MeV_error",&vis_pic_0_5MeV_error);
  
  final_tree.Branch("bc_pic_1MeV",&bc_pic_1MeV);
  final_tree.Branch("bc_pic_1MeV_error",&bc_pic_1MeV_error);
  final_tree.Branch("bc_pic_0_5MeV",&bc_pic_0_5MeV);
  final_tree.Branch("bc_pic_0_5MeV_error",&bc_pic_0_5MeV_error);
  
  final_tree.Branch("u_pic_1MeV",&u_pic_1MeV);
  final_tree.Branch("u_pic_1MeV_error",&u_pic_1MeV_error);
  final_tree.Branch("u_pic_0_5MeV",&u_pic_0_5MeV);
  final_tree.Branch("u_pic_0_5MeV_error",&u_pic_0_5MeV_error);
  
  final_tree.Branch("nb_entries",&nb_entries);
  final_tree.Branch("time_bi",&time_bi);

  TFile *file = new TFile("simu/simu_flat_70M.root", "READ");
  TTree* Result_tree = (TTree*)file->Get("Event");
  vector<int>* num_om = nullptr;
  vector<int>* number_of_kinks_per_track = nullptr;
  int number_of_electrons;
  vector<double>* evis = nullptr;
  vector<double>* evis_bcu = nullptr;
  vector<double>* evis_bc = nullptr;
  vector<double>* evis_u = nullptr;
  vector<double>* vertex_3D_start_y = nullptr;
  vector<double>* vertex_3D_start_z = nullptr;
  vector<double>* ellipse_source = nullptr;

  gROOT->cd();
  Result_tree->SetBranchStatus("*",0);
  Result_tree->SetBranchStatus("num_om",1);
  Result_tree->SetBranchAddress("num_om", &num_om);
  Result_tree->SetBranchStatus("number_of_kinks_per_track",1);
  Result_tree->SetBranchAddress("number_of_kinks_per_track", &number_of_kinks_per_track);
  Result_tree->SetBranchStatus("number_of_electrons",1);
  Result_tree->SetBranchAddress("number_of_electrons", &number_of_electrons);
  Result_tree->SetBranchStatus("evis",1);
  Result_tree->SetBranchAddress("evis", &evis);
  Result_tree->SetBranchStatus("evis_bcu",1);
  Result_tree->SetBranchAddress("evis_bcu", &evis_bcu);
  Result_tree->SetBranchStatus("evis_bc",1);
  Result_tree->SetBranchAddress("evis_bc", &evis_bc);
  Result_tree->SetBranchStatus("evis_u",1);
  Result_tree->SetBranchAddress("evis_u", &evis_u);
  Result_tree->SetBranchStatus("vertex_3D_start_y",1);
  Result_tree->SetBranchAddress("vertex_3D_start_y", &vertex_3D_start_y);
  Result_tree->SetBranchStatus("vertex_3D_start_z",1);
  Result_tree->SetBranchAddress("vertex_3D_start_z", &vertex_3D_start_z);
  Result_tree->SetBranchStatus("ellipse_source",1);
  Result_tree->SetBranchAddress("ellipse_source", &ellipse_source);
  


  std::array<TH1D*, 712> histo_evis;
  std::array<TH1D*, 712> histo_evis_bcu;
  std::array<TH1D*, 712> histo_evis_bc;
  std::array<TH1D*, 712> histo_evis_u;

  int nb_bins = 300;
  int max_charge = 2000;
  int min_charge = 0;
  for(int i=0; i<712; i++){
    if(i<520){
    histo_evis[i] = new TH1D(Form("om_%d",i),Form("om_%d",i),nb_bins,min_charge,max_charge);
    histo_evis_bcu[i] = new TH1D(Form("om_bcu_%d",i),Form("om_bcu_%d",i),nb_bins,min_charge,max_charge);
    histo_evis_bc[i] = new TH1D(Form("om_bc_%d",i),Form("om_bc_%d",i),nb_bins,min_charge,max_charge);
    histo_evis_u[i] = new TH1D(Form("om_u_%d",i),Form("om_u_%d",i),nb_bins,min_charge,max_charge);
    }
    else{
      histo_evis[i] = new TH1D(Form("om_%d",i),Form("om_%d",i),nb_bins/2,min_charge,max_charge);
      histo_evis_bcu[i] = new TH1D(Form("om_bcu_%d",i),Form("om_bcu_%d",i),nb_bins/2,min_charge,max_charge);
      histo_evis_bc[i] = new TH1D(Form("om_bc_%d",i),Form("om_bc_%d",i),nb_bins/2,min_charge,max_charge);
      histo_evis_u[i] = new TH1D(Form("om_u_%d",i),Form("om_u_%d",i),nb_bins/2,min_charge,max_charge);
    }
  }
  
  for(int entry=0; entry<Result_tree->GetEntries();entry++){
    Result_tree->GetEntry(entry);
    for(int k=0; k<number_of_electrons; k++){
      if(number_of_kinks_per_track->at(k)==0 && ellipse_source->at(k)<1){	
	  //if((abs(vertex_3D_start_y->at(k))<500 && vertex_3D_start_z->at(k)<-300)){
	  histo_evis[num_om->at(k)]->Fill(evis->at(k)*1000);
	  histo_evis_bcu[num_om->at(k)]->Fill(evis_bcu->at(k)*1000);
	  histo_evis_bc[num_om->at(k)]->Fill(evis_bc->at(k)*1000);
	  histo_evis_u[num_om->at(k)]->Fill(evis_u->at(k)*1000);
      }
    } 
  }

  TF1* f_MultipleGaus = new TF1 ("f_MultipleGaus","[0]*(7.11*TMath::Gaus(x[0],[1],[2]) + 1.84*TMath::Gaus(x[0],[1]*(1047.8/975.7),[2]*sqrt((1047.8/975.7))) + 0.44*TMath::Gaus(x[0],([1]*(1059.8/975.7)),[2]*sqrt((1059.8/975.7))))", 0, 2000);
  TF1* f_MultipleGaus_pic1 = new TF1 ("f_MultipleGaus_pic1","[0]*(1.54*TMath::Gaus(x[0],[1],[2]) + 0.44*TMath::Gaus(x[0],[1]*(553.8/481.7),[2]*sqrt((553.8/481.7))) + 0.11*TMath::Gaus(x[0],([1]*(565.8 /481.7)),[2]*sqrt((565.8 /481.7))))", 0, 2000);
  file_create->cd();
  for(int i=0; i<712; i++){
    //fit                                                                                               
    //cout<<"hist"<<i<<endl;                                                                            
    //cout<<histograms[i][0]->GetEntries()<<endl;
    nb_entries = histo_evis[i]->GetEntries();
    if(histo_evis[i]->GetEntries()>500 && histo_evis[i]->GetMean()>0){
      TFitResultPtr fitResult;      
      for (int j = 0; j < 10; j++) {
        if(j!=0){//if not first time
	  // f_MultipleGaus->SetParLimits(2, 700, 8000);
	  // f_MultipleGaus_pic1->SetParLimits(2, 300, 700);

	  //f_MultipleGaus->SetParLimits(1,histograms[i]->GetMaximumBin()*(max_charge/nb_bins)-histograms[i]->GetMaximumBin()*(max_charge/nb_bins)/5,histograms[i]->GetMaximumBin()*(max_charge/nb_bins)+histograms[i]->GetMaximumBin()*(max_charge/nb_bins)/5);
	  TSpectrum spectrum;
          if (spectrum.Search(histo_evis[i], 10, "", 0.05) == 2){
            Double_t *xPeaks = spectrum.GetPositionX();
            f_MultipleGaus->SetParameters(250,xPeaks[0],xPeaks[0]/10);
	    //f_MultipleGaus->SetParLimits(1,xPeaks[0]-xPeaks[0]/5,xPeaks[0]+xPeaks[0]/5);
	    f_MultipleGaus->SetRange(f_MultipleGaus->GetParameter(1)-1*f_MultipleGaus->GetParameter(2), f_MultipleGaus->GetParameter(1)+3*f_MultipleGaus->GetParameter(2));
	    
	    f_MultipleGaus_pic1->SetParameters(25,xPeaks[1],xPeaks[1]/10);
	    f_MultipleGaus_pic1->SetParLimits(1,xPeaks[1]-xPeaks[1]/5,xPeaks[1]+xPeaks[1]/5);
	    //f_MultipleGaus_pic1->SetParLimits(2,0,10);
	    f_MultipleGaus_pic1->SetRange(f_MultipleGaus_pic1->GetParameter(1)-3*f_MultipleGaus_pic1->GetParameter(2), f_MultipleGaus_pic1->GetParameter(1)+2*f_MultipleGaus_pic1->GetParameter(2));	    
	    //cout<<"PEAK "<<xPeaks[0]<<" "<<xPeaks[1]<<endl;
	  }
	  else{
	    Double_t *xPeaks = spectrum.GetPositionX();
            f_MultipleGaus->SetParameters(250,xPeaks[0],xPeaks[0]/10);
	    //f_MultipleGaus->SetParLimits(1,xPeaks[0]-xPeaks[0]/5,xPeaks[0]+xPeaks[0]/5);
	    f_MultipleGaus->SetRange(f_MultipleGaus->GetParameter(1)-1*f_MultipleGaus->GetParameter(2), f_MultipleGaus->GetParameter(1)+3*f_MultipleGaus->GetParameter(2));

            f_MultipleGaus_pic1->SetParameters(250,xPeaks[1],xPeaks[1]/10);
	    //f_MultipleGaus_pic1->SetParLimits(1,xPeaks[1]-xPeaks[1]/5,xPeaks[1]+xPeaks[1]/5);
	    f_MultipleGaus_pic1->SetRange(f_MultipleGaus_pic1->GetParameter(1)-1*f_MultipleGaus_pic1->GetParameter(2), f_MultipleGaus_pic1->GetParameter(1)+3*f_MultipleGaus_pic1->GetParameter(2));
	  }
	}
	else{//if first time	  
	  TSpectrum spectrum;                                                                
	  if (spectrum.Search(histo_evis[i], 10, "", 0.05) == 2){
	    Double_t *xPeaks = spectrum.GetPositionX();
	    f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);
	    f_MultipleGaus_pic1->SetParameters(25,xPeaks[1],xPeaks[1]/10);		    
	  }
	  else{
	    Double_t *xPeaks = spectrum.GetPositionX();
            f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);
 	    cout<<"histo "<<i<<"has problem"<<endl;
	  }
	}
        fitResult = histo_evis[i]->Fit(f_MultipleGaus, "RQ0");
	histo_evis[i]->Fit(f_MultipleGaus_pic1, "RQ0+");
	histo_evis_bcu[i]->Fit(f_MultipleGaus_pic1, "RQ0+");
	histo_evis_bc[i]->Fit(f_MultipleGaus_pic1, "RQ0+");
	histo_evis_u[i]->Fit(f_MultipleGaus_pic1, "RQ0+");

	histo_evis[i]->Fit(f_MultipleGaus, "RQ0+");	
	histo_evis_bcu[i]->Fit(f_MultipleGaus, "RQ0+");
	histo_evis_bc[i]->Fit(f_MultipleGaus, "RQ0+");
	histo_evis_u[i]->Fit(f_MultipleGaus, "RQ0+");
      }
      ndf = nb_bins*(f_MultipleGaus->GetParameter(1)+3*f_MultipleGaus->GetParameter(2)-(f_MultipleGaus->GetParameter(1)-1*f_MultipleGaus->GetParameter(2)))/max_charge; //500 bins de 0 a 4 MeV
    
      TCanvas* canvas = new TCanvas(Form("om_%d",i), Form("om_%d",i), 800, 600);
      if(fitResult == 1){
	//cout<<"om"<<i<<endl;	
        // canvas->cd();
        // gStyle->SetOptFit(1111);
        // histo_evis[i]->Draw();
	//canvas->SaveAs(Form("/home/granjon/Bi/histo/histo_%d.png",i));                     
	// canvas->Write();
        // canvas->Close();
	om_number=i;
        chi2 = 10000;
	bcu_pic_1MeV = 0;
	bcu_pic_1MeV_error = 0;
	bcu_pic_0_5MeV = 0;
	bcu_pic_0_5MeV_error=0;
	bc_pic_1MeV = 0;
	bc_pic_1MeV_error = 0;
	bc_pic_0_5MeV = 0;
	bc_pic_0_5MeV_error=0;
	u_pic_1MeV = 0;
	u_pic_1MeV_error = 0;
	u_pic_0_5MeV = 0;
	u_pic_0_5MeV_error=0;
	sigma_vis_pic_1MeV=0;
	sigma_bcu_pic_1MeV=0;
	vis_pic_1MeV = 0;
        vis_pic_1MeV_error = 0;
	vis_pic_0_5MeV = 0;
	vis_pic_0_5MeV_error=0;
        time_bi = time_bi;
	ndf = 1;
        final_tree.Fill();
        //cout<<"hist "<<i<<"pb fit"<<endl;                                                             
      }
      else{
	//cout<<" om bon"<<i<<endl;
        //cout<<"histo"<<i<<"chi2 "<<f_MultipleGaus->GetChisquare()<<endl;                              
	// canvas->cd();
	// gStyle->SetOptFit(1111);
	// histo_evis[i]->Draw();
	// f_MultipleGaus->Draw("same");
	// f_MultipleGaus_pic1->Draw("same");
	// gSystem->mkdir("e_vis_u");
	// canvas->SaveAs(Form("e_vis_u/histo_%d.png",i));             
	// canvas->Write();
	// canvas->Close();
	om_number=i;
	chi2 = f_MultipleGaus->GetChisquare();
	
	TF1*func = histo_evis_bcu[i]->GetFunction("f_MultipleGaus_pic1");
	TF1*func_pic = histo_evis_bcu[i]->GetFunction("f_MultipleGaus");
	if(func && func_pic){
	  TCanvas* canvas = new TCanvas(Form("om_bcu_%d",i), Form("om_bcu_%d",i), 800, 600);
	  canvas->cd();
	  gStyle->SetOptFit(1111111);
	  histo_evis_bcu[i]->Draw();
	  func->Draw("same");
	  func_pic->Draw("same");
	  canvas->Write();
	  canvas->Close();
	  bcu_pic_1MeV = histo_evis_bcu[i]->GetFunction("f_MultipleGaus")->GetParameter(1);
	  bcu_pic_1MeV_error = histo_evis_bcu[i]->GetFunction("f_MultipleGaus")->GetParError(1);
	  bcu_pic_0_5MeV =  histo_evis_bcu[i]->GetFunction("f_MultipleGaus_pic1")->GetParameter(1);
	  bcu_pic_0_5MeV_error =  histo_evis_bcu[i]->GetFunction("f_MultipleGaus_pic1")->GetParError(1);
	  sigma_bcu_pic_1MeV = histo_evis_bcu[i]->GetFunction("f_MultipleGaus")->GetParameter(2); 
	}
	TF1*func_bc = histo_evis_bc[i]->GetFunction("f_MultipleGaus_pic1");
	TF1*func_bc_pic = histo_evis_bc[i]->GetFunction("f_MultipleGaus");
	if(func_bc && func_bc_pic){
	  TCanvas* canvas = new TCanvas(Form("om_bc_%d",i), Form("om_bc_%d",i), 800, 600);
	  canvas->cd();
	  gStyle->SetOptFit(1111111);
	  histo_evis_bc[i]->Draw();
	  func_bc->Draw("same");
	  func_bc_pic->Draw("same");
	  canvas->Write();
	  canvas->Close();	  
	  bc_pic_1MeV = histo_evis_bc[i]->GetFunction("f_MultipleGaus")->GetParameter(1);
	  bc_pic_1MeV_error = histo_evis_bc[i]->GetFunction("f_MultipleGaus")->GetParError(1);
	  bc_pic_0_5MeV =  histo_evis_bc[i]->GetFunction("f_MultipleGaus_pic1")->GetParameter(1);
	  bc_pic_0_5MeV_error =  histo_evis_bc[i]->GetFunction("f_MultipleGaus_pic1")->GetParError(1);
	}
	TF1*func_u = histo_evis_u[i]->GetFunction("f_MultipleGaus_pic1");
	TF1*func_u_pic = histo_evis_u[i]->GetFunction("f_MultipleGaus");
	if(func_u && func_u_pic){
	  TCanvas* canvas = new TCanvas(Form("om_u_%d",i), Form("om_u_%d",i), 800, 600);
	  canvas->cd();
	  gStyle->SetOptFit(1111111);
	  histo_evis_u[i]->Draw();
	  func_u->Draw("same");
	  func_u_pic->Draw("same");
	  canvas->Write();
	  canvas->Close();

	  u_pic_1MeV = histo_evis_u[i]->GetFunction("f_MultipleGaus")->GetParameter(1);
	  u_pic_1MeV_error = histo_evis_u[i]->GetFunction("f_MultipleGaus")->GetParError(1);
	  u_pic_0_5MeV =  histo_evis_u[i]->GetFunction("f_MultipleGaus_pic1")->GetParameter(1);
	  u_pic_0_5MeV_error =  histo_evis_u[i]->GetFunction("f_MultipleGaus_pic1")->GetParError(1);
	}
	TF1* func_vis = histo_evis[i]->GetFunction("f_MultipleGaus_pic1");
	TF1* func_vis_pic = histo_evis[i]->GetFunction("f_MultipleGaus");
	if(func_vis && func_vis_pic){
	  TCanvas* canvas = new TCanvas(Form("om_vis_%d",i), Form("om_vis_%d",i), 800, 600);
	  canvas->cd();
	  gStyle->SetOptFit(1111111);
	  histo_evis[i]->Draw();
	  func_vis->Draw("same");
	  func_vis_pic->Draw("same");
	  canvas->Write();
	  canvas->Close();
	  sigma_vis_pic_1MeV = histo_evis[i]->GetFunction("f_MultipleGaus")->GetParameter(2);
	  vis_pic_1MeV = histo_evis[i]->GetFunction("f_MultipleGaus")->GetParameter(1);
	  vis_pic_1MeV_error = histo_evis[i]->GetFunction("f_MultipleGaus")->GetParError(1);
	  vis_pic_0_5MeV =  histo_evis[i]->GetFunction("f_MultipleGaus_pic1")->GetParameter(1);
	  vis_pic_0_5MeV_error =  histo_evis[i]->GetFunction("f_MultipleGaus_pic1")->GetParError(1);
	}
	time_bi = time_bi;
	final_tree.Fill();
      }
    }
    else{
      TCanvas* canvas = new TCanvas(Form("om_vis_%d",i), Form("om_vis_%d",i), 800, 600);
      canvas->cd();
      gStyle->SetOptFit(1111111);
      histo_evis[i]->Draw();
      canvas->Write();
      canvas->Close();
	  
      f_MultipleGaus->SetParameters(0,0,0);
      om_number=i;
      chi2 = 2000;
      ndf = 1;
      time_bi = time_bi;
      bcu_pic_1MeV = 0;
      sigma_bcu_pic_1MeV = 0;
      bcu_pic_1MeV_error = 0;
      bcu_pic_0_5MeV = 0;
      bcu_pic_0_5MeV_error=0;
      bc_pic_1MeV = 0;
      bc_pic_1MeV_error = 0;
      bc_pic_0_5MeV = 0;
      bc_pic_0_5MeV_error=0;
      u_pic_1MeV = 0;
      u_pic_1MeV_error = 0;
      u_pic_0_5MeV = 0;
      u_pic_0_5MeV_error=0;
      sigma_vis_pic_1MeV=0;
      vis_pic_1MeV = 0;
      vis_pic_1MeV_error = 0;
      vis_pic_0_5MeV = 0;
      vis_pic_0_5MeV_error=0;
      final_tree.Fill();      
    }


    	histo_evis[i]->Delete();
  }
  final_tree.Write();
  file_create->Close();

}



int main(int argc, char const *argv[]){
  std::cout << "" << '\n';
  std::cout << "START OF THE CALORIMETER FIT" << '\n';
  std::cout << "" << '\n';
  compute_simulated_parameters();
  
  return 0;
}









