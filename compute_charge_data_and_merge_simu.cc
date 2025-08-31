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
//summer 2023 cara
//int vecteur_Bi[17] = {1090,1097,1105,1112,1119,1126,1133,1149,1156,1163,1170,1175,1180,1190,1195,1200,1205};
//int vecteur_Bi[15] = {1105,1112,1119,1126,1133,1149,1156,1163,1170,1175,1180,1190,1195,1200,1205};
//int vecteur_Bi[1] = {1112};

//March 2023 cara
//int vecteur_Bi[23] ={1274,1275,1276,1277,1281,1286,1287,1292,/*1297 only 1 min,*/1303,1313,1320,/*1326too short*/ /*1327 root is truncated*/1337,1351,1359,1367,/*1374 root has problem*/1381,1389,1396,1403,1412,1413,1414,1415/*1416 too short + source move*/};
int vecteur_Bi[1] = {1000};

void fit_Bi_energy(int run_number){
  double chi2,chi2_0_5_MeV, Q_1_MeV, Q_1_MeV_error,time_bi,ndf, Q_0_5_MeV, Q_0_5_MeV_error, Center_Q_1_MeV, Center_Q_1_MeV_error, Center_Q_0_5_MeV, Center_Q_0_5_MeV_error, Edge_Q_1_MeV, Edge_Q_1_MeV_error, Edge_Q_0_5_MeV, Edge_Q_0_5_MeV_error, sigma_1_MeV, sigma_1_MeV_error,sigma_0_5_MeV, sigma_0_5_MeV_error;
  int om_number,nb_entries;
  TFile *file_create = new TFile(Form("data_runs/data_mean_%d_extracted.root",run_number), "RECREATE");
  TTree final_tree("final_tree","");
  final_tree.Branch("chi2_0_5_MeV",&chi2_0_5_MeV);
  final_tree.Branch("chi2_1_MeV",&chi2);
  final_tree.Branch("om_number",&om_number);
  final_tree.Branch("run_number",&run_number);
  final_tree.Branch("Q_1_MeV",&Q_1_MeV);
  final_tree.Branch("Q_1_MeV_error",&Q_1_MeV_error);
  final_tree.Branch("sigma_1_MeV",&sigma_1_MeV);
  final_tree.Branch("sigma_1_MeV_error",&sigma_1_MeV_error);
  final_tree.Branch("Q_0_5_MeV",&Q_0_5_MeV);
  final_tree.Branch("Q_0_5_MeV_error",&Q_0_5_MeV_error);
  final_tree.Branch("sigma_0_5_MeV",&sigma_0_5_MeV);
  final_tree.Branch("sigma_0_5_MeV_error",&sigma_0_5_MeV_error);
  final_tree.Branch("Center_Q_1_MeV",&Center_Q_1_MeV);
  final_tree.Branch("Center_Q_1_MeV_error",&Center_Q_1_MeV_error);
  final_tree.Branch("Center_Q_0_5_MeV",&Center_Q_0_5_MeV);
  final_tree.Branch("Center_Q_0_5_MeV_error",&Center_Q_0_5_MeV_error);
  final_tree.Branch("Edge_Q_1_MeV",&Edge_Q_1_MeV);
  final_tree.Branch("Edge_Q_1_MeV_error",&Edge_Q_1_MeV_error);
  final_tree.Branch("Edge_Q_0_5_MeV",&Edge_Q_0_5_MeV);
  final_tree.Branch("Edge_Q_0_5_MeV_error",&Edge_Q_0_5_MeV_error);
  final_tree.Branch("time_bi",&time_bi);
  final_tree.Branch("nb_entries",&nb_entries);
  final_tree.Branch("ndf",&ndf);


  TFile *file = new TFile(Form("data_runs/after_gain_correction_data_run_%d.root",run_number), "READ");
  TTree* Result_tree = (TTree*)file->Get("Event");
  vector<int>* num_om = nullptr;
  vector<int>* number_of_kinks_per_track = nullptr;
  vector<double>* charge_elec = nullptr;
  vector<double>* delta_y_calo_center = nullptr;
  vector<double>* delta_z_calo_center = nullptr;
  vector<double>* vertex_3D_start_y = nullptr;
  vector<double>* vertex_3D_start_z = nullptr;
  vector<double>* ellipse_source = nullptr;

  
  int number_of_electrons;  
  gROOT->cd();
  Result_tree->SetBranchStatus("*",0);
  Result_tree->SetBranchStatus("*",0);
  Result_tree->SetBranchStatus("num_om",1);
  Result_tree->SetBranchAddress("num_om", &num_om);
  Result_tree->SetBranchStatus("number_of_kinks_per_track",1);
  Result_tree->SetBranchAddress("number_of_kinks_per_track", &number_of_kinks_per_track);
  Result_tree->SetBranchStatus("number_of_electrons",1);
  Result_tree->SetBranchAddress("number_of_electrons", &number_of_electrons);
  Result_tree->SetBranchStatus("charge_elec",1);
  Result_tree->SetBranchAddress("charge_elec", &charge_elec);
  Result_tree->SetBranchStatus("delta_y_calo_center",1);
  Result_tree->SetBranchAddress("delta_y_calo_center", &delta_y_calo_center);
  Result_tree->SetBranchStatus("delta_z_calo_center",1);
  Result_tree->SetBranchAddress("delta_z_calo_center", &delta_z_calo_center);
  Result_tree->SetBranchStatus("vertex_3D_start_y",1);
  Result_tree->SetBranchAddress("vertex_3D_start_y", &vertex_3D_start_y);
  Result_tree->SetBranchStatus("vertex_3D_start_z",1);
  Result_tree->SetBranchAddress("vertex_3D_start_z", &vertex_3D_start_z);
  Result_tree->SetBranchStatus("ellipse_source",1);
  Result_tree->SetBranchAddress("ellipse_source", &ellipse_source);


  std::array<TH1D*, 712> histograms;
  std::array<TH1D*, 712> histograms_OM_center;
  std::array<TH1D*, 712> histograms_OM_edge;
  int nb_bins = 500;
  int max_charge = 60000;
  int min_charge = 0;
  for(int i=0; i<712; i++){
    if(i<520){
    histograms[i] = new TH1D(Form("om_%d",i),Form("om_%d",i),nb_bins,min_charge,max_charge);
    histograms_OM_center[i]=new TH1D(Form("om_center_%d",i),Form("om_center_%d",i),nb_bins,min_charge,max_charge);
    histograms_OM_edge[i]=new TH1D(Form("om_edge_%d",i),Form("om_edge_%d",i),nb_bins,min_charge,max_charge);
    }
    else{
      histograms[i] = new TH1D(Form("om_%d",i),Form("om_%d",i),nb_bins/2,min_charge,max_charge);
      histograms_OM_center[i]=new TH1D(Form("om_center_%d",i),Form("om_center_%d",i),nb_bins/2,min_charge,max_charge);
      histograms_OM_edge[i]=new TH1D(Form("om_edge_%d",i),Form("om_edge_%d",i),nb_bins,min_charge,max_charge);
    }
  }
  
  for(int entry=0; entry<Result_tree->GetEntries();entry++){
    Result_tree->GetEntry(entry);
    for(int k=0; k<number_of_electrons; k++){
      if(number_of_kinks_per_track->at(k)==0 && ellipse_source->at(k)<1){
	//if(!(abs(vertex_3D_start_y->at(k))<500 && vertex_3D_start_z->at(k)<-300)){
	//if((abs(vertex_3D_start_y->at(k))<500 && vertex_3D_start_z->at(k)<-300)){
	  histograms[num_om->at(k)]->Fill(charge_elec->at(k));
	  if(delta_y_calo_center->at(k)<100 && delta_z_calo_center->at(k)<100){
	    histograms_OM_center[num_om->at(k)]->Fill(charge_elec->at(k));
	  }
	  else{
	    histograms_OM_edge[num_om->at(k)]->Fill(charge_elec->at(k));
	  }
	  //}
      }
    }
  }
  
  std::ofstream log_file(Form("log/log_file_%d.txt",run_number));  
  TF1* f_MultipleGaus = new TF1 ("f_MultipleGaus","[0]*(7.08*TMath::Gaus(x[0],[1],[2]) + 1.84*TMath::Gaus(x[0],[1]*(1047.8/975.7),[2]*sqrt((1047.8/975.7))) + 0.44*TMath::Gaus(x[0],([1]*(1059.8/975.7)),[2]*sqrt((1059.8/975.7))))", 0, 60000);
  TF1* f_MultipleGaus_pic1 = new TF1 ("f_MultipleGaus_pic1","[0]*(1.54*TMath::Gaus(x[0],[1],[2]) + 0.44*TMath::Gaus(x[0],[1]*(553.8/481.7),[2]*sqrt((553.8/481.7))) + 0.11*TMath::Gaus(x[0],([1]*(565.8 /481.7)),[2]*sqrt((565.8 /481.7))))", 0, 60000);
  file_create->cd();
  for(int i=0; i<712; i++){
    //fit                                                                                               
    //cout<<"hist"<<i<<endl;                                                                            
    //cout<<histograms[i][0]->GetEntries()<<endl;
    nb_entries = histograms[i]->GetEntries();
    if(histograms[i]->GetEntries()>400 && histograms[i]->GetMean()>0){
      TFitResultPtr fitResult;      
      for (int j = 0; j < 10; j++) {
        if(j!=0){
	  f_MultipleGaus->SetParLimits(2, 800, 8000);
	  //f_MultipleGaus->SetParLimits(1,histograms[i]->GetMaximumBin()*(max_charge/nb_bins)-histograms[i]->GetMaximumBin()*(max_charge/nb_bins)/5,histograms[i]->GetMaximumBin()*(max_charge/nb_bins)+histograms[i]->GetMaximumBin()*(max_charge/nb_bins)/5);
	  TSpectrum spectrum;
          if (spectrum.Search(histograms[i], 7, "", 0.11) == 2){
            Double_t *xPeaks = spectrum.GetPositionX();
            f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);
	    f_MultipleGaus->SetParLimits(1,xPeaks[0]-xPeaks[0]/5,xPeaks[0]+xPeaks[0]/5);
	    f_MultipleGaus->SetRange(f_MultipleGaus->GetParameter(1)-1*f_MultipleGaus->GetParameter(2), f_MultipleGaus->GetParameter(1)+2*f_MultipleGaus->GetParameter(2));
	    
	    f_MultipleGaus_pic1->SetParameters(7,xPeaks[1],xPeaks[1]/10);	    
	    f_MultipleGaus_pic1->SetRange(f_MultipleGaus_pic1->GetParameter(1)-2*f_MultipleGaus_pic1->GetParameter(2), f_MultipleGaus_pic1->GetParameter(1)+3*f_MultipleGaus_pic1->GetParameter(2));	    
	    //f_MultipleGaus_pic1->SetParLimits(0,0,40);
	    f_MultipleGaus_pic1->SetParLimits(1,xPeaks[1]-xPeaks[1]/5,xPeaks[1]+xPeaks[1]/5);
	    //f_MultipleGaus_pic1->SetParLimits(2,1500,5000);

	    //cout<<"PEAK "<<xPeaks[0]<<" "<<xPeaks[1]<<endl;
	  }
	  else{
	     int sigma=7;
            if(spectrum.Search(histograms[i], sigma, "", 0.11)==1 && j==9){
              log_file<<" om "<<i<<" has only one peak "<<endl;
            }
            while(spectrum.Search(histograms[i], sigma, "", 0.11) > 2){
              sigma++;
            }
	    Double_t *xPeaks = spectrum.GetPositionX();
	     if (xPeaks[1] > xPeaks[0]) {
              std::swap(xPeaks[0], xPeaks[1]);
            }

            f_MultipleGaus->SetParameters(5,xPeaks[0],xPeaks[0]/10);
	    f_MultipleGaus->SetParLimits(1,xPeaks[0]-xPeaks[0]/5,xPeaks[0]+xPeaks[0]/5);
	    f_MultipleGaus->SetRange(f_MultipleGaus->GetParameter(1)-1*f_MultipleGaus->GetParameter(2), f_MultipleGaus->GetParameter(1)+2*f_MultipleGaus->GetParameter(2));

            f_MultipleGaus_pic1->SetParameters(25,xPeaks[1],xPeaks[1]/10);
	    f_MultipleGaus_pic1->SetParLimits(1,xPeaks[1]-xPeaks[1]/5,xPeaks[1]+xPeaks[1]/5);
	    f_MultipleGaus_pic1->SetRange(f_MultipleGaus_pic1->GetParameter(1)-2*f_MultipleGaus_pic1->GetParameter(2), f_MultipleGaus_pic1->GetParameter(1)+3*f_MultipleGaus_pic1->GetParameter(2));
	  }
	}
	else{//if first time	  
	  TSpectrum spectrum;                                                                
	  if (spectrum.Search(histograms[i], 7, "", 0.11) == 2){
	    Double_t *xPeaks = spectrum.GetPositionX();
	    f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);
	    f_MultipleGaus_pic1->SetParameters(250,xPeaks[1],xPeaks[1]/10);		    
	  }
	  else{
	    Double_t *xPeaks = spectrum.GetPositionX();
            int sigma=7;
            while(spectrum.Search(histograms[i], sigma, "", 0.11) > 2){
              sigma++;
            }
            if (xPeaks[1] > xPeaks[0]) {
              std::swap(xPeaks[0], xPeaks[1]);
            }

            f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);

	  }
	}
        fitResult = histograms[i]->Fit(f_MultipleGaus, "RQ0");
	histograms[i]->Fit(f_MultipleGaus_pic1, "RQ0+");
	histograms_OM_center[i]->Fit(f_MultipleGaus_pic1, "RQ0+");
	histograms_OM_center[i]->Fit(f_MultipleGaus, "RQ0+");
	histograms_OM_edge[i]->Fit(f_MultipleGaus, "RQ0+");
	histograms_OM_edge[i]->Fit(f_MultipleGaus_pic1, "RQ0+");
      }
      ndf = nb_bins*(f_MultipleGaus->GetParameter(1)+3*f_MultipleGaus->GetParameter(2)-(f_MultipleGaus->GetParameter(1)-1*f_MultipleGaus->GetParameter(2)))/max_charge; //500 bins de 0 a 4 MeV
      TCanvas* canvas = new TCanvas(Form("om_%d",i), Form("om_%d",i), 800, 600);
      if(fitResult == 1){
	//cout<<"om"<<i<<endl;	
        canvas->cd();
        gStyle->SetOptFit(1111);
        histograms[i]->Draw();
	//canvas->SaveAs(Form("/home/granjon/Bi/histo/histo_%d.png",i));                     
	canvas->Write();
        canvas->Close();
	om_number=i;
        chi2 = 10000;
	chi2_0_5_MeV=10000;
	Q_1_MeV = 0;
	Q_0_5_MeV = 0;
	Q_0_5_MeV_error = 0;
	Q_1_MeV_error = 0;
	Center_Q_1_MeV = 0;
	Center_Q_0_5_MeV = 0;
	Center_Q_0_5_MeV_error = 0;
	Center_Q_1_MeV_error = 0;
	run_number = run_number;
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
	// histograms[i]->Draw();
	// f_MultipleGaus->Draw("same");
	// f_MultipleGaus_pic1->Draw("same");
	//histograms[i]->ShowPeaks(5, "", 0.2);
  	// gSystem->mkdir(Form("histo_data/run_%d",run_number));
	// canvas->SaveAs(Form("histo_data/run_%d/histo_%d.png",run_number,i));             
	// canvas->Write();
	// canvas->Close();
	om_number=i;
	//cout<<histograms[i]->GetFunction("f_MultipleGaus_pic1")->GetParameter(0)<< " "<<histograms[i]->GetFunction("f_MultipleGaus_pic1")->GetParameter(1) << " "<< histograms[i]->GetFunction("f_MultipleGaus_pic1")->GetParameter(2)<<endl;
	TF1* func_ = histograms[i]->GetFunction("f_MultipleGaus");
        TF1* func_1 = histograms[i]->GetFunction("f_MultipleGaus_pic1");
        if(func_ && func_1){
	  canvas->cd();
          gStyle->SetOptFit(1111111);
	  histograms[i]->Draw();
	  func_->Draw("same");
          func_1->Draw("same");
	  canvas->Write();
          canvas->Close();	   
	  Q_1_MeV = histograms[i]->GetFunction("f_MultipleGaus")->GetParameter(1);
	  Q_0_5_MeV = histograms[i]->GetFunction("f_MultipleGaus_pic1")->GetParameter(1);
	  sigma_1_MeV = histograms[i]->GetFunction("f_MultipleGaus")->GetParameter(2);
	  sigma_1_MeV_error = histograms[i]->GetFunction("f_MultipleGaus")->GetParError(2);
	  Q_0_5_MeV_error = histograms[i]->GetFunction("f_MultipleGaus_pic1")->GetParError(1);
	  Q_1_MeV_error = histograms[i]->GetFunction("f_MultipleGaus")->GetParError(1);
	  sigma_0_5_MeV = histograms[i]->GetFunction("f_MultipleGaus_pic1")->GetParameter(2);
	  sigma_0_5_MeV_error = histograms[i]->GetFunction("f_MultipleGaus_pic1")->GetParError(2);
	  chi2 = func_->GetChisquare()/func_->GetNDF();
	  chi2_0_5_MeV = func_1->GetChisquare()/func_1->GetNDF();

	}
	TF1* func = histograms_OM_center[i]->GetFunction("f_MultipleGaus");
	TF1* func1 = histograms_OM_center[i]->GetFunction("f_MultipleGaus_pic1");
	if(func && func1){
	  TCanvas* canvas = new TCanvas(Form("om_center_%d",i), Form("om_center_%d",i), 800, 600);
	  canvas->cd();
          gStyle->SetOptFit(1111111);
	  histograms_OM_center[i]->Draw();
	  func->Draw("same");
          func1->Draw("same");
	  canvas->Write();
          canvas->Close();

	  Center_Q_1_MeV = histograms_OM_center[i]->GetFunction("f_MultipleGaus")->GetParameter(1);
	  Center_Q_0_5_MeV = histograms_OM_center[i]->GetFunction("f_MultipleGaus_pic1")->GetParameter(1);
	  Center_Q_0_5_MeV_error = histograms_OM_center[i]->GetFunction("f_MultipleGaus_pic1")->GetParError(1);
	  Center_Q_1_MeV_error = histograms_OM_center[i]->GetFunction("f_MultipleGaus")->GetParError(1);
	}
	TF1* func_e = histograms_OM_edge[i]->GetFunction("f_MultipleGaus");
	TF1* func1_e = histograms_OM_edge[i]->GetFunction("f_MultipleGaus_pic1");
	if(func_e && func1_e){
	  TCanvas* canvas = new TCanvas(Form("om_edge_%d",i), Form("om_edge_%d",i), 800, 600);
	  canvas->cd();
          gStyle->SetOptFit(1111111);
	  histograms_OM_edge[i]->Draw();
	  func_e->Draw("same");
          func1_e->Draw("same");
	  canvas->Write();
          canvas->Close();
	  
	  Edge_Q_1_MeV = histograms_OM_edge[i]->GetFunction("f_MultipleGaus")->GetParameter(1);
	  Edge_Q_0_5_MeV = histograms_OM_edge[i]->GetFunction("f_MultipleGaus_pic1")->GetParameter(1);
	  Edge_Q_0_5_MeV_error = histograms_OM_edge[i]->GetFunction("f_MultipleGaus_pic1")->GetParError(1);
	  Edge_Q_1_MeV_error = histograms_OM_edge[i]->GetFunction("f_MultipleGaus")->GetParError(1);
	}
	run_number = run_number;
	time_bi = time_bi;
	final_tree.Fill();
	//histograms[i][0]->Write();                                                           

        //cout<<"hist "<<i<<"good"<<endl;                                                               
      }
    }
    else{
      TCanvas* canvas = new TCanvas(Form("om_%d",i), Form("om_%d",i), 800, 600);
      canvas->cd();
      gStyle->SetOptFit(1111);
      histograms[i]->Draw();
      //canvas->SaveAs(Form("/home/granjon/Bi/histo/run_%d/histo_%d.png",run_number,i));
      canvas->Write();
      canvas->Close();
      f_MultipleGaus->SetParameters(0,0,0);
      om_number=i;
      chi2 = 2000;
      chi2_0_5_MeV=2000;
      ndf = 1;
      run_number = run_number;
      time_bi = time_bi;
      Q_1_MeV = 0;
      Q_0_5_MeV=0;
      Q_0_5_MeV_error = 0;
      Q_1_MeV_error = 0;
      Center_Q_1_MeV = 0;
      Center_Q_0_5_MeV = 0;
      Center_Q_0_5_MeV_error = 0;
      Center_Q_1_MeV_error = 0;

      final_tree.Fill();      
      //cout<<"hist without entry "<<i<<endl;                                                          
    }
    	histograms[i]->Delete();
  }
  final_tree.Write();
  file_create->Close();

}


void Merger_data_simu(int run_number){
  //variables for data
  TFile *f_data = TFile::Open(Form("data_runs/data_mean_%d_extracted.root",run_number));  
  TTree *t_data = (TTree*)f_data->Get("final_tree");
  int om_number_data, run_number_data, nb_entries;
  double Q_1_MeV, Q_1_MeV_error, Q_0_5_MeV, Q_0_5_MeV_error, Center_Q_1_MeV, Center_Q_1_MeV_error, Center_Q_0_5_MeV, Center_Q_0_5_MeV_error, Edge_Q_1_MeV, Edge_Q_1_MeV_error, Edge_Q_0_5_MeV, Edge_Q_0_5_MeV_error, sigma_1_MeV, sigma_1_MeV_error, sigma_0_5_MeV, sigma_0_5_MeV_error;
  t_data->SetBranchAddress("om_number", &om_number_data);
  t_data->SetBranchAddress("run_number", &run_number_data);
  t_data->SetBranchAddress("Q_1_MeV", &Q_1_MeV);
  t_data->SetBranchAddress("Q_1_MeV_error", &Q_1_MeV_error);
  t_data->SetBranchAddress("Q_0_5_MeV", &Q_0_5_MeV);
  t_data->SetBranchAddress("Q_0_5_MeV_error", &Q_0_5_MeV_error);
  t_data->SetBranchAddress("sigma_1_MeV",&sigma_1_MeV);
  t_data->SetBranchAddress("sigma_1_MeV_error",&sigma_1_MeV_error);
  t_data->SetBranchAddress("sigma_0_5_MeV",&sigma_0_5_MeV);
  t_data->SetBranchAddress("sigma_0_5_MeV_error",&sigma_0_5_MeV_error);
  t_data->SetBranchAddress("Center_Q_1_MeV", &Center_Q_1_MeV);
  t_data->SetBranchAddress("Center_Q_1_MeV_error", &Center_Q_1_MeV_error);
  t_data->SetBranchAddress("Center_Q_0_5_MeV", &Center_Q_0_5_MeV);
  t_data->SetBranchAddress("Center_Q_0_5_MeV_error", &Center_Q_0_5_MeV_error);
  t_data->SetBranchAddress("Edge_Q_1_MeV", &Edge_Q_1_MeV);
  t_data->SetBranchAddress("Edge_Q_1_MeV_error", &Edge_Q_1_MeV_error);
  t_data->SetBranchAddress("Edge_Q_0_5_MeV", &Edge_Q_0_5_MeV);
  t_data->SetBranchAddress("Edge_Q_0_5_MeV_error", &Edge_Q_0_5_MeV_error);
  t_data->SetBranchAddress("nb_entries", &nb_entries);
  
  //variables for simu
  TFile *f_simu = TFile::Open("simu/simu_flat_mean_extracted.root");
  //TFile *f_simu = TFile::Open("simu/simu_mean_extracted_without_problem_sources.root");
  TTree *t_simu = (TTree*)f_simu->Get("final_tree");
  int om_number_simu;
  double bcu_pic_1MeV, bcu_pic_1MeV_error, bcu_pic_0_5MeV, bcu_pic_0_5MeV_error, vis_pic_1MeV, vis_pic_1MeV_error, vis_pic_0_5MeV, vis_pic_0_5MeV_error, bc_pic_1MeV, bc_pic_1MeV_error, bc_pic_0_5MeV, bc_pic_0_5MeV_error, u_pic_1MeV, u_pic_1MeV_error, u_pic_0_5MeV, u_pic_0_5MeV_error;
  t_simu->SetBranchAddress("om_number", &om_number_simu);
  t_simu->SetBranchAddress("bcu_pic_1MeV", &bcu_pic_1MeV);
  t_simu->SetBranchAddress("bcu_pic_1MeV_error", &bcu_pic_1MeV_error);
  t_simu->SetBranchAddress("bcu_pic_0_5MeV", &bcu_pic_0_5MeV);
  t_simu->SetBranchAddress("bcu_pic_0_5MeV_error", &bcu_pic_0_5MeV_error);
  t_simu->SetBranchAddress("vis_pic_1MeV", &vis_pic_1MeV);
  t_simu->SetBranchAddress("vis_pic_1MeV_error", &vis_pic_1MeV_error);
  t_simu->SetBranchAddress("vis_pic_0_5MeV", &vis_pic_0_5MeV);
  t_simu->SetBranchAddress("vis_pic_0_5MeV_error", &vis_pic_0_5MeV_error);
  t_simu->SetBranchAddress("bc_pic_1MeV", &bc_pic_1MeV);
  t_simu->SetBranchAddress("bc_pic_1MeV_error", &bc_pic_1MeV_error);
  t_simu->SetBranchAddress("bc_pic_0_5MeV", &bc_pic_0_5MeV);
  t_simu->SetBranchAddress("bc_pic_0_5MeV_error", &bc_pic_0_5MeV_error);
  t_simu->SetBranchAddress("u_pic_1MeV", &u_pic_1MeV);
  t_simu->SetBranchAddress("u_pic_1MeV_error", &u_pic_1MeV_error);
  t_simu->SetBranchAddress("u_pic_0_5MeV", &u_pic_0_5MeV);
  t_simu->SetBranchAddress("u_pic_0_5MeV_error", &u_pic_0_5MeV_error);

  //Merged tree
  double a, a_error, E_vis_MC, E_vis_MC_error, E_vis_data,E_vis_data_error, corr_simu, corr_simu_error, E_vis_center_data, E_vis_center_data_error, E_vis_edge_data, E_vis_edge_data_error, fwhm, fwhm_error;
  TFile *fout = new TFile(Form("Final_merged_%d.root",run_number), "RECREATE");
  TTree *t_out = new TTree("MergedTree", "Fusion sans map (double)");
  t_out->Branch("om_number", &om_number_data);
  t_out->Branch("E_vis_MC", &E_vis_MC);
  t_out->Branch("E_vis_MC_error", &E_vis_MC_error);
  t_out->Branch("E_vis_center_data", &E_vis_center_data);
  t_out->Branch("E_vis_center_data_error", &E_vis_center_data_error);
  t_out->Branch("E_vis_edge_data", &E_vis_edge_data);
  t_out->Branch("E_vis_edge_data_error", &E_vis_edge_data_error);
  t_out->Branch("E_vis_data", &E_vis_data);
  t_out->Branch("E_vis_data_error", &E_vis_data_error);
  t_out->Branch("corr_simu", &corr_simu);
  t_out->Branch("corr_simu_error", &corr_simu_error);
  t_out->Branch("a", &a);
  t_out->Branch("a_error", &a_error);
  t_out->Branch("Q_1_MeV", &Q_1_MeV);
  t_out->Branch("Q_1_MeV_error", &Q_1_MeV_error);
  t_out->Branch("Q_0_5_MeV", &Q_0_5_MeV);
  t_out->Branch("Q_0_5_MeV_error", &Q_0_5_MeV_error);
  t_out->Branch("sigma_1_MeV",&sigma_1_MeV);
  t_out->Branch("sigma_1_MeV_error",&sigma_1_MeV_error);
  t_out->Branch("sigma_0_5_MeV",&sigma_0_5_MeV);
  t_out->Branch("sigma_0_5_MeV_error",&sigma_0_5_MeV_error);
  t_out->Branch("Center_Q_1_MeV", &Center_Q_1_MeV);
  t_out->Branch("Center_Q_1_MeV_error", &Center_Q_1_MeV_error);
  t_out->Branch("Center_Q_0_5_MeV", &Center_Q_0_5_MeV);
  t_out->Branch("Center_Q_0_5_MeV_error", &Center_Q_0_5_MeV_error);
  t_out->Branch("Edge_Q_1_MeV", &Edge_Q_1_MeV);
  t_out->Branch("Edge_Q_1_MeV_error", &Edge_Q_1_MeV_error);
  t_out->Branch("Edge_Q_0_5_MeV", &Edge_Q_0_5_MeV);
  t_out->Branch("Edge_Q_0_5_MeV_error", &Edge_Q_0_5_MeV_error);
  t_out->Branch("bcu_pic_1MeV", &bcu_pic_1MeV);
  t_out->Branch("bcu_pic_1MeV_error", &bcu_pic_1MeV_error);
  t_out->Branch("bcu_pic_0_5MeV", &bcu_pic_0_5MeV);
  t_out->Branch("bcu_pic_0_5MeV_error", &bcu_pic_0_5MeV_error);
  t_out->Branch("vis_pic_1MeV", &vis_pic_1MeV);
  t_out->Branch("vis_pic_1MeV_error", &vis_pic_1MeV_error);
  t_out->Branch("vis_pic_0_5MeV", &vis_pic_0_5MeV);
  t_out->Branch("vis_pic_0_5MeV_error", &vis_pic_0_5MeV_error);
  t_out->Branch("bc_pic_1MeV", &bc_pic_1MeV);
  t_out->Branch("bc_pic_1MeV_error", &bc_pic_1MeV_error);
  t_out->Branch("bc_pic_0_5MeV", &bc_pic_0_5MeV);
  t_out->Branch("bc_pic_0_5MeV_error", &bc_pic_0_5MeV_error);
  t_out->Branch("u_pic_1MeV", &u_pic_1MeV);
  t_out->Branch("u_pic_1MeV_error", &u_pic_1MeV_error);
  t_out->Branch("u_pic_0_5MeV", &u_pic_0_5MeV);
  t_out->Branch("u_pic_0_5MeV_error", &u_pic_0_5MeV_error);
  t_out->Branch("nb_entries", &nb_entries);
  t_out->Branch("fwhm", &fwhm);
  t_out->Branch("fwhm_error", &fwhm_error);

  std::ifstream infile("list_om_tab.txt");
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

  //Fill simulation file
  std::ofstream csv_file("simulation_calibration.csv");
  csv_file << "om_number  corr_simu\n";
  //std::ofstream csv_file("data_calibration.csv");
  //csv_file << "om_number  a \n";
  for (int i = 0; i < t_data->GetEntries(); ++i) {
    t_data->GetEntry(i);
    for (int j = 0; j < t_simu->GetEntries(); ++j) {
      t_simu->GetEntry(j);
      if(om_number_data==om_number_simu){
	if(Q_1_MeV!=0){
	  a = vis_pic_1MeV/Q_1_MeV;
	  a_error = a*sqrt(pow(vis_pic_1MeV_error/vis_pic_1MeV,2) + pow(Q_1_MeV_error/Q_1_MeV,2));
	  corr_simu = vis_pic_1MeV/bcu_pic_1MeV;
	  corr_simu_error = corr_simu*sqrt(pow(vis_pic_1MeV_error/vis_pic_1MeV,2) + pow(bcu_pic_1MeV_error/bcu_pic_1MeV,2));
	  E_vis_MC = bcu_pic_0_5MeV * corr_simu;
	  E_vis_MC_error = E_vis_MC*sqrt(pow(bcu_pic_0_5MeV_error/bcu_pic_0_5MeV,2) + pow(corr_simu_error/corr_simu,2));	  
	  E_vis_data = a * Q_0_5_MeV;
	  E_vis_data_error = E_vis_data*sqrt(pow(a_error/a,2) + pow(Q_0_5_MeV_error/Q_0_5_MeV,2));
	  //compute coeff for center of the OM
	  double a_center = vis_pic_1MeV/Center_Q_1_MeV;
	  double a_center_error = a_center*sqrt(pow(vis_pic_1MeV_error/vis_pic_1MeV,2) + pow(Center_Q_1_MeV_error/Center_Q_1_MeV,2));	 
	  E_vis_center_data = a_center * Center_Q_0_5_MeV;
	  E_vis_center_data_error = E_vis_data*sqrt(pow(a_center_error/a_center,2) + pow(Center_Q_0_5_MeV_error/Center_Q_0_5_MeV,2));
	  //compute coeff for edges of the OM
	  double a_edge = vis_pic_1MeV/Edge_Q_1_MeV;
	  double a_edge_error = a_edge*sqrt(pow(vis_pic_1MeV_error/vis_pic_1MeV,2) + pow(Edge_Q_1_MeV_error/Edge_Q_1_MeV,2));
	  E_vis_edge_data = a_edge * Edge_Q_0_5_MeV;
	  E_vis_edge_data_error = E_vis_data*sqrt(pow(a_edge_error/a_edge,2) + pow(Edge_Q_0_5_MeV_error/Edge_Q_0_5_MeV,2));
	  if (!(std::find(long_term_jump_bi_vec.begin(), long_term_jump_bi_vec.end(), om_number_data) != long_term_jump_bi_vec.end() ||
	      std::find(long_term_jump_li_vec.begin(), long_term_jump_li_vec.end(), om_number_data) != long_term_jump_li_vec.end() ||
		//std::find(Se_num_vec.begin(), Se_num_vec.end(), om_number_data) != Se_num_vec.end() ||
		std::find(b_num_vec.begin(), b_num_vec.end(), om_number_data) != b_num_vec.end())) {

	    fwhm = 100*2.355*(sigma_1_MeV/Q_1_MeV)*sqrt(vis_pic_1MeV/1000);
	    fwhm_error = fwhm * sqrt(pow(sigma_1_MeV_error/sigma_1_MeV,2)+pow(Q_1_MeV_error/Q_1_MeV,2)+pow(vis_pic_1MeV_error/(2*vis_pic_1MeV),2));}
	  else{
	    fwhm=0;
	  }

	}
	else{
	  fwhm=0;
	  a = 0;
	  a_error = 0;
	  corr_simu = 0;
	  corr_simu_error = 0;
	  E_vis_MC = 0;
	  E_vis_MC_error = 0;
	  E_vis_data = 0;
	  E_vis_data_error = 0;	  
	}
	csv_file << om_number_data << "  " << corr_simu << "\n";
	//csv_file << om_number_data << "  " << a/0.23841858 << "\n";
	t_out->Fill();
	break;
      }
    }
  }
  t_out->Write();
  fout->Close();
  //  csv_file.close();

}




int main(int argc, char const *argv[]){
  int n_run, run;
  std::vector<int> ref_run_number, ref_time, energy_run_number/*, run_number*/;
  int compteur = 0;
  string file, ref_correction, calo_correction;
  bool add = false;
  bool energy = false;

  //Passage au calo 
   std::cout << "" << '\n';
   std::cout << "START OF THE CALORIMETER FIT" << '\n';
   std::cout << "" << '\n';

     for(int run_value : vecteur_Bi){
       fit_Bi_energy(run_value);
       Merger_data_simu(run_value);
       std::cout<<"fit_BI ok "<<run_value<<std::endl;     
     }
     return 0;
}









