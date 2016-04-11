/*

Functions
@sebastian.laurien@gmx.de
load: .L anaTileTester.C+

Fit the gain of a SiPM  
TF1 * MultiPFit(TH1F *  peaks, Int_t nOfpeaks, Float_t res,Float_t sigmapeak,Float_t threshold)

peaks: input histogram
nOfpeaks: how many peaks shoud be fitted
res: nresolution between peaks (min. value 1)
sigmapeak:  sigma of the peaks (minimum)
threshold: peaks with amplitude less than threshold*highest_peak are discarded


*/

#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <time.h>
/****************************/
/******* ROOT **************/
#include <TF1.h>
#include <TH1F.h>
#include "TH2F.h"
#include <TH1D.h>
#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TChain.h>
#include "TMath.h"
#include "TLegend.h"
#include "TPad.h"
#include "TObject.h"
#include "TList.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH3F.h"
#include "TLine.h"
#include "TList.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TSpectrum.h"
#include "TText.h"
#include "TROOT.h"
#include "TVirtualFFT.h"
#include "ExtractGain.cc"



template<typename T>

std::string to_string(T val) {
  std::stringstream ss;
  ss<<val;
  return ss.str();
}


using namespace std;

Double_t MultiGaus(Double_t *x, Double_t *par){//multi peak gaussian function
  //sum of npeaks gaussians
  //peak distance is a fit parameter
  //
  Double_t xx = x[0];
  Double_t npeaks = par[0];//number of peaks
  Double_t gain = par[1];//gain
  Double_t p0 = par[2];//first peak position

  Double_t mg=0;
  for(Int_t i=0; i<npeaks; i++){

    mg += par[2*i+3]*TMath::Gaus(xx,gain*i+p0,par[2*i+4]);
  }
  return mg;
}

Double_t langau(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

Float_t FindPedestal(TH1F* hist){
  Float_t p(0.);
  TSpectrum * spectrum = new TSpectrum(1,1);//searching for peaks
  spectrum->Search(hist,10,"goff",0.5);
  Float_t * maxbin = new Float_t[1];
  maxbin = spectrum->GetPositionX();
  Float_t rms = hist->GetRMS();
  Float_t xinf = std::max(0.,(double)maxbin[0]-2*rms);
  hist->GetXaxis()->SetRangeUser(xinf,maxbin[0]+2*rms);
  p = hist->GetMean();
  return p;
}
Float_t FindZP(TH1F * hist){
   Float_t p(0.);
   Double_t maxbin = hist->GetBinContent(hist->GetMaximumBin());
    for (int i=0 ; i<hist->GetNbinsX();  i++){
        if (hist->GetBinContent(i)<maxbin*0.9){
        hist->SetBinContent(i,0);
        }
    }
  p = hist->GetMean();
  return p;
}

Double_t FindZPinMip(TH1F * hist,int scale=1){
  // make a copy, kill everything over 1000
  // find the highest bin! 
  Float_t p(0.);
  
  TH1F * temphist = (TH1F*)hist->Clone("temphist");
    for (int i=1000/scale ; i<4000/scale;  i++){
          temphist->SetBinContent(i,0);
        }
    
    p = temphist->GetMaximumBin();
  return p;

}
Double_t FindZPinMipValue(TH1F * hist,int scale =1){
  // make a copy, kill everything over 1000
  // find the highest bin! 
  Float_t p(0.);
  TH1F * temphist = (TH1F*)hist->Clone("temphist");
    for (int i=1500/scale ; i<4000/scale;  i++){
          temphist->SetBinContent(i,0);
        }
    
    p = temphist->GetBinContent(temphist->GetMaximumBin());
  return p;

}
TH1F * PedestalSuppressedMip(TH1F * mip_hist, TH1F * ped_hist, Double_t scale){
  // make a copy of the mip get FindZPinMIP normalize highest bin of ped to this and substract it! 
 TH1F * temphist = (TH1F*)mip_hist->Clone("temphist");
 temphist->Rebin(scale);
 TH1F * temphist_ped= (TH1F*)ped_hist->Clone("temphist_ped");
 TH1F* shift_ped =(TH1F*)ped_hist->Clone("shift_ped");
 temphist_ped->Rebin(scale);
 shift_ped->Rebin(scale);
 // TH1F* shift_ped = new TH1F("tmp_ped","tmp_ped",204,1,4001);
 // TH1F* shift_mip = new TH1F("tmp_mip","tmp_mip",4096,1,4096);
 // temphist_ped->Scale(FindZPinMipValue(temphist)/FindZPinMipValue(temphist_ped));
 // temphist_ped->Scale(FindZPinMipValue(temphist,1)/FindZPinMipValue(temphist_ped,1));
 Double_t Scaler= (FindZPinMipValue(temphist,scale)/FindZPinMipValue(temphist_ped,scale));
 // shift_ped->Rebin(10);
 // temphist->Rebin(10);
 //temphist_ped->Rebin(10);
 // temphist = (TH1F*)mip_hist->Clone("temphist");
 // temphist_ped= (TH1F*)ped_hist->Clone("temphist_ped");
 temphist_ped->Scale(Scaler);
 Int_t Shift = FindZPinMip(temphist,scale)-FindZPinMip(temphist_ped,scale);
 for(int i=0;i<=3000;i++){
        Float_t bincontent = temphist_ped->GetBinContent(i);
	if(bincontent>0) shift_ped->SetBinContent(Shift+i+1,bincontent);
 }
 // TH1F * temphist = (TH1F*)mip_hist->Clone("temphist");
 temphist->Add(shift_ped,-1);
 
 return temphist;
}



TF1 * FitwithLG(TH1F * hin, Int_t rebin=0){
  if(rebin >1) hin->Rebin(rebin);

  Float_t mean = hin->GetMean();
  Float_t rms = hin->GetRMS();
  Float_t rmsmin = 50.;
  Float_t rmsmax = 1000.;

  Float_t xinf(0.);
  Float_t xsup(0.);

  xinf = std::max(mean-1*rms,(Float_t)hin->GetBinLowEdge(hin->FindFirstBinAbove(0)));
  xsup = std::min((Float_t)hin->GetBinLowEdge(hin->FindLastBinAbove(0)),(Float_t)(mean+rms));


  Float_t rmssup = std::min(2*rms,rmsmax);
  if(rmssup<=rmsmin) rmssup = rmsmax;

  Float_t A = hin->Integral();
  if (xsup>3500) xsup=3500;
  TF1 *fit = new TF1("fit",langau,xinf,xsup,4);
  fit->SetNpx(100000);
  fit->SetParameter(0,rms);
  fit->SetParLimits(0,0.2*rms,3*rms);
  fit->SetParameter(1,mean);
  fit->SetParLimits(1,xinf,mean+2*rms);
  fit->SetParameter(2,A);
  fit->SetParLimits(2,A-TMath::Sqrt(A),A+10*TMath::Sqrt(A));
  fit->SetParameter(3,rms);
  fit->SetParLimits(3,rmsmin,rmssup);



  hin->Fit(fit,"RQN0");

  return fit;
}


TF1 * MultiPFit(TH1F *  peaks, Int_t nOfpeaks, Float_t res,Float_t sigmapeak,Float_t threshold){

  //nOfpeaks -  sets the max n. of peaks stored
  //res -  resolution between peaks (min. value 1)
  //sigmapeak - sigma
  // threshold - peaks with amplitude less than threshold*highest_peak are discarded

  /*************finding peaks*****************/
  TSpectrum * spectrum = new TSpectrum(nOfpeaks,res);//searching for peaks
  spectrum->Search(peaks,sigmapeak,"goff",threshold);

  Int_t nPeaks = spectrum->GetNPeaks();
  TF1 * mpeak = NULL;
  if (nPeaks > 1) {

    Float_t * x = new Float_t[nPeaks];
    Float_t * y = new Float_t[nPeaks];

    x = spectrum->GetPositionX();//getting peak positions
    y = spectrum->GetPositionY();//getting amplitude

    Int_t lookup[nPeaks];//creating a lookup table to order peaks
    Float_t dummyx[nPeaks];
    Float_t max = 20000;
    for(Int_t i=0; i<nPeaks; i++){
      dummyx[i] = x[i];
    }
    for(Int_t j=0; j<nPeaks; j++){//ordering peaks
      for(Int_t i=0; i<nPeaks; i++){
	if(dummyx[i]<max){
	  max = dummyx[i];
	  lookup[j]=i;
	}
      }
      max = 20000;
      dummyx[lookup[j]]=30000;
    }



    Int_t istart(0);
    Int_t PeakN(nPeaks);
    if(nPeaks>2){
      if((x[lookup[1]]-x[lookup[0]]) > 1.5*(x[lookup[2]]-x[lookup[1]])){
	istart = 1;
	PeakN = nPeaks-1;
      }
    }
    Float_t gain = x[lookup[istart+1]]-x[lookup[istart]]; //gain distance between first 2 peaks
    Float_t szero = gain/5; //initial sigma
    Float_t xinf = std::max((Float_t)(x[lookup[istart]]-2.*szero),(Float_t)0.);
    Float_t xsup = std::min((Float_t)(x[lookup[nPeaks-1]]+2.*szero),(Float_t)4000.);
    mpeak = new TF1("mpeak",MultiGaus,xinf,xsup,2*PeakN+3);
    mpeak->SetNpx(100000);
    mpeak->SetLineColor(kRed);
    mpeak->FixParameter(0,PeakN);//number of peaks FIXED
    mpeak->SetParameter(1,gain);//gain
    mpeak->SetParLimits(1,0.6*gain,1.2*gain);
    mpeak->SetParameter(2,x[lookup[istart]]);//pede
    mpeak->SetParLimits(2,x[lookup[istart]]-szero/2,x[lookup[istart]]+szero/2);
    for(Int_t i=0; i<PeakN; i++){
      mpeak->SetParameter(2*i+3,y[lookup[i+istart]]);
      mpeak->SetParLimits(2*i+3,0.8*y[lookup[i+istart]],1.5*y[lookup[i+istart]]);//height
      mpeak->SetParameter(2*i+4,TMath::Sqrt(i+1)*szero);
      mpeak->SetParLimits(2*i+4,TMath::Sqrt(i+1)*szero/5,1.5*TMath::Sqrt(i+1)*szero);//sigma
    }
    peaks->Fit(mpeak,"RNOQ","",xinf,xsup);
    delete[] x;
    delete[] y;
    return mpeak;
  } else {
    return mpeak;
  }
}


Float_t RoundTemperature(Float_t tempIn){
  Float_t tempOut(0.);
  Float_t res = tempIn - (Float_t)floor(tempIn);

  if(res < 0.5){
    tempOut = (Float_t)floor(tempIn);
  }else{
    tempOut = (Float_t)floor(tempIn) +1.;
  }
  return tempOut;

}
Float_t ComputeLYError(Float_t mip,Float_t p, Float_t g, Float_t miperr, Float_t perr, Float_t gerr){

  Float_t ly = (mip-p)/g;
  Float_t lymipe = miperr/g;
  Float_t lyperr = -perr/g;
  Float_t lygerr = -gerr*ly/g;
  Float_t lyerr = TMath::Sqrt(lymipe*lymipe + lyperr*lyperr + lygerr*lygerr);
  return lyerr;
}
// here comes the whole DCR part... intput should be histobgramm, threshold gate.. out comes a dcr value... the rest goes into the script.. 

Float_t * DCR_rate (TH1F * hDCR, Float_t threshold, Float_t gate, Float_t gate_err){
  Float_t * DCR = new Float_t[2];
  Double_t Noise_Nerr(0);
  Double_t  Noise_Derr(0);
  Float_t Noise_N = hDCR->IntegralAndError(hDCR->FindBin(threshold),3500,Noise_Nerr);
  Float_t Noise_D = hDCR->IntegralAndError(1,3500,Noise_Derr);
  Float_t Proba =Noise_N/Noise_D;
  cout<<Proba<<endl;
  Float_t Proba_err = sqrt(pow(Noise_Derr/Noise_D,2)+pow(Noise_Nerr/Noise_N,2));
  DCR[0] =  (1.0/(gate))*-log(1.0-Proba);
 // DCR[1] = sqrt(pow((1.0/gate)*Proba_err,2));
  DCR[1] = sqrt(pow(gate_err/gate,2)+pow((1.0/gate)*Proba_err,2));

 return DCR;
}
 
Float_t * DCR_rate (TH1F * hDCR, Float_t gate, Float_t gate_err,TF1* myfit){
  Float_t * DCR = new Float_t[2];
  Double_t Noise_Nerr(0);
  Double_t  Noise_Derr(0);
  Float_t Noise_D = hDCR->IntegralAndError(1,3500,Noise_Derr);
  Float_t Noise_N = myfit->Integral(1,3500);
  Noise_Nerr = 10;
  Float_t Proba =Noise_N/Noise_D;
  cout<<"Probat :"<<Proba<<endl;
  Float_t Proba_err = sqrt(pow(Noise_Derr/Noise_D,2)+pow(Noise_Nerr/Noise_N,2));
  DCR[0] =  (1.0/(gate))*-log(Proba);
  cout<<"DCR: "<<DCR[0]<<endl;
  DCR[1] = sqrt(pow(gate_err/gate,2)+pow((1.0/gate)*Proba_err,2));

 return DCR;
}
 
