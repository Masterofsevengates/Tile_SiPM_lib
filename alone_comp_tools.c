/*
Tool Library for analysing SiPMs and tiles from the Tiletester
@C Sebastian Laurien (sebastian.laurine@gmx.de) 
loaded in root with .L alone_comp_tools.c+

Functions:
Double_t langau(Double_t *x, Double_t *par) 
	Landau Gauss convolution. Just the function returning a double

TTree * load(TString inputfile,TString Name) 
	Extracting a TTree from the inputfile


TH1F * Histo2 (TTree * treeout, int channel, TString nameadd = "", Double_t min=1, Double_t max=4096, int bins=4096)
	Extracting a TH1F from our data format... only works with files from daqmul and related... 
	Needs to be fed with a TTree (see above). Name of the histogramm is "aaa_channel_nameadd 
	


*/
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include "TGraph.h"



/*Double_t fit_response(Double_t *x, Double_t *par){
  Double_t xx=x[0];
  Double_t result = TMath::Exp(par[0](xx-par[1])); //Gausfunction
  return result;
}

*/
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


TTree * load(TString inputfile, TString name="data"){
  TFile * datafile = new TFile((TString)inputfile);
  TTree * treeout  = (TTree*)datafile->Get(name);
  return treeout;
}

TH1F * Histo2 (TTree * treeout, int channel, TString nameadd = "", Double_t min=1, Double_t max=4096, int bins=4096){
  Int_t nEnt = treeout->GetEntries();
  TString branch = "qdcch";
  branch += channel;
  TBranch * b_qdchit = treeout->GetBranch(branch);
  Int_t qdc;
  b_qdchit->SetAddress(&qdc);
  TString name = "aaa";
  name += "_";
  name += channel;
  name+=nameadd;
  TH1F * hqdc = new TH1F(name,name,bins,min,max); //filling histo
  for(Int_t i=0; i<nEnt; i++){

    b_qdchit->GetEntry(i);
    hqdc->Fill(qdc);
    //  cout<<qdc<<endl;

  }
  return hqdc;
}

TH1F * Histo2_cut (TTree * treeout, int channel, int cut , TString nameadd = "", Double_t min=1, Double_t max=4096, int bins=4096){
  Int_t nEnt = treeout->GetEntries();
  TString branch = "qdcch";
  branch += channel;
  TBranch * b_qdchit = treeout->GetBranch(branch);
  Int_t qdc;
  b_qdchit->SetAddress(&qdc);
  TString name = "aaa";
  name += "_";
  name += channel;
  name+=nameadd;
  TH1F * hqdc = new TH1F(name,name,bins,min,max); //filling histo
  for(Int_t i=0; i<nEnt; i++){

    b_qdchit->GetEntry(i);
    if (qdc>cut) hqdc->Fill(qdc);
    //  cout<<qdc<<endl;

  }
  return hqdc;
}




 TH1F * Histo_LG (TTree * treeout, int channel, TString nameadd = "", Double_t min=1, Double_t max=4096, int bins=40960){
    channel+=16;
  Int_t nEnt = treeout->GetEntries();
  TString branch = "qdcch";
  branch += channel;
  TBranch * b_qdchit = treeout->GetBranch(branch);
  Int_t qdc;
  b_qdchit->SetAddress(&qdc);
  TString name = "aaa";
  name += "_";
  name += channel;
  name+=nameadd;

  //get the conversion! 
  // Double_t LG[100000],HG[100000];
  TH1F * hqdc = new TH1F(name,name,bins,min,max*10); //filling histo
  for(Int_t i=0; i<nEnt; i++){

    b_qdchit->GetEntry(i);
    hqdc->Fill(qdc*10);
    //  cout<<qdc<<endl;

  }
  return hqdc;
}
TH1I * Histo2I (TTree * treeout, int channel, TString nameadd = "", Double_t min=1, Double_t max=4096, int bins=4096){
  Int_t nEnt = treeout->GetEntries();
  TString branch = "qdcch";
  branch += channel;
  TBranch * b_qdchit = treeout->GetBranch(branch);
  Int_t qdc;
  b_qdchit->SetAddress(&qdc);
  TString name = "aaa";
  name += "_";
  name += channel;
  name+=nameadd;
  TH1I * hqdc = new TH1I(name,name,bins,min,max); //filling histo
  for(Int_t i=0; i<nEnt; i++){

    b_qdchit->GetEntry(i);
    hqdc->Fill(qdc);
    //  cout<<qdc<<endl;

  }
  return hqdc;
}

/*TH1I * Pedestal_stabi (int step,TTree * thetree,int channel){
  TH1I * theresult = new TH1I("result","result",2096,1,4096);
  Double_t x[1000];
  Double_t x_err[1000];
  Double_t y_err[1000];
  Double_t y[1000];
  for (int i = 0 ; i<1000; i++){
    x[i]=0;
    x_err[i]=0;
    y[i]=0;
    y_err[i]=0;
  }  

  for (int i = 0; i<999; i++){
    TH1F * histo = Histo3(thetree,channel,i*step,i*step+step,"",1,800,800);
    x[i]=getzpeak(histo);
    // histo->Fit("gaus");
    // TF1 * myfit = (TF1*) histo->GetFunction("gaus");
    //  x[i] = histo->GetFunction("gaus")->GetParameter(1);     
    x[i]=histo->GetMean();
    y[i]=i;

    TString branch = "qdcch";
    branch += channel;
    TBranch * b_qdchit = thetree->GetBranch(branch);
    Int_t qdc;
    b_qdchit->SetAddress(&qdc);
    for (int j = i*step; j<i*step+step; j++){
      if (j>thetree->GetEntries()-2) break;
      b_qdchit->GetEntry(j);
      // cout<<"loading : "<<j<<" of"<<thetree->GetEntries()<<" "<<qdc<<endl;
      theresult->Fill(qdc-x[i]+x[0]);
    }
  }
  TGraph * theplot = new TGraph(1000,y,x);
  theplot->Draw("A*");
 
  return theresult;

}
*/
/*
TH1F * Histo3 (TTree * treeout, int channel, int start, int entries, TString nameadd = "", Double_t min=1, Double_t max=4096, int bins=4096){
  Int_t nEnt = treeout->GetEntries();
  TString branch = "qdcch";
  branch += channel;
  TBranch * b_qdchit = treeout->GetBranch(branch);
  Int_t qdc;
  b_qdchit->SetAddress(&qdc);
  TString name = "aaa";
  name += "_";
  name += channel;
  name+=nameadd;
  TH1F * hqdc = new TH1F(name,name,bins,min,max); //filling histo
  for(Int_t i=start; i<entries; i++){

    b_qdchit->GetEntry(i);
    hqdc->Fill(qdc);
    //  cout<<qdc<<endl;

  }
  return hqdc;
}

TH1F * Temperature_hist(TTree * treeout){
  TBranch * b_qdchit = treeout->GetBranch("temperature");
  Float_t qdc;
  b_qdchit->SetAddress(&qdc);
  TH1F * hqdc = new TH1F("Temperature","Temperature",80,0,40);
  for(Int_t i=0; i<treeout->GetEntries(); i++){
    b_qdchit->GetEntry(i);
    hqdc->Fill(qdc);
  }
  return hqdc;
}
Double_t * Mean_Temperature(TTree * treeout){
  TH1F * temp = Temperature_hist(treeout);
  Double_t * temperature = new Double_t[2];
  temperature[0] = temp->GetMean();
  temperature[1] = 0.1 + temp->GetMeanError();
  return temperature;


}
*/
TH1F * Histo(TString inputfile, int channel, Double_t min=-0.5, Double_t max=4000.5, int bins=4001){
  //procedure from Marco to get the data from a histogramm
  //opening file and reding tree

  TFile * datafile = new TFile((TString)inputfile);
  TTree * treeout  = (TTree*)datafile->Get("data");

  Int_t nEnt = treeout->GetEntries();
  TString branch = "qdcch";
  branch += channel;
  TBranch * b_qdchit = treeout->GetBranch(branch);
  Int_t qdc;
  b_qdchit->SetAddress(&qdc);
  TString name = inputfile;
  name += "_";
  name += channel;
  TH1F * hqdc = new TH1F(name,name,bins,min,max); //filling histo
  for(Int_t i=0; i<nEnt; i++){
    
    b_qdchit->GetEntry(i);
    hqdc->Fill(qdc);
   
    //  cout<<qdc<<endl;

  }
  return hqdc;

}
//to get data out of file
Double_t * xdata(TTree * treeout, int channel){
  Int_t nEnt = treeout->GetEntries();
  Double_t * x;
  x = new Double_t[nEnt];
  TString branch = "qdcch";
  branch += channel;
  TBranch * b_qdchit = treeout->GetBranch(branch);
  Int_t qdc;
  b_qdchit->SetAddress(&qdc);
  for (int i = 0; i<nEnt; i++){
    b_qdchit->GetEntry(i);
    x[i]=qdc;
  }
  return x;
}
//doublepeakfit
int dfit(TString input1,int ch,Double_t *xl,Double_t *xh,Double_t *xlerr, Double_t *xherr){
  int chl = ch;
  int chh = ch+16;
  TF1 *gfit1 = new TF1("Gaussian","gaus");
  TF1 *gfit2 = new TF1("Gaussian","gaus");
  TH1F * hist1 = Histo( input1, chl);
  TH1F * hist2 = Histo( input1 ,chh);
  cout<<"here i am"<<endl;
  hist1->Fit(gfit1,"SN");
  hist2->Fit(gfit2,"SN");
  if (gfit1->GetParameter(1) > 4000 or gfit1->GetParameter(1)<0){
    *xl = 0;
    *xlerr = 0;
  }
  else {
    *xl=gfit1->GetParameter(1);
    *xlerr=gfit1->GetParError(1);
  }
  if (gfit1->GetParameter(1) > 4000 or gfit1->GetParameter(1)<0){
    *xh = 0;
    *xherr = 0;
  }
  else{

    *xh=gfit2->GetParameter(1);
    *xherr=gfit2->GetParError(1);
  }
  return 0;

}
//again different parameters
int dfit(TString input1,int ch,Double_t *xl,Double_t *xh){
  int chl = ch;
  int chh = ch+16;
  TF1 *gfit1 = new TF1("Gaussian","gaus");
  TF1 *gfit2 = new TF1("Gaussian","gaus");
  TH1F * hist1 = Histo( input1, chl);
  TH1F * hist2 = Histo( input1 ,chh);
  cout<<"here i am"<<endl;
  hist1->Fit(gfit1,"SN");
  hist2->Fit(gfit2,"SN");
  if (gfit1->GetParameter(1) > 4000 or gfit1->GetParameter(1)<0){
    *xl = 0;
  }
  else {
    *xl=gfit1->GetParameter(1);

  }
  if (gfit1->GetParameter(1) > 4000 or gfit1->GetParameter(1)<0){
    *xh = 0;

  }
  else{

    *xh=gfit2->GetParameter(1);

  }
  return 0;

}
int ffit(TString input1,int ch1,int ch2,Double_t *xl1 ,Double_t *xh1, Double_t *xl2, Double_t *xh2){
  Double_t x1,x2,x3,x4;
  dfit(input1,ch1,&x1,&x2);
  dfit(input1,ch2,&x3,&x4);
  *xl1 = x1;
  *xh1 = x2;
  *xl2 = x3;
  *xh2 = x4;
  return 0;
}
void xdata(TTree * treeout,int ch1, Double_t * xl1,Double_t * xh1){
  Int_t nEnt = treeout->GetEntries();
  TString branch = "qdcch";
  branch += ch1;
  TBranch * b_qdchit = treeout->GetBranch(branch);
  Int_t qdc;
  b_qdchit->SetAddress(&qdc);
  for (int i = 0; i<nEnt; i++){
    b_qdchit->GetEntry(i);
    xl1[i]=qdc;
    //cout<<x[i]<<"   "<<i<<endl;
  }
  branch = "qdcch";
  ch1 +=16;
  branch += ch1;
  b_qdchit = treeout->GetBranch(branch);
  b_qdchit->SetAddress(&qdc);
  for (int i = 0; i<nEnt; i++){
    b_qdchit->GetEntry(i);
    xh1[i]=qdc;
    //cout<<x[i]<<"   "<<i<<endl;
  }
}
Int_t get_runs(TString inputfile){
  TFile * datafile = new TFile((TString)inputfile);
  TTree * treeout = (TTree*)datafile->Get("data");
  Int_t nEnt = treeout->GetEntries();
  return nEnt;
}
void bubbleSort(Double_t *  numbers, int array_size)
{
  int i, j;
  Double_t temp;
  for (i = (array_size - 1); i > 0; i--)
    {
      for (j = 1; j <= i; j++)
	{
	  if (numbers[j-1] > numbers[j])
	    {
	      temp = numbers[j-1];
	      numbers[j-1] = numbers[j];
	      numbers[j] = temp;
	    }
	}
    }
}
Double_t maxval(Double_t *ary, int length){
  Double_t val = 0;
  for (int i=0; i<length; i++){
    if (ary[i] >val) val=ary[i];
  }
  return val;
}
Double_t minval(Double_t *ary, int length){
  Double_t val = 50000000000000;
  for (int i=0; i<length; i++){
    if (ary[i] <val) val=ary[i];
  }
  return val;
}
TH1F * Histo(TString inputfile){
  //procedure from Marco to get the data from a histogramm
  //opening file and reding tree
  TFile * datafile = new TFile((TString)inputfile);
  TTree * treeout  = (TTree*)datafile->Get("data");
  Int_t nEnt = treeout->GetEntries();
  TBranch * b_qdchit = treeout->GetBranch("qdcch");
  Int_t qdc;
  b_qdchit->SetAddress(&qdc);
  TH1F * hqdc = new TH1F("hqdc","hqdc",4001,-0.5,4000.5); //filling histo
  for(Int_t i=0; i<nEnt; i++){
    b_qdchit->GetEntry(i);
    hqdc->Fill(qdc);
  }
  return hqdc;
}
/*Double_t *  Recast(Float_t * ar, int found){//returns a sorted Doulbe from floats
  Double_t  x[found];
  for (int i =0; i<found;i++){
    x[i]=ar[i];
  }
  bubbleSort(x,found);
  return x;
}
*/
int GetSize (Double_t *par){//get the size of an array ending with 0
  int length=0;
  while(par[length]!=0)length++;
  return length+1;
}
Double_t multigaus(Double_t *x, Double_t *par){
  Double_t xx=x[0];
  Double_t result = par[0];
  int number= (GetSize(par)-2)/3;
  for (int i= 0; i<number;i++){
    result+=TMath::Gaus(xx,par[i*3+2],par[i*3+3])*par[i*3+1]; //Gausfunction
  }
  return result;
}
Double_t phpeakfunc(Double_t *x, Double_t *par){
  //par[0] = distance
  //par[1] = offset
  //par[2] = postition of the fisrt peak
  //par[3] = width first peak
  //par[4] = Gain first peak
  //par[5] = width of second peak
  //par[6] = Gain secopnd peak
  Double_t xx=x[0];
  Double_t result = par[1];
  int number= (GetSize(par)-4)/2;
  for (int i= 0; i<number;i++){
    result+=((TMath::Gaus(xx,par[2]+(i*par[0]),par[3+(i*2)]))*par[4+(i*2)]); //Gausfunction
  }
  return result;
}
TF1 * fitmultigaus(int min, int max,int par) {
  TF1 * fit = new TF1("fitmultigaus",multigaus ,min,max,par);
  return fit;
}
TF1 * phpeakfit(int min, int max, int par){
  TF1 * fit = new TF1("phpeakfit",phpeakfunc ,min,max,par);
  return fit;
}
Double_t Gain (TH1F * input){

  TSpectrum *spectrum = new TSpectrum(100);
  int pfound = spectrum->Search(input,2,"goff",0.05);
  cout<<pfound<<endl;
  Float_t *xtemp = spectrum->GetPositionX();
  Double_t x[pfound];
  for (int i =0; i<pfound;i++){
    x[i]=xtemp[i];
  }
  bubbleSort(x,pfound);
  int max = (int)(x[pfound-1]+x[1]-x[0]);
  int length = 4+pfound*2;
  TF1 *fitfunc = phpeakfit(0,max,length);
  Double_t  parameter[length];
  parameter[length-1]= 0;
  parameter[1] = 0.01;
  parameter[2] = x[0]; //position of the first peak
  parameter[0] = x[1] -x[0];
  for (int i=0; i<pfound; i++){
    parameter[3+(i*2)] = 1;//set the width
    // parameter[2+i*3] = x[i]; //mean no nedd to any more
    parameter[4+(i*2)]=1000;//set the gain
  }
  parameter[length-1]= 0;
  fitfunc->SetParameters(parameter);
  fitfunc->SetLineColor(kRed);
  fitfunc->FixParameter(length-1,0);
  fitfunc->SetParLimits(2, x[0]-50,x[0]+50);
  fitfunc->SetParLimits(0, 20,100);
  input->Fit("phpeakfit");
  fitfunc->SetNpx(10000);
  TF1 * result = input->GetFunction("phpeakfit");
  cout<<result->GetParameter(0)<<"   "<<result->GetParError(0)<<endl;
  //now we have an fit to the gain.. make it better!!
  //how many peaks should there be?
  int peaks = 1+((max-result->GetParameter(2))/result->GetParameter(0));
  length = 4+peaks*2;
  cout<<peaks<<"   "<<length<<endl;
  //deja vue
  Double_t  parameter2[length];
  parameter2[length-1]= 0;
  parameter2[1] = 0.01;
  parameter2[2] = result->GetParameter(2); //position of the first peak
  parameter2[0] = result->GetParameter(0);
  for (int i=0; i<pfound; i++){
    parameter2[3+(i*2)] = 1;//set the width
    // parameter[2+i*3] = x[i]; //mean no nedd to any more
    parameter2[4+(i*2)]=1000;//set the gain
  }
  parameter2[length-1]= 0;
  fitfunc->SetParameters(parameter2);
  fitfunc->SetLineColor(kRed);
  fitfunc->FixParameter(length-1,0);
  fitfunc->SetParLimits(2, x[0]-50,x[0]+50);
  fitfunc->SetParLimits(0, 10,150);
  fitfunc->SetParLimits(1,0,10);
  //doto here!!!
  for (int i=0; i<pfound; i++){
    fitfunc->SetParLimits(3+i*2 ,1,10000);
    fitfunc->SetParLimits(4+i*2 ,0.1,1000000000);
  }
  input->Fit("phpeakfit");
  fitfunc->SetNpx(10000);
  result = input->GetFunction("phpeakfit");
  cout<<result->GetParameter(0)<<"   "<<result->GetParError(0)<<endl;
  return result->GetParameter(0);
}
//-----------------------------------------------------------------------
//
//	Convoluted Landau and Gaussian Fitting Function
//         (using ROOT's Landau and Gauss functions)
//
//  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
//  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
//   Markus Friedl (Markus.Friedl@cern.ch)
//
//  to execute this example, do:
//  root > .x langaus.C
// or
//  root > .x langaus.C++
//
//-----------------------------------------------------------------------
Double_t langaufun(Double_t *x, Double_t *par) {
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //this one i wanna have!!!!
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
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  mpc = par[1] - mpshift * par[0];
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow) / np;
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
Double_t getped(int channel,TTree * filename){
  TH1F * ch_2 = Histo2(filename,channel);
  return ch_2->GetMean();
}
Double_t getzpeak(TH1F * ch_2){
  Double_t maxbin = ch_2->GetBinContent(ch_2->GetMaximumBin());
  for (int i=0 ; i<4000;  i++){
    if (ch_2->GetBinContent(i)<maxbin*0.5){
      ch_2->SetBinContent(i,0);
    }
  }
  return ch_2->GetMean();
}
Double_t getzpeak_error(TH1F * ch_2){
  Double_t maxbin = ch_2->GetBinContent(ch_2->GetMaximumBin());
  for (int i=0 ; i<4000;  i++){
    if (ch_2->GetBinContent(i)<maxbin*0.5){
      ch_2->SetBinContent(i,0);
    }
  }
  return ch_2->GetMeanError();
}

Double_t * Effi(int tile, Double_t voltage, Double_t Temperature, Double_t * threshold, TSQLServer * db){
  Double_t * x;
  int length = 0; 
  int i = 0;
  while (threshold[i]!=0){
    length++;
    i++;
    
  }
  x = new Double_t[length];
  Double_t zp=0.0;
  TSQLRow *row;
  TSQLResult *res;
  char sql[200];
  sprintf(sql, "select * from tile_Source_fit where (number=%d and TEMPERATURE =%d and V>%f and V <%f)",tile,22,voltage-0.1,voltage+0.1);
  res = db->Query(sql);
  if (res->GetRowCount()!=1){
    cout<<"ERROR "<<sql<<endl;
    for (int i = 0; i<length; i++){
      x[i] = 0;
    }
  }
  else{
  row = res->Next();
  zp=atof(row->GetField(8));
  TF1 *fit = new TF1("fit",langau,0,4000,4);
  fit->SetParameter(0,atof(row->GetField(4)));
  fit->SetParameter(1,atof(row->GetField(5)));
  fit->SetParameter(2,atof(row->GetField(6)));
  fit->SetParameter(3,atof(row->GetField(7)));
  fit->SetRange(0,40000);
  // get the maximum! 
  Double_t Mip=  fit->GetMaximumX();
  for (int i = 0; i<length; i++){
    x[i] = fit->Integral(zp+(threshold[i]*(Mip-zp)),40000)/fit->Integral(0,40000) ;
  }
  }
  return  x;
}


TGraph * Response_Curve(int tile, Double_t temp, TSQLServer * db){
  Double_t x[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t y[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  char sqla[2000];
  TSQLRow *row;
  TSQLResult *res;
  sprintf(sqla,"select tile_Source.V,(tile_Source.Source-tile_Pedestal.pedestal)/tile_gain.gain from (tile_Source inner join tile_Pedestal on (tile_Source.number = tile_Pedestal.number and tile_Source.TEMPERATURE=tile_Pedestal.TEMPERATURE and tile_Source.V = tile_Pedestal.V)) inner join tile_gain on (tile_Source.number = tile_gain.number and tile_Source.TEMPERATURE=tile_gain.TEMPERATURE and tile_Source.V = tile_gain.V) where tile_Source.number =%d and tile_Source.TEMPERATURE =%f and (tile_Source.Source-tile_Pedestal.pedestal)>5",tile,temp);
  //  printf(sqla,"select tile_Source.V,(tile_Source.Source-tile_Pedestal.pedestal) from (tile_Source inner join tile_Pedestal on (tile_Source.number = tile_Pedestal.number and tile_Source.TEMPERATURE=tile_Pedestal.TEMPERATURE and tile_Source.V = tile_Pedestal.V)) inner join tile_gain on (tile_Source.number = tile_gain.number and tile_Source.TEMPERATURE=tile_gain.TEMPERATURE and tile_Source.V = tile_gain.V) where tile_Source.number =%d and tile_Source.TEMPERATURE =%f",tile,temp);
  res = db->Query(sqla);
  int length = res->GetRowCount();
  if (length==0) cout<<"NoTile!"<<endl;
  for (int i = 0; i < length; i++) {
  	row = res->Next();
        x[i] =atof(row->GetField(0));
        y[i] =atof(row->GetField(1));
	//	cout<<x[i]<<" "<<y[i]<<endl;
  }
  TGraph * graph = new TGraph(length,x,y);
  return graph;
}


TGraph * Pedestal_Curve(int tile, Double_t temp, TSQLServer * db,bool norm = false){
  Double_t x[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t y[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  char sql[2000];
  TSQLRow *row;
  TSQLResult *res;
  if (norm) sprintf(sql,"select tile_Pedestal.V,tile_Pedestal.Pedestal/tile_gain.gain from tile_Pedestal inner join tile_gain on (tile_Pedestal.number = tile_gain.number and tile_Pedestal.TEMPERATURE=tile_gain.TEMPERATURE and tile_Pedestal.V = tile_gain.V) where tile_Pedestal.number =%d and tile_Pedestal.TEMPERATURE =%f;",tile,temp);
  else sprintf(sql,"select tile_Pedestal.V,tile_Pedestal.Pedestal from tile_Pedestal where tile_Pedestal.number =%d and tile_Pedestal.TEMPERATURE =%f;",tile,temp);
  res = db->Query(sql);
  int length = res->GetRowCount();
  if (length==0) cout<<"NoTile!"<<endl;
  for (int i = 0; i < length; i++) {
  	row = res->Next();
        x[i] =atof(row->GetField(0));
        y[i] =atof(row->GetField(1));
  }
  TGraph * graph = new TGraph(14,x,y);
  return graph;
}


Double_t * Px_Effi(int tile, Double_t voltage, Double_t Temperature, Double_t * threshold, TSQLServer * db){
  Double_t * x;
  int length = 1; 
  x = new Double_t[length];
  TGraph * resp = Response_Curve(tile,Temperature,db);
  TF1 * myfit = new TF1("myfit","pol2");
  resp->Fit(myfit,"Q");
  for (int i = 0; i<length; i++){
      x[i] = myfit->Eval(voltage)*threshold[i];
      //    cout<<x[i]<<" "<<voltage<<" "<<threshold[i]<<endl;
  
  }
  return x;
}

Double_t * Px_DCR(int tile, Double_t voltage, Double_t Temperature,Double_t Gate,Double_t * values, TSQLServer * db){
  Double_t * x;
  int length = 1;
  x = new Double_t[length];
  TGraph * resp = Pedestal_Curve(tile,Temperature,db,true);
  TF1 * myfit = new TF1("myfit","pol2");
  resp->Fit(myfit,"Q");
  Double_t mean_ped = myfit->Eval(voltage);
  // ok lets make a poissonian out of it.. 
  TF1 *poiss = new TF1("poiss","TMath::Poisson(x,[0])",0,200);
  poiss->SetParameter(0,1*mean_ped);// overestimate the noise! 
  Double_t Int_all = poiss->Integral(0,200);
  for (int i = 0; i<length;i ++){
    x[i] =(1.0/(Gate))*-log(1.0- (poiss->Integral(values[i],200)/poiss->Integral(0,200)));
    //   cout<<x[i]<<endl;
  }
  


  return x;
}

