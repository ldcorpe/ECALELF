#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TAxis.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <cmath>

#include <string>
#include <cstring>
#include <sstream>


//returns a pair <effSigma,median>
std::pair <Double_t,Double_t> effSigma_function(RooDataSet *dat, RooRealVar invMass)
{

	//prepare finely binned histogram to receive RooDataSet
	TH1F *hist = new TH1F("hist","hist",1000,0,200);

	//Fill Histogram with RooDataSet without cuts
	dat->fillHistogram(hist, invMass,"",0);

	TAxis *xaxis = hist->GetXaxis();
	Int_t nb = xaxis->GetNbins();

	//exclude non valid cases
	if(nb < 10) {
		cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
		return std::make_pair(0.,0);
	}
	Double_t bwid = xaxis->GetBinWidth(1);
	if(bwid == 0) {
		cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
		return std::make_pair(0.,0);
	}


	Double_t xmin = xaxis->GetXmin();
	Double_t ave = hist->GetMean();
	Double_t rms = hist->GetRMS();

	Double_t total=0.;
	for(Int_t i=0; i<nb+2; i++) {
		total+=hist->GetBinContent(i);
	}


	Int_t ierr=0;
	Int_t ismin=999;

	Double_t rlim=0.683*total;
	Int_t nrms=rms/(bwid); // Set scan size to +/- rms
	if(nrms > nb/10) nrms=nb/10;//  Could be tuned...

	Double_t widmin=9999999.;
	for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
		Int_t ibm=(ave-xmin)/bwid+1+iscan;
		Double_t x=(ibm-0.5)*bwid+xmin;
		Double_t xj=x;
		Double_t xk=x;
		Int_t jbm=ibm;
		Int_t kbm=ibm;
		Double_t bin=hist->GetBinContent(ibm);
		total=bin;
		for(Int_t j=1;j<nb;j++){
			if(jbm < nb) {
				jbm++;
				xj+=bwid;
				bin=hist->GetBinContent(jbm);
				total+=bin;
				if(total > rlim) break;
			}
			else ierr=1;
			if(kbm > 0) {
				kbm--;
				xk-=bwid;
				bin=hist->GetBinContent(kbm);
				total+=bin;
				if(total > rlim) break;
			}
			else ierr=1;
		}
		Double_t dxf=(total-rlim)*bwid/bin;
		Double_t wid=(xj-xk+bwid-dxf)*0.5;
		if(wid < widmin) {
			widmin=wid;
			ismin=iscan;
		}
	}
	if(ismin == nrms || ismin == -nrms) ierr=3;
	if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;


 //compute the median for 1-d histogram hist
	 Int_t nbins = hist->GetXaxis()->GetNbins();
	 Double_t *x = new Double_t[nbins];
	 Double_t *y = new Double_t[nbins];
	 for (Int_t i=0;i<nbins;i++) {
	 x[i] = hist->GetXaxis()->GetBinCenter(i+1);
	 y[i] = hist->GetBinContent(i+1);
	 }
	 Double_t median = TMath::Median(nbins,x,y);

	return std::make_pair (widmin,median);

}


void effSigma()
{
	RooRealVar invMass("invMass","",91,0,200);
	RooRealVar m("mass","",91.188);
	RooRealVar g("width","",2);
	RooBreitWigner bw("bw","", invMass, m, g);
	RooDataSet *dataset = bw.generate(invMass, 1e6);

	std::pair <Double_t, Double_t> t = effSigma_function(dataset,invMass);

	std::cout << t.first << " "<<t.second << std::endl;

};


