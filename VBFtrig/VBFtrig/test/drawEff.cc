#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TPad.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
#include "tdrstyle.C"

using namespace std;

void drawEff(){
string rname = "trigtree.root"; // file name

gROOT -> Reset();
setTDRStyle();

// syntax: drive( variable name >> histo name (nbin, min, max), denominator cut, numerator cut, labels, plot name, legend, root file name) ;
drive("genmet>>hden(10,0,300)","ak5gendeta>3.5","trig1==1", "; gen MET true [GeV];Efficiency;","genmetEff.png", "MET > 140 Mjj > 600 HTT > 200", rname);
}

void drive(string sden, string dencut, string numcut, string caption, string fname, string leg, string rname){
TCanvas *c1 = new TCanvas("c1","multipads",800,600);
c1 -> SetFillColor (0);

TLegend *myleg = new TLegend(0.20,0.74,0.78,0.92);
myleg -> SetFillStyle(0);

cout << " Start " << endl;

TFile* ftt = new TFile(rname.c_str(),"READ");

cout << " File open ..." << endl;

TTree* ttt = (TTree *) (ftt -> Get("vbftrig/mytree"));

cout << " Tree open ..." << endl;

float ndata = ttt->Draw(sden.c_str(), dencut.c_str() );
TH1F* hden = (TH1F*)gDirectory->Get("hden");

// prepare numerator
string snum = sden;
string sfrom = "hden";
string sto   = "hnum";
size_t start_pos = snum.find(sfrom);
if(start_pos == std::string::npos) return;
else{
	snum.replace(start_pos, sfrom.length(), sto);
}

if(dencut != ""){
	if(numcut != "") numcut += " & " + dencut;
	else numcut = dencut ;
}
else numcut = dencut;

float ndata1 = ttt->Draw(snum.c_str(), numcut.c_str());
TH1F* hnum = (TH1F*)gDirectory->Get("hnum");

cout << " Plot 1 drawn ... " << endl;

hnum->Sumw2();
hden->Sumw2();

TGraphAsymmErrors * heff = new TGraphAsymmErrors(hnum,hden);

heff -> SetTitle(caption.c_str() );
heff -> SetMarkerColor(kBlack);
heff -> SetLineColor(kBlack);
heff -> SetLineWidth(2);
heff -> SetMarkerStyle(23);
heff -> SetMarkerSize(2);
heff -> SetMinimum(0);
heff -> SetMaximum(0.25);
heff -> Draw("ap");

myleg -> AddEntry( heff,     leg.c_str(),"LEP");
myleg -> Draw();

c1 -> SaveAs(fname.c_str());
ftt -> Close();

cout << " Macro ends ... " << endl;
}
