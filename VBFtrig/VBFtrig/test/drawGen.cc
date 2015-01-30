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

void drawGen(){
string rname = "trigtree.root"; // file name

gROOT -> Reset();
setTDRStyle();

// syntax: drive( variable name >> histo name (nbin, min, max), cut, labels, plot name, root file name) ;
drive("ak5genj1pt>>hden(25,0,300)","ak5genj1pt>0", "; highest p_{T} gen AK5 Jet p_{T} [GeV];a.u.;", "jetPt1.png", rname);
}

void drive(string hden1, string dencut1, string caption, string fname, string rname ){
TCanvas *c1 = new TCanvas("c1","multipads",800,600);
c1 -> SetFillColor (0);

cout << " Start " << endl;

TFile* ftt = new TFile(rname.c_str(),"READ");

cout << " File open ..." << endl;

TTree* ttt = (TTree *) (ftt -> Get("vbftrig/mytree"));

cout << " Tree open ..." << endl;

float ndata = ttt->Draw(hden1.c_str(), dencut1.c_str() );
TH1F* hden = (TH1F*)gDirectory->Get("hden");

cout << " Plot 1 drawn ... with nevents = " << ndata << endl;

hden -> SetTitle(caption.c_str() );
hden->Sumw2();

hden -> SetMarkerColor(kBlack);
hden -> SetLineColor(kBlack);
hden -> SetLineWidth(2);
hden -> SetMarkerStyle(23);
hden -> SetMarkerSize(2);

hden -> DrawNormalized("ep");

c1 -> SaveAs(fname.c_str());
ftt -> Close();

cout << " Macro ends ... " << endl;
delete c1;
}
