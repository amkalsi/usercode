// ROOT specific
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"
// C files
#include <string>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <regex.h>
#include <stdio.h>
// Lorentz Vector in CMSSW
#include "DataFormats/Math/interface/LorentzVectorFwd.h"
// CMSSW
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
// MC
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/GenJet.h"
// for trigger
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "HLTrigger/Muon/interface/HLTMuonIsoFilter.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

using namespace edm;
using namespace std;
using namespace reco;
using namespace trigger;

class VBFtrig : public edm::EDAnalyzer {
public:
  explicit VBFtrig(const edm::ParameterSet&);
  ~VBFtrig();
private:

  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
// for trigger
  virtual void beginRun  (edm::Run const & iRun, edm::EventSetup const& iSetup);

  double GenStableMET(edm::Handle<reco::GenParticleCollection>& genParticles);
  double GenMET(edm::Handle<reco::GenParticleCollection>& genParticles);
  double GenHT(edm::Handle<reco::GenParticleCollection>& genParticles);
  std::pair <double, double> GenMjjDeta(edm::Handle<reco::GenParticleCollection>& genParticles, double, double);
  std::pair <double, double> GenPtEta(edm::Handle<reco::GenParticleCollection>& genParticles, double, int);
  std::pair <double, double> GenTauPtEta(edm::Handle<reco::GenParticleCollection>& genParticles, double);
  std::pair <double, double> GenChi0PtEta(edm::Handle<reco::GenParticleCollection>& genParticles);
  std::pair <double, double> GenChiPmPtEta(edm::Handle<reco::GenParticleCollection>& genParticles);
  void Ak4GenMjjDeta(edm::Handle<reco::GenJetCollection>& genJets, double jetptcut, double detacut);
  double Ak4GenHT(edm::Handle<reco::GenJetCollection>& genJets, double jetptcut, double etacut);

  TTree* mytree;

    double genstablemet;

    double genht;
    double genmet;
    double genmjj;
    double gendeta;

    double genmettrue;
    double ak4genmjj ;
    double ak4gendeta ;
    double ak4genj1pt ;
    double ak4genj2pt ;
    double ak4genj1eta ;
    double ak4genj2eta ;
    double ak4genj1phi ;
    double ak4genj2phi ;
    int  ak4genj1hltmatch ;
    int  ak4genj2hltmatch ;

    double genleppt ;
    double genlepeta ;
    double gentaupt ;
    double gentaueta ;

    double genchi0pt ;
    double genchi0eta ;
    double genchipmpt ;
    double genchipmeta ;


    double l1pt ;
    double l3leppt1 ;
    double l3leppt2 ;
    double l3leppt3 ;

    double calomet ;
    double calometID ;
    double ht ;
    double dijet ;
    double pfmet ;

bool isSignal ;
edm::InputTag triggerResultsTag_;
edm::InputTag genJetTag_;
edm::InputTag genMETTag_;

int leptonPid ;

double jetEtaHT ;
double jetptcut ;
double detacut ;
double lepptcut ;
double tauptcut ;
// for trigger
string processName;// = "HLT";
HLTConfigProvider hltConfig_;
unsigned int ntrig;


int trig1;
char ctrig1[128];
string strig1 ;
string sl1name ;
string sl3lep1 ;
string sl3lep2 ;
string sl3lep3 ;

string scalomet ;
string scalometID ;
string sht ;
string sdijet ;
string spfmet ;

};

VBFtrig::VBFtrig(const edm::ParameterSet& iConfig):
isSignal(iConfig.getUntrackedParameter<bool>("isSignal", false  )),
triggerResultsTag_(iConfig.getParameter<edm::InputTag>("triggerResults")) ,
genJetTag_(iConfig.getParameter<edm::InputTag>("genJetTag")),
genMETTag_(iConfig.getParameter<edm::InputTag>("genMETTag")),
leptonPid(iConfig.getUntrackedParameter<int>("leptonPid", 0 )),
jetEtaHT(iConfig.getUntrackedParameter<double>("jetEtaHT", 0 )),
jetptcut(iConfig.getUntrackedParameter<double>("jetPtCut", 0 )),
detacut(iConfig.getUntrackedParameter<double>("dEtaCut", 0 )),
lepptcut(iConfig.getUntrackedParameter<double>("lepPtCut", 0 )),
tauptcut(iConfig.getUntrackedParameter<double>("tauPtCut", 0 )),
strig1(iConfig.getUntrackedParameter<string>("triggerName1", "none" )),
sl1name(iConfig.getUntrackedParameter<string>("triggerL1Name", "none" )),
sl3lep1(iConfig.getUntrackedParameter<string>("triggerL3LeptonName1", "none" )),
sl3lep2(iConfig.getUntrackedParameter<string>("triggerL3LeptonName2", "none" )),
sl3lep3(iConfig.getUntrackedParameter<string>("triggerL3LeptonName3", "none" )),
scalomet(iConfig.getUntrackedParameter<string>("triggerCaloMETName", "none" )),
scalometID(iConfig.getUntrackedParameter<string>("triggerCaloMETIDName", "none" )),
sht(iConfig.getUntrackedParameter<string>("triggerHTName", "none" )),
sdijet(iConfig.getUntrackedParameter<string>("triggerDiJetName", "none" )),
spfmet(iConfig.getUntrackedParameter<string>("triggerPFMETName", "none" ))
{
processName = "TEST";

  edm::Service<TFileService> fs;
  mytree = fs->make<TTree>("mytree", "TTree to store Ntuple");

    mytree->Branch("trig1",    &trig1,  "trig1/I");
    mytree->Branch("ctrig1",    &ctrig1,  "ctrig1/C");

 mytree->Branch("l1pt",&l1pt,"l1pt/D");
 mytree->Branch("l3leppt1",&l3leppt1,"l3leppt1/D");
 mytree->Branch("l3leppt2",&l3leppt2,"l3leppt2/D");
 mytree->Branch("l3leppt3",&l3leppt3,"l3leppt3/D");

 mytree->Branch("calomet",&calomet,"calomet/D");
 mytree->Branch("calometID",&calometID,"calometID/D");
 mytree->Branch("ht",&ht,"ht/D");
 mytree->Branch("dijet",&dijet,"dijet/D");
 mytree->Branch("pfmet",&pfmet,"pfmet/D");

if(isSignal == true){
 mytree->Branch("genstablemet",&genstablemet,"genstablemet/D");
 mytree->Branch("genleppt",&genleppt,"genleppt/D");
 mytree->Branch("genlepeta",&genlepeta,"genlepeta/D");
 mytree->Branch("gentaupt",&gentaupt,"gentaupt/D");
 mytree->Branch("gentaueta",&gentaueta,"gentaueta/D");

 mytree->Branch("genchi0pt",&genchi0pt,"genchi0pt/D");
 mytree->Branch("genchi0eta",&genchi0eta,"genchi0eta/D");
 mytree->Branch("genchipmpt",&genchipmpt,"genjchipmpt/D");
 mytree->Branch("genchipmeta",&genchipmeta,"genjchipmeta/D");

 mytree->Branch("genht",&genht,"genht/D");
 mytree->Branch("genmet",&genmet,"genmet/D");
 mytree->Branch("genmjj",&genmjj,"genmjj/D");
 mytree->Branch("gendeta",&gendeta,"gendeta/D");
 mytree->Branch("genmettrue",&genmettrue,"genmettrue/D");

 mytree->Branch("ak4genmjj",&ak4genmjj,"ak4genmjj/D");
 mytree->Branch("ak4gendeta",&ak4gendeta,"ak4gendeta/D");
 mytree->Branch("ak4genj1pt",&ak4genj1pt,"ak4genj1pt/D");
 mytree->Branch("ak4genj2pt",&ak4genj2pt,"ak4genj2pt/D");
 mytree->Branch("ak4genj1eta",&ak4genj1eta,"ak4genj1eta/D");
 mytree->Branch("ak4genj2eta",&ak4genj2eta,"ak4genj2eta/D");
 mytree->Branch("ak4genj1phi",&ak4genj1phi,"ak4genj1phi/D");
 mytree->Branch("ak4genj2phi",&ak4genj2phi,"ak4genj2phi/D");
 mytree->Branch("ak4genj1hltmatch",&ak4genj1hltmatch,"ak4genj1hltmatch/I");
 mytree->Branch("ak4genj2hltmatch",&ak4genj2hltmatch,"ak4genj2hltmatch/I");
}

}

VBFtrig::~VBFtrig()
{
}

// for trigger
void VBFtrig::beginRun  (edm::Run const & iRun, edm::EventSetup const& iSetup){
        bool changed(true);
        hltConfig_.init(iRun,iSetup,processName,changed);
        ntrig = hltConfig_.size();
        hltConfig_.dump("Triggers");
}

double VBFtrig::GenStableMET(edm::Handle<reco::GenParticleCollection>& genParticles){
        double tempmet = 0;
        vector <unsigned int > vit;
        for( unsigned int i = 0; i < genParticles->size(); i++ ) {
                int id = abs(genParticles->at(i).pdgId() ); // child's abs id
                if(id != 1000022 && id != 1000023 && id != 12 && id != 14 && id != 16  ){
                        if( genParticles->at(i).status() == 1 ) vit.push_back( i );
                }
        }

        double tpx = 0;
        double tpy = 0;
        double tpz = 0;
        double ten = 0;
        for( unsigned int i = 0; i < vit.size(); i++){
                tpx +=  genParticles->at(vit.at(i) ) . px();
                tpy +=  genParticles->at(vit.at(i) ) . py();
                tpz +=  genParticles->at(vit.at(i) ) . pz();
                ten +=  genParticles->at(vit.at(i) ) . energy();
        }
        TLorentzVector lmet (tpx, tpy, tpz, ten);
        tempmet = lmet.Pt();
        return tempmet;
}

double VBFtrig::GenHT(edm::Handle<reco::GenParticleCollection>& genParticles){
        double tempht = 0;
        for( unsigned int i = 0; i < genParticles->size(); i++ ) {
                int id = abs(genParticles->at(i).pdgId() ); // child's abs id
                if(id >= 1 && id <= 4 ){
                        if( genParticles->at(i).status() == 1 ){
				if( genParticles->at(i).pt() > 10 && fabs(genParticles->at(i).eta() ) < 3.0 ){
					tempht += genParticles->at(i).pt() ;
				}
			}
                }
        }
        return tempht;
}

pair <double, double> VBFtrig::GenChiPmPtEta(edm::Handle<reco::GenParticleCollection>& genParticles){
        pair <double, double> mypteta;
        mypteta.first  = 0;
        mypteta.second = 0;
        int  uit = -1;
        for( unsigned int i = 0; i < genParticles->size(); i++ ) {
                int id = abs(genParticles->at(i).pdgId() ); // child's abs id
                if(id == 1000024 ){ 
                        uit = (int) i ;
			break;
                }
        }
        double tpx = 0;
        double tpy = 0;
        double tpz = 0;
        double ten = 0;
        if( uit != -1){
                tpx =  genParticles->at(uit ) . px();
                tpy =  genParticles->at(uit ) . py();
                tpz =  genParticles->at(uit ) . pz();
                ten =  genParticles->at(uit ) . energy();
        	TLorentzVector lchi (tpx, tpy, tpz, ten);
	        mypteta.first  = lchi.Pt();
        	mypteta.second = lchi.Eta();
        }
        return mypteta;
}

pair <double, double> VBFtrig::GenChi0PtEta(edm::Handle<reco::GenParticleCollection>& genParticles){
        pair <double, double> mypteta;
        mypteta.first  = 0;
        mypteta.second = 0;
        int  uit = -1;
        for( unsigned int i = 0; i < genParticles->size(); i++ ) {
                int id = abs(genParticles->at(i).pdgId() ); // child's abs id
                if(id == 1000022 ){
                        uit = (int) i ;
                        break;
                }
        }
        double tpx = 0;
        double tpy = 0;
        double tpz = 0;
        double ten = 0;
        if( uit != -1){
                tpx =  genParticles->at(uit ) . px();
                tpy =  genParticles->at(uit ) . py();
                tpz =  genParticles->at(uit ) . pz();
                ten =  genParticles->at(uit ) . energy();
                TLorentzVector lchi (tpx, tpy, tpz, ten);
                mypteta.first  = lchi.Pt();
                mypteta.second = lchi.Eta();
        }
        return mypteta;
}

double VBFtrig::GenMET(edm::Handle<reco::GenParticleCollection>& genParticles){
        double tempmet = 0;
        vector <unsigned int > vit;
        for( unsigned int i = 0; i < genParticles->size(); i++ ) {
                int id = abs(genParticles->at(i).pdgId() ); // child's abs id
                if(id == 1000022 || id == 12 || id == 14 || id == 16 ){ 
                        if( genParticles->at(i).status() == 1 ) vit.push_back( i );
                }
        }

        double tpx = 0;
        double tpy = 0;
        double tpz = 0;
        double ten = 0;
        for( unsigned int i = 0; i < vit.size(); i++){
                tpx +=  genParticles->at(vit.at(i) ) . px();
                tpy +=  genParticles->at(vit.at(i) ) . py();
                tpz +=  genParticles->at(vit.at(i) ) . pz();
                ten +=  genParticles->at(vit.at(i) ) . energy();
        }
        TLorentzVector lmet (tpx, tpy, tpz, ten);
        tempmet = lmet.Pt();
        return tempmet;
}

pair <double, double> VBFtrig::GenMjjDeta(edm::Handle<reco::GenParticleCollection>& genParticles, double ptcut, double detacut ){
        pair <double, double> mymjjdeta;
        mymjjdeta.first  = 0;
        mymjjdeta.second = 0;
        double tempmjj = 0;
        double tempdeta = 0;
        vector <unsigned int > vit;
        for( unsigned int i = 0; i < genParticles->size(); i++ ) {
                int id = abs(genParticles->at(i).pdgId() ); // child's abs id
                if(id >= 1 && id <= 4  ){ // light flavor quarks only (no gluon)
                        if( fabs( genParticles->at( i ) . eta() ) < 5.0){
                                if( genParticles->at( i ) . pt() > ptcut )
                                        vit.push_back( i );
                        }
                }
        }
        for( unsigned int j = 0; j < vit.size(); j++){
                for( unsigned int k = j; k < vit.size(); k++){
                        if( j < k && deltaR(genParticles->at( vit.at(j) ), genParticles->at(vit.at(k)) ) > 0.5){
                                double tpx =  (genParticles->at(vit.at(j) ) . px() + genParticles->at(vit.at(k) ) . px());
                                double tpy =  (genParticles->at(vit.at(j) ) . py() + genParticles->at(vit.at(k) ) . py());
                                double tpz =  (genParticles->at(vit.at(j) ) . pz() + genParticles->at(vit.at(k) ) . pz());
                                double ten =  (genParticles->at(vit.at(j) ) . energy() + genParticles->at(vit.at(k) ) . energy());
                                TLorentzVector ljj (tpx, tpy, tpz, ten);
                                tempmjj = ljj.M();
                                tempdeta = fabs( genParticles->at(vit.at(j) ) . eta() - genParticles->at(vit.at(k) ) . eta());
                                if( genParticles->at(vit.at(j) ) . eta() * genParticles->at(vit.at(k) ) . eta() < 0)
                                if(tempdeta > detacut && tempmjj > mymjjdeta.first){
                                        mymjjdeta.first = tempmjj;
                                        mymjjdeta.second = tempdeta;
                                }
                        }
                }
        }
        return mymjjdeta;
}

pair <double, double> VBFtrig::GenPtEta(edm::Handle<reco::GenParticleCollection>& genParticles, double ptcut, int thispid){
        pair <double, double> mymupteta;
        mymupteta.first  = 0;
        mymupteta.second = 0;
        vector <unsigned int > vit;
        for( unsigned int i = 0; i < genParticles->size(); i++ ) {
                int id = abs(genParticles->at(i).pdgId() ); // child's abs id
                if(id == thispid  && genParticles->at(i).status() == 1 ){ // stable
                        if( fabs( genParticles->at( i ) . eta() ) < 2.5){
                                if( genParticles->at( i ) . pt() > ptcut ){
                                        vit.push_back( i );
                                }
                        }
                }
        }
        for( unsigned int j = 0; j < vit.size(); j++){
                const Candidate * mom = genParticles->at( vit.at(j) ).mother();
		if( mom != NULL){
			int momid = abs( mom -> pdgId()) ;
			if ( abs( mom -> pdgId()) != 15 ){ 
				while ( momid == thispid ){
                			const Candidate * granma = mom -> mother();
					mom = granma;
					momid = abs( mom -> pdgId() ) ;
					if ( momid != thispid ) break;
				}
			}
                        if( momid == 15  ){ // tau ancestor
				if( genParticles->at( vit.at(j) ) . pt() > mymupteta.first ){
                                        mymupteta.first  =  genParticles->at( vit.at(j) ) . pt();
                                        mymupteta.second =  genParticles->at( vit.at(j) ) . eta() ;
                                }
                        }
                }
        }
        return mymupteta;
}

pair <double, double> VBFtrig::GenTauPtEta(edm::Handle<reco::GenParticleCollection>& genParticles, double ptcut){
        pair <double, double> mytaupteta;
        mytaupteta.first  = 0;
        mytaupteta.second = 0;
        for( unsigned int i = 0; i < genParticles->size(); i++ ) {
                int id = abs(genParticles->at(i).pdgId() ); // child's abs id
                if(id == 15 ){ 
                        if( fabs( genParticles->at( i ) . eta() ) < 2.5){
                                if( genParticles->at( i ) . pt() > 10 ){
                                        if( genParticles->at( i ) . pt() > mytaupteta.first ){
                                                mytaupteta.first =  genParticles->at( i ) . pt() ;
                                                mytaupteta.second =  genParticles->at( i ) .  eta() ;
                                        }
                                }
                        }
                }
        }
        return mytaupteta;
}


void VBFtrig::Ak4GenMjjDeta(edm::Handle<reco::GenJetCollection>& genJets, double jetptcut, double detacut){
  for( unsigned int i = 0; i < genJets -> size(); i++){
	if( fabs(genJets -> at(i).eta()) < 5.0 ){
		if( genJets -> at(i).pt() > jetptcut ){
  			for( unsigned int j = i; j < genJets -> size(); j++){
				if( i < j){
					if( fabs(genJets -> at(j).eta()) < 5.0 ){
						if( genJets -> at(j).pt() > jetptcut ){
							if( genJets -> at(i).eta() *  genJets -> at(j).eta()  < 0 ){
								double deleta = fabs(  genJets -> at(i).eta() - genJets -> at(j).eta());
								double dijms = (  genJets -> at(i).p4() + genJets -> at(j).p4()).mass();
								if( deleta > detacut && dijms > ak4genmjj) {
									ak4genmjj = dijms ;
									ak4gendeta = deleta;
									ak4genj1eta = genJets -> at(i).eta() ;
									ak4genj1phi = genJets -> at(i).phi() ;
									ak4genj1pt = genJets -> at(i).pt() ;
									ak4genj2eta = genJets -> at(j).eta() ;
									ak4genj2phi = genJets -> at(j).phi() ;
									ak4genj2pt = genJets -> at(j).pt() ;
								}
							}
						}
					}
				}
			}
		}
	}
  }
}

double VBFtrig::Ak4GenHT(edm::Handle<reco::GenJetCollection>& genJets, double jetptcut, double etacut){
	double tempht = 0;
	for( unsigned int i = 0; i < genJets -> size(); i++){
        	if( fabs(genJets -> at(i).eta()) < etacut ){
                	if( genJets -> at(i).pt() > jetptcut ){
                		tempht += genJets -> at(i).pt() ;
			}
		}
	}
	return tempht;
}

void
VBFtrig::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

if(isSignal){
genstablemet = -1;

genht = -1;
genmet = -1;
genmjj = -1;
gendeta = 99;
genmettrue = -1;

genlepeta = 99;
genleppt = -1;
gentaueta = 99;
gentaupt = -1;

genchi0eta = 99;
genchi0pt = -1;
genchipmeta = 99;
genchipmpt = -1;

ak4genmjj = -1;
ak4gendeta = -1;
ak4genj1pt = -1;
ak4genj2pt = -1;
ak4genj1eta = 99;
ak4genj2eta = 99;
ak4genj1phi = 9;
ak4genj2phi = 9;
ak4genj1hltmatch = 0;
ak4genj2hltmatch = 0;
}

l1pt = -1;
l3leppt1 = -1;
l3leppt2 = -1;
l3leppt3 = -1;

calomet = -1;
calometID = -1;
ht = -1;
dijet = -1;
pfmet = -1;

bool fl1 = false;
bool fl3lep1 = false;
bool fl3lep2 = false;
bool fl3lep3 = false;

bool fcalomet = false;
bool fcalometID = false;
bool fht = false;
bool fdijet = false;
bool fpfmet = false;

trig1 = -1;
snprintf(ctrig1, 127, "none");

// first fill gen jets for later use
if(isSignal){
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel(genJetTag_, genJets);
  if( genJets.isValid() ){
	Ak4GenMjjDeta( genJets, jetptcut, detacut);
	genht = Ak4GenHT( genJets, jetptcut, jetEtaHT);
  }
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  if(genParticles.isValid() ){
        genstablemet = GenStableMET(genParticles);
        //genht = GenHT(genParticles);
        genmet = GenMET(genParticles);
        pair <double, double> pmjjdeta = GenMjjDeta(genParticles, jetptcut, detacut);
        genmjj  = pmjjdeta.first;
        gendeta  = pmjjdeta.second;
        pair <double, double> pmupteta = GenPtEta(genParticles, lepptcut, leptonPid);
        genleppt  = pmupteta.first;
        genlepeta = pmupteta.second;
        pair <double, double> ptaupteta = GenTauPtEta(genParticles, tauptcut);
        gentaupt  = ptaupteta.first;
        gentaueta = ptaupteta.second;
        pair <double, double> pchi0pteta = GenChi0PtEta(genParticles);
        genchi0pt  = pchi0pteta.first;
        genchi0eta = pchi0pteta.second;
        pair <double, double> pchipmpteta = GenChiPmPtEta(genParticles);
        genchipmpt  = pchipmpteta.first;
        genchipmeta = pchipmpteta.second;
  }
  edm::Handle< vector< reco::GenMET> > genMETs;
  iEvent.getByLabel(genMETTag_, genMETs);
  if( genMETs.isValid() ){
	genmettrue = genMETs -> at(0).et();
  }
}

  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  edm::Handle<trigger::TriggerEvent>           triggerEventHandle_;
  iEvent.getByLabel("hltTriggerSummaryAOD", triggerEventHandle_);
  if(triggerResultsHandle_.isValid() == true ){
  	const edm::TriggerNames vTrigNames = iEvent.triggerNames(*triggerResultsHandle_);
  	const unsigned int ntrigs = triggerResultsHandle_->size();

        for(unsigned int u = 0; u < hltConfig_.triggerNames().size(); u++){
        	const unsigned int triggerIndex(hltConfig_.triggerIndex(hltConfig_.triggerNames().at(u)));
	        if( (ntrigs > triggerIndex) ){
                	const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));
	                string moduleLabel = "";
			if (hltConfig_.triggerNames().at(u).find(strig1) != std::string::npos ){ // mu path
                	        if(triggerResultsHandle_->accept(triggerIndex) ){ 
					trig1 = 1;
                        		snprintf(ctrig1, 127, "%s", hltConfig_.triggerNames().at(u).c_str() );
				}
				fl1 = false;
				fl3lep1 = false;
				fl3lep2 = false;
				fl3lep3 = false;

				fcalomet = false;
				fcalometID = false;
				fht = false;
				fdijet = false;
				fpfmet = false;
				for(unsigned int l = 0; l < moduleLabels.size(); l++){
					moduleLabel = moduleLabels.at(l);
					if(moduleLabel.find(sl1name) != std::string::npos) fl1 = true;
					if(moduleLabel.find(sl3lep1) != std::string::npos) fl3lep1 = true;
					if(moduleLabel.find(sl3lep2) != std::string::npos) fl3lep2 = true;
					if(moduleLabel.find(sl3lep3) != std::string::npos) fl3lep3 = true;

					if(moduleLabel.find(scalomet) != std::string::npos) fcalomet = true;
					if(moduleLabel.find(scalometID) != std::string::npos) fcalometID = true;
					if(moduleLabel.find(sht) != std::string::npos) fht = true;
					if(moduleLabel.find(sdijet) != std::string::npos) fdijet = true;
					if(moduleLabel.find(spfmet) != std::string::npos) fpfmet = true;
				}
				if(fl1 == true ){
					const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(sl1name,"",processName)));
				        if (filterIndex < triggerEventHandle_->sizeFilters() ){
					        const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
			        		const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
				        	const size_type nI(VIDS.size());
					        const size_type nK(KEYS.size());
					        assert(nI==nK);
        					const size_type n2(max(nI,nK));
				        	const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
					        for (size_type i=0; i!=n2; ++i) {
				        	        const TriggerObject& TO(TOC[KEYS[i]]);
							if( l1pt < TO.pt() ) l1pt = TO.pt() ;
        					}
					} // filer size
				}
				if(fl3lep1 == true ){
					const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(sl3lep1,"",processName)));
				        if (filterIndex < triggerEventHandle_->sizeFilters() ){
					        const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
			        		const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
				        	const size_type nI(VIDS.size());
					        const size_type nK(KEYS.size());
					        assert(nI==nK);
        					const size_type n2(max(nI,nK));
				        	const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
					        for (size_type i=0; i!=n2; ++i) {
				        	        const TriggerObject& TO(TOC[KEYS[i]]);
							if( l3leppt1 < TO.pt() ) l3leppt1 = TO.pt() ;
        					}
					} // filer size
				}

				if(fl3lep2 == true ){
					const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(sl3lep2,"",processName)));
				        if (filterIndex < triggerEventHandle_->sizeFilters() ){
					        const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
			        		const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
				        	const size_type nI(VIDS.size());
					        const size_type nK(KEYS.size());
					        assert(nI==nK);
        					const size_type n2(max(nI,nK));
				        	const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
					        for (size_type i=0; i!=n2; ++i) {
				        	        const TriggerObject& TO(TOC[KEYS[i]]);
							if( l3leppt2 < TO.pt() ) l3leppt2 = TO.pt() ;
        					}
					} // filer size
				}
				if(fl3lep3 == true ){
					const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(sl3lep3,"",processName)));
				        if (filterIndex < triggerEventHandle_->sizeFilters() ){
					        const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
			        		const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
				        	const size_type nI(VIDS.size());
					        const size_type nK(KEYS.size());
					        assert(nI==nK);
        					const size_type n2(max(nI,nK));
				        	const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
					        for (size_type i=0; i!=n2; ++i) {
				        	        const TriggerObject& TO(TOC[KEYS[i]]);
							if( l3leppt3 < TO.pt() ) l3leppt3 = TO.pt() ;
        					}
					} // filer size
				}

				if(fcalomet == true ){
					const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(scalomet,"",processName)));
				        if (filterIndex < triggerEventHandle_->sizeFilters() ){
					        const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
			        		const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
				        	const size_type nI(VIDS.size());
					        const size_type nK(KEYS.size());
					        assert(nI==nK);
        					const size_type n2(max(nI,nK));
				        	const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
					        for (size_type i=0; i!=n2; ++i) {
				        	        const TriggerObject& TO(TOC[KEYS[i]]);
							if( calomet < TO.pt() ) calomet = TO.pt() ;
        					}
					} // filer size
				}
				if(fcalometID == true ){
					const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(scalometID,"",processName)));
				        if (filterIndex < triggerEventHandle_->sizeFilters() ){
					        const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
			        		const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
				        	const size_type nI(VIDS.size());
					        const size_type nK(KEYS.size());
					        assert(nI==nK);
        					const size_type n2(max(nI,nK));
				        	const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
					        for (size_type i=0; i!=n2; ++i) {
				        	        const TriggerObject& TO(TOC[KEYS[i]]);
							if( calometID < TO.pt() ) calometID = TO.pt() ;
        					}
					} // filer size
				}
				if(fht == true ){
					const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(sht,"",processName)));
				        if (filterIndex < triggerEventHandle_->sizeFilters() ){
					        const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
			        		const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
				        	const size_type nI(VIDS.size());
					        const size_type nK(KEYS.size());
					        assert(nI==nK);
        					const size_type n2(max(nI,nK));
				        	const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
					        for (size_type i=0; i!=n2; ++i) {
				        	        const TriggerObject& TO(TOC[KEYS[i]]);
							if( ht < TO.pt() ) ht = TO.pt() ;
        					}
					} // filer size
				}
				if(fdijet == true ){
					const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(sdijet,"",processName)));
				        if (filterIndex < triggerEventHandle_->sizeFilters() ){
					        const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
			        		const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
				        	const size_type nI(VIDS.size());
					        const size_type nK(KEYS.size());
					        assert(nI==nK);
        					const size_type n2(max(nI,nK));
				        	const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
					        for (size_type i=0; i < n2 ; i++ ) {
				        	        const TriggerObject& TO1(TOC[KEYS[i]]);
							if( TO1.pt() > jetptcut ){
if(isSignal){
// now match genjets to hlt jets
	if( ak4genj1hltmatch == 0 && ak4genj1eta != 99 && ak4genj1phi != 9){
		double detajj = ak4genj1eta - TO1.eta();
		double dphijj = ak4genj1phi - TO1.phi();
		double drjj = sqrt( detajj*detajj + dphijj*dphijj);
		if( drjj < 0.5) ak4genj1hltmatch = 1;
	}
	if( ak4genj2hltmatch == 0 && ak4genj2eta != 99 && ak4genj2phi != 9){
		double detajj = ak4genj2eta - TO1.eta();
		double dphijj = ak4genj2phi - TO1.phi();
		double drjj = sqrt( detajj*detajj + dphijj*dphijj);
		if( drjj < 0.5) ak4genj2hltmatch = 1;
	}
// 
}
						        	for (size_type j=i; j < n2; j++) {
									if( i < j ){
						        	        	const TriggerObject& TO2(TOC[KEYS[j]]);
										if( TO2.pt() > jetptcut ){
											double tmpdeta = fabs( TO1.eta() - TO2.eta() );
											double tmpopposite = TO1.eta() * TO2.eta() ;
											if( tmpdeta > detacut && tmpopposite < 0){
												TLorentzVector j1 ( TO1.px(),  TO1.py(),  TO1.pz(),  TO1.energy());
												TLorentzVector j2 ( TO2.px(),  TO2.py(),  TO2.pz(),  TO2.energy());
												double tmpmass = ( j1 + j2 ).M();
												if( dijet < tmpmass ) dijet = tmpmass ;
        										}
        									}
        								}
        							}
        						}
        					}
					} // filer size
				}
				if(fpfmet == true ){
					const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(spfmet,"",processName)));
				        if (filterIndex < triggerEventHandle_->sizeFilters() ){
					        const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
			        		const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
				        	const size_type nI(VIDS.size());
					        const size_type nK(KEYS.size());
					        assert(nI==nK);
        					const size_type n2(max(nI,nK));
				        	const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
					        for (size_type i=0; i!=n2; ++i) {
				        	        const TriggerObject& TO(TOC[KEYS[i]]);
							if( pfmet < TO.pt() ) pfmet = TO.pt() ;
        					}
					} // filer size
				}
			}
		}
	}// trig path
  }

  mytree->Fill(); // fill even if no tri-lepton
}// ana ends

void 
VBFtrig::beginJob()
{
}
void 
VBFtrig::endJob() 
{
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(VBFtrig);
