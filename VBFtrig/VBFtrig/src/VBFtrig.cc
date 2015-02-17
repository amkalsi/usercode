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
  std::pair <double, double> GenMuPtEta(edm::Handle<reco::GenParticleCollection>& genParticles, double);
  std::pair <double, double> GenTauPtEta(edm::Handle<reco::GenParticleCollection>& genParticles, double);
  std::pair <double, double> GenChi0PtEta(edm::Handle<reco::GenParticleCollection>& genParticles);
  std::pair <double, double> GenChiPmPtEta(edm::Handle<reco::GenParticleCollection>& genParticles);
  void Ak5GenMjjDeta(edm::Handle<reco::GenJetCollection>& genJets, double jetptcut, double detacut);
  double Ak5GenHT(edm::Handle<reco::GenJetCollection>& genJets, double jetptcut, double etacut);

  TTree* mytree;

    double genstablemet;

    double genht;
    double genmet;
    double genmjj;
    double gendeta;

    double genmettrue;
    double ak5genmjj ;
    double ak5gendeta ;
    double ak5genj1pt ;
    double ak5genj2pt ;
    double ak5genj1eta ;
    double ak5genj2eta ;
    double ak5genj1phi ;
    double ak5genj2phi ;
    int  ak5genj1hltmatch ;
    int  ak5genj2hltmatch ;

    double genmupt ;
    double genmueta ;
    double gentaupt ;
    double gentaueta ;

    double genchi0pt ;
    double genchi0eta ;
    double genchipmpt ;
    double genchipmeta ;


    double l1mupt ;
    double l3mupt ;
    double l3isomupt ;
    double calomet ;
    double calometID ;
    double ht ;
    double dijet ;
    double pfmet ;

bool isSignal ;
edm::InputTag triggerResultsTag_;
edm::InputTag genJetTag_;
edm::InputTag genMETTag_;

double jetEtaHT ;
double jetptcut ;
double detacut ;
double muptcut ;
double tauptcut ;
// for trigger
string processName;// = "HLT";
HLTConfigProvider hltConfig_;
unsigned int ntrig;


int trig1;
char ctrig1[128];
string strig1 ;
string sl1mu ;
string sl3mu ;
string sl3isomu ;
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
jetEtaHT(iConfig.getUntrackedParameter<double>("jetEtaHT", 0 )),
jetptcut(iConfig.getUntrackedParameter<double>("jetPtCut", 0 )),
detacut(iConfig.getUntrackedParameter<double>("dEtaCut", 0 )),
muptcut(iConfig.getUntrackedParameter<double>("muPtCut", 0 )),
tauptcut(iConfig.getUntrackedParameter<double>("tauPtCut", 0 )),
strig1(iConfig.getUntrackedParameter<string>("triggerName1", "none" )),
sl1mu(iConfig.getUntrackedParameter<string>("triggerL1MuName", "none" )),
sl3mu(iConfig.getUntrackedParameter<string>("triggerL3MuName", "none" )),
sl3isomu(iConfig.getUntrackedParameter<string>("triggerL3IsoMuName", "none" )),
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

 mytree->Branch("l1mupt",&l1mupt,"l1mupt/D");
 mytree->Branch("l3mupt",&l3mupt,"l3mupt/D");
 mytree->Branch("l3isomupt",&l3isomupt,"l3isomupt/D");
 mytree->Branch("calomet",&calomet,"calomet/D");
 mytree->Branch("calometID",&calometID,"calometID/D");
 mytree->Branch("ht",&ht,"ht/D");
 mytree->Branch("dijet",&dijet,"dijet/D");
 mytree->Branch("pfmet",&pfmet,"pfmet/D");

if(isSignal == true){
 mytree->Branch("genstablemet",&genstablemet,"genstablemet/D");
 mytree->Branch("genmupt",&genmupt,"genmupt/D");
 mytree->Branch("genmueta",&genmueta,"genmueta/D");
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

 mytree->Branch("ak5genmjj",&ak5genmjj,"ak5genmjj/D");
 mytree->Branch("ak5gendeta",&ak5gendeta,"ak5gendeta/D");
 mytree->Branch("ak5genj1pt",&ak5genj1pt,"ak5genj1pt/D");
 mytree->Branch("ak5genj2pt",&ak5genj2pt,"ak5genj2pt/D");
 mytree->Branch("ak5genj1eta",&ak5genj1eta,"ak5genj1eta/D");
 mytree->Branch("ak5genj2eta",&ak5genj2eta,"ak5genj2eta/D");
 mytree->Branch("ak5genj1phi",&ak5genj1phi,"ak5genj1phi/D");
 mytree->Branch("ak5genj2phi",&ak5genj2phi,"ak5genj2phi/D");
 mytree->Branch("ak5genj1hltmatch",&ak5genj1hltmatch,"ak5genj1hltmatch/I");
 mytree->Branch("ak5genj2hltmatch",&ak5genj2hltmatch,"ak5genj2hltmatch/I");
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

pair <double, double> VBFtrig::GenMuPtEta(edm::Handle<reco::GenParticleCollection>& genParticles, double ptcut){
        pair <double, double> mymupteta;
        mymupteta.first  = 0;
        mymupteta.second = 0;
        vector <unsigned int > vit;
        for( unsigned int i = 0; i < genParticles->size(); i++ ) {
                int id = abs(genParticles->at(i).pdgId() ); // child's abs id
                if(id == 13  && genParticles->at(i).status() == 1 ){ // stable
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
				while ( momid == 13 ){
                			const Candidate * granma = mom -> mother();
					mom = granma;
					momid = abs( mom -> pdgId() ) ;
					if ( momid != 13 ) break;
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


void VBFtrig::Ak5GenMjjDeta(edm::Handle<reco::GenJetCollection>& genJets, double jetptcut, double detacut){
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
								if( deleta > detacut && dijms > ak5genmjj) {
									ak5genmjj = dijms ;
									ak5gendeta = deleta;
									ak5genj1eta = genJets -> at(i).eta() ;
									ak5genj1phi = genJets -> at(i).phi() ;
									ak5genj1pt = genJets -> at(i).pt() ;
									ak5genj2eta = genJets -> at(j).eta() ;
									ak5genj2phi = genJets -> at(j).phi() ;
									ak5genj2pt = genJets -> at(j).pt() ;
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

double VBFtrig::Ak5GenHT(edm::Handle<reco::GenJetCollection>& genJets, double jetptcut, double etacut){
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

genmueta = 99;
genmupt = -1;
gentaueta = 99;
gentaupt = -1;

genchi0eta = 99;
genchi0pt = -1;
genchipmeta = 99;
genchipmpt = -1;

ak5genmjj = -1;
ak5gendeta = -1;
ak5genj1pt = -1;
ak5genj2pt = -1;
ak5genj1eta = 99;
ak5genj2eta = 99;
ak5genj1phi = 9;
ak5genj2phi = 9;
ak5genj1hltmatch = 0;
ak5genj2hltmatch = 0;
}

l1mupt = -1;
l3mupt = -1;
l3isomupt = -1;
calomet = -1;
calometID = -1;
ht = -1;
dijet = -1;
pfmet = -1;

bool fl1mu = false;
bool fl3mu = false;
bool fl3isomu = false;
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
	Ak5GenMjjDeta( genJets, jetptcut, detacut);
	genht = Ak5GenHT( genJets, jetptcut, jetEtaHT);
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
        pair <double, double> pmupteta = GenMuPtEta(genParticles, muptcut);
        genmupt  = pmupteta.first;
        genmueta = pmupteta.second;
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
				fl1mu = false;
				fl3mu = false;
				fl3isomu = false;
				fcalomet = false;
				fcalometID = false;
				fht = false;
				fdijet = false;
				fpfmet = false;
				for(unsigned int l = 0; l < moduleLabels.size(); l++){
					moduleLabel = moduleLabels.at(l);
					if(moduleLabel.find(sl1mu) != std::string::npos) fl1mu = true;
					if(moduleLabel.find(sl3mu) != std::string::npos) fl3mu = true;
					if(moduleLabel.find(sl3isomu) != std::string::npos) fl3isomu = true;
					if(moduleLabel.find(scalomet) != std::string::npos) fcalomet = true;
					if(moduleLabel.find(scalometID) != std::string::npos) fcalometID = true;
					if(moduleLabel.find(sht) != std::string::npos) fht = true;
					if(moduleLabel.find(sdijet) != std::string::npos) fdijet = true;
					if(moduleLabel.find(spfmet) != std::string::npos) fpfmet = true;
				}
				if(fl1mu == true ){
					const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(sl1mu,"",processName)));
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
							if( l1mupt < TO.pt() ) l1mupt = TO.pt() ;
        					}
					} // filer size
				}
				if(fl3mu == true ){
					const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(sl3mu,"",processName)));
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
							if( l3mupt < TO.pt() ) l3mupt = TO.pt() ;
        					}
					} // filer size
				}
				if(fl3isomu == true ){
					const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(sl3isomu,"",processName)));
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
							if( l3isomupt < TO.pt() ) l3isomupt = TO.pt() ;
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
	if( ak5genj1hltmatch == 0 && ak5genj1eta != 99 && ak5genj1phi != 9){
		double detajj = ak5genj1eta - TO1.eta();
		double dphijj = ak5genj1phi - TO1.phi();
		double drjj = sqrt( detajj*detajj + dphijj*dphijj);
		if( drjj < 0.5) ak5genj1hltmatch = 1;
	}
	if( ak5genj2hltmatch == 0 && ak5genj2eta != 99 && ak5genj2phi != 9){
		double detajj = ak5genj2eta - TO1.eta();
		double dphijj = ak5genj2phi - TO1.phi();
		double drjj = sqrt( detajj*detajj + dphijj*dphijj);
		if( drjj < 0.5) ak5genj2hltmatch = 1;
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
