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
//#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "DataFormats/Common/interface/TriggerResults.h"
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
  std::pair <double, double> GenMjjDeta(edm::Handle<reco::GenParticleCollection>& genParticles, double, double);
  std::pair <double, double> GenMuPtEta(edm::Handle<reco::GenParticleCollection>& genParticles, double);
  std::pair <double, double> GenTauPtEta(edm::Handle<reco::GenParticleCollection>& genParticles, double);
  std::pair <double, double> GenChi0PtEta(edm::Handle<reco::GenParticleCollection>& genParticles);
  std::pair <double, double> GenChiPmPtEta(edm::Handle<reco::GenParticleCollection>& genParticles);
  void Ak5GenMjjDeta(edm::Handle<reco::GenJetCollection>& genJets, double jetptcut, double detacut);

  TTree* mytree;

    double genstablemet;

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

    double genmupt ;
    double genmueta ;
    double gentaupt ;
    double gentaueta ;

    double genchi0pt ;
    double genchi0eta ;
    double genchipmpt ;
    double genchipmeta ;


double jetptcut ;
double detacut ;
double muptcut ;
double tauptcut ;
// for trigger
string processName;// = "HLT";
HLTConfigProvider hltConfig_;
unsigned int ntrig;
edm::InputTag triggerResultsTag_;

int trig1;
char ctrig1[128];
string strig1 ;
};

VBFtrig::VBFtrig(const edm::ParameterSet& iConfig):
jetptcut(iConfig.getUntrackedParameter<double>("jetPtCut", 0 )),
detacut(iConfig.getUntrackedParameter<double>("dEtaCut", 0 )),
muptcut(iConfig.getUntrackedParameter<double>("muPtCut", 0 )),
tauptcut(iConfig.getUntrackedParameter<double>("tauPtCut", 0 )),
strig1(iConfig.getUntrackedParameter<string>("triggerName1", "none" ))
{
triggerResultsTag_ = iConfig.getParameter<edm::InputTag>("triggerResults") ;
processName = "TEST";

  edm::Service<TFileService> fs;
  mytree = fs->make<TTree>("mytree", "TTree to store Ntuple");

    mytree->Branch("trig1",    &trig1,  "trig1/I");
    mytree->Branch("ctrig1",    &ctrig1,  "ctrig1/C");

 mytree->Branch("genstablemet",&genstablemet,"genstablemet/D");
 mytree->Branch("genmupt",&genmupt,"genmupt/D");
 mytree->Branch("genmueta",&genmueta,"genmueta/D");
 mytree->Branch("gentaupt",&gentaupt,"gentaupt/D");
 mytree->Branch("gentaueta",&gentaueta,"gentaueta/D");

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

 mytree->Branch("genchi0pt",&genchi0pt,"genchi0pt/D");
 mytree->Branch("genchi0eta",&genchi0eta,"genchi0eta/D");
 mytree->Branch("genchipmpt",&genchipmpt,"genjchipmpt/D");
 mytree->Branch("genchipmeta",&genchipmeta,"genjchipmeta/D");

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
									ak5genj1pt = genJets -> at(i).pt() ;
									ak5genj2eta = genJets -> at(j).eta() ;
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

void
VBFtrig::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
genstablemet = -1;

genmet = -1;
genmjj = -1;
gendeta = 99;
genmettrue = -1;
ak5genmjj = -1;
ak5gendeta = -1;

genmueta = 99;
genmupt = -1;
gentaueta = 99;
gentaupt = -1;

ak5genj1pt = -1;
ak5genj2pt = -1;
ak5genj1eta = 99;
ak5genj2eta = 99;

genchi0eta = 99;
genchi0pt = -1;
genchipmeta = 99;
genchipmpt = -1;

trig1 = -1;
snprintf(ctrig1, 127, "none");

  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if(triggerResultsHandle_.isValid() == true ){
  	const edm::TriggerNames vTrigNames = iEvent.triggerNames(*triggerResultsHandle_);
  	const unsigned int ntrigs = triggerResultsHandle_->size();
  	for (unsigned int itr=0; itr<ntrigs; itr++){
    		TString trigName = vTrigNames . triggerName(itr);
    		if (triggerResultsHandle_->accept(itr)){
			if(trigName.Contains(strig1.c_str())){ 
				trig1 = 1;
                        	snprintf(ctrig1, 127, "%s", trigName.Data() );
			}
    		}
  	}
  }

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  if(genParticles.isValid() ){
        genstablemet = GenStableMET(genParticles);
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
  iEvent.getByLabel("genMetTrue", genMETs);
  if( genMETs.isValid() ){
	genmettrue = genMETs -> at(0).et();
  }
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel("ak5GenJets", genJets);
  if( genJets.isValid() ){
	Ak5GenMjjDeta( genJets, jetptcut, detacut);
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
