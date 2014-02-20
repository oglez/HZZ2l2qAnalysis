/// @file
/// Code to create the UserData for electrons for the HZZ2l2q framework.
/// <PRE>
/// Written by Francesco Fabozzi
///            O. Gonzalez (20/X/2013) added the regression energy part
/// </PRE>

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
// #include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"

#include <vector>
#include <TMath.h>

using namespace edm;
using namespace std;
using namespace reco;

class Higgs2l2bElectronUserData : public edm::EDProducer {
public:
  Higgs2l2bElectronUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  InputTag src_, rho_, primaryVertices_;

  bool regressionEnergy_;   // To know whether we should use regression energy
  //  const float R03;
};

Higgs2l2bElectronUserData::Higgs2l2bElectronUserData( const ParameterSet & cfg ):
  src_( cfg.getParameter<InputTag>("src") ),
  rho_( cfg.getParameter<edm::InputTag>("rho")),
  primaryVertices_(cfg.getParameter<InputTag>("primaryVertices"))//,
						    //  R03(0.3)
  ,regressionEnergy_(cfg.getUntrackedParameter<bool>("applyRegressionEnergy",false))
{
  produces<std::vector<pat::Electron> >();
}

void Higgs2l2bElectronUserData::produce( Event & evt, const EventSetup & ) {

  Handle<vector<pat::Electron>  > electrons;
  evt.getByLabel(src_, electrons);

  Handle<double> rhoHandle;
  evt.getByLabel(rho_,rhoHandle);

  double rho = *rhoHandle; 
  //  float PUEnergyInCone = (TMath::Pi()) * R03 * R03 * rho;  

  Handle<reco::VertexCollection> primaryVertices;  // Collection of primary Vertices
  evt.getByLabel(primaryVertices_, primaryVertices);
  const reco::Vertex &pv = (*primaryVertices)[0];

  Handle<reco::BeamSpot> bsHandle;
  evt.getByLabel("offlineBeamSpot", bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();
    
  Handle<reco::ConversionCollection> conversions;
  evt.getByLabel("allConversions", conversions);
  
  auto_ptr<vector<pat::Electron> > electronColl( new vector<pat::Electron> (*electrons) );
  for (unsigned int i = 0; i< electronColl->size();++i){
    pat::Electron & el = (*electronColl)[i];

    // The first part is to compute the energy regression and apply it to the electron

    float uncorrE = -1;
    if (regressionEnergy_) {
      uncorrE = el.energy();
      float e=el.ecalRegressionEnergy();
      if (e>0 && uncorrE>0) {
	uncorrE=e/uncorrE;
	math::XYZTLorentzVector newP4(el.px()*uncorrE,el.py()*uncorrE,el.pz()*uncorrE,e);
	el.setP4(newP4);
      }
    }
    el.addUserFloat("uncorrEnergy",uncorrE);

    // Doing trigger matches

    const pat::TriggerObjectStandAloneCollection elHLTMatches = el.triggerObjectMatches();
    float elHLTBit =-1 ;
    unsigned int elHLTSize = elHLTMatches.size();
    elHLTSize>0 ? elHLTBit = 1 : elHLTBit = 0;  
    el.addUserFloat("elHLTBit", elHLTBit);

    //    float absCombIsoPUCorrected = -1.0;
    //    if( el.isEB() )
    //      // pedestal subtraction for barrel 
    //      absCombIsoPUCorrected = el.dr03TkSumPt() + max(0., el.dr03EcalRecHitSumEt() - 1.) + el.dr03HcalTowerSumEt() - PUEnergyInCone;
    //    else
    //      absCombIsoPUCorrected = el.dr03TkSumPt() + el.dr03EcalRecHitSumEt() + el.dr03HcalTowerSumEt() - PUEnergyInCone;
    //
    //    el.addUserFloat("absCombIsoPUCorrected", absCombIsoPUCorrected);

    float dzVtx(-1000.0);
    float dxyVtx(-1000.0);
    float missHits(-1000.0);
    if( el.gsfTrack().isNonnull() ) {
      dzVtx = el.gsfTrack()->dz(pv.position());
      dxyVtx = el.gsfTrack()->dxy(pv.position());
      missHits = el.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    }
    el.addUserFloat("dz", dzVtx);
    el.addUserFloat("dxy", dxyVtx);
    el.addUserFloat("mHits", missHits);

    // conversion rejection match
    bool hasMatchConv = ConversionTools::hasMatchedConversion(el, conversions, beamspot.position());
    el.addUserFloat("hasMatchConv", float(hasMatchConv));

    // check electron ID
    // kinematic variables
    bool isEBEle           = el.isEB() ? true : false;
    float ptEle            = el.pt();
    float etaEle           = el.superCluster()->eta();
    // id variables
    float dEtaInEle        = el.deltaEtaSuperClusterTrackAtVtx();
    float dPhiInEle        = el.deltaPhiSuperClusterTrackAtVtx();
    float sigmaIEtaIEtaEle = el.sigmaIetaIeta();
    float hoeEle           = el.hadronicOverEm();
    float ooemoopEle       = (1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy());

    // pf isolation variables
    float iso_chEle = el.chargedHadronIso();
    float iso_nhEle = el.neutralHadronIso();
    float iso_emEle = el.photonIso();

    // detector isolation variables
    float trackIsoEle = el.dr03TkSumPt();
    float ecalIsoEle = el.dr03EcalRecHitSumEt();
    float hcalIsoEle = el.dr03HcalTowerSumEt();

//    cout<<"MIERDA "//OLD <<std::setprecision(6)
//        <<isEBEle<<" "<<ptEle<<" "<<etaEle<<" "<<dEtaInEle<<" "<<dPhiInEle<<" "<<sigmaIEtaIEtaEle<<" "<<hoeEle<<" "
//        <<ooemoopEle<<" "<<dxyVtx<<" "<<dzVtx<<" "<<iso_chEle<<" "<<iso_emEle<<" "<<iso_nhEle<<" "<<hasMatchConv
//        <<" "<<missHits<<" "<<rho<<endl;

    bool passVetoID = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::VETO, isEBEle, ptEle, etaEle,
						  dEtaInEle, dPhiInEle, sigmaIEtaIEtaEle, hoeEle,
						  ooemoopEle, dxyVtx, dzVtx, iso_chEle, iso_emEle, iso_nhEle, 
						  hasMatchConv, missHits, rho);
    
    bool passLooseID = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE, isEBEle, ptEle, etaEle,
						  dEtaInEle, dPhiInEle, sigmaIEtaIEtaEle, hoeEle,
						  ooemoopEle, dxyVtx, dzVtx, iso_chEle, iso_emEle, iso_nhEle, 
						  hasMatchConv, missHits, rho);
    
    bool passMediumID = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, isEBEle, ptEle, etaEle,
						  dEtaInEle, dPhiInEle, sigmaIEtaIEtaEle, hoeEle,
						  ooemoopEle, dxyVtx, dzVtx, iso_chEle, iso_emEle, iso_nhEle, 
						  hasMatchConv, missHits, rho);
    
    bool passTightID = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, isEBEle, ptEle, etaEle,
						  dEtaInEle, dPhiInEle, sigmaIEtaIEtaEle, hoeEle,
						  ooemoopEle, dxyVtx, dzVtx, iso_chEle, iso_emEle, iso_nhEle, 
						  hasMatchConv, missHits, rho);

    bool passTriggerTight = EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERTIGHT, 
								 isEBEle, ptEle, dEtaInEle, dPhiInEle, 
								 sigmaIEtaIEtaEle, hoeEle,
								 trackIsoEle, ecalIsoEle, hcalIsoEle);

    //    cout << "passveto = " << passVetoID <<  " passloose = " << passLooseID << " passmedium = " << passMediumID << " passtight = " << passTightID << " passtriggertight = " << passTriggerTight << endl;

    // embed maximum passed criterion in a code     
    int cutIDCode(0);    
    if(passTightID)
      cutIDCode = 4;
    else {
      if(passMediumID)
	cutIDCode = 3;
      else {
	if(passLooseID)
	  cutIDCode = 2;
	else {
	  if(passVetoID)
	    cutIDCode = 1;
	}
      }
    }

    //    cout << "CUT ID CODE = " << cutIDCode << endl;
    
    //el.addUserFloat("passTriggerTight", float(passTriggerTight));
    el.addUserInt("passTriggerTight", passTriggerTight);
    //el.addUserFloat("cutIDCode", float(cutIDCode));
    el.addUserInt("cutIDCode",cutIDCode);

    int isIsolated=0;  // Marks passing ID and isolation!
    if (cutIDCode>1 && passTriggerTight>0) isIsolated=1;
    el.addUserInt("isIsolated",isIsolated);

    // Some feedback for testing:

    //cout<<"ELECTRON: "<<cutIDCode<<" "<<el.pt()<<" "<<endl;


    
  }

  evt.put(electronColl);

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( Higgs2l2bElectronUserData );


