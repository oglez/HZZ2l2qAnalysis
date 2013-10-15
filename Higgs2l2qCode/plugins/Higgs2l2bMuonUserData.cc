/// @file
/// Addition to the muon candidate of the required userdata

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include <vector>
#include <TMath.h>

using namespace edm;
using namespace std;
using namespace reco;

class Higgs2l2bMuonUserData : public edm::EDProducer {
public:
  Higgs2l2bMuonUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  InputTag src_, primaryVertices_;
  //  const float R03;
};

Higgs2l2bMuonUserData::Higgs2l2bMuonUserData( const ParameterSet & cfg ):
  src_( cfg.getParameter<InputTag>("src") ),
  //  rho_( cfg.getParameter<edm::InputTag>("rho")),
  primaryVertices_(cfg.getParameter<InputTag>("primaryVertices"))//,
  //  R03(0.3)
{
  produces<std::vector<pat::Muon> >();
}

void Higgs2l2bMuonUserData::produce( Event & evt, const EventSetup & ) {

  Handle<vector<pat::Muon>  > muons;
  evt.getByLabel(src_,muons);

  //  Handle<double> rhoHandle;
  //  evt.getByLabel(rho_,rhoHandle);

  //  double rho = *rhoHandle; 
  //  float PUEnergyInCone = (TMath::Pi()) * R03 * R03 * rho;  

  Handle<reco::VertexCollection> primaryVertices;  // Collection of primary Vertices
  evt.getByLabel(primaryVertices_, primaryVertices);
  const reco::Vertex &pv = (*primaryVertices)[0];

  auto_ptr<vector<pat::Muon> > muonColl( new vector<pat::Muon> (*muons) );
  for (unsigned int i = 0; i< muonColl->size();++i){
    pat::Muon & m = (*muonColl)[i];

    const pat::TriggerObjectStandAloneCollection muHLTMatches =  m.triggerObjectMatches();
    float muHLTBit =-1 ;
    unsigned int muHLTSize = muHLTMatches.size();
    muHLTSize>0 ? muHLTBit = 1 : muHLTBit = 0;  
    m.addUserFloat("muHLTBit", muHLTBit);


    //    float absCombIsoPUCorrected = m.trackIso() + m.caloIso() - PUEnergyInCone;
    //    m.addUserFloat("absCombIsoPUCorrected", absCombIsoPUCorrected);

    float dzVtx(-1000.0);
    if( m.innerTrack().isNonnull() )
      dzVtx = m.innerTrack()->dz(pv.position());  // Same as muonBestTrack
    m.addUserFloat("dzVtx", dzVtx);

    //cout<<"MIERDA: "<<dzVtx<<" "<<m.muonBestTrack()->dz(pv.position())<<endl;

    // The isolation, with delta-beta correction, is computed using the following
    // recipe.
    // It should be noted that we checked this is the same as the following
    // that was used in the python: 
    //       (max(0., (neutralHadronIso + photonIso - 0.5*puChargedHadronIso) ) 
    //               + chargedHadronIso) < 0.12 * pt
    //
    int isolated=0;
    float iso04 = m.pfIsolationR04().sumChargedHadronPt
      + max(0.,m.pfIsolationR04().sumNeutralHadronEt+m.pfIsolationR04().sumPhotonEt-0.5*m.pfIsolationR04().sumPUPt);
    if (iso04<0.12*m.pt()) isolated=1;   // IT IS ISOLATED (and also passed ID)

    m.addUserInt("isIsolated",isolated);

    // For checking values and muons:
    // cout<<"MUON: "<<i<<" "<<m.pt()<<" "<<m.eta()<<" "<<isolated<<endl;
    //    cout<<"MUON      "<<m.neutralHadronIso()<<" "<<m.photonIso()<<" "<<m.puChargedHadronIso()<<" "<<m.chargedHadronIso()<<" "<<endl;
    //    cout<<"MUON      "<<m.pfIsolationR04().sumNeutralHadronEt<<" "<<m.pfIsolationR04().sumPhotonEt<<" "<<m.pfIsolationR04().sumPUPt<<" "<<m.pfIsolationR04().sumChargedHadronPt<<endl;

  }

  evt.put( muonColl);

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( Higgs2l2bMuonUserData );


