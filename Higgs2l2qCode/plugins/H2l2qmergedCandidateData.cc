/// @file
///
/// The code computes some information associated to the final candidates
/// when there is a merged jet and stores it as a PAT collection.
///
/// <PRE>
/// Written by Oscar Gonzalez: 30/IX/2013   Copied from the dijet one
///             (note: //OGL on useful stuff we might need to change)
/// </PRE>
///
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "PhysicsTools/CandUtils/interface/CenterOfMassBooster.h"
#include "PhysicsTools/CandUtils/interface/Booster.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/Track.h"

// PF candidates
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Helicity angles 
#include "HZZ2l2qAnalysis/Higgs2l2qCode/interface/Helicity.h"
// KinFit
#include "HZZ2l2qAnalysis/Higgs2l2qCode/interface/JetKinFitter.h"
// Likelihood discriminant
#include "HZZ2l2qAnalysis/Higgs2l2qCode/interface/HelicityLikelihoodDiscriminant.h"

#include <Math/VectorUtil.h>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

class H2l2qmergedCandidateData : public edm::EDProducer {
public:
  H2l2qmergedCandidateData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  void helicityAngles(const reco::Candidate *, const reco::Candidate *);
  //OGL  void helicityAnglesRefit(const reco::Candidate *, TLorentzVector, TLorentzVector);
  //OGL  int runJetKinFit(TLorentzVector &, TLorentzVector &, 
  //OGL  		  const TLorentzVector &, TLorentzVector &, TLorentzVector &,
  //OGL                   float &, float &);
  void getLDVariables(float, float, float, float, float, float, 
  		      float &, float &, float &);  

  InputTag higgsTag;//OGL, gensTag;
  InputTag prunedJetTag;// pruned jets for pruned mass, subjet kinematics and subjet btags
  //OGL, PFCandTag, vtxTag;
  //OGL  double deltaZCut_;
  //OGL  PFJetIDSelectionFunctor jetIDLoose;
//OGL   pat::strbitset ret; 
  // angular variables with unrefitted jet quantities
  double costhetaNT1, costhetaNT2, costhetastarNT, phiNT, phiNT1;
  // angular variables with refitted jet quantities
//OGL   double costhetaNT1Refit, costhetaNT2Refit, costhetastarNTRefit, phiNTRefit, phiNT1Refit;
//OGL   double zNominalMass_;
//OGL   JetKinFitter kinFitter_;
};

H2l2qmergedCandidateData::H2l2qmergedCandidateData( const ParameterSet & cfg ):
  higgsTag( cfg.getParameter<InputTag>( "higgs" ) ),
  prunedJetTag( cfg.getParameter<InputTag>( "prunedjets" ) )
  //OGL  ,gensTag( cfg.getParameter<edm::InputTag>("gensTag"))
//OGL  ,PFCandTag( cfg.getParameter<edm::InputTag>("PFCandidates"))
//OGL  ,vtxTag(cfg.getParameter<InputTag>("primaryVertices"))
//OGL  ,deltaZCut_( cfg.getParameter<double>("dzCut"))
//OGL  ,jetIDLoose( PFJetIDSelectionFunctor::FIRSTDATA,
//OGL	      PFJetIDSelectionFunctor::LOOSE )
//OGL  ,zNominalMass_( 91.1876 )
//OGL  ,kinFitter_( zNominalMass_, 0.0 )
{
  //OGL   ret = jetIDLoose.getBitTemplate();
  produces<vector<pat::CompositeCandidate> >("h").setBranchAlias( "h" );
}

void H2l2qmergedCandidateData::produce( Event & evt, const EventSetup & ) {

  Handle<std::vector<reco::CompositeCandidate> > higgsH;
  evt.getByLabel(higgsTag, higgsH);

  Handle<std::vector<pat::Jet>> prunedH;
  evt.getByLabel(prunedJetTag, prunedH);


//OGL   Handle<GenParticleCollection> gensH;
//OGL   evt.getByLabel(gensTag, gensH);
//OGL   GenParticleCollection gens = *gensH;
//OGL 
//OGL   // get PFCandidates
//OGL   Handle<PFCandidateCollection> pfCandidates;
//OGL   evt.getByLabel(PFCandTag, pfCandidates);
//OGL 
//OGL   // get Primary vtx
//OGL   //  Handle<reco::VertexCollection> primaryVertices;  // Collection of primary Vertices
//OGL   //  evt.getByLabel(vtxTag, primaryVertices);
//OGL   //  const reco::Vertex &pv = (*primaryVertices)[0];

  auto_ptr<vector<pat::CompositeCandidate> > higgsColl( new vector<pat::CompositeCandidate> () );

  float zzdPhi, zzdEta, zzdr, lldPhi, lldEta,lldr, jjdPhi, jjdEta, jjdr; 
//OGL   float neutralEmEnergy, chargedEmEnergy, chargedHadronEnergy, energy;
//OGL   float jminid, jmaxid;
//OGL   float j0LooseID, j1LooseID;
//OGL   bool  jminbmatch, jmincmatch, jmaxbmatch, jmaxcmatch;
//OGL   //  float  lminpt, lmineta, lminphi, lmaxpt, lmaxeta, lmaxphi, jminpt, jmineta, jminphi, jmaxpt,  jmaxeta,  jmaxphi;
//OGL   float  jmineta, jminphi, jmaxeta, jmaxphi;
//OGL   float j1RefitPt, j2RefitPt;
//OGL   float j1RefitEta, j2RefitEta;
//OGL   float j1RefitPhi, j2RefitPhi;
//OGL   float j1RefitE, j2RefitE;
//OGL   float ZjjRefitMass;
//OGL   float HZZRefitMass;
//OGL   float KFchiSquare, KFchiSquareProb;
//OGL   TLorentzVector j1corr;
//OGL   TLorentzVector j2corr;
//OGL   TLorentzVector HZZKinFit4mom, ZLL4mom, Zjj4mom; //initialized to (0, 0, 0 ,0)
  float helyLD;
  float ldSig, ldBkg;
//OGL   float helyLDRefit;
//OGL   float ldSigRefit, ldBkgRefit;
//OGL   //  float trkMetX, trkMetY, trkMet, trkPlusNeuMet;
//OGL   //  float neutralContributionX, neutralContributionY;
//OGL   //  float trkCorrectedMetX, trkCorrectedMetY;

  for (unsigned int i = 0; i< higgsH->size();++i){
    const reco::CompositeCandidate & H = (*higgsH)[i];
    edm::Ref<std::vector<reco::CompositeCandidate> > hRef(higgsH, i);
    pat::CompositeCandidate h(H);

    const reco::Candidate * Zll = h.daughter(0);     // The leptonic Z

    const Candidate * zDauRefl0 = Zll->daughter(0);
    const Candidate * zDauRefl1 = Zll->daughter(1);

    // We have only one hadronic jet... that is the daughter (CA8 jet) itself

    const Candidate * zDauRefca8 = H.daughter(1);
    const pat::Jet & jca8 = dynamic_cast<const pat::Jet &>(*(zDauRefca8));

    // find the appropriate pruned jet.
    // Due to the pruning they are not exactly identical, so we use Dr matching
    long unsigned indexTmp=0;
    int index = -1;
    double dRmin=99.0;
    for(pat::JetCollection::const_iterator pj = prunedH->begin();
	pj != prunedH->end(); ++pj, ++indexTmp ){
      
      edm::Ptr<pat::Jet> pjTMPPtr( prunedH , indexTmp );
      double dRtmp=deltaR(jca8,*pjTMPPtr);
      if(dRtmp<dRmin && dRtmp<0.7 ){//matching failed if greater than jet radius
	dRmin=dRtmp;
	index=indexTmp;
      }
    }//end loop on pruned jet collection

    //OLDstd::cout << index <<std::endl;
    edm::Ptr<pat::Jet> prunedJet(prunedH,index);
    
    //OLDstd::cout<< prunedJet->mass()<<std::endl;
    //OLDstd::cout<< prunedJet->numberOfDaughters()<<std::endl;
    
    if(prunedJet->mass()>5){  //IGNORE SUBSTRUCTURE FOR SINGLE SUBJET
      // compute helicity angles with unrefitted jet quantities
      helicityAngles( Zll, &(*prunedJet) );
      
      //OLDstd::cout << "Helicity angles: "<<
      //OLDcosthetaNT1 <<" , "<< costhetaNT2 <<" , "<< costhetastarNT <<" , "<<
      //OLDphiNT       <<" , "<< phiNT1      << std::endl;
      
      // dPhi, dEta, dr between H and Zs daughters
      zzdPhi = fabs( deltaPhi(Zll->phi(), jca8.phi()) ) ;
      zzdEta = fabs( (Zll->eta() - jca8.eta()) );
      zzdr = deltaR(Zll->eta(), Zll->phi(), jca8.eta(), jca8.phi() );
      
      const Candidate * subjet1ref = prunedJet->daughter(0);
      const pat::Jet & subjet1 = dynamic_cast<const pat::Jet &>(*(subjet1ref));
      const Candidate * subjet2ref = prunedJet->daughter(1);
      const pat::Jet & subjet2 = dynamic_cast<const pat::Jet &>(*(subjet2ref));
      
      //OLDstd::cout << "Subjet 1 btags: " << std::endl;
      //OLDstd::vector<std::pair<std::string, float> > btaginfo = subjet1.getPairDiscri();
      //OLDfor(std::vector<std::pair<std::string, float> >::const_iterator tagger = btaginfo.begin(); tagger!=btaginfo.end();tagger++)
      //OLDstd::cout << tagger->first << " : " << tagger->second << std::endl;
      //OLD}
      //OLDelse{
      //OLDstd::cout << "pruned mass smaller than 5GeV, don't look at substructure"<<std::endl;
    }
 
//OGL     
//OGL     lldPhi = fabs(deltaPhi(zDauRefl0->phi(),zDauRefl1->phi() ) ) ;
//OGL     lldEta = fabs(zDauRefl0->eta() - zDauRefl1->eta());
//OGL     lldr = deltaR(zDauRefl0->eta(), zDauRefl0->phi(), zDauRefl1->eta(), zDauRefl1->phi() );
//OGL       
//OGL     jjdPhi = fabs(deltaPhi(zDauRefj0->phi(),zDauRefj1->phi() ) ) ;
//OGL     jjdEta = fabs(zDauRefj0->eta() - zDauRefj1->eta());
//OGL     jjdr = deltaR(zDauRefj0->eta(), zDauRefj0->phi(), zDauRefj1->eta(), zDauRefj1->phi() );
//OGL 
//OGL     // store jetID for Z->jj daughters
//OGL     j0LooseID = (float) jetIDLoose( j0, ret );
//OGL     j1LooseID = (float) jetIDLoose( j1, ret );
//OGL    
//OGL     if(j0.pt() < j1.pt() ){      
//OGL       neutralEmEnergy = j0.neutralEmEnergy();
//OGL       chargedEmEnergy = j0.chargedEmEnergy() ;
//OGL       chargedHadronEnergy =j0.chargedHadronEnergy() ;
//OGL       energy = j0.energy() ;
//OGL       if (neutralEmEnergy/energy < 1.00 && chargedEmEnergy/energy < 1.00 && chargedHadronEnergy/energy > 0) jminid = true;
//OGL       else jminid = false;
//OGL   
//OGL       neutralEmEnergy = j1.neutralEmEnergy();
//OGL       chargedEmEnergy = j1.chargedEmEnergy() ;
//OGL       chargedHadronEnergy =j1.chargedHadronEnergy() ;
//OGL       energy = j1.energy() ;
//OGL       if (neutralEmEnergy/energy < 1.00 && chargedEmEnergy/energy < 1.00 && chargedHadronEnergy/energy > 0) jmaxid = true;
//OGL       else jmaxid = false;
//OGL     }
//OGL     else{
//OGL       neutralEmEnergy = j1.neutralEmEnergy();
//OGL       chargedEmEnergy = j1.chargedEmEnergy() ;
//OGL       chargedHadronEnergy = j1.chargedHadronEnergy() ;
//OGL       energy = j1.energy() ;
//OGL       if (neutralEmEnergy/energy < 1.00 && chargedEmEnergy/energy < 1.00 && chargedHadronEnergy/energy > 0) jminid = true;
//OGL       else jminid = false;
//OGL 
//OGL       neutralEmEnergy = j0.neutralEmEnergy();
//OGL       chargedEmEnergy = j0.chargedEmEnergy() ;
//OGL       chargedHadronEnergy = j0.chargedHadronEnergy() ;
//OGL       energy = j0.energy() ;
//OGL       if (neutralEmEnergy/energy < 1.00 && chargedEmEnergy/energy < 1.00 && chargedHadronEnergy/energy > 0) jmaxid = true;
//OGL       else jmaxid = false;
//OGL     }
//OGL 
//OGL     // gen info    
//OGL     //    if (zDauRefl0->pt() < zDauRefl1->pt()) {
//OGL     //      lminpt  = zDauRefl0->pt();
//OGL     //      lmineta = zDauRefl0->eta();
//OGL     //      lminphi = zDauRefl0->phi();
//OGL     //      lmaxpt  = zDauRefl1->pt();
//OGL     //      lmaxeta = zDauRefl1->eta();
//OGL     //      lmaxphi = zDauRefl1->phi();
//OGL     //    }
//OGL     //    
//OGL     //    else {
//OGL     //      lminpt  = zDauRefl1->pt();
//OGL     //      lmineta = zDauRefl1->eta();
//OGL     //      lminphi = zDauRefl1->phi();
//OGL     //      lmaxpt  = zDauRefl0->pt();
//OGL     //      lmaxeta = zDauRefl0->eta();
//OGL     //      lmaxphi = zDauRefl0->phi();
//OGL     //    }
//OGL     
//OGL     if (zDauRefj0->pt() < zDauRefj1->pt()) {
//OGL       //      jminpt  = zDauRefj0->pt();
//OGL       jmineta = zDauRefj0->eta();
//OGL       jminphi = zDauRefj0->phi();
//OGL       //      jmaxpt  = zDauRefj1->pt();
//OGL       jmaxeta = zDauRefj1->eta();
//OGL       jmaxphi = zDauRefj1->phi();
//OGL     }
//OGL     
//OGL     else {
//OGL       //      jminpt  = zDauRefj1->pt();
//OGL       jmineta = zDauRefj1->eta();
//OGL       jminphi = zDauRefj1->phi();
//OGL       //      jmaxpt  = zDauRefj0->pt();
//OGL       jmaxeta = zDauRefj0->eta();
//OGL       jmaxphi = zDauRefj0->phi();
//OGL     }
//OGL     
//OGL     jminbmatch = false;
//OGL     jmincmatch = false;
//OGL     jmaxbmatch = false;
//OGL     jmaxcmatch = false;
//OGL     
//OGL     for (size_t k = 0; k < gens.size(); k++) {
//OGL       if (abs(gens[k].pdgId()) == 5 && gens[k].status() != 3 && deltaR(jmineta, jminphi, gens[k].eta(), gens[k].phi()) < 0.3) jminbmatch = true; 
//OGL       if (abs(gens[k].pdgId()) == 4 && gens[k].status() != 3 && deltaR(jmineta, jminphi, gens[k].eta(), gens[k].phi()) < 0.3) jmincmatch = true; 
//OGL       if (abs(gens[k].pdgId()) == 5 && gens[k].status() != 3 && deltaR(jmaxeta, jmaxphi, gens[k].eta(), gens[k].phi()) < 0.3) jmaxbmatch = true; 
//OGL       if (abs(gens[k].pdgId()) == 4 && gens[k].status() != 3 && deltaR(jmaxeta, jmaxphi, gens[k].eta(), gens[k].phi()) < 0.3) jmaxcmatch = true; 
//OGL     }
//OGL     
//OGL 
//OGL     // prepare input for KinFit
//OGL     double j1en = j0.energy(); 
//OGL     double j1pt = j0.pt();
//OGL     double j1eta = j0.eta(); 
//OGL     double j1phi = j0.phi();
//OGL     double j2en = j1.energy(); 
//OGL     double j2pt = j1.pt(); 
//OGL     double j2eta = j1.eta(); 
//OGL     double j2phi = j1.phi();
//OGL     j1corr.SetPtEtaPhiE(j1pt,j1eta,j1phi,j1en);
//OGL     j2corr.SetPtEtaPhiE(j2pt,j2eta,j2phi,j2en);
//OGL     ZLL4mom.SetPtEtaPhiM(Zll->pt(), Zll->eta(), Zll->phi(), Zll->mass());
//OGL 
//OGL     // run KinFit
//OGL     int kinfitstatus = runJetKinFit(j1corr, j2corr, ZLL4mom, Zjj4mom, HZZKinFit4mom, KFchiSquare, KFchiSquareProb);
//OGL 
//OGL     if (kinfitstatus==0) {
//OGL       j1RefitPt = j1corr.Pt();
//OGL       j2RefitPt = j2corr.Pt();
//OGL       j1RefitEta = j1corr.Eta();
//OGL       j2RefitEta = j2corr.Eta();
//OGL       j1RefitPhi = j1corr.Phi(); 
//OGL       j2RefitPhi = j2corr.Phi();
//OGL       j1RefitE = j1corr.E(); 
//OGL       j2RefitE = j2corr.E();
//OGL       ZjjRefitMass = Zjj4mom.M();
//OGL       HZZRefitMass = HZZKinFit4mom.M();
//OGL     } else {
//OGL       //kinematic fit failed
//OGL       j1RefitPt = 0;
//OGL       j2RefitPt = 0;
//OGL       j1RefitEta = 0;
//OGL       j2RefitEta = 0;
//OGL       j1RefitPhi = 0; 
//OGL       j2RefitPhi = 0;
//OGL       j1RefitE = 0; 
//OGL       j2RefitE = 0;
//OGL       ZjjRefitMass = 0;
//OGL       HZZRefitMass = 0;
//OGL       KFchiSquare = -1. ; 
//OGL       KFchiSquareProb = -1.;
//OGL     }
//OGL 
//OGL     // compute helicity angles with refitted jet quantities
//OGL     helicityAnglesRefit( Zll, j1corr, j2corr);
//OGL 
//OGL 
    // Get LD variable
    getLDVariables( costhetaNT1, costhetaNT2, costhetastarNT,
		    phiNT, phiNT1, H.mass(), 
		    ldSig, ldBkg, helyLD );

//OGL     if (kinfitstatus==0)
//OGL       getLDVariables( costhetaNT1Refit, costhetaNT2Refit, costhetastarNTRefit,
//OGL 		      phiNTRefit, phiNT1Refit, HZZRefitMass, 
//OGL 		      ldSigRefit, ldBkgRefit, helyLDRefit );
//OGL     else{
//OGL       ldSigRefit = -100.0;
//OGL       ldBkgRefit = -100.0;
//OGL       helyLDRefit = -100.0;
//OGL     }
//OGL 
//OGL 
//OGL     h.addUserFloat("zzdPhi", zzdPhi);
//OGL     h.addUserFloat("zzdEta", zzdEta);
//OGL     h.addUserFloat("zzdr", zzdr);
//OGL     h.addUserFloat("lldPhi", lldPhi);
//OGL     h.addUserFloat("lldEta", lldEta);
//OGL     h.addUserFloat("lldr", lldr);
//OGL     h.addUserFloat("jjdPhi", jjdPhi);
//OGL     h.addUserFloat("jjdEta", jjdEta);
//OGL     h.addUserFloat("jjdr", jjdr);
//OGL     h.addUserFloat("jminbmatch",jminbmatch );
//OGL     h.addUserFloat("jmincmatch",jmincmatch );
//OGL     h.addUserFloat("jmaxbmatch",jmaxbmatch );
//OGL     h.addUserFloat("jmaxcmatch",jmaxcmatch );
//OGL     h.addUserFloat("jminid",jminid);
//OGL     h.addUserFloat("jmaxid",jmaxid);
//OGL     h.addUserFloat("jet1LooseID",j0LooseID);
//OGL     h.addUserFloat("jet2LooseID",j1LooseID);
    h.addUserFloat("costhetaNT1",costhetaNT1);
    h.addUserFloat("costhetaNT2",costhetaNT2);
    h.addUserFloat("phiNT",phiNT);
    h.addUserFloat("phiNT1",phiNT1);
    h.addUserFloat("costhetastarNT",costhetastarNT);
//OGL     h.addUserFloat("costhetaNT1Refit",costhetaNT1Refit);
//OGL     h.addUserFloat("costhetaNT2Refit",costhetaNT2Refit);
//OGL     h.addUserFloat("phiNTRefit",phiNTRefit);
//OGL     h.addUserFloat("phiNT1Refit",phiNT1Refit);
//OGL     h.addUserFloat("costhetastarNTRefit",costhetastarNTRefit);
//OGL     h.addUserFloat("j1RefitPt", j1RefitPt);
//OGL     h.addUserFloat("j2RefitPt", j2RefitPt);
//OGL     h.addUserFloat("j1RefitEta", j1RefitEta);
//OGL     h.addUserFloat("j2RefitEta", j2RefitEta);
//OGL     h.addUserFloat("j1RefitPhi", j1RefitPhi);
//OGL     h.addUserFloat("j2RefitPhi", j2RefitPhi);
//OGL     h.addUserFloat("j1RefitE", j1RefitE);
//OGL     h.addUserFloat("j2RefitE", j2RefitE);
//OGL     h.addUserFloat("ZjjRefitMass", ZjjRefitMass);
//OGL     h.addUserFloat("HZZRefitMass", HZZRefitMass);
//OGL     h.addUserFloat("KFchiSquare", KFchiSquare);
//OGL     h.addUserFloat("KFchiSquareProb", KFchiSquareProb);
    h.addUserFloat("helyLD", helyLD);
    h.addUserFloat("ldSig", ldSig);
    h.addUserFloat("ldBkg", ldBkg);
//OGL     h.addUserFloat("helyLDRefit", helyLDRefit);
//OGL     h.addUserFloat("ldSigRefit", ldSigRefit);
//OGL     h.addUserFloat("ldBkgRefit", ldBkgRefit);
//OGL     // new met variables
//OGL     //    h.addUserFloat("trkMetX", trkMetX);
//OGL     //    h.addUserFloat("trkMetY", trkMetY);
//OGL     //    h.addUserFloat("trkCorrectedMetX", trkCorrectedMetX);
//OGL     //    h.addUserFloat("trkCorrectedMetY", trkCorrectedMetY);
//OGL     //    h.addUserFloat("trkMet", trkMet);
//OGL     //    h.addUserFloat("trkPlusNeuMet", trkPlusNeuMet);

    higgsColl->push_back(h);
  }
  
  evt.put( higgsColl, "h");
}

void H2l2qmergedCandidateData::getLDVariables( float costhetaNT1, float costhetaNT2, float costhetastarNT,
					       float phiNT, float phiNT1, float Hmass, 
					       float & ldSig, float & ldBkg, float & helyLD ){
  ldSig = -100.0; ldBkg=-100.0; helyLD=-100.0;
  HelicityLikelihoodDiscriminant LD_;
  HelicityLikelihoodDiscriminant::HelicityAngles myha;
  myha.helCosTheta1    = costhetaNT1;
  myha.helCosTheta2    = costhetaNT2;
  myha.helCosThetaStar = costhetastarNT;
  myha.helPhi          = phiNT;
  myha.helPhi1         = phiNT1;
  myha.mzz             = Hmass;
  LD_.setMeasurables(myha);
  ldSig = LD_.getSignalProbability();
  ldBkg = LD_.getBkgdProbability();
  helyLD = ldSig / (ldSig + ldBkg);
}

void H2l2qmergedCandidateData::helicityAngles (const reco::Candidate *Zll, const reco::Candidate *Zjj) {
  // prepare for helicity angles computation
  costhetaNT1 = -8.80; costhetaNT2 = -8.80; costhetastarNT = -8.80;
  phiNT       = -8.80; phiNT1      = -8.80;
  TLorentzVector p4lept1(0.0,0.0,0.0,0.0);
  TLorentzVector p4lept2(0.0,0.0,0.0,0.0);
  TLorentzVector p4jet1(0.0,0.0,0.0,0.0);
  TLorentzVector p4jet2(0.0,0.0,0.0,0.0);
  //set as lepton #1 the negative one
  int lM, lP;
  if(Zll->daughter(0)->charge()<0.0) 
    lM = 0;
  else   
    lM = 1;
  lP = 1-lM;
   
  p4lept1.SetPxPyPzE(Zll->daughter(lM)->p4().x(),Zll->daughter(lM)->p4().y(),Zll->daughter(lM)->p4().z(),Zll->daughter(lM)->p4().e());
  p4lept2.SetPxPyPzE(Zll->daughter(lP)->p4().x(),Zll->daughter(lP)->p4().y(),Zll->daughter(lP)->p4().z(),Zll->daughter(lP)->p4().e());
  
  p4jet1.SetPxPyPzE(Zjj->daughter(0)->p4().x(),Zjj->daughter(0)->p4().y(),Zjj->daughter(0)->p4().z(),Zjj->daughter(0)->p4().e());
  p4jet2.SetPxPyPzE(Zjj->daughter(1)->p4().x(),Zjj->daughter(1)->p4().y(),Zjj->daughter(1)->p4().z(),Zjj->daughter(1)->p4().e());
  //compute helicity angles
  Helicity myAngles;
  myAngles.calculateAngles(p4lept1, p4lept2, p4jet1, p4jet2, costhetaNT1, costhetaNT2, costhetastarNT, phiNT, phiNT1);
}


//OGL void H2l2qmergedCandidateData::helicityAnglesRefit(const reco::Candidate *Zll, TLorentzVector p4jet1, TLorentzVector p4jet2) 
//OGL {  
//OGL   // prepare for helicity angles computation
//OGL   costhetaNT1Refit = -8.80; costhetaNT2Refit = -8.80; costhetastarNTRefit = -8.80;
//OGL   phiNTRefit       = -8.80; phiNT1Refit      = -8.80;
//OGL   TLorentzVector p4lept1(0.0,0.0,0.0,0.0);
//OGL   TLorentzVector p4lept2(0.0,0.0,0.0,0.0);
//OGL   //set as lepton #1 the negative one
//OGL   int lM, lP;
//OGL   if(Zll->daughter(0)->charge()<0.0) 
//OGL     lM = 0;
//OGL   else   
//OGL     lM = 1;
//OGL   lP = 1-lM;
//OGL   
//OGL   p4lept1.SetPxPyPzE(Zll->daughter(lM)->p4().x(),Zll->daughter(lM)->p4().y(),Zll->daughter(lM)->p4().z(),Zll->daughter(lM)->p4().e());
//OGL   p4lept2.SetPxPyPzE(Zll->daughter(lP)->p4().x(),Zll->daughter(lP)->p4().y(),Zll->daughter(lP)->p4().z(),Zll->daughter(lP)->p4().e());
//OGL   
//OGL   //compute helicity angles
//OGL   Helicity myAngles;
//OGL   myAngles.calculateAngles(p4lept1, p4lept2, p4jet1, p4jet2, costhetaNT1Refit, costhetaNT2Refit, costhetastarNTRefit, phiNTRefit, phiNT1Refit);
//OGL }

//OGL int H2l2qmergedCandidateData::runJetKinFit(TLorentzVector &j1,TLorentzVector &j2,
//OGL 				    const TLorentzVector &ZLL, TLorentzVector & Zjj, 
//OGL 				    TLorentzVector &XZZ, float & chiSquare, 
//OGL 				    float & chiSquareProb) {
//OGL   
//OGL   int status=0;
//OGL 
//OGL   //pass the two four momenta and initialize the Kinfit object
//OGL   kinFitter_.setJet4Mom(j1,j2);
//OGL 
//OGL   //ask the kinfit object to correct the four momenta
//OGL   status=kinFitter_.Refit();
//OGL   if(status==0){
//OGL     j1=kinFitter_.getCorrJets().at(0);
//OGL     j2=kinFitter_.getCorrJets().at(1);
//OGL     chiSquare     = kinFitter_.chiSquare();
//OGL     chiSquareProb = kinFitter_.chiSquareProb();
//OGL   }
//OGL 
//OGL   //update also 4-mom of XZZ and of ZJJ
//OGL   Zjj = j1+j2;
//OGL   XZZ = ZLL+Zjj;
//OGL 
//OGL   return status;
//OGL }


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( H2l2qmergedCandidateData );


