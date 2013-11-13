///
/// Code modified by Oscar Gonzalez (8/X/2013) after inheriting it from Francesco.
/// Modified to add the PU MVA information to the jets... as I best can, since Yun Ju did
/// not give me enough feedback about the needed information.

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
//#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
//#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "DataFormats/Candidate/interface/Candidate.h"
//#include "DataFormats/CaloTowers/interface/CaloTower.h"
//#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometrySurface/interface/Line.h"

#include "TLorentzVector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include <memory>
#include <Riostream.h>
#include <string>
#include <vector>

using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;

typedef std::vector< pat::Jet > PFJetCollectionAB;
typedef std::vector< reco::PFCandidate > PFCandCollectionAB;

class PFJetUserData : public edm::EDProducer{

public:
  PFJetUserData(const edm::ParameterSet&);
  ~PFJetUserData(){
    //
  }

private:
  virtual void beginJob();
  virtual void endJob();
  void produce( edm::Event &, const edm::EventSetup &);

  // Utilities for beta/beta*
  /// Computes the beta and betastar variables for the given PAT jet.
  void computeBeta (const pat::Jet &jjet,float *beta,float *betastar);
  /// Reads the vertices from the event record and process them to be used
  /// for the calculation of beta and betastar. Note that some selection is
  /// performed on what is a vertex.
  /// It also selects which is the main vertex of the analysis... taken to
  /// be the first that is valid and not fake (if it passes the cuts).
  //  void readVertices (const edm::Event &iEvent);
  // changing the interface
  void readVertices ( edm::Handle<reco::VertexCollection> recVtxs );
  // check if the jet satisfies the taggeability criteria

  void isTaggableJet( const pat::Jet & jet, const reco::Vertex & primVertex,
		      const edm::ESHandle<TransientTrackBuilder> & builder,
		      int & nTracks,
		      int & nChargedTracks, int & nChargedTracksSV,
		      bool & isTaggable, bool & isTaggableSV ); 

  /// Gets the corrections to be applied to the jet for smearing (called if enabled).
  void computeSmearing (float pt, float eta, float genPt, std::vector<float> *val) const;


    
  //data members
  edm::InputTag   jetLabel_;
  //  bool is2012Data_;
  //  edm::InputTag   qgMap_;

  //Z variable of the reconstructed vertices.
  std::vector<float> _verticesZ;  
  //Index of the main vertex (-1 if no main vertex).  
  int _mainVertex;  
  // Variables for statistics
  int _nValidVertices;   // Number of valid vertices.
  int _nSelectedVertices;  // Number of selected vertices.
  int _nEventsWithValidVtx;  // Number of events with a "valid" vertex.
  int _nEventsWithMainVtx;  // Number of events with a main vertex.

  bool applySmearing_;   // Aply the smearing to the jets.

  bool verbose_;
};

PFJetUserData::PFJetUserData(const edm::ParameterSet &pSet) :
  applySmearing_(pSet.getUntrackedParameter<bool>("applySmearing",false))
{
  jetLabel_ =pSet.getUntrackedParameter<edm::InputTag>("JetInputCollection");
  //  is2012Data_ =pSet.getUntrackedParameter<bool>("is2012Data");
  //  qgMap_ =pSet.getUntrackedParameter<edm::InputTag>("qgMap");
  verbose_=pSet.getUntrackedParameter<bool>("Verbosity");

  // for beta/beta* 
  _verticesZ.clear();
  _mainVertex=-1;
  _nValidVertices=0;
  _nSelectedVertices=0;
  _nEventsWithValidVtx=0;
  _nEventsWithMainVtx=0;

  // issue the produce<>
 produces< std::vector< pat::Jet > >();
 
}

void PFJetUserData::beginJob(){
  //
}

void PFJetUserData::endJob(){
}

void PFJetUserData::produce(edm::Event &iEvt,  const edm::EventSetup &iSetup){
  if(verbose_)std::cout<<"Processing run "<<iEvt.id().run()<<", event "<<iEvt.id().event()<<std::endl;

  //TransientTrackBuilder
  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

  // Reading the vertices
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvt.getByLabel("offlinePrimaryVertices",recVtxs);
  // primary vtx
  const reco::Vertex & primVertex = (*recVtxs)[0];

 //pick from the event the input jet collection 
  Handle< PFJetCollectionAB > jetColl;
  //edm::Handle<edm::View<pat::Jet> > jetColl;
  iEvt.getByLabel(jetLabel_, jetColl);

  //  edm::Handle<edm::View<pat::Jet> > jets;
  //  iEvt.getByLabel("selectedPatJetsAK5",jets);

  /*
    Handle< std::vector<pat::Muon> > muColl;
    iEvt.getByLabel(muSrc_, muColl);
    std::auto_ptr< std::vector<pat::Muon> > muCollNew(new std::vector<pat::Muon>(*muColl))
    for(unsigned int i=0; i<muColl->size();++i ){
    //  const pat::Muon& myMu=(*imu);
    pat::Muon& myMu=(*muColl)[i];
    myMu.addUserFloat("dummyFloat", 9999.0);
    }
  */

  //  Handle< ValueMap<float> > qglMap;
  //  if(is2012Data_){
  //    iEvt.getByLabel(qgMap_,qglMap);
  //  }

  int nTotInJets=0;
  nTotInJets += jetColl->size();

  if(nTotInJets>0){
    if(!jetColl->at(0).isPFJet() ){
      throw cms::Exception("Bad Input") <<"ERROR in PFJetUserData::produce ! Jet collection is NOT made by PF Jets. "<<std::endl;
    }
  }

  if(verbose_)std::cout<<"In the Event are present "<<nTotInJets<<" PF jets."<<std::endl;

   //create output collection with PF Candidates in the jets
  std::auto_ptr< std::vector< pat::Jet > > outputPFJets(new std::vector< pat::Jet >(*jetColl));
  //std::auto_ptr< std::vector< pat::Jet > > outputPFJets(new std::vector< pat::Jet >);

  // Adding the information we need from the PU tagger (official)
  Handle<ValueMap<float> > puJetIdMVA;
  iEvt.getByLabel("puJetMva","full53xDiscriminant",puJetIdMVA);

  Handle<ValueMap<int> > puJetIdFlag;
  iEvt.getByLabel("puJetMva","full53xId",puJetIdFlag);

  //OLD  Handle<ValueMap<StoredPileupJetIdentifier> > jets;
  //OLD  iEvt.getByLabel("puJetId",jets);

//OLD   // Basic check:
//OLD   cout<<"MIERDA: "<<jetColl->size()<<" "
//OLD       <<outputPFJets->size()<<" "<<puJetIdMVA->size()<<" "<<puJetIdFlag->size()<<endl;
//OLD 
//OLD   for (unsigned int i=0; i<jetColl->size(); ++i ) {
//OLD     edm::RefToBase<pat::Jet> jetRef(edm::Ref<PFJetCollectionAB>(jetColl,i));
//OLD 
//OLD     cout<<"   JET "<<" "<<jetColl->at(i).pt()<<endl;
//OLD 
//OLD     const pat::Jet & patjet = jetColl->at(i);
//OLD     float mva   = (*puJetIdMVA)[jetRef];
//OLD     //int    idflag = (*puJetIdFlag)[i];
//OLD     cout << "MIERDA jet " << i << " pt " << patjet.pt() << " eta " << patjet.eta() << " PU JetID MVA " << mva <<endl; 
//OLD   }

  // Initial setup for beta/beta* variables
  _verticesZ.clear();
  _mainVertex=-1;
  // read the vertices! 
  //  readVertices(iEvt);
  // changinf the interface: giving the vtx collection as input
  readVertices(recVtxs);

  uint njetInColl(0);

  for(PFJetCollectionAB::iterator ijet=outputPFJets->begin(); 
      ijet!=outputPFJets->end();++ijet,++njetInColl ){

    // computing beta/beta* variables
    //    pat::Jet *jjet = &(*ijet);
    const pat::Jet & patjjet = *ijet;
    // We ask for the variables to store
    float beta=-1;
    float betastar=-1;
    //    computeBeta(*jjet,&beta,&betastar);
    computeBeta(patjjet,&beta,&betastar);
      
    // get taggeability
    // auxiliary variables, made external for monitoring
    int nChargedTracks_Aux = 0;
    int nChargedTracksSV_Aux = 0;
    int nTracks_Aux = 0;
    
    bool isTaggable(false), isTaggableSV(false);
    isTaggableJet( patjjet, primVertex, builder, nTracks_Aux,
		   nChargedTracks_Aux, nChargedTracksSV_Aux,
		   isTaggable, isTaggableSV );
    

    edm::RefToBase<pat::Jet> jetRef(edm::Ref<PFJetCollectionAB>(jetColl,njetInColl));
    //    edm::RefToBase<reco::Jet> jetRef = (*jetColl)[njetInColl];

    //    edm::RefToBase<reco::Jet> jetRef = patjjet.originalObjectRef()

    //    float qgl = -100.0;
    //    if(is2012Data_)
    //      qgl = (*qglMap)[jetRef];
    //
    //    cout << "(eta, qgl) for jet " << njetInColl << " = " << ijet->eta() << ", " << qgl << endl;
    
    if(verbose_){
      std::cout<<"From PFJetUserData::produce: PF Jet substructure (q/g separation) -> "<<
	"pt= "<<ijet->pt()<<"  eta= "<<ijet->eta()<<"   phi= "<<ijet->phi()<<flush;
    }
    
    // Processing the constituents!
    
    std::vector<reco::PFCandidatePtr> pfCandidates = ijet->getPFConstituents();
    
    float sumPt_cands=0.;
    float sumPt2_cands=0.;
    float rms_cands=0.;
    
    //loop on jet constituents       
    for (vector<reco::PFCandidatePtr>::const_iterator jt = pfCandidates.begin(); jt != pfCandidates.end(); ++jt) {
      
      //PFCandidate::ParticleType id = (*jt)->particleId();
      // Convert particle momentum to normal TLorentzVector, wrong type :(
      math::XYZTLorentzVectorD const& p4t = (*jt)->p4();
      TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
      TLorentzVector jetp4;
      //	 jetp4.SetPtEtaPhiE(ijet->pt(), ijet->eta(), ijet->phi(), ijet->energy());
      jetp4.SetPtEtaPhiE(ijet->pt(), ijet->eta(), ijet->phi(), ijet->energy());
      
      sumPt_cands += p4.Pt();
      sumPt2_cands += (p4.Pt()*p4.Pt());
      //float deltaR = ijet->p4().DeltaR(p4);
      float deltaR = jetp4.DeltaR(p4);
      rms_cands += (p4.Pt()*p4.Pt()*deltaR*deltaR);
      
    }//end loop on PF cands
    
    float nChrgdMult=float(ijet->chargedHadronMultiplicity());
    float nNeutrMult=float(ijet->neutralHadronMultiplicity());
    float ptDJet = sqrt( sumPt2_cands )/sumPt_cands;
    float rmsCandJet = rms_cands/(sumPt_cands*sumPt_cands);
    
    if(verbose_){
      std::cout<<"  NchgdHadrMult="<<nChrgdMult<<"  NneutrHadrMult= "<<nNeutrMult<<"   ptD= "<<ptDJet<<"   rmsJET= "<<rmsCandJet<<std::endl;
    }
    
    // beta/beta* variables
    ijet->addUserFloat("puBeta",beta);
    ijet->addUserFloat("puBetaStar",betastar);
    // taggability variables
    ijet->addUserInt("isTaggable",int(isTaggable));
    ijet->addUserInt("isTaggableSV",int(isTaggableSV));
    ijet->addUserInt("nChTks",nChargedTracks_Aux);
    ijet->addUserInt("nChTksSV",nChargedTracksSV_Aux);
    ijet->addUserInt("nTks",nTracks_Aux);
    // QG discriminant
    //    ijet->addUserFloat("qgLike",qgl);
    //
    ijet->addUserFloat("nChrgdHadrMult", float(ijet->chargedHadronMultiplicity()));
    ijet->addUserFloat("nNeutrHadrMult", float(ijet->neutralHadronMultiplicity()));
    ijet->addUserFloat("nPhotMult", float(ijet->photonMultiplicity()));
    ijet->addUserFloat("nElecMult", float(ijet->electronMultiplicity()));
    ijet->addUserFloat("nMuonMult", float(ijet->muonMultiplicity()));
    ijet->addUserFloat("nHFHadrMult", float(ijet->HFHadronMultiplicity ()));
    ijet->addUserFloat("nHFEMMult", float(ijet->HFEMMultiplicity()));
    ijet->addUserFloat("pfChrgdHadrEnergy", float(ijet->chargedHadronEnergy()));
    ijet->addUserFloat("pfNeutrHadrEnergy", float(ijet->neutralHadronEnergy()));
    ijet->addUserFloat("pfPhotEnergy", float(ijet->photonEnergy()));
    ijet->addUserFloat("pfElecEnergy", float(ijet->electronEnergy()));
    ijet->addUserFloat("pfHFHadrEnergy", float(ijet->HFHadronEnergy()));
    ijet->addUserFloat("pfHFEMEnergy", float(ijet->HFEMEnergy()));
    ijet->addUserFloat("ptDJet", ptDJet);
    ijet->addUserFloat("RMSJet", rmsCandJet);

    // MVA PU information
    ijet->addUserFloat("puJetIdMVA",(*puJetIdMVA)[jetRef]);
    ijet->addUserInt("puJetIdFlag",(*puJetIdFlag)[jetRef]);

    // Testing MC: if requested

    float originalpt=-1;

    //OLD    cout<<"JETS "<<applySmearing_<<" "<<originalpt<<" "<<patjjet.pt()<<" "<<patjjet.eta()<<" "<<patjjet.phi()<<" "<<ijet->genJet()<<endl;

    if (applySmearing_) {
      const reco::GenJet *genJet = ijet->genJet();
      if (genJet!=NULL) {
	std::vector<float> corr;
	originalpt=patjjet.pt();
	computeSmearing(originalpt,patjjet.eta(),patjjet.genJet()->pt(),&corr);

	float newJERSF=corr[0]/originalpt;
	
	// We do not apply the smearing if the new momentum                                                                                  
        // is too small                                                                                                                      
        if (newJERSF>0.1) {
	  math::XYZTLorentzVector newP4(patjjet.px()*newJERSF,patjjet.py()*newJERSF,patjjet.pz()*newJERSF,patjjet.energy()*newJERSF);
	  
	  ijet->setP4(newP4);
	}
      }
    }
    //OLD    cout<<"       "<<originalpt<<" "<<patjjet.pt()<<" "<<patjjet.eta()<<" "<<patjjet.phi()<<" "<<patjjet.genJet();
    //OLD     if (patjjet.genJet()!=NULL) {
    //OLD       cout<<" "<<patjjet.genJet()->pt()<<" "<<patjjet.genJet()->eta()<<" "<<patjjet.genJet()->phi();
    //OLD     }
    //OLD     else cout<<"UNMATCHED";
    //OLD     cout<<endl;
    
    ijet->addUserFloat("unsmearedPt",originalpt);

    //OLD njetInColl++;

    // Printing info on jets

//    cout<<"JET "<<originalpt<<" "<<patjjet.pt()<<" "<<patjjet.eta()<<" "<<patjjet.phi()<<" "<<patjjet.genJet()<<endl;
//    cout<<"    "<< int(isTaggable)<<" "<<int(isTaggableSV)<<" "<<patjjet.bDiscriminator("jetBProbabilityBJetTags")<<" "<<patjjet.bDiscriminator("combinedSecondaryVertexMVABJetTags")<<endl;
//

  }//end loop on jets
  
  //  for(std::vector<pat::Muon>::const_iterator imu=muColl->begin(); imu!=muColl->end();++imu ){
  if(verbose_) std::cout<<"Ended run "<<iEvt.id().run()<<", event "<<iEvt.id().event()<<", processed "<< outputPFJets->size()<<" jets."<<std::endl;
  iEvt.put(outputPFJets);
  if(verbose_)std::cout<<"UserJetCollection added to edmEvent."<<std::endl;

}//end produce


//-----------------------------------------------------------------------
void PFJetUserData::computeBeta (const pat::Jet &jjet,float *beta,float *betastar)
// Computes the beta and betastar variables for the given PAT jet.
{
  // We set the values to 0 only if there will be information
  // associated to them... i.e. if there are vertices.
  if (_verticesZ.size()>0) {
    *betastar=0;
    if (_mainVertex!=-1) *beta=0;
  }

  float totalpt=0;  // Scalar sum of the pt of the charged PF constituent of the jet

  // We loop over the charged particles in the jet 

  for (std::vector<reco::PFCandidatePtr>::const_iterator xpart = jjet.getPFConstituents().begin();
         xpart!=jjet.getPFConstituents().end();++xpart) {

    reco::PFCandidate::ParticleType typ = (*xpart)->particleId();
    if (typ!=reco::PFCandidate::h
	&& typ!=reco::PFCandidate::e
	&& typ!=reco::PFCandidate::mu) continue;  // Neutral, ignored

    // We check the distance wrt to the vertices:

    float zpfo = (*xpart)->vz();

    totalpt += (*xpart)->pt();

    int minvtx=-1;
    float mindist=0.2;

    {int ivtx=0;
    for (std::vector<float>::const_iterator zvtx = _verticesZ.begin();
	 zvtx!=_verticesZ.end();++zvtx,++ivtx) {

      float d = fabs(zpfo-*zvtx);
      
      if (d<mindist) {
	minvtx=ivtx;
	mindist=d;
      }
    }}

    // Depending if the clostest vertex (if any) is the main or other
    // we fill beta or betastar
    if (minvtx==_mainVertex) {
      *beta += (*xpart)->pt();
    }
    else if (minvtx>=0) {
      *betastar += (*xpart)->pt();
    }
  }

  // We normalize the variables properly:

  if (totalpt>0) {
    if (*beta>0) *beta /= totalpt;
    if (*betastar>0) *betastar /= totalpt;
  }  

  //TEST  std::cout<<"JET: "<<jjet.pt()<<" "<<jjet.rapidity()<<" "<<jjet.phi()<<" "<<*beta<<" "<<*betastar<<std::endl;
}

//-----------------------------------------------------------------------
//void PFJetUserData::readVertices (const edm::Event &iEvent)
void PFJetUserData::readVertices (edm::Handle<reco::VertexCollection> recVtxs)
// Reads the vertices from the event record and process them to be used
// for the calculation of beta and betastar. Note that some selection is
// performed on what is a vertex.
// It also selects which is the main vertex of the analysis... taken to
// be the first that is valid and not fake (if it passes the cuts).
{
  // Reading the vertices
  //  edm::Handle<reco::VertexCollection> recVtxs;
  //  iEvent.getByLabel("offlinePrimaryVertices",recVtxs);

  // Processing the information

  int formain=0;

  for (reco::VertexCollection::const_iterator xvtx = recVtxs->begin();
       xvtx!=recVtxs->end();++xvtx) {

    if (!xvtx->isValid() || xvtx->isFake()) continue;  // Only valid and not fake vertices

    if (formain==0) ++_nEventsWithValidVtx;
    
    ++formain;
    /// We are selecting the first valid and not fake vertex as the
    /// main one... there is no main vertex if it fails the selection.

    // Cuts to select the vertices:

    if (xvtx->ndof()<4) continue;  // Number of dof

    if (fabs(xvtx->z())>24) continue;  // Z of the vertex

    {float rho = sqrt(xvtx->x()*xvtx->x()+xvtx->y()*xvtx->y());
      if (rho>2) continue;  // Rho of the vertex
    }
    
    // Here we have a good vertex:

    if (formain==1) {  // It is the first that is valid... main vertex.
      _mainVertex=_verticesZ.size();
      ++_nEventsWithMainVtx;
    }
    _verticesZ.push_back(xvtx->z());
  }

  // Counting the valid and not fake vertices.
  _nValidVertices+=formain;

  // Counting the selected vertices.
  _nSelectedVertices += _verticesZ.size();

}


void PFJetUserData::isTaggableJet( const pat::Jet & jet, const reco::Vertex & primVertex,
				   const edm::ESHandle<TransientTrackBuilder> & builder,
				   int & nTracks,
				   int & nChargedTracks, int & nChargedTracksSV,
				   bool & isTaggable, bool & isTaggableSV ) 
{

  bool passBaseSel = false;
  if ( jet.isPFJet() ) {

    passBaseSel = (
		   jet.pt() > 10.0 &&
		   TMath::Abs(jet.eta()) < 2.4 &&
		   jet.neutralHadronEnergyFraction() < 0.99 &&
		   jet.neutralEmEnergyFraction() < 0.99 &&
		   jet.nConstituents() > 1 &&
		   jet.chargedHadronEnergyFraction() > 0.0 &&
		   jet.chargedMultiplicity() > 0.0 &&
		   jet.chargedEmEnergyFraction() < 0.99
		   );

  } else if( jet.isCaloJet() ) {
    passBaseSel = false;
  }

  nChargedTracks = 0;
  nChargedTracksSV = 0;
  nTracks = 0;
  std::vector<Measurement1D> ipValErr;
  std::vector<Measurement1D> aipValErr;
  reco::TrackBase::TrackQuality quality=reco::TrackBase::qualityByName("highPurity");

  GlobalPoint Pv_point = GlobalPoint(primVertex.x(),
				     primVertex.y(),
				     primVertex.z());
  math::XYZPointD pv (primVertex.x(), primVertex.y(), primVertex.z());
  GlobalVector direction( jet.momentum().x(),
			  jet.momentum().y(),
			  jet.momentum().z() );

  const reco::TrackRefVector & tracks = jet.associatedTracks();

  reco::TrackRefVector::const_iterator lastTrack = tracks.end();
  for (  reco::TrackRefVector::const_iterator trackRef = tracks.begin();
	 trackRef != lastTrack; ++trackRef ) {
    ++nTracks;

    const reco::Track * track = trackRef->get();


    if( track->charge() &&
	track->hitPattern().numberOfValidPixelHits() >= 2 &&
	track->hitPattern().numberOfValidHits() >= 8 &&
	//track->hitPattern().numberOfValidTrackerHits() >= 8 &&
	track->pt() > 1.0 &&
	track->normalizedChi2() < 5 &&
	TMath::Abs( track->dxy(pv) ) < 0.2 && //?
	TMath::Abs( track->dz(pv) ) < 17 
	) {


      /// just for charged
      if ( reco::deltaR(jet, *track) < 0.5 ) {

	double distancetojet = 999.;
	double decayLen =  999.;

	if ( builder.isValid() ) {
	  // have new track within jet (after all cuts but distancetojet )
	  // now compute distancetojet 
	  reco::TransientTrack transientTrack = builder->build( track );
	  TrajectoryStateOnSurface stateAtOrigin = transientTrack.impactPointState(); 
	  if ( stateAtOrigin.isValid() ) {

	    // track axis
	    Line::PositionType posTrack(stateAtOrigin.globalPosition());
	    Line::DirectionType dirTrack((stateAtOrigin.globalMomentum()).unit());
	    Line trackLine(posTrack, dirTrack);

	    // jet axis
	    GlobalVector jetVector = direction.unit();
	    Line::DirectionType dirJet( jetVector );    
	    Line::PositionType posJet( Pv_point );
	    Line jetLine( posJet, dirJet );

	    /// distace from track closest approach point to the jet axis
	    distancetojet = ( jetLine.distance(trackLine) ).mag();


	    //  dtojet -> Fill(distancetojet);
	    TrajectoryStateOnSurface closest =
	      IPTools::closestApproachToJet(stateAtOrigin,
					    primVertex, direction,
					    transientTrack.field());
	    if ( closest.isValid() ) {
	      decayLen = (closest.globalPosition() - (Pv_point)).mag();
	    }
		
	  }
	  if ( distancetojet < 0.07 && decayLen < 5 ) {
	    ++nChargedTracks;
	    ipValErr.push_back( IPTools::signedImpactParameter3D(transientTrack, direction, primVertex ).second );
	    aipValErr.push_back( IPTools::signedImpactParameter3D(transientTrack, direction, primVertex ).second );
	  }
	  if ( track->quality(quality) &&
	       reco::deltaR(jet.momentum(), track->momentum()) < 0.3 && 
	       distancetojet < 0.2 ) {
	    ++nChargedTracksSV;
	  }
	  // if no valid builder
	} else {
	  ++nChargedTracks;
	  ++nChargedTracksSV;
	}
      }
    }
  }
  ///
  isTaggable = ( passBaseSel && nChargedTracks >= 1 && nTracks >= 2  );
  isTaggableSV = ( passBaseSel && nChargedTracksSV >= 1 && nTracks >= 2  );
}

void PFJetUserData::computeSmearing (float pt, float eta, float genPt, std::vector<float> *val) const
// Gets the corrections to be applied to the jet for smearing (called if enabled).
{
  // Reading the parameters to compute the correction:
  eta=fabs(eta);
  double ptSF(1.0), ptSF_err(0.06);
  if(eta<0.5) { ptSF=1.052; ptSF_err=sqrt(pow(0.012,2)+pow(0.5*(0.062+0.061),2)); }
  else if(eta>=0.5 && eta<1.1) { ptSF=1.057; ptSF_err=sqrt(pow(0.012,2)+pow(0.5*(0.056+0.055),2)); }
  else if(eta>=1.1 && eta<1.7) { ptSF=1.096; ptSF_err=sqrt(pow(0.017,2)+pow(0.5*(0.063+0.062),2)); }
  else if(eta>=1.7 && eta<2.3) { ptSF=1.134; ptSF_err=sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2)); }
  else if(eta>=2.3 && eta<5.0) { ptSF=1.288; ptSF_err=sqrt(pow(0.127,2)+pow(0.5*(0.155+0.153),2)); }
      
  // Setting the values:

  val->clear();
  val->push_back(TMath::Max(0.,(genPt+ptSF*(pt-genPt))));  // Central: index 0
  val->push_back(TMath::Max(0.,(genPt+(ptSF+ptSF_err)*(pt-genPt))));   // Up error: index 1
  val->push_back(TMath::Max(0.,(genPt+(ptSF-ptSF_err)*(pt-genPt))));   // Down error: index 2
}

// ========= MODULE DEF ==============
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFJetUserData);
