/// @file
/// Modified version of the TotalKinematicFilter, with more options
/// and creating variables to be used in the analysis.
/// Filter is disabled by default.
///
/// <PRE>
/// Written by Oscar Gonzalez (28/XI/2013)
/// </PRE>

#include "../interface/CiematKinematicsTest.h"

using namespace edm;

CiematKinematicsTest::CiematKinematicsTest(const edm::ParameterSet& iPSet):  
  edm::EDProducer()
  ,src_(iPSet.getParameter<edm::InputTag>("src"))
					  //  ,tolerance_(iPSet.getParameter<double>("tolerance")),
					  //  ,verbose_(iPSet.getUntrackedParameter<bool>("verbose",false))
{ 

  produces<float>("maxMCKinematicTest").setBranchAlias("maxMCKinematicTest");  
   
}

CiematKinematicsTest::~CiematKinematicsTest() {}

void CiematKinematicsTest::produce(edm::Event& iEvent,const edm::EventSetup& iSetup)
{ 

  float nEcms = 0.;
  unsigned int nInit = 0;

  std::vector<float> p4tot(4,0.);
  unsigned int nPart = 0;

  // Gather information on the reco::GenParticle collection
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(src_, genParticles );
  
  for (reco::GenParticleCollection::const_iterator iter=genParticles->begin();iter!=genParticles->end();++iter){
    if ( nInit < 3 && (*iter).status() == 3 && (*iter).pdgId() == 2212 ) {
      nInit++;
      nEcms += (*iter).energy();
    }
    if ( (*iter).status() == 1) { 
      nPart++;

      p4tot[1] += (*iter).px();
      p4tot[2] += (*iter).py();
      p4tot[3] += (*iter).pz();
      p4tot[0] += std::sqrt( (*iter).px()*(*iter).px() + 
                           (*iter).py()*(*iter).py() + 
                           (*iter).pz()*(*iter).pz() + 
                           (*iter).mass()*(*iter).mass()) ; 
    }
  }

  std::auto_ptr<float> max(new float(0));

  (*max)=std::abs(p4tot[0]-nEcms);
  for (int i=1;i<4;++i) if (std::abs(p4tot[i]) > (*max)) (*max)=std::abs(p4tot[i]);

  // Now we have te maximum... if it is larger than 50 GeV we report with a WARNING

  if ((*max)>50) {
    std::cout << "WARNING: MC-Kinematic test failed for a lot: "<<(*max)
	      <<" for event: "<<iEvent.id().run()<<" "
	      << iEvent.id().event()<<std::endl;
  }

  // Saving the element
  iEvent.put(max,"maxMCKinematicTest");
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(CiematKinematicsTest);
//=======================================================================
