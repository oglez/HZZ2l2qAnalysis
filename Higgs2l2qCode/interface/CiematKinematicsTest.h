#ifndef CiematKinematicsTest_h
#define CiematKinematicsTest_h

// Filter to select events with deviation from nominal total kinematics smaller than a tolerance parameter 

// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Particle.h"

class CiematKinematicsTest : public edm::EDProducer 
{
  public:
  explicit CiematKinematicsTest(const edm::ParameterSet&);
  ~CiematKinematicsTest();

  private:
  virtual void produce (edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------

  edm::InputTag src_;
  //  double tolerance_;
  //  bool verbose_;

};

#endif
