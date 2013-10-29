///@file
/// Module to read the candidates from several data modules and make a print
/// out of some basic information of the candidates... to make simple debugging.
/// <PRE>
/// Written by O. Gonzalez (14/X/2013)
/// </PRE>

#include "../interface/Higgs2l2qDebugModule.h"


// CMS classes

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//OLD #include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/MET.h"



// C++ classes

#include <memory>

#include <vector>
using std::vector;

#include <iostream>
using std::endl;
using std::cout;
using std::cerr;

#include <string>
using std::string;

//-----------------------------------------------------------------------
Higgs2l2qDebugModule::Higgs2l2qDebugModule (const edm::ParameterSet& iConfig) :
  edm::EDAnalyzer()
  ,collections_(iConfig.getParameter< vector<string> >("collections"))
// Constructor of the class
{
}

//-----------------------------------------------------------------------
Higgs2l2qDebugModule::~Higgs2l2qDebugModule (void)
// Destructor of the class
{
  // The collection handlers must be deleted
  collections_.clear();
}

//-----------------------------------------------------------------------
void Higgs2l2qDebugModule::beginJob (void)
// Method runs for the EDAnalyzer at the beginning of the job.
{
  cout<<"Start-up of Higgs2l2qDebugModule: "<<collections_.size()<<endl;

}

//-----------------------------------------------------------------------
void Higgs2l2qDebugModule::endJob (void)
// Method run for the EDAnalyzer at the end of the job.
{
  cout<<"--------------------------------------------------------------------------"<<endl;
  cout<<"Report from Higgs2l2qDebugModule: "<<endl;
  cout<<"--------------------------------------------------------------------------"<<endl;
}

//-----------------------------------------------------------------------
void Higgs2l2qDebugModule::analyze (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// Method run for each event in the analysis.
{
  // We open each collection and process

  for (vector<string>::const_iterator xcol = collections_.begin(); 
       xcol!=collections_.end();++xcol) {

    cout<<"CANDIDATE-DEBUG: candidates from "<<(*xcol);//<<endl;

    edm::Handle<vector<reco::CompositeCandidate> > cand_H;
    iEvent.getByLabel((*xcol),cand_H);

    cout<<" size: "<<cand_H->size()<<endl;

    for (vector<reco::CompositeCandidate>::const_iterator xcand = cand_H->begin();
	 xcand!=cand_H->end();++xcand) {
      // Printing out the information for the candidate
      printInfoCandidate(*xcand);
    }
  }

  // Check on MET:

  edm::Handle<pat::METCollection> metHandle;
  iEvent.getByLabel("patMETsPFJetsAK5", metHandle);
  pat::METCollection met_h = *metHandle;

  cout<<"MET: "<<met_h.front().et()<<" "<<met_h.front().sumEt()<<" "<<met_h.front().mEtSig()<<endl;
  // met significance
  TMatrixD metmat = met_h.front().getSignificanceMatrix();
  //  metmat.Print();
  if( (metmat < 1.0e10) && (metmat > -1.0e10) ) {
    //    float determ_mat = metmat.Determinant();
    //    cout << "determinant = " << determ_mat << endl;
    cout<<"     significance: "<< met_h.front().significance()<<endl;
  } 

}

//-----------------------------------------------------------------------
void Higgs2l2qDebugModule::printInfoCandidate (const reco::CompositeCandidate &cand) 
// Prints the information for a candidate to identify it kinematically.
{
  cout<<"Candidate: "<<cand.pt()<<" "<<cand.eta()<<" "<<cand.phi()<<" "<<cand.mass()<<" "<<cand.numberOfDaughters()<<endl;
  cout<<"         Daugh-0: "<<cand.daughter(0)->pt()<<" "<<cand.daughter(0)->eta()<<" "<<cand.daughter(0)->phi()<<" "<<cand.daughter(0)->mass()<<endl;
  cout<<"         Daugh-1: "<<cand.daughter(1)->pt()<<" "<<cand.daughter(1)->eta()<<" "<<cand.daughter(1)->phi()<<" "<<cand.daughter(1)->mass()<<endl;
}

//-----------------------------------------------------------------------
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Higgs2l2qDebugModule);
//=======================================================================
