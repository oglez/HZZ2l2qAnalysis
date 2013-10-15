///@file
/// Module to read the candidates from several data modules and make a print
/// out of some basic information of the candidates... to make simple debugging.
/// <PRE>
/// Written by O. Gonzalez (14/X/2013)
/// </PRE>                

#ifndef __HIGGS2L2Q_DEBUGMODULE__H_
#define __HIGGS2L2Q_DEBUGMODULE__H_

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"

// CMS-based classes

class Event;
class ParameterSet;
class EventSetup;
//OLD #include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

//-----------------------------------------------------------------------
/// This class is the module to read candidates from several data
/// modules and make a print out of some basic information of the
/// candidates... to make simple debugging.

class Higgs2l2qDebugModule : public edm::EDAnalyzer 
{
  // Internal variables

  std::vector<std::string> collections_;    ///< collections of candidates to process.

  // Internal methods

  /// Method runs for the EDAnalyzer at the beginning of the job.
  virtual void beginJob (void);

  /// Method run for each event in the analysis.
  virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);

  /// Method run for the EDAnalyzer at the end of the job.
  virtual void endJob (void);

  /// Prints the information for a candidate to identify it kinematically. 
  void printInfoCandidate (const reco::CompositeCandidate &cand);


public:

  /// Destructor of the class.
  virtual ~Higgs2l2qDebugModule (void);

  /// Constructor of the class.
  explicit Higgs2l2qDebugModule (const edm::ParameterSet &iConfig);

};
#endif

//====================================================================
