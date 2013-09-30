#include <vector>

#include "TLorentzVector.h"
#include "TMatrixD.h"

//#include "TFitConstraintM.h"
//#include "TFitParticleEtEtaPhi.h"
//#include "TKinFitter.h"
#include "HZZ2l2qAnalysis/Higgs2l2qCode/interface/DiJetKinFitter.h"
using namespace std;

class JetKinFitter{

public:
  JetKinFitter();
  JetKinFitter(float m, float merr );//mass constraint and its error
 ~JetKinFitter();
  void setJet4Mom(TLorentzVector j1in,TLorentzVector j2in);
  vector<TLorentzVector> getCorrJets();
  int Refit();
  double chiSquare() { return chiSquare_;};
  double chiSquareProb() { return chiSquareProb_;} ;

private:
  float M_,Merr_;
  TLorentzVector j1in_, j2in_,j1out_, j2out_;
  vector<TLorentzVector> corrjets_;
  //  Double_t ErrEt(Float_t Et, Float_t Eta);
  //  Double_t ErrEta(Float_t Et, Float_t Eta);
  //  Double_t ErrPhi(Float_t Et, Float_t Eta);
  double chiSquare_, chiSquareProb_;
  // void init();
};
