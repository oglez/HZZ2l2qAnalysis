// ------------------------------------------------
// 
// DiJetKinFitter - kinematic fit for X->jj
//
// ------------------------------------------------


#ifndef DiJetKinFitter_h
#define DiJetKinFitter_h

#include "../interface/TKinFitter.h"
#include "TLorentzVector.h"



class DiJetKinFitter : public TKinFitter {

 public:
  
  DiJetKinFitter( const TString&  name, const TString& title, double mass = 91.1876 );

  void set_mass( double mass ) { mass_ = mass; };

  double ErrEt(double Et, double Eta); 
  double ErrEta(double Et, double Eta); 
  double ErrPhi(double Et, double Eta); 

  std::pair<TLorentzVector,TLorentzVector> fit( TLorentzVector jet1, TLorentzVector jet2 );

 private:

  double mass_;
  TString name_;
  TString title_;

};



#endif
