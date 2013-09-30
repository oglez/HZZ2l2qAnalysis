#include "../interface/JetKinFitter.h"
#include <Riostream.h>
using namespace std;


//
JetKinFitter::JetKinFitter(){
  //
  M_=91.1876;//gev
  Merr_=0.0;
  j1in_.SetXYZM(0.0, 0.0, 0.0, 0.0);
  j2in_.SetXYZM(0.0, 0.0, 0.0, 0.0);
  j1out_.SetXYZM(0.0, 0.0, 0.0, 0.0);
  j2out_.SetXYZM(0.0, 0.0, 0.0, 0.0);

}

JetKinFitter::JetKinFitter(float m, float merr ){
 M_=m;//gev
  Merr_=merr;
  j1in_.SetXYZM(0.0, 0.0, 0.0, 0.0);
  j2in_.SetXYZM(0.0, 0.0, 0.0, 0.0);
  j1out_.SetXYZM(0.0, 0.0, 0.0, 0.0);
  j2out_.SetXYZM(0.0, 0.0, 0.0, 0.0);
}
JetKinFitter::~JetKinFitter(){
  //
}
void JetKinFitter::setJet4Mom(TLorentzVector j1in,TLorentzVector j2in){
  j1in_=j1in;
  j2in_=j2in;
  //cout<<"JetKinFitter::setJet4Mom, changing inputs j1: "<< j1in_.Pt()<<" , "<<j1in_.Eta()<<" , "<<j1in_.Phi()<<" ) "<<"   J2 ("<<j2in_.Pt()<<" , "<<j2in_.Eta()<<" , "<<j2in_.Phi()<<" ) ] "<<endl;
  corrjets_.clear();
}
vector<TLorentzVector> JetKinFitter::getCorrJets(){
  //cout<<"JetKinFitter::getCorrJets(), returning jets j1: "<< corrjets_.at(0).Pt()<<" , "<<corrjets_.at(0).Eta()<<" , "<<corrjets_.at(0).Phi()<<" ) "<<"   J2 ("<<corrjets_.at(1).Pt()<<" , "<<corrjets_.at(1).Eta()<<" , "<<corrjets_.at(1).Phi()<<" ) "<<endl;
  return corrjets_;
}
int JetKinFitter::Refit(){
 
  DiJetKinFitter* fitter_jets = new DiJetKinFitter( "fitter_jets", "fitter_jets" );
  //perform the fit
  std::pair<TLorentzVector,TLorentzVector> jets_kinfit = fitter_jets->fit(j1in_, j2in_);

  corrjets_.push_back(jets_kinfit.first);
  corrjets_.push_back(jets_kinfit.second);
 
    //    TLorentzVector Zjj = j1in_ + j2in_;
    // TLorentzVector Zjj_kinfit_jets = j1out_ + j2out_;

  // store chi2/ndf and prob
  chiSquare_ = fitter_jets->getS()/fitter_jets->getNDF();
  chiSquareProb_ = TMath::Prob(fitter_jets->getS(), fitter_jets->getNDF());

  delete fitter_jets;

    return 0;
}


// error functions for jets:
/*
Double_t JetKinFitter::ErrEt(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 5.6;
    b = 1.25;
    c = 0.033;
  }
  else{
    a = 4.8;
    b = 0.89;
    c = 0.043;
  }
  InvPerr2 = (a * a) + (b * b) * Et + (c * c) * Et * Et;
  return InvPerr2;
}



Double_t JetKinFitter::ErrEta(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 1.215;
    b = 0.037;
    c = 7.941 * 0.0001;
  }
  else{
    a = 1.773;
    b = 0.034;
    c = 3.56 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}

Double_t JetKinFitter::ErrPhi(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 6.65;
    b = 0.04;
    c = 8.49 * 0.00001;
  }
  else{
    a = 2.908;
    b = 0.021;
    c = 2.59 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}

*/
