#include "TROOT.h"
#include "TChain.h"
#include "TCut.h"
#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"
#include "TMath.h"

#include <iostream>
#include <algorithm> 
#include <exception>
#include <iterator>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include<boost/tokenizer.hpp>
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

using namespace std;

void setGraphics(TH1F *histo){

  histo->SetFillColor(kAzure+7);
  histo->SetLineColor(kBlue+1);
}


TH1F * makeHist(string variable, string title, int nBins, float min, float max, TChain & Events,float scale=1.){

  TH1F * h = new TH1F(title.c_str(),title.c_str(), nBins, min, max);
  Events.Project(title.c_str(), variable.c_str());
  h->Scale(scale);
  setGraphics(h);
  return h;
}

/* add two histos with two different weights*/

TH1F * addHistos(TH1F*h1,TH1F*h2, float w1, float w2, string title){
  TH1F * h = new TH1F(title.c_str(),title.c_str(), h1->GetNbinsX(),h1->GetXaxis()->GetXmin() ,h1->GetXaxis()->GetXmax());
  h->Add(h1,h2,w1,w2);
  //h->Write();
  setGraphics(h);
  return h;
}


/* create a single channel histo adding dau1 and dau2 variable for VBF and GF*/

void singleChannel(string channel, string variable1,string variable2, string title, int nBins, float min, float max, TChain & Events, float scale){

  TH1F * h1 = makeHist((channel+"Higgs"+variable1).c_str(), (title+"1").c_str(), nBins, min,  max, Events,scale);
  TH1F * h2 = makeHist((channel+"Higgs"+variable2).c_str(), (title+"2").c_str(), nBins, min,  max, Events,scale);
  TH1F * h = addHistos(h1,h2,1,1,(title+"Histo").c_str());
  h->Write();
  delete h1; 
  delete h2; 
  delete h; 

}

void histo(string variable, string title, int nBins, float min, float max, TChain & Events,  TCut & cutMu, TCut & cutEl, float scale,  bool w=false){

  TH1F * hmu = new TH1F((title+"mu").c_str(),(title+"mu").c_str(), nBins, min, max);
  TH1F * hel = new TH1F((title+"el").c_str(),(title+"el").c_str() , nBins, min, max);

 /* muon channel + electron channel for VBF events */

  Events.Project((title+"mu").c_str(), ("muHiggs"+variable).c_str(), cutMu);
  Events.Project((title+"el").c_str(), ("elHiggs"+variable).c_str(), cutEl);
  hmu->Sumw2();
  hmu->Add(hel);
  hmu->Scale(scale);
  hmu->SetTitle(title.c_str());

  float error = sqrt(hmu->GetEntries())*scale;
  setGraphics(hmu);
  if (w==true) hmu->Write();
  
  delete hmu;
  delete hel;

}

void createHistos(string variable, string title, int nBins, float min, float max, TChain & Events, float scale, bool b = false ){

  /* ****************** */
  /* SELECTION FOR H350 */
  /* ****************** */

    TCut baseSelMu = "(muHiggsLeptDau1Pt>20 && muHiggsLeptDau2Pt>20  && TMath::Abs(muHiggsLeptDau1dB)<0.02 && TMath::Abs(muHiggsLeptDau2dB)<0.02 && (muHiggsLeptDau1Eta < 2.1 || muHiggsLeptDau2Eta < 2.1) && muHiggsJetDau1Pt>30 && muHiggsJetDau2Pt>30 )"; 
    TCut baseSelEl = "elHiggsLeptDau1Pt>20 && elHiggsLeptDau2Pt>20 && elHiggsJetDau1Pt>30 && elHiggsJetDau2Pt>30 && (elHiggsEleDau1VBTF80CombID==7 || elHiggsEleDau2VBTF80CombID==7)";

    TCut massSelMu = baseSelMu && "TMath::Abs(muHiggszllMass - 91)< 10 && TMath::Abs(muHiggszjjMass - 91)< 15";
    TCut massSelEl = baseSelEl && "TMath::Abs(elHiggszllMass - 91)< 10 && TMath::Abs(elHiggszjjMass - 91)< 15";

    TCut btagSelMu = massSelMu && "((muHiggsJet1CSVMVA>0.5 && muHiggsJet2JbProb>0.9 )||( muHiggsJet2CSVMVA>0.5 && muHiggsJet1JbProb>0.9))";
    TCut btagSelEl = massSelEl && "((elHiggsJet1CSVMVA>0.5 && elHiggsJet2JbProb>0.9 )||(elHiggsJet2CSVMVA>0.5 && elHiggsJet1JbProb>0.9))";

    TCut ptllSelMu = btagSelMu && " muHiggszllPt > 90";
    TCut ptllSelEl = btagSelEl && " elHiggszllPt > 90";
    
    TCut metSelMu = ptllSelMu && " met < 35";
    TCut metSelEl = ptllSelEl && " met < 35";

    TCut drjjSelMu = metSelMu && " muHiggsjjdr < 1.7";
    TCut drjjSelEl = metSelEl && " elHiggsjjdr < 1.7";
    
    TCut hmassSelMu = drjjSelMu && " muHiggsMass < 385 &&  muHiggsMass > 315 ";
    TCut hmassSelEl = drjjSelEl && " elHiggsMass < 385 &&  elHiggsMass > 315 ";


    histo(variable,(title+"BaseSel").c_str(), nBins, min, max, Events, baseSelMu, baseSelEl,  scale, true);
    histo(variable,(title+"MassSel").c_str() , nBins, min, max, Events, massSelMu, massSelEl,  scale, true);
    histo(variable,(title+"BtagSel").c_str() , nBins, min, max, Events,  btagSelMu, btagSelEl,  scale, true);
    histo(variable,(title+"ptllSel").c_str() , nBins, min, max, Events,  ptllSelMu, ptllSelEl, scale, true);
    histo(variable,(title+"metSel").c_str() , nBins, min, max, Events,  metSelMu, metSelEl, scale, true );
    histo(variable,(title+"drjjSel").c_str() , nBins, min, max, Events,  drjjSelMu, drjjSelEl,  scale, true);
    histo(variable,(title+"hmassSel").c_str() , nBins, min, max, Events, hmassSelMu, hmassSelEl, scale, true);
}




void bkgPlots(string s, TChain & Events, float sigma=1., float nevtin=1000.) {

  TFile * output_file = TFile::Open( (s+".root").c_str(), "RECREATE");

  float L= 1000;
  float scaleFact = L*sigma/nevtin ;

  createHistos("Mass", "HMass", 200, 0,1000, Events, scaleFact, true);
  createHistos("Pt", "HPt", 200, 0,1000, Events,  scaleFact);
  createHistos("Eta", "HEta", 13, -6.5,6.5, Events,  scaleFact);
  createHistos("Phi", "HPhi", 7, -3.5,3.5, Events,  scaleFact);
  createHistos("zllMass", "ZllMass", 200, 0,200, Events,  scaleFact);
  createHistos("zjjMass", "ZjjMass", 200, 0,200, Events,  scaleFact);
  createHistos("zllPt", "ZllPt", 200, 0,200, Events,  scaleFact);
  createHistos("zjjPt", "ZjjPt", 200, 0,200, Events,  scaleFact);
  createHistos("zllEta", "ZllEta", 13, -6.5,6.5, Events,  scaleFact);
  createHistos("zjjEta", "ZjjEta", 13, -6-5,6.5, Events,  scaleFact);
  createHistos("zllPhi", "ZllPhi", 13, -6.5,6.5, Events,  scaleFact);
  createHistos("zjjPhi", "ZjjPhi", 13, -6-5,6.5, Events,  scaleFact);
  createHistos("Jet1CSVMVA", "Jet1CSVMVA", 100, -100,100, Events,  scaleFact);
  createHistos("Jet2CSVMVA", "Jet2CSVMVA", 100, -100,100, Events,  scaleFact);
  createHistos("Jet1JbProb", "Jet1JbProb", 100, -100,100, Events,  scaleFact);
  createHistos("Jet2JbProb", "Jet2JbProb", 100, -100,100, Events,  scaleFact);
  createHistos("jjdr", "jjdr", 100, -100,100, Events,  scaleFact);
    
  TH1F * h = makeHist("met", "met", 100, -100,100, Events, scaleFact);
  h->Write();
  delete h;
  /*  singleChannel("mu", "LeptDau1Pt","LeptDau2Pt", "muPt", 100, 0, 200,  Events, scaleFact);
  singleChannel("mu", "LeptDau1Eta","LeptDau2Eta", "muEta", 100, -6., 6.,  Events, scaleFact);
  singleChannel("mu", "LeptDau1Phi","LeptDau2Phi", "muPhi", 100, -3., 3.,  Events, scaleFact);
  singleChannel("el", "LeptDau1Pt","LeptDau2Pt", "elPt", 100, 0, 200,  Events, scaleFact);
  singleChannel("el", "LeptDau1Eta","LeptDau2Eta", "elEta", 100, -6., 6.,  Events, scaleFact);
  singleChannel("el", "LeptDau1Phi","LeptDau2Phi", "elPhi", 100, -3., 3.,  Events, scaleFact);
  singleChannel("mu", "JetDau1Pt","JetDau2Pt", "muChanJetPt", 100, 0, 200,  Events, scaleFact);
  singleChannel("mu", "JetDau1Eta","JetDau2Eta", "muChanJetEta", 100, -6., 6.,  Events, scaleFact);
  singleChannel("mu", "JetDau1Phi","JetDau2Phi", "muChanJetPhi", 100, -3., 3.,  Events, scaleFact);
  singleChannel("el", "JetDau1Pt","JetDau2Pt", "elChanJetPt", 100, 0, 200,  Events, scaleFact);
  singleChannel("el", "JetDau1Eta","JetDau2Eta", "elChanJetEta", 100, -6., 6.,  Events, scaleFact);
  singleChannel("el", "JetDau1Phi","JetDau2Phi", "elChanJetPhi", 100, -3., 3.,  Events, scaleFact);
  */

  output_file->Close();
 }


int main() {


  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  gROOT->SetStyle("Plain");
  TChain EventsTT("Events"); 
  EventsTT.Add( "rfio:/castor/cern.ch/user/d/decosa/Higgs/TT/edmntp/HTTEdmNtuples.root");
  TChain EventsWZ("Events"); 
  EventsWZ.Add( "rfio:/castor/cern.ch/user/d/decosa/Higgs/WZ/edmntp/HWZEdmNtuples.root");
  TChain EventsZZ("Events"); 
  EventsZZ.Add( "rfio:/castor/cern.ch/user/d/decosa/Higgs/ZZ/edmntp/HZZEdmNtuples.root");


  bkgPlots("TT", EventsTT, 16.5,940000);
  bkgPlots("ZZ", EventsZZ, 7.76,1000000);
  bkgPlots("WZ", EventsWZ, 18.2,2194752);




}
