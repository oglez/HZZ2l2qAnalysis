#include "TROOT.h"
#include "TChain.h"
#include "TCut.h"
#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"

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
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"


#if defined(__CINT__) && !defined(__MAKECINT__)
namespace edm {
    typedef
    edm::Wrapper<vector<double> >
    Wrapper<vector<double,allocator<double> > >;
}
#endif




using namespace std;


/*   *************         */
/*   Value cuts for H350   */
/*   *************         */


void setGraphics(TH1F *histo){

  histo->SetFillColor(kAzure+7);
  histo->SetLineColor(kBlue+1);
}


TH1F * makeHist(string variable, string title, int nBins, float min, float max, TChain & Events,float scale){
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

void singleChannel(string channel, string variable1,string variable2, string title, int nBins, float min, float max, TChain & EventsVBF, TChain & EventsGF, float VBF, float GF){

  TH1F * h1VBF = makeHist((channel+"Higgs"+variable1).c_str(), (title+"VBF1").c_str(), nBins, min,  max, EventsVBF,1);
  TH1F * h2VBF = makeHist((channel+"Higgs"+variable2).c_str(), (title+"VBF2").c_str(), nBins, min,  max, EventsVBF,1);
  TH1F * hVBF = addHistos(h1VBF,h2VBF,1,1,(title + "VBF").c_str());  
  TH1F * h1GF = makeHist((channel+"Higgs"+variable1).c_str(), (title+"GF1").c_str(), nBins, min,  max, EventsGF,1);
  TH1F * h2GF = makeHist((channel+"Higgs"+variable2).c_str(), (title+"GF2").c_str(), nBins, min,  max, EventsGF,1);
  TH1F * hGF = addHistos(h1GF,h2GF,1,1,(title+"GF").c_str());
  
  TH1F * h = addHistos(hVBF, hGF, VBF, GF, (title+"Hist").c_str());
  
  //setGraphics(h);
  h->Write();
  delete h1VBF; 
  delete h2VBF; 
  delete hVBF; 
  delete h1GF; 
  delete h2GF; 
  delete hGF; 
  delete h;
}

/* create an histo combining electron and muon channel for that variable (VBF and GF)*/

TH1F * combHistos(string variable, string title, int nBins, float min, float max, TChain & EventsVBF, TChain & EventsGF,float VBF, float GF){
  TH1F * hVBF = makeHist(variable, title + "VBF", nBins, min, max, EventsVBF, VBF);
  TH1F * hGF = makeHist(variable, title + "GF", nBins, min, max, EventsGF, GF);
  TH1F * h = new TH1F((title).c_str(),title.c_str(), nBins, min, max);

  h->Add(hVBF,hGF);
  

  delete hVBF;
  delete hGF;
  setGraphics(h);
  return h;

}

void histo(string variable, string title, int nBins, float min, float max, TChain & EventsVBF, TChain & EventsGF,float VBF, float GF, TCut & cutMu, TCut & cutEl,  bool w=false){

  TH1F * h = new TH1F((title).c_str(),title.c_str(), nBins, min, max);
  TH1F * hVBFmu = new TH1F((title+"VBFmu").c_str(),(title+"VBFmu").c_str(), nBins, min, max);
  TH1F * hVBFel = new TH1F((title+"VBFel").c_str(),(title+"VBFel").c_str() , nBins, min, max);
  TH1F * hGFmu = new TH1F((title + "GFmu").c_str(),(title+"GFmu").c_str() , nBins, min, max);
  TH1F * hGFel = new TH1F((title+"GFel").c_str(), (title+"GFel").c_str(), nBins, min, max);

  EventsVBF.Project((title+"VBFmu").c_str(), ("muHiggs"+variable).c_str(), cutMu);
  EventsVBF.Project((title+"VBFel").c_str(), ("elHiggs"+variable).c_str(), cutEl);
  hVBFmu->Sumw2();
  hVBFmu->Add(hVBFel);
  
  EventsGF.Project((title+"GFmu").c_str(), ("muHiggs"+variable).c_str(),cutMu);
  EventsGF.Project((title+"GFel").c_str(), ("elHiggs"+variable).c_str(),cutEl );
  hGFmu->Sumw2();
  hGFmu->Add(hGFel);
  
  float VBFerror = sqrt(hVBFmu->GetEntries())*VBF;
  float GFerror = sqrt(hGFmu->GetEntries())*GF;
  hVBFmu->Scale(VBF);
  hGFmu->Scale(GF);
  h->Add(hVBFmu,hGFmu);
  float error = sqrt(pow(VBFerror,2)+pow(GFerror,2));
  setGraphics(h);
  if(w==true)h->Write();
  
  delete hVBFmu;
  delete hVBFel;
  delete hGFmu;
  delete hGFel;
  delete h;
}

void createHistos(string variable, string title, int nBins, float min, float max, TChain & EventsVBF, TChain & EventsGF,float VBF, float GF ){

  /* ****************** */
  /* SELECTION FOR H350 */
  /* ****************** */

    TCut baseSelMu = "(muHiggsLeptDau1Pt>20 && muHiggsLeptDau2Pt>20  && abs(muHiggsLeptDau1dB)<0.02 && abs(muHiggsLeptDau2dB)<0.02 && (muHiggsLeptDau1Eta < 2.1 || muHiggsLeptDau2Eta < 2.1) && muHiggsJetDau1Pt>30 && muHiggsJetDau2Pt>30 )"; 
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


    histo(variable,(title+"BaseSel").c_str(), nBins, min, max, EventsVBF, EventsGF, VBF,GF, baseSelMu, baseSelEl, true);
    histo(variable,(title+"MassSel").c_str() , nBins, min, max, EventsVBF, EventsGF, VBF,GF, massSelMu, massSelEl);
    histo(variable,(title+"BtagSel").c_str() , nBins, min, max, EventsVBF, EventsGF, VBF,GF,  btagSelMu, btagSelEl);
    histo(variable,(title+"ptllSel").c_str() , nBins, min, max, EventsVBF, EventsGF, VBF,GF,  ptllSelMu, ptllSelEl);
    histo(variable,(title+"metSel").c_str() , nBins, min, max, EventsVBF, EventsGF, VBF,GF,  metSelMu, metSelEl);
    histo(variable,(title+"drjjSel").c_str() , nBins, min, max, EventsVBF, EventsGF, VBF,GF,  drjjSelMu, drjjSelEl, true);
    histo(variable,(title+"hmassSel").c_str() , nBins, min, max, EventsVBF, EventsGF, VBF,GF, hmassSelMu, hmassSelEl );
}






int main() {

  //#include "DataFormats/FWLite/interface/Handle.h"

  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  gSystem->Load("libDataFormatsFWLite.so");
  gSystem->Load("libDataFormatsPatCandidates.so");


  TFile * output_file = TFile::Open("h350.root", "RECREATE");

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" ");
  const int max_line_len = 1024;
  gROOT->SetStyle("Plain");
  
  char str[max_line_len];
  ifstream inFile("test.txt");

  bool hFlag = true;
  map<string,vector<double> > data;

  
  while ( inFile ) {
    inFile.getline(str, max_line_len);
    //cout << str << endl;
    string s(str);
    tokenizer tok(s, sep);
    tokenizer::iterator i = tok.begin();
    if(i == tok.end() || (*i)[0] == '#') continue;
    string key = *i; ++i;
    if(hFlag && key != "hMass") {
      cerr << "first datacard line must be hMass" << endl;
      exit(100);
    }
    hFlag = false;
    vector<double> & v = data[key];
    for(; i!=tok.end();++i){
      stringstream ss(stringstream::in | stringstream::out);
      ss << *i;
      double x;
      ss>>x;
      v.push_back(x);
    }
    if(v.size() != data["hMass"].size()) {
      cerr << "all data lines must have the same number of entries" << endl;
      exit (101);
    }
  }


  

  cout<<"vbf events"<<data["VBFevents"][0]<<endl;
  TChain EventsVBF("Events"); 
  TChain EventsGF("Events"); 


  //  for(int k; k < v.size(); ++k){

  //TDirectory * dir = output_file->mkdir("hmassPlots");

  //cout<<"size "<<data["hMass"].size()<<endl; 

  for(unsigned int k = 2; k <3; ++k){
    //  for(unsigned int k = 0; k < data["hMass"].size(); ++k){
  
    //cout<< "MASS " << data["hMass"][k]<< endl;
  //    cout<<"mass"<<mass<<endl;
    
    EventsVBF.Add("rfio:/castor/cern.ch/user/d/decosa/Higgs/h350/edmntp/H350VBFEdmNtuples.root");
    EventsGF.Add("rfio:/castor/cern.ch/user/d/decosa/Higgs/h350/edmntp/H350GFEdmNtuples.root");

    float VBF = data["VBFxsec"][k]*1000/data["VBFevents"][k] ;
    float GF = data["GFxsec"][k]*1000/data["GFevents"][k];

    cout<<"VBF "<<VBF <<endl;
    cout<<"GF "<<GF <<endl;


    createHistos("Mass", "HMass", 200, 0,1000, EventsVBF, EventsGF, VBF, GF);
    createHistos("Pt", "HPt", 200, 0,1000, EventsVBF, EventsGF, VBF,GF);
    createHistos("Eta", "HEta", 13, -6.5,6.5, EventsVBF, EventsGF, VBF, GF);
    createHistos("Phi", "HPhi", 7, -3.5,3.5, EventsVBF, EventsGF, VBF, GF);
    createHistos("zllMass", "ZllMass", 200, 0,200, EventsVBF, EventsGF, VBF, GF);
    createHistos("zjjMass", "ZjjMass", 200, 0,200, EventsVBF, EventsGF, VBF, GF);
    createHistos("zllPt", "ZllPt", 200, 0,200, EventsVBF, EventsGF, VBF, GF);
    createHistos("zjjPt", "ZjjPt", 200, 0,200, EventsVBF, EventsGF, VBF, GF);
    createHistos("zllEta", "ZllEta", 13, -6.5,6.5, EventsVBF, EventsGF, VBF, GF);
    createHistos("zjjEta", "ZjjEta", 13, -6-5,6.5, EventsVBF, EventsGF, VBF, GF);
    createHistos("zllPhi", "ZllPhi", 13, -6.5,6.5, EventsVBF, EventsGF, VBF, GF);
    createHistos("zjjPhi", "ZjjPhi", 13, -6-5,6.5, EventsVBF, EventsGF, VBF, GF);
    createHistos("Jet1CSVMVA", "Jet1CSVMVA", 100, -100,100, EventsVBF, EventsGF, VBF, GF);
    createHistos("Jet2CSVMVA", "Jet2CSVMVA", 100, -100,100, EventsVBF, EventsGF, VBF, GF);
    createHistos("Jet1JbProb", "Jet1JbProb", 100, -100,100, EventsVBF, EventsGF, VBF, GF);
    createHistos("Jet2JbProb", "Jet2JbProb", 100, -100,100, EventsVBF, EventsGF, VBF, GF);
    createHistos("jjdr", "jjdr", 100, -100,100, EventsVBF, EventsGF, VBF, GF);
       
    TH1F * h = combHistos("met", "met", 100, -100,100, EventsVBF, EventsGF, VBF, GF);
    h->Write();
    delete h;
    

    /*    singleChannel("mu", "LeptDau1Pt","LeptDau2Pt", "muPt", 100, 0, 200,  EventsVBF, EventsGF, VBF, GF);
    singleChannel("mu", "LeptDau1Eta","LeptDau2Eta", "muEta", 100, -6., 6.,  EventsVBF, EventsGF, VBF, GF);
    singleChannel("mu", "LeptDau1Phi","LeptDau2Phi", "muPhi", 100, -3., 3.,  EventsVBF, EventsGF, VBF, GF);
    singleChannel("el", "LeptDau1Pt","LeptDau2Pt", "elPt", 100, 0, 200,  EventsVBF, EventsGF, VBF, GF);
    singleChannel("el", "LeptDau1Eta","LeptDau2Eta", "elEta", 100, -6., 6.,  EventsVBF, EventsGF, VBF, GF);
    singleChannel("el", "LeptDau1Phi","LeptDau2Phi", "elPhi", 100, -3., 3.,  EventsVBF, EventsGF, VBF, GF);
    singleChannel("mu", "JetDau1Pt","JetDau2Pt", "muChanJetPt", 100, 0, 200,  EventsVBF, EventsGF, VBF, GF);
    singleChannel("mu", "JetDau1Eta","JetDau2Eta", "muChanJetEta", 100, -6., 6.,  EventsVBF, EventsGF, VBF, GF);
    singleChannel("mu", "JetDau1Phi","JetDau2Phi", "muChanJetPhi", 100, -3., 3.,  EventsVBF, EventsGF, VBF, GF);
    singleChannel("el", "JetDau1Pt","JetDau2Pt", "elChanJetPt", 100, 0, 200,  EventsVBF, EventsGF, VBF, GF);
    singleChannel("el", "JetDau1Eta","JetDau2Eta", "elChanJetEta", 100, -6., 6.,  EventsVBF, EventsGF, VBF, GF);
    singleChannel("el", "JetDau1Phi","JetDau2Phi", "elChanJetPhi", 100, -3., 3.,  EventsVBF, EventsGF, VBF, GF);
    */



  }
  /*
  for(map<string,vector<double> >::const_iterator k = data.begin(); k != data.end(); ++k) {
  cout << k->first << ":";
  
  for(vector<double>::const_iterator j = k->second.begin(); j != k->second.end(); ++j) {
  cout << " " << *j;
  }
  cout << endl;
  }
  */

  output_file->Close();
 }
