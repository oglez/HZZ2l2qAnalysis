#include <Riostream.h>
#include <vector>
#include <string>


//PDFs 
#include "../interface/RooSpinZero5DV2.h"
#include "../interface/RooBkgd2L2JV2.h"

using namespace std;
//using namespace RooFit;


class HelicityLikelihoodDiscriminant{

 public:
  HelicityLikelihoodDiscriminant();
  ~HelicityLikelihoodDiscriminant();

  struct HelicityAngles{   
    float helCosTheta1;
    float helCosTheta2;
    float helCosThetaStar;
    float helPhi;
    float helPhi1;
    float mzz;
  };

  void init();
  void setParamFile(string myfilename);
  void setParametersFromFile();
  double getSignalProbability();
  double getBkgdProbability();
  void setMeasurables(vector<double> myvars);
  void setMeasurables(double newmzz,double newcostheta1,double newcostheta2,double newcosthetastar, double newphi, double newphistar1);
  void setMeasurables(HelicityAngles ha);


 private:

  //DATA MEMBERS
  string parfilename_;
  int NVARS;//this should be set in the constructor equal to the number of measurables below


  // ---------------- measurables ---------------------- 
  double mzz;	
  double costheta1;
  double costheta2;
  double costhetastar;
  double phi;
  double phistar1;
  
  //---- parameters ---------
  //
  //for background PDFs
//===================== cos theta 1 ====================                                 
  double bkgd_h1_acca0_;
  double bkgd_h1_acca2_;
  double bkgd_h1_acca4_;
  
  //===================== cos theta 2 ===================                                   
  double bkgd_h2_acca0_;
  double bkgd_h2_acca2_;
  double bkgd_h2_acca4_;
  double bkgd_h2_g_;
  double bkgd_h2_cutOff_;
  
  //===================== cos theta * ==================       
  double bkgd_hs_para2_;
  double bkgd_hs_para4_;

  //===================== phi =========================
  double bkgd_p_acca0_;
  double bkgd_p_acca1_;
  double bkgd_p_acca2_;
  
  //===================== phi1* =======================    
  double bkgd_p1_acca0_;
  double bkgd_p1_acca1_;
  double bkgd_p1_acca2_;  

  //for signal PDFs ===================================
  double sig_fppVal_;
  double sig_fmmVal_;
  double sig_fpmVal_;
  double sig_fp0Val_;
  double sig_f0mVal_;
  
  double sig_phippVal_;
  double sig_phimmVal_;
  double sig_phipmVal_;
  double sig_phip0Val_;
  double sig_phi0mVal_;
  
  double sig_fz1Val_;
  double sig_fz2Val_;
  
  double sig_R1Val_;
  double sig_R2Val_;
  
  double sig_para2_;
  double sig_para4_;
  double sig_acca0_;
  double sig_acca1_;
  double sig_acca2_;
  double sig_a2_;
  double sig_a4_;
  double sig_b2_;
  double sig_b4_;
  double sig_N_;
  double sig_g_;
  double sig_cutOff_;
   //
  //------ end parameters ----------------


  //--------------- PDFs -----------------------
  //
  //backgrounds
  RooBkgd2L2JV2 *background;
  //signal
  RooSpinZero5DV2 *signal;
  //
  //--------------- end PDFs ----------------





  //  void initMyRooRealVar(RooRealVar &rrv, const char* newname, const char* newtitle, double newmin, double newmax);
  // void initMyRooRealVar(RooRealVar &rrv, const char* newname, const char* newtitle, double newvalue, double newmin, double newmax);
  void initBkgdPDFs();
  void initSignalPDFs();
  void setBkgdParameters();
  void setSignalParameters();

};
