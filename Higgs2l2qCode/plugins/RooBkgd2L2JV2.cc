#include <Riostream.h>
#include <math.h>
#include "../interface/RooBkgd2L2JV2.h"

using namespace std;

//ClassImp(RooBkgd2L2JV2) 

RooBkgd2L2JV2::RooBkgd2L2JV2():
  name_("background"),
  title_("backgroundtitle"),
  h1(0.0),
  h2(0.0),
  hs(0.0),
  phi(0.0),
  phi1(0.0),
  mZZ(0.0),
  h1param0(0.0),
  h1param2(0.0),
  h1param4(0.0),
  h2param0(0.0),
  h2param2(0.0),
  h2param4(0.0),
  h2g(0.0),
  h2cutOff(0.0),
  hsparam2(0.0),
  hsparam4(0.0),
  pacca0(0.0),
  pacca1(0.0),
  pacca2(0.0),
  p1acca0(0.0),
  p1acca1(0.0),
  p1acca2(0.0)
{
}

RooBkgd2L2JV2::RooBkgd2L2JV2(const char *name, const char *title,
			     double mycostheta1, double mycostheta2,
			     double mycosthetas, double myphi, 
			     double myphi1, double mymZZ,
			     vector<double>allPars):
  name_(name),
  title_(title)
{
  setVars(mycostheta1,mycostheta2,mycosthetas,myphi,myphi1,mymZZ);
  setPars(allPars);
}

RooBkgd2L2JV2::RooBkgd2L2JV2(const char *name, const char *title, 
	     double _h1,	  
	     double _h2,	  
	     double _hs, 	  
	     double _phi,	  
	     double _phi1,	  
	     double _mZZ,	  
	     double _h1param0, 
	     double _h1param2, 
	     double _h1param4, 
	     double _h2param0, 
             double _h2param2, 
             double _h2param4, 
             double _h2g,      	  
             double _h2cutOff, 
	     double _hsparam2, 
             double _hsparam4, 
             double _pacca0,   
             double _pacca1,   
             double _pacca2,   
             double _p1acca0,  
	     double _p1acca1,  
             double _p1acca2):                              
  name_(name),
  title_(title),
  h1(_h1),
  h2(_h2),
  hs(_hs),
  phi(_phi),
  phi1(_phi1),
  mZZ(_mZZ),
  h1param0(_h1param0),
  h1param2(_h1param2),
  h1param4(_h1param4),
  h2param0(_h2param0),
  h2param2(_h2param2),
  h2param4(_h2param4),
  h2g(_h2g),
  h2cutOff(_h2cutOff),
  hsparam2(_hsparam2),
  hsparam4(_hsparam4),
  pacca0(_pacca0),
  pacca1(_pacca1),
  pacca2(_pacca2),
  p1acca0(_p1acca0),
  p1acca1(_p1acca1),
  p1acca2(_p1acca2)
{ 
} 


 RooBkgd2L2JV2::RooBkgd2L2JV2(const RooBkgd2L2JV2& other, const char* name) :  
   name_(name),
   title_("dummytitle"),
   h1(other.h1),
   h2(other.h2),
   hs(other.hs),
   phi(other.phi),
   phi1(other.phi1),
   mZZ(other.mZZ),
   h1param0(other.h1param0),
   h1param2(other.h1param2),
   h1param4(other.h1param4),
   h2param0(other.h2param0),
   h2param2(other.h2param2),
   h2param4(other.h2param4),
   h2g(other.h2g),
   h2cutOff(other.h2cutOff),
   hsparam2(other.hsparam2),
   hsparam4(other.hsparam4),
   pacca0(other.pacca0),
   pacca1(other.pacca1),
   pacca2(other.pacca2),
   p1acca0(other.p1acca0),
   p1acca1(other.p1acca1),
   p1acca2(other.p1acca2)
{ 
} 



 double RooBkgd2L2JV2::evaluate() const 
 { 
   double h1Func = (1.0+h1param0*h1*h1)*(1.0+h1param2*h1*h1+h1param4*h1*h1*h1*h1);
   //   std::cout << "h1Func: " << h1Func << endl;
   //double h2Func = (1.0+h2param0*h2*h2)*(1.0+h2param2*h2*h2+h2param4*h2*h2*h2*h2)
   //  /(1 + exp((-h2cutOff + h2)/h2g));
   double h2Func = (1.0+h2param2*h2*h2-h2param2*h2*h2*h2*h2)/(1 + exp((-h2cutOff + h2)/h2g));
   //   std::cout <<"h2Func: " << h2Func <<endl;
   double hsFunc =  (1.0+hsparam2*hs*hs+hsparam4*hs*hs*hs*hs);
   //   std::cout <<"hsFunc: " << hsFunc <<endl;
   double phiFunc = (pacca1*cos(phi)+pacca2*cos(2.0*phi)+pacca0);
   //   std::cout <<"phiFunc: " << phiFunc <<endl;
   double phi1Func= (p1acca1*cos(phi1)+p1acca2*cos(2.0*phi1)+p1acca0);
   //   std::cout <<"phi1Func: " << phi1Func <<endl;

   double h1_Norm = 2*(1 + (h1param0+h1param2)/3 + (h1param4+h1param0*h1param2)/5
		       + h1param0*h1param4/7);
   //   std::cout <<"h1Norm: " << h1_Norm <<endl;
   //double h2_Norm = (1 + (h2param0+h2param2)/3 + (h2param4+h2param0*h2param2)/5
   //	     + h2param0*h2param4/7);
   double h2_Norm = 1+h2param2*(1/3-1/5);
   //   std::cout <<"h2Norm: " << h2_Norm<<endl;
   double hs_Norm = 2*(1 + hsparam2/3 + hsparam4/5);
   //   std::cout <<"hsNorm: " << hs_Norm <<endl;
   double phi_Norm = 2*3.1415*pacca0;
   //   std::cout <<"phiNorm: " << phi_Norm <<endl;
   double phi1_Norm= 2*3.1415*p1acca0;
   //   std::cout <<"phi1Norm: " << phi1_Norm <<endl;

   if(h2Func<0){ 
     //     std::cout << "JUST SO YOU KNOW: h2FUNC was less than zero, specifically: " << h2Func << " at h2: " << h2 << " and mzz: " << mZZ << std::endl;
     //h2Func=.0001;
   }
   if(hsFunc<0){
     // std::cout << "JUST SO YOU KNOW: hsFUNC was less than zero, specifically: " << hsFunc << " at hs: " << hs << " and mzz: " << mZZ << std::endl;
     //hsFunc=.0001;
   }  
   if(h1Func<0){
     //std::cout << "JUST SO YOU KNOW: h1FUNC was less than zero, specifically: " << h1Func << std::endl;
     //h1Func=.0001;
   }
   if(phiFunc<0){
     // std::cout << "JUST SO YOU KNOW: phiFUNC was less than zero, specifically: " << phiFunc << std::endl;
     //phiFunc=.0001;
   }
   if(phi1Func<0){
     //std::cout << "JUST SO YOU KNOW: phi1FUNC was less than zero, specifically: " << phi1Func << std::endl;
     // cout << "phi1: " << phi1 << "hs: " << hs << endl;
     //phi1Func=.0001;
   }   
   //std::cout << "Bkg Prod: " << h1Func*h2Func*hsFunc*phiFunc*phi1Func/(h1_Norm*h2_Norm*hs_Norm*phi_Norm*phi1_Norm) << endl;
   return(h1Func*h2Func*hsFunc*phiFunc*phi1Func/(h1_Norm*h2_Norm*hs_Norm*phi_Norm*phi1_Norm));
 } 



void RooBkgd2L2JV2::setVars(vector<double> allVars){
  if(int(allVars.size())<6){
    std::cout<<"ERROR from RooBkgd2L2JV2::setVars. Size of vector with measurables is less than 6: "<<allVars.size()<<std::endl;
  }

  h1=allVars.at(0);
  h2=allVars.at(1);
  hs=allVars.at(2);
  phi=allVars.at(3);
  phi1=allVars.at(4);
  mZZ=allVars.at(5);
}


void RooBkgd2L2JV2::setVars(double newcostheta1, double newcostheta2,
			    double newcosthetas, double newphi,
			    double newphi1, double newmZZ){
  h1=newcostheta1;
  h2=newcostheta2;
  hs=newcosthetas;
  phi=newphi;
  phi1=newphi1;
  mZZ=newmZZ;
}

void RooBkgd2L2JV2::setPars(vector<double> allPars){
  if(int(allPars.size())<20){
    std::cout<<"ERROR in RooBkgd2L2JV2::setPars! Size of input params less than minimum: " 
	     << allPars.size() << std::endl;
  }

  h1param0=allPars.at(0); 
  h1param2=allPars.at(1);
  h1param4=allPars.at(2);
  h2param0=allPars.at(3);
  h2param2=allPars.at(4);
  h2param4=allPars.at(5); 
  h2g=allPars.at(6);
  h2cutOff=allPars.at(7);
  hsparam2=allPars.at(8);
  hsparam4=allPars.at(9);
  pacca0=allPars.at(13);   
  pacca1=allPars.at(14);
  pacca2=allPars.at(15);
  p1acca0=allPars.at(17);
  p1acca1=allPars.at(18);
  p1acca2=allPars.at(19);
	    
}  


void RooBkgd2L2JV2::SetParameters(){

  h1param0_0= 0.575536;
  h1param0_1=0.000852046;
  h1param0_2=0.;

  h1param2_0=0.382343;
  h1param2_1=-0.000626345;
  h1param2_2=0.;

  h1param4_0= -1.07047;
  h1param4_1= 0.000670441;
  h1param4_2=0.;

  h2param0_0= 2.93756;
  h2param0_1= -0.0239946;
  h2param0_2= 3.61955e-05;

  h2param2_0= -1.01298;
  h2param2_1= .0047841;
  h2param2_2= 0.;

  h2param4_0= -0.357522;
  h2param4_1= -0.00109411;
  h2param4_2=0.;

  if(mZZ>=250){
    h2g_0= .0594631;
    h2g_1= -.0000696859;
    h2g_2=0.0;
  }else{
    h2g_0= .488962;
    h2g_1= -.00177202;
    h2g_2= 0.0;
  }

  h2cutOff_0= .692987;
  h2cutOff_1= .000206678;
  h2cutOff_2=0.;

  hsparam2_0= 1.84553;
  hsparam2_1= -0.0100222;
  hsparam2_2=0.;

  hsparam4_0=-7.39537;
  hsparam4_1= 0.0403784;
  hsparam4_2=0.;

  pacca0_0= 14.6416;
  pacca0_1= 0.010085;
  pacca0_2=0.;

  pacca1_0= 0.107773;
  pacca1_1= -0.000254654;
  pacca1_2=0.;

  pacca2_0= -1.78714;
  pacca2_1= 0.0037279;
  pacca2_2=0.;

  p1acca0_0 = 100.;
  p1acca0_1 = 0.;
  p1acca0_2 = 0.;

  p1acca1_0= 4.39047;
  p1acca1_1= -0.0254559;
  p1acca1_2= 3.56541e-05;

  p1acca2_0= 11.6347;
  p1acca2_1= -0.0750686;
  p1acca2_2= 8.39119e-05;

  p1acca0 = p1acca0_0 + p1acca0_1*mZZ + p1acca0_2*mZZ*mZZ;
  p1acca1 = p1acca1_0 + p1acca1_1*mZZ + p1acca1_2*mZZ*mZZ;
  p1acca2 = p1acca2_0 + p1acca2_1*mZZ + p1acca2_2*mZZ*mZZ;

  pacca0 = pacca0_0 + pacca0_1*mZZ + pacca0_2*mZZ*mZZ; 
  pacca1 = pacca1_0 + pacca1_1*mZZ + pacca1_2*mZZ*mZZ; 
  pacca2 = pacca2_0 + pacca2_1*mZZ + pacca2_2*mZZ*mZZ; 

  hsparam2 = hsparam2_0 + hsparam2_1*mZZ + hsparam2_2*mZZ*mZZ;
  hsparam4 = hsparam4_0 + hsparam4_1*mZZ + hsparam4_2*mZZ*mZZ;

  h2cutOff = h2cutOff_0 + h2cutOff_1*mZZ + h2cutOff_2*mZZ*mZZ;

  if(mZZ>=700)
    h2g = 0.01068297;
  else
    h2g = h2g_0 + h2g_1*mZZ + h2g_2*mZZ*mZZ;

  h2param0 = h2param0_0 + h2param0_1*mZZ + h2param0_2*mZZ*mZZ;
  h2param2 = h2param2_0 + h2param2_1*mZZ + h2param2_2*mZZ*mZZ;
  h2param4 = h2param4_0 + h2param4_1*mZZ + h2param4_2*mZZ*mZZ;

  h1param0 = h1param0_0 + h1param0_1*mZZ + h1param0_2*mZZ*mZZ;
  h1param2 = h1param2_0 + h1param2_1*mZZ + h1param2_2*mZZ*mZZ;
  h1param4 = h1param4_0 + h1param4_1*mZZ + h1param4_2*mZZ*mZZ;

}
