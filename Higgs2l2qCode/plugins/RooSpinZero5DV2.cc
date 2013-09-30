#include <Riostream.h>
#include <math.h>
#include "../interface/RooSpinZero5DV2.h"

using namespace std;


RooSpinZero5DV2::RooSpinZero5DV2():
  name_("Higgs"),
  title_("HiggsTitle"),
  h1(0.0),
  h2(0.0),
  hs(0.0),
  Phi(0.0),
  Phi1(0.0),
  mZZ(0.0),
  fppVal(0.0),
  fmmVal(0.0),
  fpmVal(0.0),
  fp0Val(0.0),
  f0mVal(0.0),
  phippVal(0.0),
  phimmVal(0.0),
  phipmVal(0.0),
  phip0Val(0.0),
  phi0mVal(0.0),
  fz1Val(0.0),
  fz2Val(0.0),
  R1Val(0.0),
  R2Val(0.0),
  para2(0.0),
  para4(0.0),
  acca0(0.0),
  acca1(0.0),
  acca2(0.0),
  a2(0.0),
  a4(0.0),
  cutOff(0.0),
  g(0.0),
  b2(0.0),
  b4(0.0),
  N(0.0)
{
}



RooSpinZero5DV2::RooSpinZero5DV2(const char *name, const char *title, vector<double> allVars, vector<double> allPars):
  name_(name),
  title_(title){
  if(!(allVars.size()==6)){
    std::cout << "ERROR RooSpinZero5DV2::RooSpinZero5DV2, number of measurables passed is wrong :-/" << std::endl;
    setVars(allVars.at(0),allVars.at(1),
	    allVars.at(2),allVars.at(3),
	    allVars.at(4),allVars.at(5));
  }
  setParams(allPars);
}


RooSpinZero5DV2::RooSpinZero5DV2(const char *name, const char *title, 
				  double _h1,
				  double _h2,
                                  double _hs,
				 double _Phi,
				  double _Phi1,
			          double _mZZ,
				  double _fppVal,
				  double _fmmVal,
				  double _fpmVal,
				  double _fp0Val,
				  double _f0mVal,
				  double _phippVal,
				  double _phimmVal,
				  double _phipmVal,
				  double _phip0Val,
				  double _phi0mVal,
				  double _fz1Val,
				  double _fz2Val,
				  double _R1Val,
				  double _R2Val,
				  double _para2,
				  double _para4,
				  double _acca0,
				  double _acca1,
				  double _acca2,
				  double _a2,
				  double _a4,
				  double _cutOff,
				  double _g,
                                  double _b2,
				  double _b4,
				  double _N):
   name_(name),
   title_(title),
   h1(_h1),
   h2(_h2),
   hs(_hs),
   Phi(_Phi),
   Phi1(_Phi1),
   mZZ(_mZZ),
   fppVal(_fppVal),
   fmmVal(_fmmVal),
   fpmVal(_fpmVal),
   fp0Val(_fp0Val),
   f0mVal(_f0mVal),
   phippVal(_phippVal),
   phimmVal(_phimmVal),
   phipmVal(_phipmVal),
   phip0Val(_phip0Val),
   phi0mVal(_phi0mVal),
   fz1Val(_fz1Val),
   fz2Val(_fz2Val),
   R1Val(_R1Val),
   R2Val(_R2Val),
   para2(_para2),
   para4(_para4),
   acca0(_acca0),
   acca1(_acca1),
   acca2(_acca2),
   a2(_a2),
   a4(_a4),
   cutOff(_cutOff),
   g(_g),
   b2(_b2),
   b4(_b4),
   N(_N)
{ 
} 


 RooSpinZero5DV2::RooSpinZero5DV2(const RooSpinZero5DV2& other, const char* name) :  
   name_(name),
   h1(other.h1),
   h2(other.h2),
   hs(other.hs),
   Phi(other.Phi),
   Phi1(other.Phi1),
   mZZ(other.mZZ),
   fppVal(other.fppVal),
   fmmVal(other.fmmVal),
   fpmVal(other.fpmVal),
   fp0Val(other.fp0Val),
   f0mVal(other.f0mVal),
   phippVal(other.phippVal),
   phimmVal(other.phimmVal),
   phipmVal(other.phipmVal),
   phip0Val(other.phip0Val),
   phi0mVal(other.phi0mVal),
   fz1Val(other.fz1Val),
   fz2Val(other.fz2Val),
   R1Val(other.R1Val),
   R2Val(other.R2Val),
   para2(other.para2),
   para4(other.para4),
   acca0(other.acca0),
   acca1(other.acca1),
   acca2(other.acca2),
   a2(other.a2),
   a4(other.a4),
   cutOff(other.cutOff),
   g(other.g),
   b2(other.b2),
   b4(other.b4),
   N(other.N)
 {
 } 



double RooSpinZero5DV2::evaluate() const
 { 

   


   double shs = sqrt(1-hs*hs);
   double sh1 = sqrt(1-h1*h1);
   double sh2 = sqrt(1-h2*h2);
   
   if ((1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal) < 0) return 1e-9;
   
   double term1Coeff = (2.-2.*fz1Val+fz2Val-6.*(2.-4.*fz1Val-fz2Val)*pow(hs,2)+3.*(6.-10.*fz1Val-5.*fz2Val)*pow(hs,4));
   double term1A = 4.*(1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal)*pow(sh1,2)*pow(sh2,2);
   double term1B = (fppVal+fmmVal)*((1.+h1*h1)*(1.+h2*h2)+4.*R1Val*R2Val*h1*h2);
   double term1C = -2.*(fppVal-fmmVal)*(R1Val*h1*(1.+h2*h2)+R2Val*h2*(1.+h1*h1));
   double term1D = 4.*sqrt(fppVal*(1-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*(R1Val-h1)*(R2Val-h2)*sh1*sh2*cos(Phi+phippVal);
   double term1E = 4.*sqrt(fmmVal*(1-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*(R1Val+h1)*(R2Val+h2)*sh1*sh2*cos(Phi-phimmVal);
   double term1F = 2.*sqrt(fppVal*fmmVal)*pow(sh1,2)*pow(sh2,2)*cos(2.*Phi+phippVal-phimmVal);
   double term1 = term1Coeff*(term1A+term1B+term1C+term1D+term1E+term1F);
   
   double term2Coeff = 8.*(fz1Val+fz2Val+3.*(2.-3.*fz1Val-2.*fz2Val)*pow(hs,2)-(6.-10.*fz1Val-5.*fz2Val)*pow(hs,4));
   double term2A = (fp0Val+f0mVal)*(1.-h1*h1*h2*h2) - (fp0Val-f0mVal)*(R1Val*h1*pow(sh2,2)+R2Val*h2*pow(sh1,2));
   double term2B = 2.*sqrt(fp0Val*f0mVal)*sh1*sh2*(R1Val*R2Val-h1*h2)*cos(Phi+phip0Val-phi0mVal);
   double term2 = term2Coeff*(term2A+term2B);
   
   double term3Coeff = -8.*(fz1Val-fz2Val+(6.-10.*fz1Val-5.*fz2Val)*pow(hs,2))*pow(shs,2)*sh1*sh2*cos(Phi + 2.*Phi1);
   double term3A = (fp0Val+f0mVal)*(R1Val*R2Val+h1*h2)-(fp0Val-f0mVal)*(R1Val*h2+R2Val*h1)+2.*sqrt(fp0Val*f0mVal)*sh1*sh2*cos(Phi+phip0Val-phi0mVal);
   double term3 = term3Coeff*term3A;
   
   double term4Coeff = 6.-2.*fz1Val-5.*fz2Val-6.*(2.-2.*fz1Val-3.*fz2Val)*pow(hs,2)+(6.-10.*fz1Val-5.*fz2Val)*pow(hs,4);
   double term4A = fpmVal*((1.+h1*h1)*(1.+h2*h2)-4.*R1Val*R2Val*h1*h2);
   double term4 = term4Coeff*term4A;
   
   double term5 = pow(shs,4)*(6.-10.*fz1Val-5.*fz2Val)*fpmVal*pow(sh1,2)*pow(sh2,2)*cos(2.*Phi+4.*Phi1);
   
   double term6Coeff = (-1.)*sqrt(6.)*(2.-2.*fz1Val-3.*fz2Val-(6.-10.*fz1Val-5.*fz2Val)*pow(hs,2))*pow(shs,2);
   double term6A = 2.*sqrt(fpmVal*(1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*sh1*sh2*((R1Val-h1)*(R2Val+h2)*cos(Phi+2.*Phi1-phipmVal)+(R1Val+h1)*(R2Val-h2)*cos(Phi+2.*Phi1+phipmVal));
   double term6B = sqrt(fpmVal*fmmVal)*(sh1*sh1*(1.+2.*R2Val*h2+h2*h2)*cos(2.*Phi1-phipmVal+phimmVal)+sh2*sh2*(1.+2.*R1Val*h1+h1*h1)*cos(2.*Phi+2.*Phi1+phipmVal-phimmVal));
   double term6C = sqrt(fppVal*fpmVal)*(sh1*sh1*(1.-2.*R2Val*h2+h2*h2)*cos(2.*Phi1+phipmVal-phippVal)+sh2*sh2*(1.-2.*R1Val*h1+h1*h1)*cos(2.*Phi+2.*Phi1-phipmVal+phippVal));
   double term6 = term6Coeff*(term6A+term6B+term6C);
   
   /// new mixing terms
   double term7Coeff = -4.*sqrt(3.)*(2.-4.*fz1Val-fz2Val-(6.-10.*fz1Val-5.*fz2Val)*pow(hs,2))*hs*shs;
   double term7A = sqrt(fmmVal*f0mVal)*(sh1*(R1Val+h1)*(1.+2.*R2Val*h2+h2*h2)*cos(Phi1+phimmVal-phi0mVal)-sh2*(R2Val+h2)*(1.+2.*R1Val*h1+h1*h1)*cos(Phi+Phi1-phimmVal+phi0mVal));
   double term7B = sqrt(fmmVal*fp0Val)*(sh1*sh1*sh2*(R2Val+h2)*cos(Phi-Phi1-phimmVal+phip0Val)-sh2*sh2*sh1*(R1Val+h1)*cos(2.*Phi+Phi1-phimmVal+phip0Val));
   double term7C = -sqrt(fppVal*f0mVal)*(sh1*sh1*sh2*(R2Val-h2)*cos(Phi-Phi1+phippVal-phi0mVal)-sh2*sh2*sh1*(R1Val-h1)*cos(2.*Phi+Phi1+phippVal-phi0mVal));
   double term7D = -sqrt(fppVal*fp0Val)*(sh1*(R1Val-h1)*(1.-2.*R2Val*h2+h2*h2)*cos(Phi1+phip0Val-phippVal)-sh2*(R2Val-h2)*(1.-2.*R1Val*h1+h1*h1)*cos(Phi+Phi1+phippVal-phip0Val));
   double term7E = -2.*sqrt(f0mVal*(1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*(sh1*sh2*sh2*(R1Val+h1)*cos(Phi1+phi0mVal)-sh2*sh1*sh1*(R2Val+h2)*cos(Phi+Phi1-phi0mVal));
   double term7F = 2.*sqrt(fp0Val*(1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*(sh1*sh2*sh2*(R1Val-h1)*cos(Phi1-phip0Val)-sh2*sh1*sh1*(R2Val-h2)*cos(Phi+Phi1+phip0Val));
   double term7 = term7Coeff*(term7A+term7B+term7C+term7D+term7E+term7F);
   
   double term8Coeff = 2.*sqrt(2.)*hs*shs*(6.-6.*fz1Val-9.*fz2Val-(6.-10.*fz1Val-5.*fz2Val)*pow(hs,2));
   double term8A = sqrt(fpmVal*f0mVal)*(sh1*(R1Val-h1)*(1.+2.*R2Val*h2+h2*h2)*cos(Phi1-phipmVal+phi0mVal)-sh2*(R2Val-h2)*(1.+2.*R1Val*h1+h1*h1)*cos(Phi+Phi1+phipmVal-phi0mVal));
   double term8B = sqrt(fpmVal*fp0Val)*(sh2*(R2Val+h2)*(1.-2.*R1Val*h1+h1*h1)*cos(Phi+Phi1-phipmVal+phip0Val)-sh1*(R1Val+h1)*(1.-2.*R2Val*h2+h2*h2)*cos(Phi1+phipmVal-phip0Val));
   double term8 = term8Coeff*(term8A+term8B);
   
   double term9Coeff = -2.*sqrt(2.)*hs*pow(shs,3)*(6.-10.*fz1Val-5.*fz2Val);
   double term9A = sqrt(fpmVal*f0mVal)*(sh1*sh1*sh2*(R2Val+h2)*cos(Phi+3.*Phi1-phipmVal+phi0mVal)-sh2*sh2*sh1*(R1Val+h1)*cos(2.*Phi+3.*Phi1+phipmVal-phi0mVal));
   double term9B = sqrt(fpmVal*fp0Val)*(sh1*sh2*sh2*(R1Val-h1)*cos(2.*Phi+3.*Phi1-phipmVal+phip0Val)-sh2*sh1*sh1*(R2Val-h2)*cos(Phi+3.*Phi1+phipmVal-phip0Val));
   double term9 = term9Coeff*(term9A+term9B);
   
   // signs of the interference terms are flipped!!!!
   if (true){
     term7 *= (-1.);
     term8 *= (-1.);
     term9 *= (-1.);
   }
   
   double sum = term1+term2+term3+term4+term5+term6+term7+term8+term9;

   double hs_accp = (1.0+para2*pow(hs,2)+para4*pow(hs,4));
   double h1_accp = (1+b2*pow(h1,2)+b4*pow(h1,4));
   double p1_accp = (acca0+acca1*cos(Phi1)+acca2*cos(2.0*Phi1));
   double h2_accp = (1.0+a2*pow(h2,2)+a4*pow(h2,4))/(1.0+exp((h2-cutOff)/g));
   double hs_Norm = 2*(1 + para2/3 + para4/5 );
   double h1_Norm = 2*(1+b2/3+b4/5);
   double h2_Norm = 1+a2/3+a4/5;
   double p1_Norm = 2*3.1415*acca0;

   double accp = hs_accp*h1_accp*h2_accp*p1_accp/(hs_Norm*h1_Norm*h2_Norm*p1_Norm);

   //double ResidNorm = .323*exp(-.007*mZZ);
   double ResidNorm = .04;

   if(sum*accp<=0){sum=accp=.0001;}

   return sum*accp*ResidNorm ;
 } 

void RooSpinZero5DV2::SetAcceptanceParameters(){

  para2_0= -0.819042;
  para2_1= 0.0016744;
  para2_2=0.;

  para4_0= 1.0512;
  para4_1= -0.00317063;
  para4_2=0.;

  acca0_0=990.;
  acca0_1=0.;
  acca0_2=0.;

  acca1_0= -45.5558;
  acca1_1= 0.0945721;
  acca1_2=0.;

  acca2_0= -173.089;
  acca2_1= 0.172187;
  acca2_2=0.;

  a2_0= -2.36442;
  a2_1= 0.00556876;
  a2_2=0.;

  a4_0= 1.50929;
  a4_1= -0.00584568;
  a4_2=0.;

  cutOff_0= 1.31613;
  cutOff_1= -0.000331295;
  cutOff_2=0.;

  g_0= 0.00503465;
  g_1= -5.96656e-06;
  g_2= 0.;

  b2_0=.13;
  b2_1=0.;
  b2_2=0.;

  b4_0= -1.02326;
  b4_1= 0.00136598;
  b4_2=0.;

  double gamma=mZZ*mZZ/(2*91.1876*91.1876)-1;
  fmmVal=fppVal=1/(gamma*gamma+2);

  para2 = para2_0+para2_1*mZZ+para2_2*mZZ*mZZ;
  para4 = para4_0+para4_1*mZZ+para4_2*mZZ*mZZ;
  acca0 = acca0_0+acca0_1*mZZ+acca0_2*mZZ*mZZ;
  acca1 = acca1_0+acca1_1*mZZ+acca1_2*mZZ*mZZ;
  acca2 = acca2_0+acca2_1*mZZ+acca2_2*mZZ*mZZ;
  a2 = a2_0+a2_1*mZZ+a2_2*mZZ*mZZ;
  a4 = a4_0+a4_1*mZZ+a4_2*mZZ*mZZ;
  cutOff = cutOff_0+cutOff_1*mZZ+cutOff_2*mZZ*mZZ;
  b2 = b2_0+b2_1*mZZ+a2_2*mZZ*mZZ;
  b4 = b4_0+b4_1*mZZ+b4_2*mZZ*mZZ;

}

void RooSpinZero5DV2::setVars(double newcostheta1, double newcostheta2,
			      double newcosthetas,double newphi,
		       double newphi1, double newmZZ){
  h1=newcostheta1;
  h2=newcostheta2;
  hs=newcosthetas;
  Phi=newphi;
  Phi1=newphi1;
  mZZ=newmZZ;

}


void RooSpinZero5DV2::setVars(vector<double> allVars){
  if(int(allVars.size())<6){
    std::cout<<"ERROR from RooSpinZero5DV2::setVars. Size of vector with measurables is less than 6: "<<allVars.size()<<std::endl;
  }

  h1=allVars.at(0);
  h2=allVars.at(1);
  hs=allVars.at(2);
  Phi=allVars.at(3);
  Phi1=allVars.at(4);
  mZZ=allVars.at(5);
}


void RooSpinZero5DV2::setParams(vector<double> allPars){

  if(int(allPars.size())<29){
    std::cout << "ERROR from RooSpinZero5DV2::setParams. Size of vector with params is less the minimum: "
	      << allPars.size() << std::endl;
  }

  fppVal=allPars.at(0);
  fmmVal=allPars.at(1);
  fpmVal=allPars.at(2);
  fp0Val=allPars.at(3);
  f0mVal=allPars.at(4);
  phippVal=allPars.at(5);
  phimmVal=allPars.at(6);
  phipmVal=allPars.at(7);
  phip0Val=allPars.at(8);
  phi0mVal=allPars.at(9);
  fz1Val=allPars.at(10);
  fz2Val=allPars.at(11);
  R1Val=allPars.at(12);
  R2Val=allPars.at(13);
  para2=allPars.at(14);
  para4=allPars.at(15);
  acca0=allPars.at(18);
  acca1=allPars.at(19);
  acca2=allPars.at(20);
  a2=allPars.at(22);
  a4=allPars.at(23);
  cutOff=allPars.at(24);
  g=allPars.at(25);
  b2=allPars.at(26);
  b4=allPars.at(27);
  N=allPars.at(28);
}
