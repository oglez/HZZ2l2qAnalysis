#ifndef _HELICITY_H_
#define _HELICITY_H_
#include "TLorentzVector.h"
#include <iostream>

class Helicity {
   public:
      Helicity() {};
      ~Helicity() {};

      //      void calculateAngles(TLorentzVector, TLorentzVector, TLorentzVector,
      //           TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector,
      //           double&, double&, double&, double&, double&, double&,
      //           double&, double&, double&, bool&);

      void calculateAngles(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, 
			   double&, double&, double&, double&, double&);
};


#endif
