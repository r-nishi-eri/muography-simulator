#ifndef _GENERATOR_
#define _GENERATOR_ 1

#include "globals.hh"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <string>

using namespace std;

class Generator{
public:
 
  double getGamma(double, double);
  double getDifferentialFlux(double, double);
  
  // FILE MODE INJECTION
  Generator(double radius_xy, double radius_z, double emin, double emax, 
	    int parcode, int parcharge, 
	    string fnameflux, int numcz, int nummom, double coszstart, double dcosz);
  // ISOTROPIC INJECTION
  Generator(double radius_xy, double radius_z, double e, int parcode, int parcharge);
  // COSINE INJECTION
  Generator(double radius_xy, double radius_z, double emin, double emax, int parcode, int parcharge, double gammaindex, double nindex);

  int shoot(double, double);
  double getRandomizedTheta(){return RandomizedTheta;};
  double getRandomizedEnergy(){return RandomizedEnergy;};
  double getNumberOfInjection() {return num_shoot;};
  double getElapsedTime();
  double getTotalIntensity();
  int getParCode() {return parCode;};
  int getParCharge() {return parCharge;};
  void setDifluxAt1GeV(double flux){diFluxAt1GeV = flux;};
  void setTotalFlux(double flux);

protected:
  void read_flux_file();
  void calc_table_iso();
  void calc_table_cos();
  void calc_table_file();

  // contents of energy spectrum file
  static const int num_cz_MAX  = 100;
  static const int num_mom_MAX = 200;
  int num_cz, num_mom;
  double mom     [num_cz_MAX][num_mom_MAX];
  double diflux  [num_cz_MAX][num_mom_MAX];
  //double mu_mins [num_cz][num_mom];

  // tables for MC generation
  static const int numThetaDiv = 100;
  double acceptances[numThetaDiv];
  double thetaMins  [numThetaDiv];
  double thetaMaxs  [numThetaDiv];
  double sumTable;
  double skyOblateness;
  double radiusXY, radiusZ;
  double sum_acceptance;
  static const int numEnergyDiv = 100;
  double energyMin, energyMax;
  double energies[numEnergyDiv];
  double f [numThetaDiv][numEnergyDiv];
  double g [numThetaDiv];
  double sum_g;
  
  // Shooting
  double RandomizedTheta;
  double RandomizedEnergy;
  int num_shoot;

  string modeInjection;
  int parCode;
  int parCharge;

  string fNameFlux;
  int numFluxFile;
  double coszMax;
  double coszMin;
  double dCosz;

  // n index for cos^{n} type injection.
  double nIndex;
  // gamma index for E^{-gamma} type injection.
  double gammaIndex;
  // differential flux at 1 GeV
  double diFluxAt1GeV; // unit is m^{-2} sr^{-1} sec^{-1} GeV^{gamma-1}
  // acceptance * integrated flux
  double acceptance_integrated_flux;

};


#endif
