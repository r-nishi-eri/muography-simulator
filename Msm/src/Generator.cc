#include "Generator.hh"
#include <iostream>
using namespace std;


// Constructor for FILE MODE
Generator::Generator(double radius_xy, double radius_z, double emin, double emax, 
		     int parcode, int parcharge, 
		     string fnameflux, int numcz, int nummom, double coszstart, double dcosz)
  : skyOblateness (radius_z / radius_xy),
    radiusXY (radius_xy),
    radiusZ (radius_z),
    energyMin (emin),
    energyMax (emax),
    parCode(parcode),
    parCharge(parcharge),
    fNameFlux(fnameflux),
    num_cz(numcz),
    num_mom(nummom),
    coszMax(coszstart),
    dCosz(dcosz),
    coszMin(coszstart - double(numcz) * dcosz)
{
  modeInjection = "file";
  read_flux_file();
  calc_table_file();
  num_shoot = 0;

}

// Constructor for ISOTROPIC INJECTION
Generator::Generator(double radius_xy, double radius_z, double e, int parcode, int parcharge)
    : skyOblateness (radius_z / radius_xy),
      radiusXY (radius_xy),
      radiusZ  (radius_z),
      energyMin (e),
      parCode(parcode),
      parCharge(parcharge)      
{
  modeInjection = "iso";
  calc_table_iso();
  num_shoot = 0;
}

// Constructor for COS^{n} INJECTION
Generator::Generator(double radius_xy, double radius_z, double emin, double emax, int parcode, int parcharge, double gammaindex, double nindex)
  : skyOblateness(radius_z / radius_xy),
    radiusXY (radius_xy),
    radiusZ  (radius_z),
    energyMin (emin),
    energyMax (emax),
    parCode(parcode),
    parCharge(parcharge),
    gammaIndex(gammaindex),
    nIndex (nindex)
{
  modeInjection = "cos";
  calc_table_cos();    
  num_shoot = 0;
}


void Generator::read_flux_file(){

  for (int icz = 0; icz < num_cz; icz++){
    char czfilename[256];
    sprintf(czfilename, "%s%02d.flx", fNameFlux.c_str(), icz+1);
    ifstream fin(czfilename);
    if (fin.eof()){
      cerr << "Cannot find " << czfilename << endl;
      exit(1);
    }

    for(int iLine = 0; iLine < num_mom_MAX; iLine++){
      mom   [icz][iLine] = 0;
      diflux[icz][iLine] = 0;
    }
    
    double p, flx;
    for(int iLine = 0; iLine < num_mom; iLine++){
      fin >> p >> flx;
      mom    [icz][iLine] = p;
      diflux [icz][iLine] = flx;
      cout << p << " " << flx << endl;
    }
  }

  return;
}


void Generator::calc_table_iso()
{
  ///////////////////////////////////////////
  // acceptance table (area * solid angle) //
  ///////////////////////////////////////////
  sum_acceptance = 0;
  for (int iThetaDiv = 0; iThetaDiv < numThetaDiv; iThetaDiv++){
    double thetaMin = M_PI  * double(iThetaDiv)   / double(numThetaDiv);
    double thetaMax = M_PI  * double(iThetaDiv+1) / double(numThetaDiv); 
    double thetaCenter = (thetaMin + thetaMax) / 2.0;
    double sinTheta = sin(thetaCenter);
    double cosTheta = cos(thetaCenter);
    // acceptance unit: m^2 sr^1
    double acceptance = 2.0 * M_PI * M_PI * radiusXY * radiusZ * sinTheta * (thetaMax - thetaMin) * sqrt(sinTheta*sinTheta+cosTheta*cosTheta/skyOblateness/skyOblateness); 
    thetaMins[iThetaDiv]   = thetaMin;
    thetaMaxs[iThetaDiv]   = thetaMax;
    sum_acceptance += acceptance;
    acceptances[iThetaDiv] = sum_acceptance;
  }

}

void Generator::calc_table_cos()
{
  ///////////////////////////////////////////
  // acceptance table (area * solid angle) //
  ///////////////////////////////////////////
  sum_acceptance = 0;
  for (int iThetaDiv = 0; iThetaDiv < numThetaDiv; iThetaDiv++){
    double thetaMin = M_PI / 2.0 * double(iThetaDiv)   / double(numThetaDiv);
    double thetaMax = M_PI / 2.0 * double(iThetaDiv+1) / double(numThetaDiv); 
    double thetaCenter = (thetaMin + thetaMax) / 2.0;
    double sinTheta = sin(thetaCenter);
    double cosTheta = cos(thetaCenter);
    // acceptance unit: m^2 sr^1 
    double acceptance = 2.0 * M_PI * M_PI * pow(cosTheta, nIndex) * radiusXY * radiusZ * sinTheta * (thetaMax - thetaMin) * sqrt(sinTheta*sinTheta+cosTheta*cosTheta/skyOblateness/skyOblateness);
    thetaMins[iThetaDiv]   = thetaMin;
    thetaMaxs[iThetaDiv]   = thetaMax;
    sum_acceptance += acceptance;
    acceptances[iThetaDiv] = sum_acceptance;
  }
  
}

void Generator::calc_table_file()
{
  ///////////////////////////////////////////
  // acceptance table (area * solid angle) //
  ///////////////////////////////////////////
  double thetaMinLimit = acos(coszMax);
  double thetaMaxLimit = acos(coszMin);
  for (int iThetaDiv = 0; iThetaDiv < numThetaDiv; iThetaDiv++){
    double thetaMin = (thetaMinLimit * double(numThetaDiv - iThetaDiv) + thetaMaxLimit * double(iThetaDiv)) / double(numThetaDiv);
    double thetaMax = (thetaMinLimit * double(numThetaDiv - iThetaDiv - 1) + thetaMaxLimit * double(iThetaDiv + 1)) / double (numThetaDiv);
    double thetaCenter = (thetaMin + thetaMax) / 2.0;
    double sinTheta = sin(thetaCenter);
    double cosTheta = cos(thetaCenter);
    // acceptance unit: m^2 sr^1
    double acceptance = 2.0 * M_PI * M_PI * radiusXY * radiusZ * sinTheta * (thetaMax - thetaMin) * sqrt(sinTheta*sinTheta+cosTheta*cosTheta/skyOblateness/skyOblateness); 
    acceptances[iThetaDiv] = acceptance;
    thetaMins[iThetaDiv]   = thetaMin;
    thetaMaxs[iThetaDiv]   = thetaMax;
  }

  ///////////////////////////
  // integrated flux table //
  ///////////////////////////
  for (int iEnergyDiv = 0; iEnergyDiv < numEnergyDiv; iEnergyDiv++){
    double ene = exp(log(energyMin) + double(iEnergyDiv) / double(numEnergyDiv) * (log(energyMax)-log(energyMin)));
    energies[iEnergyDiv] = ene;
  }
  sum_g = 0;
  sum_acceptance = 0;
  for (int iThetaDiv = 0; iThetaDiv < numThetaDiv; iThetaDiv++){
    double theta = (thetaMins[iThetaDiv] + thetaMaxs[iThetaDiv]) * 0.5;
    f[iThetaDiv][0] = 0.0;
    for(int iEnergyDiv = 1; iEnergyDiv < numEnergyDiv; iEnergyDiv++){
      f[iThetaDiv][iEnergyDiv] = f[iThetaDiv][iEnergyDiv-1] + getDifferentialFlux(energies[iEnergyDiv-1], theta) * (energies[iEnergyDiv]-energies[iEnergyDiv-1]); 
    }
    double dg = acceptances[iThetaDiv] * f[iThetaDiv][numEnergyDiv-1];
    sum_acceptance += acceptances[iThetaDiv];
    sum_g += dg;
    g[iThetaDiv] = sum_g;
    cerr << iThetaDiv << " " 
         << dg << " " 
         << sum_g << endl;
  }
}


int Generator::shoot(double rnd1, double rnd2){

  //////////////////
  // Choose Theta //
  //////////////////
  int iThetaChosen;
  if (modeInjection == "iso" || modeInjection == "cos"){ // ISO or COS INJECTION
    double arnd = rnd1 * sum_acceptance;
    for(int iTable = 0; iTable < numThetaDiv; iTable++){
      if (arnd <= acceptances[iTable]){
	double amin = (iTable >= 1 ? acceptances[iTable-1] : 0);
	double amax = acceptances[iTable];
	RandomizedTheta = (thetaMins[iTable] * (amax - arnd) + thetaMaxs[iTable] * (arnd - amin)) / (amax - amin);
	break;
      }
    }
  }
  else if (modeInjection == "file"){ // FILE INJECTION
    // Choose Theta
    double grnd = rnd1 * sum_g; 
    for(int iTable = 0; iTable < numThetaDiv; iTable++){
      if (grnd <= g[iTable]){
	double gmin = (iTable >= 1 ? g[iTable-1]: 0);
	double gmax = g[iTable];
	RandomizedTheta = (thetaMins[iTable] * (gmax - grnd) + thetaMaxs[iTable] * (grnd - gmin)) / (gmax - gmin);
	iThetaChosen = iTable;
	break;
      }
    }
  }
  

  ///////////////////
  // Choose energy //
  ///////////////////
  if (modeInjection == "iso") RandomizedEnergy = energyMin;
  else if (modeInjection == "cos") {
    RandomizedEnergy = pow((1.0-rnd2)*pow(energyMin,-gammaIndex+1.0) + rnd2*pow(energyMax,-gammaIndex+1.0), 1.0/(-gammaIndex+1.0));
  }
  else if (modeInjection == "file") {
    double frnd = rnd2 * f[iThetaChosen][numEnergyDiv-1]; 
    for(int iTable = 1; iTable < numEnergyDiv; iTable++){
      if (frnd <= f[iThetaChosen][iTable]){
	double fmin = f[iThetaChosen][iTable-1];
	double fmax = f[iThetaChosen][iTable  ];
	double u = (frnd - fmin) / (fmax - fmin);
	double gamma = getGamma(energies[iTable-1], RandomizedTheta);
	RandomizedEnergy = pow((1.0-u)*pow(energies[iTable-1],-gamma+1.0) + u*pow(energies[iTable],-gamma+1.0), 1.0/(-gamma+1.0));
	break;
      }
    }    
  }
    
  // counter # particles
  num_shoot++;  
  return 0; 
}


void Generator::setTotalFlux(double flux){
  if (modeInjection == "iso"){
    diFluxAt1GeV = flux;
  }else if (modeInjection == "cos"){
    double integral = 1.0 / double (gammaIndex - 1.0) * 
      (pow(energyMin, -gammaIndex+1.0) - pow(energyMax, -gammaIndex+1.0));
    diFluxAt1GeV = flux / integral;
  }
}



double Generator::getElapsedTime(){
  if (modeInjection == "file"){
    return double(num_shoot) / sum_g;
  }else if (modeInjection == "iso"){
    return double(num_shoot) / acceptances[numThetaDiv-1] / diFluxAt1GeV;
  }else if (modeInjection == "cos"){
    double intFlux = diFluxAt1GeV / double(gammaIndex-1.0) * 
      (pow(energyMin, -gammaIndex+1.0) - pow(energyMax, -gammaIndex+1.0));
    return double(num_shoot) / acceptances[numThetaDiv-1] / intFlux; 
  }
}

double Generator::getTotalIntensity(){
  if (modeInjection == "file"){
    return sum_g;
  }else if (modeInjection == "iso"){
    return acceptances[numThetaDiv-1] * diFluxAt1GeV;
  }else if (modeInjection == "cos"){
    double intFlux = diFluxAt1GeV / double(gammaIndex-1.0) * 
      (pow(energyMin, -gammaIndex+1.0) - pow(energyMax, -gammaIndex+1.0));
    return acceptances[numThetaDiv-1] * intFlux;
  }
}



double Generator::getGamma(double p , double theta){
  
  if (modeInjection == "file"){

    // check theta range
    double cosz = cos(theta);
    if (cosz <= coszMin || coszMax <= cosz) return 2.0;
    int icz = int ((coszMax - cosz) / dCosz);
    // check p range
    if (p < mom[0][0] || mom[0][num_mom-1] < p) return 2.0;
    int ip; for (ip = 0; p < mom[0][ip] || mom[0][ip+1] < p; ip++);
    
    double mom_low   = mom[0][ip];
    double mom_high  = mom[0][ip+1];
    double flux_low  = diflux[icz][ip];
    double flux_high = diflux[icz][ip+1];
    if (flux_low > 0 && flux_high > 0){
      return -(log(flux_low)-log(flux_high)) / (log(mom_low)-log(mom_high));
    }else{
      return 0;
    }
  }

  if (modeInjection == "cos") return gammaIndex;
  else return 2.0;

}




double Generator::getDifferentialFlux(double p, double theta)
{
  // arguments: p (unit: GeV/c or GeV), theta (rad)

  // check theta range
  double cosz = cos(theta);
  if (cosz <= coszMin || coszMax <= cosz) return 0;
  
  int icz;
  if (cosz <= coszMin + dCosz / 2.0){
    icz = num_cz - 1;
  }else if (cosz >= coszMax - dCosz / 2.0){
    icz = 1;
  }else{
    icz = int(((coszMax - dCosz / 2.0) - cosz) / dCosz) + 1;
  }

  // check p range
  if (p < mom[0][0] || mom[0][num_mom-1] < p) return 0;
  int ip; 
  for (ip = 0; p < mom[0][ip] || mom[0][ip+1] < p; ip++);
  double logp = log(p);
  
  double coszmin = coszMax - double(icz+1) * dCosz;
  double coszmax = coszMax - double(icz)   * dCosz;
  double logpmin = log(mom[0][ip  ]);
  double logpmax = log(mom[0][ip+1]);
  
  double flux_coszmin_pmin = diflux[icz  ][ip  ];
  double flux_coszmin_pmax = diflux[icz  ][ip+1];
  double flux_coszmax_pmin = diflux[icz-1][ip  ];
  double flux_coszmax_pmax = diflux[icz-1][ip+1];
  if (flux_coszmin_pmin <= 0 || 
      flux_coszmin_pmax <= 0 ||
      flux_coszmax_pmin <= 0 || 
      flux_coszmax_pmax <= 0    ) return 0; // This is for avoiding log(0) error.

  double log_flux_coszmin_pmin = log(flux_coszmin_pmin);
  double log_flux_coszmin_pmax = log(flux_coszmin_pmax);
  double log_flux_coszmax_pmin = log(flux_coszmax_pmin);
  double log_flux_coszmax_pmax = log(flux_coszmax_pmax);
  
  // first interpolation
  double log_flux_pmin = (log_flux_coszmin_pmin * (coszmax - cosz) + log_flux_coszmax_pmin * (cosz - coszmin)) / (coszmax - coszmin);
  double log_flux_pmax = (log_flux_coszmin_pmax * (coszmax - cosz) + log_flux_coszmax_pmax * (cosz - coszmin)) / (coszmax - coszmin);
  // second interpolation
  double log_flux = (log_flux_pmin * (logpmax-logp) + log_flux_pmax * (logp-logpmin)) / (logpmax - logpmin);
  double flux = exp(log_flux);
  return flux;
}






