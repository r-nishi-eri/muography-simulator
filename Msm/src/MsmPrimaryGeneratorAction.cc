#include "MsmPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "MsmDetectorConstruction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Threading.hh" 




// GLOBAL VARIABLE
G4ThreadLocal G4int primaryParCode;
G4ThreadLocal G4int primaryParCharge;
G4ThreadLocal G4double primaryKE;
G4ThreadLocal G4double primaryX;
G4ThreadLocal G4double primaryY;
G4ThreadLocal G4double primaryZ;
G4ThreadLocal G4double primaryVX;
G4ThreadLocal G4double primaryVY;
G4ThreadLocal G4double primaryVZ;





MsmPrimaryGeneratorAction::MsmPrimaryGeneratorAction(string str_injectionfile, string str_injectionlog)
  : G4VUserPrimaryGeneratorAction(),
    fParticleGun(0),
    numInjection(0),
    modPrint(1000)
{

  G4RunManager* runmanager = G4RunManager::GetRunManager();
  MsmDetectorConstruction* detconst = (MsmDetectorConstruction*)runmanager->GetUserDetectorConstruction();
  
  // OPEN (RECREATE) LOG FILE
  std::stringstream ss; ss << str_injectionlog << "-" << int(G4Threading::G4GetThreadId()) << ".log";
  filenameLog = ss.str();
  std::ofstream fout(ss.str().c_str(), std::ios::out);

  // TEST
  std::ofstream fout_test ("test",std::ios::out);
  fout_test << ss.str() << endl;
  
  // READ INJECTION PARAMETER FILE
  std::ifstream fin(str_injectionfile.c_str());
  if (!fin){
    std::cerr << "file " << str_injectionfile.c_str() << " could not be opened." << std::endl;
    exit(1);
  }

  std::string buf;
  while(fin && std::getline(fin, buf)) {
    
    char str_command[256];
    char str_value[256];
    sscanf(buf.c_str(), "%s %s\n", str_command, str_value);

    if (!strcmp(str_command,"xInjCenter")){
      fInjCenterX = atof(str_value) + detconst->GetOffsetX();
      fout << "xInjCenter: " << fInjCenterX << std::endl;
    }else if (!strcmp(str_command,"yInjCenter")){
      fInjCenterY = atof(str_value) + detconst->GetOffsetY();
      fout << "yInjCenter: " << fInjCenterY << std::endl;
    }else if (!strcmp(str_command,"zInjCenter")){
      fInjCenterZ = atof(str_value) + detconst->GetOffsetZ();
      fout << "zInjCenter: " << fInjCenterZ << std::endl;
    }else if (!strcmp(str_command,"radiusInj")){
      fInjRadiusXY = atof(str_value);
      fInjRadiusZ  = atof(str_value);
      fout << "radiusInj: " << fInjRadiusXY << std::endl;
    }else if (!strcmp(str_command,"radiusInjXY")){
      fInjRadiusXY = atof(str_value);
      fout << "radiusInjXY: " << fInjRadiusXY << std::endl;
    }else if (!strcmp(str_command,"radiusInjZ")){
      fInjRadiusZ = atof(str_value);
      fout << "radiusInjZ: " << fInjRadiusZ << std::endl;      
    }else if (!strcmp(str_command, "ISO")){
      // Isotropic Injection
      fout << "Isotopic Injection: " << buf << std::endl;
      char com[256];
      double energy, totalflux;
      int parcode, parcharge;
      sscanf(buf.c_str(),
	     "%s %lf %d %d %lf\n",
	     com,
	     &energy,
	     &parcode,
	     &parcharge,
	     &totalflux);
      fout << "(XYRadius, ZRadius) = ( " << fInjRadiusXY << " , " << fInjRadiusZ << ")" << std::endl
	   << "Energy = " << energy << " GeV" << std::endl
	   << "Particle Code = " << parcode << std::endl 
	   << "Particle Charge = " << parcharge << std::endl;
      Generator *gen = new Generator(fInjRadiusXY, fInjRadiusZ, energy, parcode, parcharge);
      gen->setTotalFlux(totalflux);
      fout << "Total Flux was set to " << totalflux << "(unit)" << std::endl;
      generators.push_back(gen);

    }else if (!strcmp(str_command, "COS")){
      // Cosine dependence Injection
      char com[256];
      double emin, emax, gammaindex, nindex, totalflux;
      int parcode, parcharge;
      sscanf(buf.c_str(),
	     "%s %lf %lf %d %d %lf %lf %lf", 
	     com,
	     &emin,
	     &emax,
	     &parcode,
	     &parcharge,
	     &gammaindex,
	     &nindex,
	     &totalflux);
      fout << "Cosine Injection: " << std::endl
	   << "(XYRadius, ZRadius) = ( " << fInjRadiusXY << " , " << fInjRadiusZ << ")" << std::endl
	   << "Energy Min = " << emin << " GeV" << std::endl
	   << "Energy Max = " << emax << " GeV" << std::endl
	   << "Particle Code = " << parcode << std::endl 
	   << "Particle Charge = " << parcharge << std::endl
	   << "Gamma Index = " << gammaindex << std::endl
	   << "n Index = " << nindex << std::endl;
      Generator *gen = new Generator(fInjRadiusXY, fInjRadiusZ, emin, emax, parcode, parcharge, gammaindex, nindex);
      gen->setTotalFlux(totalflux);
      fout << "Total Flux was set to " << totalflux << "(unit)" << std::endl;
      generators.push_back(gen);

    }else if (!strcmp(str_command, "FILE")){
      // Injection using .flx files
      char com[256];
      double emin, emax, coszmax, dcosz;
      int parcode, parcharge, numcz, nummom;
      char fname_before_suffix[256];	
      sscanf(buf.c_str(),
	     "%s %lf %lf %d %d %s %d %d %lf %lf",
	     com,
	     &emin,
	     &emax,
	     &parcode,
	     &parcharge,
	     fname_before_suffix,
	     &numcz,
	     &nummom,
	     &coszmax,
	     &dcosz);
      fout << "File Injection" << std::endl
	   << "(XYRadius, ZRadius) = ( " << fInjRadiusXY << " , " << fInjRadiusZ << ")" << std::endl
	   << "Energy Min = " << emin << " GeV" << std::endl
	   << "Energy Max * " << emax << " GeV" << std::endl
	   << "Particle Code = " << parcode << std::endl
	   << "Particle Charge = " << parcharge << std::endl
	   << "Filename Before Suffix = " << fname_before_suffix << std::endl
	   << "# of Files = " << numcz << std::endl
	   << "# of Momenta = " << nummom << std::endl
	   << "Max(CosTheta) = " << coszmax << std::endl
	   << "d(CosTheta) = " << dcosz << std::endl;
      Generator *gen = new Generator(fInjRadiusXY, fInjRadiusZ, emin, emax, parcode, parcharge, fname_before_suffix, numcz, nummom, coszmax, dcosz);
      generators.push_back(gen);
    }
  }

  // Prepare total intensities
  totalIntensity = 0;
  for (int iGen = 0 ; iGen < generators.size(); iGen++){
    double intens = generators[iGen]->getTotalIntensity();
    Intensities.push_back(intens);
    totalIntensity += intens;
    fout << iGen << "-th generator: " << intens << " (particles/sec)." << std::endl;
  }
  fout << "total: " << totalIntensity << " (particles/sec)." << std::endl;

  fParticleGun = new G4ParticleGun(1);

  fout.close();

}

void MsmPrimaryGeneratorAction::PrintGeneratorInformation(){
  ofstream fout(filenameLog.c_str(), ios::app);
  fout << "Exposure (sec): " <<  generators[0]->getElapsedTime() << endl;
  fout.close();
}



MsmPrimaryGeneratorAction::~MsmPrimaryGeneratorAction()
{
  delete fParticleGun;
  PrintGeneratorInformation();
  for (int iGen = 0 ; iGen < generators.size(); iGen++){
    delete generators[iGen];
  }
}






///////////////////////////////
// Description on fInjMethod //
///////////////////////////////

void MsmPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double worldXHalfLength = 0;
  G4double worldYHalfLength = 0;
  G4LogicalVolume* worlLV
    = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox = 0;
  if ( worlLV ) worldBox = dynamic_cast< G4Box*>(worlLV->GetSolid()); 
  if ( worldBox ) {
    worldXHalfLength = worldBox->GetXHalfLength();  
    worldYHalfLength = worldBox->GetYHalfLength();
  }
  else  {
    G4cerr << "World volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }

  // Choose Generators;
  Generator *genChosen;
  if (generators.size() == 1) genChosen = generators[0];
  else{
    double rndGen = G4RandFlat::shoot(0.0, totalIntensity);
    double sumintens = 0;
    for(int iGen = 0 ; iGen < generators.size(); iGen++){
      sumintens += Intensities[iGen] ;
      if (rndGen <= sumintens){
	genChosen = generators[iGen];
	break;
      }
    }
  }

  // Shoot
  double rnd1 = G4RandFlat::shoot(0.0,1.0);
  double rnd2 = G4RandFlat::shoot(0.0,1.0);
  genChosen->shoot(rnd1, rnd2);

  // Get Direction
  double theta = genChosen->getRandomizedTheta();
  double phi   = G4RandFlat::shoot(-M_PI, +M_PI);
  double cosTheta = cos(theta);
  double sinTheta = sin(theta); 
  double cosPhi   = cos(phi);           
  double sinPhi   = sin(phi);
  // Get Energy
  double injection_energy = genChosen->getRandomizedEnergy();

  // calc injection direction (vx,vy,vz)
  double vx = -sinTheta * cosPhi;
  double vy = -sinTheta * sinPhi;
  double vz = -cosTheta;

  // calc injection point (xinj,yinj,zinj)
  double radiusSkyXY = fInjRadiusXY;
  double radiusSkyZ  = fInjRadiusZ;
  double x0, y0, z0, x, y, z, A, B, C, det;
  do{
    do{
      x0 = G4RandFlat::shoot(-radiusSkyXY, +radiusSkyXY);
      y0 = G4RandFlat::shoot(-radiusSkyXY, +radiusSkyXY);
      z0 = 0;
    } while( x0 * x0 + y0 * y0 > radiusSkyXY * radiusSkyXY ); 
    x = cosPhi * cosTheta * x0 - sinPhi * y0;
    y = sinPhi * cosTheta * x0 + cosPhi * y0;
    z = -sinTheta * x0;
    A = pow(vx / radiusSkyXY, 2) + pow(vy / radiusSkyXY, 2) + pow(vz / radiusSkyZ, 2);
    B = 2.0 * (x*vx/radiusSkyXY/radiusSkyXY + y*vy/radiusSkyXY/radiusSkyXY + z*vz/radiusSkyZ/radiusSkyZ);
    C = pow(x / radiusSkyXY, 2) + pow(y / radiusSkyXY, 2) + pow(z/ radiusSkyZ, 2) - 1.0;  
    det = B * B - A * C * 4.0; 
  } while(det < 0);
  double t = (-B-sqrt(B*B-4*A*C)) / (A*2);
  double xinj = x + vx * t + fInjCenterX;
  double yinj = y + vy * t + fInjCenterY;
  double zinj = z + vz * t + fInjCenterZ;
 

  // If the injection point is outside the terrain, Shoot in Geant4 !!!
  G4RunManager* runmanager = G4RunManager::GetRunManager();
  MsmDetectorConstruction* detconst = (MsmDetectorConstruction*)runmanager->GetUserDetectorConstruction();
  G4bool flgIsInTerrain = detconst->IsInTerrain(xinj,yinj,zinj);
  if (!flgIsInTerrain){
   
    // set position
    fParticleGun->SetParticlePosition(G4ThreeVector((xinj+fInjCenterX) * m, (yinj+fInjCenterY) * m, (zinj+fInjCenterZ) * m));
    // set direction
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));
    // set energy
    fParticleGun->SetParticleEnergy(injection_energy * GeV);

    // set particle & charge
    int parcode = genChosen->getParCode();
    int parcharge = genChosen->getParCharge();
    G4ParticleDefinition* particleDefinition;
    if (parcode == 11){
      particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
      parcharge = -1;
    }else if (parcode == -11){
      particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e+");
      parcharge = +1;
    }else if (parcode == 13){
      particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
      parcharge = -1;
    }else if (parcode == -13){
      particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("mu+");
      parcharge = +1;
    }else if (parcode == 15){
      particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("tau-");
      parcharge = -1;
    }else if (parcode == -15){
      particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("tau+");
      parcharge = +1;
    }else if (parcode == 22){ // photon
      particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
      parcharge = 0;
    }else if (parcode == 2212){
      particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("proton");
      parcharge = +1;
    }else if (parcode == 2112){
      particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
      parcharge = 0;
    }else{      
      // geantino
      particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
      parcharge = 0;
    }

    // Set Primary Info to Global Thread Local variables
    // This is for Sensitive Detector
    primaryParCode = parcode;
    primaryParCharge = parcharge;
    primaryKE = injection_energy * GeV;
    primaryX = (xinj + fInjCenterX) * m;
    primaryY = (yinj + fInjCenterY) * m;
    primaryZ = (zinj + fInjCenterZ) * m;
    primaryVX = vx;
    primaryVY = vy;
    primaryVZ = vz;

    // Start Tracking in Geant4
    fParticleGun->SetParticleDefinition(particleDefinition);
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  numInjection++;
  if (numInjection % modPrint == 0) {
    PrintGeneratorInformation();
  } 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

