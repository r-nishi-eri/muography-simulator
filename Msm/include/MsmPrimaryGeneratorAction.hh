#ifndef MsmPrimaryGeneratorAction_h
#define MsmPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include "Generator.hh"


#include <vector>

class G4ParticleGun;
class G4Event;


extern G4ThreadLocal G4int primaryParCode;
extern G4ThreadLocal G4int primaryParCharge;
extern G4ThreadLocal G4double primaryKE;
extern G4ThreadLocal G4double primaryX;
extern G4ThreadLocal G4double primaryY;
extern G4ThreadLocal G4double primaryZ;
extern G4ThreadLocal G4double primaryVX;
extern G4ThreadLocal G4double primaryVY;
extern G4ThreadLocal G4double primaryVZ;






class MsmPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  MsmPrimaryGeneratorAction(string str_injectionfile, string str_injectionlog);    
  virtual ~MsmPrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event* event);
  
  // set methods
  void SetRandomFlag(G4bool value);

  void PrintGeneratorInformation();

private:
  G4ParticleGun*  fParticleGun; // G4 particle gun

  // configuration of injection
  G4double fInjCenterX, fInjCenterY, fInjCenterZ;
  G4int fInjMethod;
  G4double fInjRadiusXY, fInjRadiusZ;
  G4int fParCode, fParCharge;
  G4double fParEnergy;
  G4int fInjNumber;

  string filenameLog;
  
  vector <Generator*> generators;
  vector <double>     Intensities;
  double totalIntensity;

  int numInjection;
  int modPrint;

};


#endif


