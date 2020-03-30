#ifndef MsmCalorimeterSD_h
#define MsmCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "MsmPrimaryGeneratorAction.hh"
#include "G4Threading.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;


class MsmCalorimeterSD : public G4VSensitiveDetector
{
  public:
  MsmCalorimeterSD(const G4String& name, const G4String& outdir);
 
  virtual ~MsmCalorimeterSD();
  
  // methods from base class
  virtual void   Initialize(G4HCofThisEvent* hitCollection);
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
  virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

private:
  G4String fnameOut;
  G4String fOutDir;
};



#endif

