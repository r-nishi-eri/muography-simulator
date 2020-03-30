#ifndef MsmDetectorConstruction_h
#define MsmDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <string>

class G4Box;
class G4VPhysicalVolume;
class G4UniformMagField;
class G4GenericMessenger;


class MsmDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    MsmDetectorConstruction(std::string, std::string);
    virtual ~MsmDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // set methods
    //
    void SetMagField(G4double fieldValue);
     
  G4double GetOffsetX(){return fOffsetX;};
  G4double GetOffsetY(){return fOffsetY;};
  G4double GetOffsetZ(){return fOffsetZ;};

  G4bool IsInTerrain(G4double, G4double, G4double);
  
  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    void ReadDemFile(std::string);
    void ReadParamFile(std::string);
  
    // data members
    //
    G4GenericMessenger*  fMessenger; // messenger 
    G4UniformMagField*   fMagField;  // magnetic field

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
  
  // DEM info
  G4double fMeshX, fMeshY;
  G4int fNumX, fNumY;
  G4double **fDemZ;
  G4double fDemZMax, fDemZMin;
  G4double fDemSizeX, fDemSizeY, fDemSizeZ;
  G4double fDemZCenter;
  G4double fOffsetX, fOffsetY, fOffsetZ;

  // Parameter detectorfile info
  G4double fXDet, fYDet, fZDet;
  G4double fRadiusDet, fLngDet, fLatDet;
  G4int fIdOutput;
  std::string fOutDir;
  G4int fIsSphere;
  G4double fTubHeight;

};



#endif

