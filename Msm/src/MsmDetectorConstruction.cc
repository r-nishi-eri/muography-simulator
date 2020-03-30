#include "MsmDetectorConstruction.hh"
#include "MsmCalorimeterSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4SDManager.hh"
#include "G4Region.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4GenericMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <stdio.h>
#include <stdlib.h>




MsmDetectorConstruction::MsmDetectorConstruction(std::string str_demfile, std::string str_detectorfile)
 : G4VUserDetectorConstruction(),
   fMessenger(0),
   fMagField(0),
   fIsSphere(true),
   fCheckOverlaps(true)
{
  // Define /Msm/det commands using generic messenger class
  fMessenger 
    = new G4GenericMessenger(this, "/Msm/det/", "Detector construction control");

  // Define /Msm/det/setMagField command
  G4GenericMessenger::Command& setMagFieldCmd
    = fMessenger->DeclareMethodWithUnit("setMagField",
					"tesla",
					&MsmDetectorConstruction::SetMagField,
					"Define magnetic field value (in X direction)");

  // Read DEMFILE
  ReadDemFile(str_demfile);

  // Read PARAMFILE
  ReadParamFile(str_detectorfile);
}




MsmDetectorConstruction::~MsmDetectorConstruction()
{ 
  delete fMagField;
  delete fMessenger;
}





// read DEM file
void MsmDetectorConstruction::ReadDemFile(std::string str_demfile){

  std::cout << "start reading file: " << str_demfile << std::endl;
  std::ifstream fin(str_demfile.c_str());
  if (!fin) {
    std::cerr << "filename " << str_demfile.c_str() << " could not be opened." << std::endl;
    exit(1);
  }
  fin >> fMeshX >> fMeshY;
  fin >> fNumX >> fNumY;
  fDemZ = new double* [fNumX];
  for (int ix = 0; ix < fNumX; ix++) fDemZ[ix] = new double [fNumY];
  double demZMax = -9999;
  double demZMin = 0;
  for (int iy = 0; iy < fNumY; iy++){
    for (int ix = 0; ix < fNumX; ix++){
      fin >> fDemZ[ix][iy];
      demZMax = (fDemZ[ix][iy] > demZMax ? fDemZ[ix][iy] : demZMax);
    }
  } 
  fDemZMax = demZMax;
  fDemZMin = demZMin;

  fDemSizeX = fMeshX * double(fNumX);
  fDemSizeY = fMeshY * double(fNumY);
  fDemSizeZ = fDemZMax-fDemZMin;
  fDemZCenter = (fDemZMax + fDemZMin) / 2.0;
  fOffsetX = - fDemSizeX / 2.0;
  fOffsetY = - fDemSizeY / 2.0;
  fOffsetZ = - fDemSizeZ / 2.0;
  std::cout << fMeshX << " " << fMeshY << std::endl;

}

G4bool MsmDetectorConstruction::IsInTerrain(G4double xWorld, G4double yWorld, G4double zWorld){
  // The unit is meter.
  G4double xMap = xWorld - fOffsetX;
  G4double yMap = yWorld - fOffsetY;
  G4double zMap = zWorld - fOffsetZ;
  if (xMap < 0 || 
      xMap > fMeshX * fNumX || 
      yMap < 0 || 
      yMap > fMeshY * fNumY){
    // Outside DEM defitinion
    if (zMap <= fDemZMin) return true;
    else return false;
  }
  G4int iDemX = xMap / G4double(fMeshX);
  G4int iDemY = yMap / G4double(fMeshY);
  G4double demz = fDemZ[iDemX][iDemY];
  if (zMap <= demz) return true;
  return false;
}

// read paramater file
void MsmDetectorConstruction::ReadParamFile(std::string str_detectorfile){

  std::cout << "start reading file: " << str_detectorfile << std::endl;
  std::ifstream fin(str_detectorfile.c_str());
  if (!fin){
    std::cerr << "file " << str_detectorfile.c_str() << " could not be opened." << std::endl;
    exit(1);
  }

  std::string buf;
  while(fin && std::getline(fin, buf)){
    
    char str_command[256];
    char str_value[256];
    sscanf(buf.c_str(), "%s %s", str_command, str_value);
    std::cout << buf << std::endl;
    std::cout << str_command << " " << str_value << std::endl;

    if (!strcmp(str_command,"xDet")){
      fXDet = atof(str_value);
      std::cout << "xDet: " << fXDet << std::endl;
    }else if (!strcmp(str_command,"yDet")){
      fYDet = atof(str_value);
      std::cout << "yDet: " << fYDet << std::endl;
    }else if (!strcmp(str_command,"zDet")){
      fZDet = atof(str_value);
      std::cout << "zDet: " << fZDet << std::endl;
    }else if (!strcmp(str_command,"radiusDet")){
      fRadiusDet = atof(str_value);
      std::cout << "radiusDet: " << fRadiusDet << std::endl;
    }else if (!strcmp(str_command,"idOutput")){
      fIdOutput = atoi(str_value);
      std::cout << "idOutput: " << fIdOutput << std::endl;
    }else if (!strcmp(str_command, "outDir")){
      fOutDir = str_value;
    }else if (!strcmp(str_command, "tub")){
      fIsSphere = false; // Tub Mode
      fTubHeight = atof(str_value); // Tub Height in meter (not in half length)
    }
  }

  // make copy of paramfile
  char outFilename [256];
  sprintf(outFilename, "detectorfile%d.save", fIdOutput);
  std::ifstream fin_copy(str_detectorfile.c_str());
  std::ofstream fout_copy(outFilename);
  fout_copy << fin_copy.rdbuf() << std::flush;

}


G4VPhysicalVolume* MsmDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  // Define volumes
  return DefineVolumes();
}




void MsmDetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;
  nistManager->FindOrBuildMaterial("G4_Fe", fromIsotopes);
  nistManager->FindOrBuildMaterial("G4_Pb", fromIsotopes);
  nistManager->FindOrBuildMaterial("G4_PHOTO_EMULSION", fromIsotopes);
  nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE", fromIsotopes);
  nistManager->FindOrBuildMaterial("G4_AIR",fromIsotopes);

  // define standard rock with density = 2.0 
  G4double density = 2.00*g/cm3;
  G4double a = 22.0*g/mole;
  G4Material* matSR = new G4Material("StandardRock", 11.0, a, density);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}




G4VPhysicalVolume* MsmDetectorConstruction::DefineVolumes()
{
  //
  // Get Materials
  //
  G4Material* airMaterial = G4Material::GetMaterial("G4_AIR");
  G4Material* terrainMaterial = G4Material::GetMaterial("StandardRock");
  G4Material* detectorMaterial = G4Material::GetMaterial("G4_AIR");

  if ( ! terrainMaterial || ! airMaterial || ! detectorMaterial){
    G4cerr << "Cannot retrieve materials already defined. " << G4endl;
    G4cerr << "Exiting application " << G4endl;
    exit(1);
  }  

  //
  // Visual Attributes
  // 
  G4VisAttributes* terrainVisAtt = new G4VisAttributes(G4Colour(0,1,1)); // cyan
  terrainVisAtt->SetForceSolid(true);
  G4VisAttributes* detectorVisAtt = new G4VisAttributes(G4Colour(1,0,0)); // red
  detectorVisAtt->SetForceSolid(true);

  //     
  // World
  //
  G4double worldSizeX = fDemSizeX * 3.0 * m;
  G4double worldSizeY = fDemSizeY * 3.0 * m;
  G4double worldSizeZ = fDemSizeZ * 3.0 * m;

  G4VSolid* worldS 
    = new G4Box("World",           // its name
                worldSizeX,
		worldSizeY,
		worldSizeZ); // its size
                         
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 airMaterial,  // its material
                 "World");         // its name
                                   
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  worldLV->SetVisAttributes (G4VisAttributes::Invisible);

  /////////////
  // Terrain //
  /////////////
  int cntDemVoxel = 0;
  for (int iDemX = 0; iDemX < fNumX; iDemX ++){
    for (int iDemY = 0; iDemY < fNumY; iDemY ++){
      char str_terrain_name[25];
      cntDemVoxel++;
      sprintf(str_terrain_name, "terrain%d", cntDemVoxel);

      G4VSolid* terrainS = new G4Box(str_terrain_name,
				     fMeshX/2*m , fMeshY/2*m, (fDemZ[iDemX][iDemY]-fDemZMin)/2*m );
      G4LogicalVolume* terrainLV = 
	new G4LogicalVolume(terrainS,
			    terrainMaterial,
			    "Terrain");
      terrainLV->SetVisAttributes(terrainVisAtt);
      
      new G4PVPlacement(0,
			G4ThreeVector((fMeshX*(double(iDemX)+0.5)+fOffsetX)*m,
				      (fMeshY*(double(iDemY)+0.5)+fOffsetY)*m,
				      ((fDemZ[iDemX][iDemY]-fDemZMin)/2+fOffsetZ)*m),
			terrainLV,
			str_terrain_name,
			worldLV,
			false,
			0,
			false);
    }
  }

  //////////////
  // Detector //
  //////////////

  G4VSolid* detectorS;
  if (fIsSphere){
    // spherical detector mode
    detectorS = new G4Orb("detector", fRadiusDet * m);
  }else{
    // tube detector mode
    detectorS = new G4Tubs("detector", 
			   (fRadiusDet       ) * m, 
			   (fRadiusDet + 0.01) * m,  // .01 m thickness
			   fTubHeight/2.0     * m, 
			   0,
			   M_PI * 2.0);
  }
  G4LogicalVolume* detectorLV = 
    new G4LogicalVolume(detectorS,
			detectorMaterial,
			"detectorLV");
  detectorLV->SetVisAttributes(detectorVisAtt);
  new G4PVPlacement(0,
		    G4ThreeVector((fXDet+fOffsetX)*m,
				  (fYDet+fOffsetY)*m,
				  (fZDet+fOffsetZ)*m),
		    detectorLV,
		    "detector",
		    worldLV,
		    false,
		    0,
		    fCheckOverlaps);

  return worldPV;
}



void MsmDetectorConstruction::ConstructSDandField()
{
  std::stringstream sss;
  sss << int(fIdOutput);
  MsmCalorimeterSD* detectorSD = new MsmCalorimeterSD("detectorSD" + sss.str(), fOutDir);
  SetSensitiveDetector("detectorLV", detectorSD);
}

void MsmDetectorConstruction::SetMagField(G4double fieldValue)
{
  // Apply a global uniform magnetic field along X axis
  G4FieldManager* fieldManager
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  // Delete the existing magnetic field
  if ( fMagField )  delete fMagField; 

  if ( fieldValue != 0. ) {
    // create a new one if not null
    fMagField 
      = new G4UniformMagField(G4ThreeVector(fieldValue, 0., 0.));
      
    fieldManager->SetDetectorField(fMagField);
    fieldManager->CreateChordFinder(fMagField);
  } 
  else {
    fMagField = 0;
    fieldManager->SetDetectorField(fMagField);
  }
}
