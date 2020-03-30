#include "MsmCalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include <iostream>

#include "G4MTRunManager.hh"

using namespace CLHEP;




MsmCalorimeterSD::MsmCalorimeterSD(const G4String& name, const G4String& outdir)
 : G4VSensitiveDetector(name)
{
  G4int threadid = G4Threading::G4GetThreadId();
  G4String strthreadid = G4UIcommand::ConvertToString(threadid);;
  if (outdir == "") outdir == ".";
  fnameOut = outdir + "/" + name + "-" + strthreadid  + ".dat";
  std::ofstream fout (fnameOut.c_str(), std::ios::app);
  fout.close();}



MsmCalorimeterSD::~MsmCalorimeterSD() 
{ 
}



void MsmCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  
}



G4bool MsmCalorimeterSD::ProcessHits(G4Step* step, 
                                     G4TouchableHistory*)
{  


  // energy deposit
  G4double edep = step->GetTotalEnergyDeposit();
  // step length
  G4double stepLength = step->GetStepLength();

  // Detect particles on the boundary
  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  if ( preStepPoint->GetStepStatus() != fGeomBoundary ) return false;

  // Get track data
  G4Track* track = step->GetTrack();

  // Get Event ID
  G4int eventid = G4MTRunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  // Get Track ID
  G4int trackid = track->GetTrackID();
    
  std::ofstream fout(fnameOut.c_str(), std::ios::app);
  fout << track->GetDefinition()->GetPDGCharge() << " "
       << track->GetDefinition()->GetPDGEncoding() << " " 
       << track->GetPosition().x() / m  << " " 
       << track->GetPosition().y() / m << " " 
       << track->GetPosition().z() / m << " " 
       << track->GetMomentumDirection().x() << " " 
       << track->GetMomentumDirection().y() << " " 
       << track->GetMomentumDirection().z() << " " 
       << track->GetKineticEnergy() / GeV << " " 
       << primaryParCharge << " " 
       << primaryParCode << " " 
       << primaryX / m << " " 
       << primaryY / m << " " 
       << primaryZ / m << " " 
       << primaryVX << " " 
       << primaryVY << " " 
       << primaryVZ << " " 
       << primaryKE / GeV << " " 
       << eventid << " " 
       << trackid << std::endl;
  fout.close();

  // Stop calculation inside the sensitive detector.
  track->SetTrackStatus(fKillTrackAndSecondaries);

  return true;
}



void MsmCalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 

  }
}


