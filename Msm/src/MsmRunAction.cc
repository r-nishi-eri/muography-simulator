#include "MsmRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"



MsmRunAction::MsmRunAction()
 : G4UserRunAction()
{ 
  // Define /Msm/run commands using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/Msm/run/", "Run control");
  
  // Define /Msm/run/setOutputFileName command
  G4GenericMessenger::Command& setOutputFileName 
    = fMessenger->DeclareProperty("setOuptutFileName",
   				  fOutputFileName,
   				  "output file name");  
  
  // default filename
  fOutputFileName = "hits.root";
  
}



MsmRunAction::~MsmRunAction()
{

}



void MsmRunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
}



void MsmRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nofEvents = aRun->GetNumberOfEvent();
  if ( nofEvents == 0 ) return;

}

