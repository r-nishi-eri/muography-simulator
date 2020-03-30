#include "MsmEventAction.hh"
#include "MsmCalorimeterSD.hh"
#include "MsmRunAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include <iostream>
#include <fstream>



MsmEventAction::MsmEventAction()
 : G4UserEventAction(),
   fPrintModulo(1)
{
}



MsmEventAction::~MsmEventAction()
{

}


void MsmEventAction::BeginOfEventAction(const G4Event* event)
{  
  G4int eventID = event->GetEventID();
  if ( eventID % fPrintModulo == 0 )  { 

  }
}



void MsmEventAction::EndOfEventAction(const G4Event* event)
{  
  G4int eventID = event->GetEventID();
  if ( eventID % fPrintModulo == 0) {

  }  
  
 
}  


