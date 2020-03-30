#include "MsmSteppingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"


MsmSteppingAction::MsmSteppingAction()
  : G4UserSteppingAction()
{
}



MsmSteppingAction::~MsmSteppingAction()
{ 
}



void MsmSteppingAction::UserSteppingAction(const G4Step* step)
{
}


