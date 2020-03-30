#include "MsmActionInitialization.hh"
#include "MsmPrimaryGeneratorAction.hh"
#include "MsmRunAction.hh"
#include "MsmEventAction.hh"
#include "MsmSteppingAction.hh"


MsmActionInitialization::MsmActionInitialization(std::string str_injectionfile, std::string str_dest_injectionlog)
: G4VUserActionInitialization()
{
  fNameInjection = str_injectionfile;
  fNameInjectionLog = str_dest_injectionlog;  
}

MsmActionInitialization::~MsmActionInitialization()
{}

void MsmActionInitialization::BuildForMaster() const
{
  SetUserAction(new MsmRunAction);
}

void MsmActionInitialization::Build() const
{
  SetUserAction(new MsmPrimaryGeneratorAction(fNameInjection,fNameInjectionLog));
  SetUserAction(new MsmRunAction);
  MsmEventAction* eventAction = new MsmEventAction;
  SetUserAction(eventAction);
}
