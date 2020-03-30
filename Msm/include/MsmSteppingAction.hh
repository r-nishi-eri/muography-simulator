#ifndef MsmSteppingAction_h
#define MsmSteppingAction_h 1

#include "G4UserSteppingAction.hh"



class MsmSteppingAction : public G4UserSteppingAction
{
public:
  MsmSteppingAction();
  virtual ~MsmSteppingAction();

  virtual void UserSteppingAction(const G4Step* step);

};



#endif
