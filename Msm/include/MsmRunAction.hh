#ifndef MsmRunAction_h
#define MsmRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4GenericMessenger.hh"
#include "globals.hh"


class G4Run;

class MsmRunAction : public G4UserRunAction
{
private:

  G4String fOutputFileName;

  G4GenericMessenger *fMessenger;

  public:
    MsmRunAction();
    virtual ~MsmRunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

};



#endif

