#ifndef MsmEventAction_h
#define MsmEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4GenericMessenger;

class MsmEventAction : public G4UserEventAction
{
public:
  MsmEventAction();
  virtual ~MsmEventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);
                     
  // set methods
  void SetPrintModulo(G4int value);
    
private:
  // methods
  G4int  fPrintModulo;
};

// inline functions

inline void MsmEventAction::SetPrintModulo(G4int value) {
  fPrintModulo = value;
}
                     


#endif

    
