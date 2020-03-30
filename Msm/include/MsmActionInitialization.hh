#ifndef MsmActionInitialization_h
#define MsmActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

#include <string>

class MsmActionInitialization : public G4VUserActionInitialization
{
public:
  MsmActionInitialization(std::string,std::string);
  virtual ~MsmActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;

  std::string fNameInjection;
  std::string fNameInjectionLog;
};

#endif
