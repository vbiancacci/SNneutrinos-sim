
#ifndef SNneutrinosActionInitialization_h
#define SNneutrinosActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class SNneutrinosActionInitialization : public G4VUserActionInitialization
{
 public:
  SNneutrinosActionInitialization();
  ~SNneutrinosActionInitialization();

  void BuildForMaster() const override;
  void Build() const override;
};

#endif