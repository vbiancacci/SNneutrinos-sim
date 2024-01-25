#include "SNneutrinosActionInitialization.hh"
#include "SNneutrinosEventAction.hh"
#include "SNneutrinosPrimaryGeneratorAction.hh"
#include "SNneutrinosRunAction.hh"
#include "SNneutrinosSteppingAction.hh"


SNneutrinosActionInitialization::SNneutrinosActionInitialization()
  : G4VUserActionInitialization()
{}

SNneutrinosActionInitialization::~SNneutrinosActionInitialization() {}


void SNneutrinosActionInitialization::BuildForMaster() const
{
  SetUserAction(new SNneutrinosRunAction());
}

void SNneutrinosActionInitialization::Build() const
{
  SNneutrinosPrimaryGeneratorAction* primary =  new SNneutrinosPrimaryGeneratorAction();
  SetUserAction(primary);
  SetUserAction(new SNneutrinosRunAction(primary));
  
  SNneutrinosEventAction* event = new SNneutrinosEventAction();
  SetUserAction(event);
  SetUserAction(new SNneutrinosSteppingAction(event));
}
