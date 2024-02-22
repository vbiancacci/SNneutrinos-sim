
#ifndef SNneutrinosRunAction_h
#define SNneutrinosRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

class SNneutrinosPrimaryGeneratorAction;
class SNneutrinosRun;

class G4Run;

class SNneutrinosRunAction : public G4UserRunAction
{
 public:
  SNneutrinosRunAction(SNneutrinosPrimaryGeneratorAction* = nullptr);
  ~SNneutrinosRunAction();

  G4Run* GenerateRun() override;
  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

 private:
  SNneutrinosRun* fRun;
  SNneutrinosPrimaryGeneratorAction* fPrimary;
  G4String FileName = "results/new/2MeV/SNneutrinos.root";
};
#endif
