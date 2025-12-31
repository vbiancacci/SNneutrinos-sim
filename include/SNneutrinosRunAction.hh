
#ifndef SNneutrinosRunAction_h
#define SNneutrinosRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"
#include "SNneutrinosRunMessenger.hh"

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
  void SetOutputFile(const G4String& filename);
  G4String GetOutputFile() const { return fOutputFileName; }

 private:
  SNneutrinosRun* fRun;
  SNneutrinosPrimaryGeneratorAction* fPrimary;
  G4String fOutputFileName;
  SNneutrinosRunMessenger* fMessenger;
  //G4String FileName = "results/L1000/e+/alternative/1MeV/SNneutrinos.root";
};
#endif
