#ifndef SNNEUTRINOSPRIMARYGENERATORMESSENGER_H
#define SNNEUTRINOSPRIMARYGENERATORMESSENGER_H

#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

class SNneutrinosPrimaryGeneratorAction;

class SNneutrinosPrimaryGeneratorMessenger : public G4UImessenger {
public:
    SNneutrinosPrimaryGeneratorMessenger(SNneutrinosPrimaryGeneratorAction*);
    ~SNneutrinosPrimaryGeneratorMessenger();

    void SetNewValue(G4UIcommand*, G4String) override;

private:
    SNneutrinosPrimaryGeneratorAction* fAction;
    G4UIcmdWithADoubleAndUnit* fEnergyCmd;
};

#endif // SNNEUTRINOSPRIMARYGENERATORMESSENGER_H