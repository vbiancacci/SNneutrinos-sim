#include "SNneutrinosPrimaryGeneratorMessenger.hh"
#include "SNneutrinosPrimaryGeneratorAction.hh"

SNneutrinosPrimaryGeneratorMessenger::SNneutrinosPrimaryGeneratorMessenger(SNneutrinosPrimaryGeneratorAction* action)
 : fAction(action)
{
    fEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gun/energy", this);
    fEnergyCmd->SetGuidance("Set primary particle energy");
    fEnergyCmd->SetParameterName("Energy", true);
    fEnergyCmd->SetUnitCategory("Energy");
}

SNneutrinosPrimaryGeneratorMessenger::~SNneutrinosPrimaryGeneratorMessenger() {
    delete fEnergyCmd;
}

void SNneutrinosPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
    if (cmd == fEnergyCmd) {
        G4double fEnergy = fEnergyCmd->GetNewDoubleValue(value);
        fAction->SetEnergy(fEnergy);
        G4cout << "Energy set to " << fEnergy/GeV << " GeV in RunConfig." << G4endl;
        G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 Messenger received energy: " << value << " MeV" << G4endl;
    }
}
