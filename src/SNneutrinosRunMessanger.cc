#include "SNneutrinosRunMessenger.hh"
#include "SNneutrinosRunAction.hh"


SNneutrinosRunMessenger::SNneutrinosRunMessenger(SNneutrinosRunAction* runAction)
 : fRunAction(runAction)
{
    fFileNameCmd = new G4UIcmdWithAString("/myRun/outputFile", this);
    fFileNameCmd->SetGuidance("Set full path for output file");
    fFileNameCmd->SetParameterName("filename", false);
}
SNneutrinosRunMessenger::~SNneutrinosRunMessenger()
{
    delete fFileNameCmd;
}

void SNneutrinosRunMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if(command == fFileNameCmd) {
        fRunAction->SetOutputFile(newValue);
    }
}