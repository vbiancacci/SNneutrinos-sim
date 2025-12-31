
#ifndef SNneutrinosRunMessenger_h
#define SNneutrinosRunMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"

class SNneutrinosRunAction;

class SNneutrinosRunMessenger : public G4UImessenger {
public:
    SNneutrinosRunMessenger(SNneutrinosRunAction* runAction);
    ~SNneutrinosRunMessenger() override;

    void SetNewValue(G4UIcommand* command, G4String newValue) override;

private:
    SNneutrinosRunAction* fRunAction;
    G4UIcmdWithAString* fFileNameCmd;
};

#endif // SNneutrinosRunMessenger_h