#ifndef SNneutrinosEventAction_h
#define SNneutrinosEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"


class SNneutrinosEventAction : public G4UserEventAction
{
 public:
  SNneutrinosEventAction();
  ~SNneutrinosEventAction();

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  void AddPMTDetectionEvent(void) {fPMT += 1;}
  void AddCerenkovEvent(G4int n ) { fCerenkov += n; }
  
 private:
  G4int event_ID;
  G4int fPMT;
  G4int fCerenkov;
  G4double efficiency;
};

#endif
