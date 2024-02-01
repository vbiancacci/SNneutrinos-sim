#ifndef SNneutrinosSteppingAction_h
#define SNneutrinosSteppingAction_h 1

#include "SNneutrinosEventAction.hh"
#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include "TGraph.h"
#include "G4SystemOfUnits.hh"

class SNneutrinosSteppingAction : public G4UserSteppingAction
{
 public:
  SNneutrinosSteppingAction(SNneutrinosEventAction*);
  ~SNneutrinosSteppingAction();

  void UserSteppingAction(const G4Step*) override;
  
  G4int GetNumberOfBounces();

  G4double PMT_QE(std::string volume, G4double energy){

    G4String nameFile= "PMT_efficiency.dat";

    G4String pathString("../materials/");
    G4String pathFile = pathString + nameFile;

    auto Graph = new TGraph(pathFile.data());

    static const G4double LambdaE = 1239.84193 * nanometer ;
    
    G4double qe;
    auto em = Graph->Eval(LambdaE/energy/nanometer);
    qe=(em >= 0 ? em : 0);
    
    return qe;
  }

 private:
  SNneutrinosEventAction* fEventAction;
  G4int fCounterBounce;
  inline void ResetBounceCounter() {
    fCounterBounce      = 0;
  }

};



#endif
