#ifndef SNneutrinosRun_h
#define SNneutrinosRun_h 1

#include "G4Run.hh"

class G4ParticleDefinition;


class SNneutrinosRun : public G4Run
{
 public:
  SNneutrinosRun();
  ~SNneutrinosRun();


  void AddWaterDetection(void) {fWaterDetection += 1;}
  
  void AddAbsorption(void) {fAbsorptionCounter += 1;}
  void AddBoundary(void) {fBoundaryCounter += 1;}
     

  void SetPrimary( G4ParticleDefinition* particle, G4double energy);


  void AddCerenkovRun(G4double n)
  {
    fCerenkovCounter += n;
    fCerenkov2 += n * n;
  };

  void AddPMTDetectionRun(G4double n)
  {
    fPMTDetectionCounter += n;
    fPMT2 += n * n;
  };

  
  void Merge(const G4Run*) override;
  void EndOfRun();

 private:
  G4ParticleDefinition* fParticle;

  G4double fCerenkovAll;
  G4double fCerenkovCounter;
  G4double fCerenkov2;

  G4double fPMTAll;
  G4double fPMTDetectionCounter;
  G4double fPMT2;

  G4double fAbsorptionCounter;
 
  G4double fBoundaryCounter;

  G4double fEnergy;
  G4int fWaterDetection;

  

};

#endif
