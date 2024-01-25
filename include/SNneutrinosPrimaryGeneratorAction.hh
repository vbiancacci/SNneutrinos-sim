
#ifndef SNneutrinosPrimaryGeneratorAction_h
#define SNneutrinosPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "SNneutrinosDetectorConstruction.hh"
#include <random>
#include "TRandom3.h"

class G4Event;
class SNneutrinosPrimaryGeneratorMessenger;


class SNneutrinosPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  SNneutrinosPrimaryGeneratorAction();
  ~SNneutrinosPrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*) override;

  G4ParticleGun* GetParticleGun() { return fParticleGun; }

 private:
  G4ParticleGun* fParticleGun = nullptr;
  SNneutrinosDetectorConstruction* detectordb;
  SNneutrinosPrimaryGeneratorMessenger* fGunMessenger;
  std::random_device rd;
  std::ranlux24      generator;
  TRandom3* xGenerator = new TRandom3(0);
  TRandom3* yGenerator = new TRandom3(0);
  TRandom3* zGenerator = new TRandom3(0);
  TRandom3* thetaGenerator = new TRandom3(0);
  TRandom3* phiGenerator = new TRandom3(0);

};



#endif
