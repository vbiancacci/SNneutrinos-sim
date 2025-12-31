
#ifndef SNneutrinosPrimaryGeneratorAction_h
#define SNneutrinosPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "SNneutrinosDetectorConstruction.hh"
#include "SNneutrinosPrimaryGeneratorMessenger.hh"
#include <random>
#include <iostream>
#include "TRandom3.h"
class G4Event;
class SNneutrinosPrimaryGeneratorMessenger;


class SNneutrinosPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  SNneutrinosPrimaryGeneratorAction();
  ~SNneutrinosPrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*) override;
  std::vector<double> ReadEnergiesFromGrabmayer(const std::string& );
  std::vector<double> ReadEnergiesFromANNRIGd(const std::string& );
  
  void SetEnergy(G4double energy);
  G4ParticleGun* GetParticleGun() { return fParticleGun; }

 private:
  G4ParticleGun* fParticleGun = nullptr;
  G4PrimaryParticle* particle = nullptr;
  std::vector<double> energies;
  G4double fEnergy;
  SNneutrinosDetectorConstruction* detectordb;
  SNneutrinosPrimaryGeneratorMessenger* fMessenger;
  TRandom3* xGenerator = new TRandom3(0);
  TRandom3* yGenerator = new TRandom3(0);
  TRandom3* zGenerator = new TRandom3(0);
  TRandom3* thetaGenerator = new TRandom3(0);
  TRandom3* phiGenerator = new TRandom3(0);
  //G4String NCaptureModel = "../../NeutronCaptureModels/156Gd-5keV-cascades.txt";
  G4String NCaptureModel = "../../NeutronCaptureModels/ANNRI-Gd155.txt";
};



#endif
