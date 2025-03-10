
#include "SNneutrinosPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4AnalysisManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VSolid.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "TFile.h"
#include "TTree.h"
#include "G4SubtractionSolid.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"



SNneutrinosPrimaryGeneratorAction::SNneutrinosPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
  , fParticleGun(nullptr)
{
  G4int n_particle = 1;
  fParticleGun     = new G4ParticleGun(n_particle);
  //energies = {1000}; //in keV
  //energies = ReadEnergiesFromANNRIGd(NCaptureModel);
  //energies = ReadEnergiesFromGrabmayer(NCaptureModel);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* def_particle = particleTable->FindParticle("e+"); //e+

  fParticleGun->SetParticleDefinition(def_particle);
  //fParticleGun->SetParticleEnergy(energies[0]*keV);
  
}


SNneutrinosPrimaryGeneratorAction::~SNneutrinosPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void SNneutrinosPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    G4LogicalVolume* tank_lvol = G4LogicalVolumeStore::GetInstance()->GetVolume("Tank_log");
    G4LogicalVolume* pmt_lvol = G4LogicalVolumeStore::GetInstance()->GetVolume("PMT_log");
    G4LogicalVolume* water_only_lvol = G4LogicalVolumeStore::GetInstance()->GetVolume("Water_only_log");

    G4Tubs* tank_solid = nullptr;
    G4Tubs* pmt_solid = nullptr;
    G4SubtractionSolid* water_only_solid = nullptr;
    auto* lvStore = G4LogicalVolumeStore::GetInstance();
    if (tank_lvol)    tank_solid = dynamic_cast<G4Tubs*>(tank_lvol->GetSolid());
    if (pmt_lvol)    pmt_solid = dynamic_cast<G4Tubs*>(pmt_lvol->GetSolid());
    if (water_only_lvol)    water_only_solid = dynamic_cast<G4SubtractionSolid*>(water_only_lvol->GetSolid());
    
    
    G4ThreeVector lo, hi;
    tank_solid->BoundingLimits(lo, hi);  
    G4double TankRadius = hi[0];
    G4double TankHeight = hi[2]*2;
    G4ThreeVector point;
    G4int maxtries=10000, itry=1;
    G4bool validPosition = false;


    do {
      point.set(lo[0] + xGenerator->Uniform(0,1)*(hi[0]-lo[0]),
              lo[1] + yGenerator->Uniform(0,1)*(hi[1]-lo[1]),
              lo[2] + zGenerator->Uniform(0,1)*(hi[2]-lo[2]));
                     
      validPosition = (water_only_solid->Inside(point) == kInside) && 
      (pmt_solid->Inside(point) != kInside);
    } while (!validPosition && ++itry < maxtries);

    if (itry == maxtries)
    G4cerr << "Unable to find a point inside your volume!" << G4endl; 
    //G4cout << "r " << pos_x*pos_x+pos_y*pos_y << " " << 550*550*cm *cm << G4endl;
    //G4cout << "h " << pos_z << " " << 650 *cm << G4endl;
    
    //fParticleGun->SetParticlePosition(point);
    //fParticleGun->SetParticlePosition(G4ThreeVector(0,500*cm,0));

    
    G4double px, py, pz;
    //G4double theta, phi;

   
    G4double theta = CLHEP::twopi/2. * thetaGenerator->Uniform(0,1);
    G4double phi   = CLHEP::twopi * phiGenerator->Uniform(0,1);
    pz = -1* std::cos(theta);
    px = -std::sin(theta) * cos(phi);
    py = -std::sin(theta) * sin(phi);
    G4ThreeVector momentumDir(px, py, pz);
    fParticleGun->SetParticleMomentumDirection(momentumDir);
    fParticleGun->SetParticlePosition(point);
    
  
    G4double energy = 60*MeV; //46.1328 *MeV; //46.1328
    //G4cout << "energy " << energy << G4endl;
    fParticleGun->SetParticleEnergy(energy);

    fParticleGun->GeneratePrimaryVertex(anEvent);
  }


std::vector<double> SNneutrinosPrimaryGeneratorAction::ReadEnergiesFromGrabmayer(const std::string& filename) {
    
    std::ifstream infile(filename);
    std::string line;

    int totalRows = 10000001; //done on terminal (wc -l 158Gd-5keV-cascades.txt)
    std::uniform_int_distribution<> dis(0, totalRows);
    std::random_device rd;
    std::mt19937 gen(rd());
    int randomIndex = dis(gen);
    //int randomIndex = std::rand() % totalRows;
    G4cout << "random index " << randomIndex << G4endl;
    for (int i = 0; i <= randomIndex; ++i) {
        std::getline(infile, line);
    }
    G4cout << "line " << line << G4endl;


    std::istringstream iss(line);
    double energy;
    while (iss >> energy) {
        energies.push_back(energy);
        //G4cout << "single energy " << energy << G4endl;
    }

    energies.erase(energies.begin(), energies.begin() + 4);  //only the last entries are energies

    return energies;
}


std::vector<double> SNneutrinosPrimaryGeneratorAction::ReadEnergiesFromANNRIGd(const std::string& filename) {
    double e;
    double energy;

    std::ifstream infile(filename);
    std::string line;

    int totalRows = 100000; 
    std::uniform_int_distribution<> dis(0, totalRows);
    std::random_device rd;
    std::mt19937 gen(rd());
    int randomIndex = dis(gen)*15;
    //int randomIndex = (std::rand() % totalRows ) * 15; 
    G4cout << "random index " << randomIndex << G4endl;

    for (int i = 0; i < randomIndex; ++i) {
        std::getline(infile, line);
    }

    int j=0;

    for (int j=0; j<1; j++){ //15
      std::getline(infile, line);
      G4cout << line << G4endl;
      std::istringstream iss(line);
      while (iss >> e) {
        energy=e;
      }
      if (energy==0)
        continue;
      else{
        energies.push_back(energy*1000);
        continue;
       // G4cout << "single energy " << energy << G4endl;
      }
    }
    return energies;
}
