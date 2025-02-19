
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

SNneutrinosPrimaryGeneratorAction::SNneutrinosPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
  , fParticleGun(nullptr)
{
  G4int n_particle = 1;
  fParticleGun     = new G4ParticleGun(n_particle);
  energies = {1000}; //in keV
  //energies = ReadEnergiesFromANNRIGd(NCaptureModel);
  //energies = ReadEnergiesFromGrabmayer(NCaptureModel);
  
  // create a messenger for this class
  // default kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* def_particle = particleTable->FindParticle("e+"); //e+

  fParticleGun->SetParticleDefinition(def_particle);
  fParticleGun->SetParticleEnergy(energies[0]*keV);
  //fParticleGun->SetParticleTime(0.0 * ns);
  
}


SNneutrinosPrimaryGeneratorAction::~SNneutrinosPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void SNneutrinosPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    /*
     //
    // Simulate only inside the water tank.
    // The volume is divided into three subvolume: top, middle, bottom.
    // The middle one contains the outer cryostat
    //

    //G4double WaterTankHeight = detectordb->GetWaterTankHeight(); //(650 + 0.8) * cm;
    //G4cout << WaterTankHeight << " !!!!!!!!!!!!!!!!!!!" << G4endl;
    G4double WaterTankHeight = (650 + 0.8) * cm;
    G4double WaterTankRadius = (550 + 0.6) * cm;
    G4double CryostatHeight = 350 * cm;
    G4double CryostatRadius = 350 * cm; 
    //G4double Offset          = (200 - 100 - (850 - 650)) * cm;
  
    G4double Vol_top    = CLHEP::pi * WaterTankRadius * WaterTankRadius * (WaterTankHeight-CryostatHeight);
    G4double Vol_middle = CLHEP::pi * (WaterTankRadius - CryostatRadius) * (WaterTankRadius - CryostatRadius)  * (CryostatHeight);
    G4double Vol_bottom = Vol_top; 

    G4double Prob_middle = Vol_middle / (Vol_bottom + Vol_middle + Vol_top);
    G4double Prob_top = (1 - Prob_middle) / 2.;
    G4double Prob_bottom = Prob_top;
   
    std::discrete_distribution<> distribution_2({ Prob_middle, Prob_top, Prob_bottom });
    generator.seed(rd());
    std::uniform_real_distribution<> rndm(0.0, 1.0);  // azimuth angle

    G4int    where = distribution_2(generator);

    G4double px, py, pz, pos_x, pos_y, pos_z;
    G4double theta = CLHEP::twopi * rndm(generator);
    G4double phi   = 2 * CLHEP::twopi *rndm(generator);
    //G4double theta = CLHEP::twopi/2. * Generator->Uniform(0,1);
    //G4double phi   = CLHEP::twopi * Generator->Uniform(0,1);
    //G4cout << theta << G4endl;
    //G4cout << phi << G4endl;
    
  if(where == 0)
    {
      G4double pos_phi    = CLHEP::twopi * rndm(generator);

      //pos_x = ( (WaterTankRadius - CryostatRadius) * Generator->Uniform(0,1) + CryostatRadius) * cos(pos_phi);
      //pos_y = ( (WaterTankRadius - CryostatRadius) * Generator->Uniform(0,1) + CryostatRadius) * sin(pos_phi);
      //pos_z = CryostatHeight * (1 - 2 * Generator->Uniform(0,1));

      pos_x = ( (WaterTankRadius - CryostatRadius) * rndm(generator) + CryostatRadius) * cos(pos_phi);
      pos_y = ( (WaterTankRadius - CryostatRadius) * rndm(generator) + CryostatRadius) * sin(pos_phi);
      pos_z = CryostatHeight * (1 - 2 * rndm(generator));
    }
  
  if(where > 0)
    {
      G4double pos_phi = CLHEP::twopi * rndm(generator);
      G4double pos_height = (WaterTankHeight - CryostatHeight) * rndm(generator) + CryostatHeight;
      G4double pos_rad = WaterTankRadius * rndm(generator);
      if(where == 2) //bottom side
        pos_height *=-1;

      pos_x = pos_rad * cos(pos_phi);
      pos_y = pos_rad * sin(pos_phi);
      pos_z = pos_height;
    }


 
    //G4double pos_z = (WaterTankHeight - CryostatHeight) * G4UniformRand() + CryostatHeight;
    //G4double pos_z = ((650-350.) * G4UniformRand() + 350.01) *cm;
    */
    //G4LogicalVolume* lvol = G4PhysicalVolumeStore::GetInstance()->GetVolume("World_phys")->GetLogicalVolume();
    
    
    G4LogicalVolume* envelope_lvol = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope_log");
    G4LogicalVolume* water_lvol = G4LogicalVolumeStore::GetInstance()->GetVolume("Water_log");
    G4Tubs* envelope_solid = nullptr;
    G4Polycone* water_solid = nullptr;
    auto* lvStore = G4LogicalVolumeStore::GetInstance();
    if (envelope_lvol)    envelope_solid = dynamic_cast<G4Tubs*>(envelope_lvol->GetSolid());
    if (water_lvol)    water_solid = dynamic_cast<G4Polycone*>(water_lvol->GetSolid());
    //auto* pvol_w = G4PhysicalVolumeStore::GetInstance()->GetVolume("World_phys");
    //G4VSolid solid_w = pvol_w->GetLogicalVolume()->GetSolid();
    G4ThreeVector lo, hi;
    envelope_solid->BoundingLimits(lo, hi);
    
    G4ThreeVector point;
    G4int maxtries=10000, itry=1;
    do {
    point.set(lo[0] + xGenerator->Uniform(0,1)*(hi[0]-lo[0]),
            lo[1] + yGenerator->Uniform(0,1)*(hi[1]-lo[1]),
            lo[2] + zGenerator->Uniform(0,1)*(hi[2]-lo[2]));
    } while (!water_solid->Inside(point) && ++itry < maxtries);

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
    
    //std::uniform_real_distribution<double> rndm_energy (1.3, 60.0);
    //G4double energy = rndm_energy(generator)*MeV;
    //energy = energy * MeV;
    G4double energy = 60*MeV; //46.1328 *MeV; //46.1328
    //G4cout << "energy " << energy << G4endl;
    fParticleGun->SetParticleEnergy(energy);

    //G4double theMass = particleTable->FindParticle("e+")->GetPDGMass();
    //G4double totMomentum = std::sqrt(energy * energy + 2 * theMass * energy);    

    fParticleGun->GeneratePrimaryVertex(anEvent);
  
    /*
    G4PrimaryVertex* vertex = new G4PrimaryVertex(point,0);


    for (const auto& energy : energies) {
      //G4cout << "energy " << energy << G4endl;
      G4PrimaryParticle* particle = new G4PrimaryParticle(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
        particle->SetTotalEnergy(energy*keV);
        theta = CLHEP::twopi/2. * thetaGenerator->Uniform(0,1);
        phi   = CLHEP::twopi * phiGenerator->Uniform(0,1);
        pz = -1* std::cos(theta);
        px = -std::sin(theta) * cos(phi);
        py = -std::sin(theta) * sin(phi);
        G4ThreeVector momentumDir(px, py, pz);
        particle->SetMomentumDirection(momentumDir);
        vertex->SetPrimary(particle);
        fParticleGun->SetParticleEnergy(energy);
    }

    


    anEvent->AddPrimaryVertex(vertex);

*/
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
