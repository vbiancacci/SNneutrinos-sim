#include "SNneutrinosRunAction.hh"
#include "SNneutrinosPrimaryGeneratorAction.hh"
#include "SNneutrinosRun.hh"
#include "G4ParticleDefinition.hh"
#include "G4Run.hh"
#include "G4AnalysisManager.hh"


SNneutrinosRunAction::SNneutrinosRunAction(SNneutrinosPrimaryGeneratorAction* prim)
  : G4UserRunAction()
  , fRun(nullptr)
  , fPrimary(prim)
{


}

SNneutrinosRunAction::~SNneutrinosRunAction()
{}


G4Run* SNneutrinosRunAction::GenerateRun()
{
  fRun = new SNneutrinosRun();
  return fRun;
}

void SNneutrinosRunAction::BeginOfRunAction(const G4Run* aRun)
{
  if(fPrimary)
  {
    G4ParticleDefinition* particle =
      fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
    fRun->SetPrimary(particle, energy);
  }

      G4cout << "### Run " << aRun->GetRunID() << " starts." << G4endl;

    G4AnalysisManager *analysis = G4AnalysisManager::Instance();
    analysis->OpenFile("test/10MeV/SNneutrinos.root"); 

    //ntupla
    analysis->CreateNtuple("Score", "Score");
    analysis->CreateNtupleDColumn("Event_ID");                    //0
    analysis->CreateNtupleDColumn("e_posx_in_cm");               //1
    analysis->CreateNtupleDColumn("e_posy_in_cm");               //2
    analysis->CreateNtupleDColumn("e_posz_in_cm");               //3
    analysis->CreateNtupleDColumn("e_momx");                     //4
    analysis->CreateNtupleDColumn("e_momy");                     //5
    analysis->CreateNtupleDColumn("e_momz");                     //6
    analysis->CreateNtupleDColumn("photon_posx_in_cm");                 //7
    analysis->CreateNtupleDColumn("photon_posy_in_cm");                 //8
    analysis->CreateNtupleDColumn("photon_posz_in_cm");                 //9
    analysis->CreateNtupleDColumn("photon_kinetic_energy_in_eV"); //10
    analysis->CreateNtupleDColumn("PMT_ID");                      //11
    analysis->CreateNtupleDColumn("PMT_efficiency");              //12

    analysis->FinishNtuple(0);

}


void SNneutrinosRunAction::EndOfRunAction(const G4Run*)
{
  if(isMaster)
    fRun->EndOfRun();

  G4AnalysisManager *analysis = G4AnalysisManager::Instance();
  analysis->Write();
  analysis->CloseFile("test/10MeV/SNneutrinos.root"); 
}

