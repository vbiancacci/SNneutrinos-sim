#include "SNneutrinosSteppingAction.hh"
#include "SNneutrinosRun.hh"
#include "G4Event.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4AnalysisManager.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4SystemOfUnits.hh"




SNneutrinosSteppingAction::SNneutrinosSteppingAction(SNneutrinosEventAction* event) : G4UserSteppingAction(), fEventAction(event) {
  ResetBounceCounter();
}

SNneutrinosSteppingAction::~SNneutrinosSteppingAction() {}

G4int SNneutrinosSteppingAction::GetNumberOfBounces() { return fCounterBounce; }

void SNneutrinosSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();

    SNneutrinosRun* run = static_cast<SNneutrinosRun*>( 
      G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  
  static G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();

  G4int event_ID;



  G4Track* track = step->GetTrack();
  
  const G4ParticleDefinition* particleDef = track->GetDynamicParticle()->GetParticleDefinition();

  auto secTrack = step->GetSecondaryInCurrentStep();
  G4double secSize = (*secTrack).size();

  G4StepPoint* endPoint   = step->GetPostStepPoint();
  G4StepPoint* startPoint = step->GetPreStepPoint();
  
  G4String startVolumeName = startPoint->GetTouchable()->GetVolume() -> GetLogicalVolume() -> GetName();

  G4String endVolumeName = endPoint->GetTouchable()->GetVolume() -> GetLogicalVolume() -> GetName();



  if (endVolumeName == "World_log" || startVolumeName == "World_log"){
    //auto voll =  startPoint->GetTouchable()->GetVolume() -> GetLogicalVolume() -> GetName();
    //auto test =  step->GetTrack()->GetNextVolume()->GetName();
    const G4Event* evttt = G4RunManager::GetRunManager()->GetCurrentEvent();
     G4int event_ID = evttt->GetEventID();

    //G4cout << "event " << event_ID << " particle " << particleDef->GetParticleName() << " kin energy " << (track->GetDynamicParticle()->GetKineticEnergy())/eV << G4endl;

    //G4cout << "step " << track->GetCurrentStepNumber() << G4endl;
    //G4cout << "vol " << voll << G4endl;
    //G4cout << "test " << test << G4endl;
    ResetBounceCounter();
    track->SetTrackStatus(fStopAndKill);
  }


  if (track->GetVolume() and track->GetNextVolume()){
   
    G4int PMT_ID;

    G4ThreeVector photon_pos (0.,0.,0.);
    G4ThreeVector vertex_pos (0., 0., 0.);

    G4double kin_energy = 0.;

    G4StepPoint* aStepPoint = 0;

    auto proc_man =
      track->GetDynamicParticle()->GetParticleDefinition()->GetProcessManager();
    G4ProcessVector* proc_vec = proc_man->GetPostStepProcessVector(typeDoIt);
    G4int n_proc              = proc_vec->entries();

    G4int n_scint = 0;
    G4int n_cer   = 0;
    G4int n_other = 0;
    for(G4int i = 0; i < n_proc; ++i)
    {
      G4String proc_name = (*proc_vec)[i]->GetProcessName();
      if(proc_name.compare("Cerenkov") == 0)
      {
        auto cer = (G4Cerenkov*) (*proc_vec)[i];
        n_cer    = cer->GetNumPhotons();
        fEventAction->AddCerenkovEvent(n_cer);
      }
      else if(proc_name.compare("Scintillation") == 0)
      {
        auto scint = (G4Scintillation*) (*proc_vec)[i];
        n_scint    = scint->GetNumPhotons();
      }
      else{
        n_other+=1;
      }
    }
    // if (n_cer>0 || n_scint>0){
    //   G4cout << "In this step, " << n_cer << " Cerenkov and " << n_scint
    //           << " scintillation photons were produced." << n_other << " other transition radiation like Rayleigh, Absorption, Mie, WLS, Boundary, Transportation." << G4endl;
    // }
       
    if(particleDef == opticalphoton){

      if(track->GetParentID() > 0){ 
        G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();

        //G4cout << procName << " " ;
        if(procName.compare("OpAbsorption") == 0) 
          run->AddAbsorption();
       
        //if(procName.compare("OpWLS") == 0) 
        //  G4cout << "!!!!! WLS" << G4endl;

        aStepPoint = endPoint;
      

        if (startVolumeName == "Water_log"){
          
            if(endVolumeName == "Water_log" ){;//|| endVolumeName == "Cryostat_log" || endVolumeName == "Tank_log"){
                track->SetTrackStatus(fStopAndKill);
                run->AddWaterDetection(); //check the absorption
                ResetBounceCounter();
            }
          
            if(endVolumeName == "PMT_log"){
                track->SetTrackStatus(fStopAndKill);
                ResetBounceCounter();
                kin_energy = track->GetDynamicParticle()->GetKineticEnergy();

                if (G4UniformRand()*100 < PMT_QE("PMT",kin_energy/eV)*collection_efficiency){

                  fEventAction->AddPMTDetectionEvent();
      
                  const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
                  //if (evt)
                  event_ID = evt->GetEventID();
                  analysis->FillNtupleDColumn(0, event_ID);
        
                  PMT_ID = endPoint->GetTouchable()->GetCopyNumber();
                  analysis->FillNtupleDColumn(11, PMT_ID);
                  
                  photon_pos = aStepPoint->GetPosition();
                  analysis->FillNtupleDColumn(7, photon_pos[0]/cm);  //pos_x
                  analysis->FillNtupleDColumn(8, photon_pos[1]/cm);  //pos_y
                  analysis->FillNtupleDColumn(9, photon_pos[2]/cm);  //pos_z
                  
                  analysis->FillNtupleDColumn(10, kin_energy/eV);

                  analysis->AddNtupleRow(0);
              }
            }
            //else G4cout << "end Volume " << endVolumeName << G4endl;
        }
        else {
            //G4cout << "startVolume " << startVolumeName << G4endl;     
        }
      
      // for boundary scattering, process name in 'transportation'.
      // Need to check differently:
      if(endPoint->GetStepStatus() == fGeomBoundary)
      {
        G4OpBoundaryProcessStatus theStatus = Undefined;
        G4ProcessManager* opManager         = opticalphoton->GetProcessManager();
        G4int n_proc = opManager->GetPostStepProcessVector(typeDoIt)->entries();
        G4ProcessVector* postStepDoItVector =
          opManager->GetPostStepProcessVector(typeDoIt);
        for(G4int i = 0; i < n_proc; ++i)
        {
          G4VProcess* currentProcess = (*postStepDoItVector)[i];

          G4OpBoundaryProcess* opProc =
            dynamic_cast<G4OpBoundaryProcess*>(currentProcess);
          if(opProc)
            theStatus = opProc->GetStatus();
        }
        if(theStatus != Undefined && theStatus != NotAtBoundary &&
          theStatus != StepTooSmall)
        {
          run->AddBoundary();
        }
        //Kill the track if it's number of bounces exceeded the limit (optical photon trapped in the foil)
        if (theStatus==SpikeReflection){
          fCounterBounce++;
        }
        
        if(theStatus==TotalInternalReflection){
          //G4cout << "here !!!!!!!!!!!!!!!!! "  << GetNumberOfBounces() << G4endl;
          G4int fBounceLimit = 300;
          if(fBounceLimit > 0 && fCounterBounce >= fBounceLimit){
            track->SetTrackStatus(fStopAndKill);
            ResetBounceCounter();
            //G4cout << "\n Bounce Limit Exceeded" << G4endl;
          }
        }
      }
    }
    else{
      G4cout << particleDef << G4endl;
    }
  }
  }
}
