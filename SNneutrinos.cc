#include "SNneutrinosDetectorConstruction.hh"
#include "SNneutrinosActionInitialization.hh"
#include "SNneutrinosPrimaryGeneratorAction.hh"
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4TrajectoryDrawByParticleID.hh"
namespace
{
  void PrintUsage()
  {
    G4cerr << " Usage: " << G4endl;

    G4cerr << " SNneutrinos  [-m macro ] [-u UIsession] [-t nThreads] [-r seed] "
           << G4endl;

    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

int main(int argc, char** argv)
{
  // Evaluate arguments
  //
  if(argc > 9)
  {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4String session;
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif

  G4long myseed = 345354;
  for(G4int i = 1; i < argc; i = i + 2)
  {
    if(G4String(argv[i]) == "-m")
      macro = argv[i + 1];
    else if(G4String(argv[i]) == "-u")
      session = argv[i + 1];
    else if(G4String(argv[i]) == "-r")
      myseed = atoi(argv[i + 1]);
#ifdef G4MULTITHREADED
    else if(G4String(argv[i]) == "-t")
    {
      nThreads = G4UIcommand::ConvertToInt(argv[i + 1]);
    }
#endif
    else
    {
      PrintUsage();
      return 1;
    }
  }

  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if(macro.size() == 0)
  {
    ui = new G4UIExecutive(argc, argv);
  }

  // Construct the default run manager
  auto runManager = G4RunManagerFactory::CreateRunManager();
#ifdef G4MULTITHREADED
  if(nThreads > 0)
    runManager->SetNumberOfThreads(nThreads); //default number of thread = 2
#endif

  // Seed the random number generator manually
  G4Random::setTheSeed(myseed);


  // Detector construction

  runManager->SetUserInitialization(new SNneutrinosDetectorConstruction());

  // Physics list
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  //G4NeutronHPphysics* neutronHP = new G4NeutronHPphysics();
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  physicsList->RegisterPhysics(opticalPhysics);
  //physicsList->RegisterPhysics(“neutronHP”);
  runManager->SetUserInitialization(physicsList);

  runManager->SetUserInitialization(new SNneutrinosActionInitialization());


  G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  G4TrajectoryDrawByParticleID* model =
  new G4TrajectoryDrawByParticleID;
  // Configure model
  model->SetDefault("cyan");
  model->Set("opticalphoton", "red");
  model->Set("gamma", "red");
  visManager->RegisterModel(model);


  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if(macro.size())
  {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macro);
  }
  else  // Define UI session for interactive mode
  {
    UImanager->ApplyCommand("/control/execute vis.mac");
    if(ui->IsGUI())
      UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
  }


  delete visManager;
  delete runManager;

  return 0;
}

