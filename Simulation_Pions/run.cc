#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "G4RunManager.hh"
#include "G4eMultipleScattering.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"
#include "G4PhysListFactory.hh"
#include "EventAction.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BIC.hh"

#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char** argv)
{
  if(argc < 1) cout<<"WTF???"<<endl;
  G4PhysListFactory factory;

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(time(NULL));

  G4RunManager* runManager = new G4RunManager;
  G4String physName = "QGSP_BERT";
  //G4String physName = "QGSP_BIC";
  G4VModularPhysicsList* phys = factory.GetReferencePhysList(physName);

  runManager->SetUserInitialization(phys);
  runManager->SetUserInitialization(new DetectorConstruction);
  runManager->SetUserAction(new PrimaryGeneratorAction);
  runManager->SetUserAction(new RunAction);
  EventAction* eventAction = new EventAction;
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(new SteppingAction);

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  runManager->Initialize();

  G4UImanager* UI = G4UImanager::GetUIpointer();

  UI->ExecuteMacroFile(argv[1]);

  delete visManager;
  delete runManager;

  return 0;
}
