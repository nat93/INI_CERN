#include "DetectorSD.hh"
#include "RunAction.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "Constants.hh"
#include <assert.h>


#include <string>

using namespace std;

DetectorSD::DetectorSD(G4String name): G4VSensitiveDetector(name)
{
    runAction = (RunAction*) G4RunManager::GetRunManager()->GetUserRunAction();
    fEventAction = (EventAction*) G4EventManager::GetEventManager()->GetUserEventAction();
}

DetectorSD::~DetectorSD() {}

void DetectorSD::Initialize(G4HCofThisEvent*)
{
    detEnergy = 0;
}

G4bool DetectorSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    G4double edep = step->GetTotalEnergyDeposit()/keV;
    detEnergy += edep;    

    return true;
}

void DetectorSD::EndOfEvent(G4HCofThisEvent*)
{
    //printf("\n--> Total deposit energy in %s %10.10f [keV]",GetName().c_str(),detEnergy);

    string detName = GetName().c_str();
    if(detName.compare("plastic_scintillator1") == 0)
    {
        runAction->detEnergy[0] = detEnergy;
    }
    else if(detName.compare("plastic_scintillator2") == 0)
    {
        runAction->detEnergy[1] = detEnergy;
    }
    else
    {
        G4cout<<"ERROR:: The current name ("<<detName<<") of the SD detector is wrong!"<<G4endl;
        assert(0);
    }
}

