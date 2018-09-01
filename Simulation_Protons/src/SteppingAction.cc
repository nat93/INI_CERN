#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "Constants.hh"

using namespace std;

SteppingAction::SteppingAction() : G4UserSteppingAction()
{
    runAction = (RunAction*) G4RunManager::GetRunManager()->GetUserRunAction();
}


SteppingAction::~SteppingAction()
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    //G4int TrackID = step->GetTrack()->GetTrackID();
    G4int ParentID = step->GetTrack()->GetParentID();

    G4StepPoint* startPoint = step->GetPreStepPoint();
    G4StepPoint* endPoint = step->GetPostStepPoint();

    G4String volumeName = startPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
    G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();
    uint length = Constants::processName.length();
    if(procName.length() >= length)
    {
        if(procName.compare(procName.size()-length,length,Constants::processName) == 0)
        {
            //printf("\n--> I.N.I. %s INIvolume: %s INItrackID: %d INIparentID: %d",procName.c_str(),volumeName.c_str(),TrackID,ParentID);
            if(volumeName.compare("crystal") == 0 && ParentID == 0)
            {
                runAction->INI++;
            }
        }
    }
}
