#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "Constants.hh"

EventAction::EventAction() : G4UserEventAction()
{
    runAction = (RunAction*) G4RunManager::GetRunManager()->GetUserRunAction();
}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event* event)
{
    G4int eventID = event->GetEventID();

    runAction->detEnergy[0] = 0.0;
    runAction->detEnergy[1] = 0.0;
    runAction->INI          = 0;

    printf("--> BeginOfEvent %12d\n",eventID);
}

void EventAction::EndOfEventAction(const G4Event*)
{
    runAction->tree->Fill();
    //printf("\n--> EndOfEvent\n");
}
