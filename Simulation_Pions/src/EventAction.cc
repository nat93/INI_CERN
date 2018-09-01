#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "Constants.hh"

#include "iostream"

using namespace std;

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

    //printf("--> BeginOfEvent %12d\n",eventID);

    if(eventID%1000 == 0)
    {
        printf("\r--> Begin of event: %10d", eventID);
        fflush(stdout);
    }
}

void EventAction::EndOfEventAction(const G4Event*)
{
    runAction->tree->Fill();
    //printf("\n--> EndOfEvent\n");
}
