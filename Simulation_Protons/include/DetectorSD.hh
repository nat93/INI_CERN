#ifndef DetectorSD_h
#define DetectorSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "EventAction.hh"
#include "G4Event.hh"

class EventAction;
class G4Step;
class RunAction;

class DetectorSD: public G4VSensitiveDetector 
{
  public:
    DetectorSD(G4String);
    ~DetectorSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);
    
  private:
    RunAction* runAction;
    EventAction*  fEventAction;

    G4double detEnergy;
};

#endif
