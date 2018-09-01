#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

using namespace std;

class G4Run;

class RunAction: public G4UserRunAction
{
public:
    RunAction();
    ~RunAction();

public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

    TFile* file;
    TTree* tree;

    Double_t    detEnergy[2];
    Int_t       INI;
};

#endif
