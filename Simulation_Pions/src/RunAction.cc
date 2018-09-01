#include <string>

#include "RunAction.hh"
#include "G4Run.hh"
#include "Randomize.hh"
#include "Constants.hh"

RunAction::RunAction() {}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*)
{
    G4String fileName = "output_"; fileName += std::to_string((int)Constants::_small_plastic_scint_gap_X); fileName += ".root";

    file = new TFile(fileName.data(),"recreate");
    tree = new TTree("Tree","A Root Tree");

    tree->Branch("detEnergy",   detEnergy,  "detEnergy[2]/D");
    tree->Branch("INI",         &INI,        "INI/I");
}

void RunAction::EndOfRunAction(const G4Run* )
{
    file->Write();
    file->Close();
}

