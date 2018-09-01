#include "PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4RandomDirection.hh"
#include "G4Proton.hh"
#include "Constants.hh"
#include "TMath.h"
#include "stdio.h"
#include "assert.h"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    particleGun = new G4ParticleGun(1);
    particleGun->SetParticleDefinition(G4Proton::Proton());
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    /*
    G4double x = Constants::_prim_particle_pos_X;
    G4double y = Constants::_prim_particle_pos_Y;
    G4double z = Constants::_prim_particle_pos_Z;

    particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,1.0));

    particleGun->GeneratePrimaryVertex(event);
    */

    G4double x = 0.0, y = 0.0, z = 0.0;
    G4bool trkIsOk = false;
    G4int nMax = 1000000;
    for(Int_t i = 0; i < nMax; i++)
    {
        x = G4RandGauss::shoot(Constants::_prim_particle_pos_X, Constants::_prim_particle_pos_X_sigma);
        y = G4RandGauss::shoot(Constants::_prim_particle_pos_Y, Constants::_prim_particle_pos_Y_sigma);
        z = Constants::_prim_particle_pos_Z;

        // STF101
        // if(Constants::_cr_type == 1 && x <= 0.45*mm && x >= -0.45*mm && TMath::Abs(y) <= 2.0*mm)
        //
        // STFLHC
        if(Constants::_cr_type == 1 && x <= 0.2*mm && x >= -0.2*mm && TMath::Abs(y) <= 2.0*mm)
        {
            trkIsOk = true;
            i = nMax;
        }
        // QMPLHC
        if(Constants::_cr_type == 2 && x <= 0.0*mm && x >= -0.5*mm && TMath::Abs(y) <= 2.0*mm)
        {
            trkIsOk = true;
            i = nMax;
        }

        if(trkIsOk)
        {
            particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
            particleGun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,1.0));

            particleGun->GeneratePrimaryVertex(event);
        }

    }
    if(!trkIsOk)
    {
        G4cout<<" ERROR --> trkIsOK == false, nMax ="<<nMax<<G4endl<<"           acseptance is too small"<<G4endl;
        assert(0);
    }

}
