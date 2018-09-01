#include "PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4RandomDirection.hh"
#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "Constants.hh"
#include "TMath.h"
#include "stdio.h"
#include "assert.h"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    particleGun = new G4ParticleGun(1);
    particleGun->SetParticleDefinition(G4PionPlus::PionPlus());
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

        G4double theta, phi, theta_x, theta_y;

        theta_x = G4RandGauss::shoot(0.0,Constants::_prim_particle_ang_X_sigma);
        theta_y = G4RandGauss::shoot(0.0,Constants::_prim_particle_ang_Y_sigma);
        theta = std::sqrt (theta_x*theta_x + theta_y*theta_y);
        if (theta != 0.0)
        {
            phi = std::acos(theta_x/theta);
            if ( theta_y < 0.0) phi = -phi;
        }
        else
        {
            phi = 0.0;
        }
        G4double px, py, pz;
        px = -std::sin(theta) * std::cos(phi);
        py = -std::sin(theta) * std::sin(phi);
        pz =  std::cos(theta);
        G4ThreeVector direction(px,py,pz);

        G4double dx = (Constants::_crystal_qmp_pos_Z - Constants::_prim_particle_pos_Z)*TMath::Tan(theta_x);
        G4double dy = (Constants::_crystal_qmp_pos_Z - Constants::_prim_particle_pos_Z)*TMath::Tan(theta_y);

        G4double xx = x + dx;
        G4double yy = y + dy;

        // STFLHC
        if(Constants::_cr_type == 1 && xx <= 0.15*mm && xx >= -0.15*mm && TMath::Abs(yy) <= 4.0*mm)
        {
            trkIsOk = true;
            i = nMax;
        }
        // QMPLHC
        if(Constants::_cr_type == 2 && xx <= -0.1*mm && xx >= -0.5*mm && TMath::Abs(yy) <= 2.0*mm)
        {
            trkIsOk = true;
            i = nMax;
        }

        if(trkIsOk)
        {
            particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
            particleGun->SetParticleMomentumDirection(direction);
            //particleGun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,1.0));
            particleGun->GeneratePrimaryVertex(event);
        }

    }
    if(!trkIsOk)
    {
        G4cout<<" ERROR --> trkIsOK == false, nMax ="<<nMax<<G4endl<<"           acseptance is too small"<<G4endl;
        assert(0);
    }

}
