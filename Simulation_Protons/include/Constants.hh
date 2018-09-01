#ifndef Constants_hh
#define Constants_hh

// GEANT
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "globals.hh"
#include "string"

namespace Constants
{
    // Physics process
    static const G4String processName = "Inelastic";
    // World ranges
    static const G4double _world_thick = 700.0*mm;
    static const G4double _world_width = 10.0*cm;
    static const G4double _world_height = 10.0*cm;
    // Primary particle position
    static const G4double _prim_particle_pos_Z = -349.0*mm;
    static const G4double _prim_particle_pos_X = 0.0*mm;
    static const G4double _prim_particle_pos_X_sigma = 1.0*mm;
    static const G4double _prim_particle_pos_Y = 0.0*mm;
    static const G4double _prim_particle_pos_Y_sigma = 1.0*mm;
    // Plastic scintillator ranges
    // Small S1 & S2
    static const G4double _small_plastic_scint_thick = 5.0*mm;
    static const G4double _small_plastic_scint_width = 10.0*mm;
    static const G4double _small_plastic_scint_height = 25.0*mm;    
    // Plastic scintillator position
    // Small S1 & S2
    static const G4double _small_plastic_scint_pos_Y = 0.0*mm;
    static const G4double _small_plastic_scint_pos_Z = 245.0*mm;
    static const G4double _small_plastic_scint_gap_X = 10.0*mm;
    // Crystal ranges STFLHC
    static const G4double _crystal_stf_thick = 4.0*mm;
    static const G4double _crystal_stf_width = 0.55*mm;
    static const G4double _crystal_stf_height = 55.0*mm;
    // Crystal ranges STF101
    // static const G4double _crystal_stf_thick = 2.0*mm;
    // static const G4double _crystal_stf_width = 1.0*mm;
    // static const G4double _crystal_stf_height = 55.0*mm;
    // Crystal position STF
    static const G4double _crystal_stf_pos_X = 0.0*mm;
    static const G4double _crystal_stf_pos_Y = 0.0*mm;
    static const G4double _crystal_stf_pos_Z = 0.0*mm;
    // Crystal ranges QMPLHC
    static const G4double _crystal_qmp_thick = 4.0*mm;
    static const G4double _crystal_qmp_width = 5.0*mm;
    static const G4double _crystal_qmp_height = 25.0*mm;
    // Crystal position QMP
    static const G4double _crystal_qmp_pos_X = 0.0*mm;
    static const G4double _crystal_qmp_pos_Y = 0.0*mm;
    static const G4double _crystal_qmp_pos_Z = 0.0*mm;
    // Foil ranges
    static const G4double _foil_thick = 0.03*mm;
    static const G4double _foil_width = 38.0*mm;
    static const G4double _foil_height = 38.0*mm;
    // Tracker plate position
    static const G4double _si_pos_X = 0.0*mm;
    static const G4double _si_pos_Y = 0.0*mm;
    static const G4double _gap_si_cr = 332.0*mm;
    // Tracker plate ranges
    static const G4double _si_thick = 0.640*mm;
    static const G4double _si_width = 38.0*mm;
    static const G4double _si_height = 38.0*mm;
    // switch on/off tracker
    static const G4bool _sw_tracker = false;//true;
    // switch on/off holder
    static const G4bool _sw_holder = false;//true;
    // switch STF(1)/QMP(2) crystal
    static const G4int _cr_type = 1; // 1 - STF; 2 - QMP
}

#endif
