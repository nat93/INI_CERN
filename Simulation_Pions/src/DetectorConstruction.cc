#include "DetectorConstruction.hh"
#include "DetectorSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4PhysicalConstants.hh"
#include "Constants.hh"

DetectorConstruction::DetectorConstruction() {}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4NistManager* nistManager = G4NistManager::Instance();

    nistManager->FindOrBuildMaterial("G4_AIR");
    nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    nistManager->FindOrBuildMaterial("G4_Si");
    nistManager->FindOrBuildMaterial("G4_Galactic");
    nistManager->FindOrBuildMaterial("G4_Ti");
    nistManager->FindOrBuildMaterial("G4_Al");


    //G4Material* AIR = G4Material::GetMaterial("G4_AIR");
    //G4Material* PLASTIC = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4Material* SILICON = G4Material::GetMaterial("G4_Si");
    G4Material* GALACTIC = G4Material::GetMaterial("G4_Galactic");
    G4Material* TITANIUM = G4Material::GetMaterial("G4_Ti");
    G4Material* ALUMINIUM = G4Material::GetMaterial("G4_Al");

    G4Element* elC = new G4Element("Carbon","C",6.,12.01*g/mole);
    G4Element* elH = new G4Element("Hydrogen","H2",1.,1.01*g/mole);
    G4Material* POLYSTYRENE = new G4Material("Polystyrene",1.032*g/cm3,2);
    POLYSTYRENE->AddElement(elC,8);
    POLYSTYRENE->AddElement(elH,8);

    // WORLD
    G4Box* world_box = new G4Box("world", Constants::_world_width/2, Constants::_world_height/2, Constants::_world_thick/2);
    //G4LogicalVolume* world_log = new G4LogicalVolume(world_box, AIR, "world");
    G4LogicalVolume* world_log = new G4LogicalVolume(world_box, GALACTIC, "world");
    G4VPhysicalVolume* world_phys = new G4PVPlacement(0, G4ThreeVector(), world_log, "world", 0, false, 0);

    if(Constants::_sw_tracker)
    {
        // SILICON TRACKER PLATE 2
        G4Box* si_box = new G4Box("tracker_plate2", Constants::_si_width/2, Constants::_si_height/2, Constants::_si_thick/2);
        G4LogicalVolume* si2_log = new G4LogicalVolume(si_box, SILICON, "tracker_plate2");
        new G4PVPlacement(0, G4ThreeVector(Constants::_si_pos_X, Constants::_si_pos_Y,Constants::_crystal_stf_pos_Z-Constants::_gap_si_cr), si2_log, "tracker_plate2", world_log, false, 0);
        // ALUMINIUM FOIL
        G4Box* fl_box = new G4Box("aluminium_foil", Constants::_foil_width/2, Constants::_foil_height/2, Constants::_foil_thick/2);
        // FOIL 1 FOR Si2
        G4LogicalVolume* fl21_log = new G4LogicalVolume(fl_box, ALUMINIUM, "foil_plate1_for_tracker_plate2");
        new G4PVPlacement(0, G4ThreeVector(Constants::_si_pos_X, Constants::_si_pos_Y,Constants::_crystal_stf_pos_Z-Constants::_gap_si_cr - (Constants::_si_thick/2 + Constants::_foil_thick/2)),
                          fl21_log, "foil_plate1_for_tracker_plate2", world_log, false, 0);
        // FOIL 2 FOR Si2
        G4LogicalVolume* fl22_log = new G4LogicalVolume(fl_box, ALUMINIUM, "foil_plate2_for_tracker_plate2");
        new G4PVPlacement(0, G4ThreeVector(Constants::_si_pos_X, Constants::_si_pos_Y,Constants::_crystal_stf_pos_Z-Constants::_gap_si_cr + (Constants::_si_thick/2 + Constants::_foil_thick/2)),
                          fl22_log, "foil_plate2_for_tracker_plate2", world_log, false, 0);

        si2_log->SetVisAttributes(G4VisAttributes(G4Color::Red()));
        fl21_log->SetVisAttributes(G4VisAttributes(G4Color::Yellow()));
        fl22_log->SetVisAttributes(G4VisAttributes(G4Color::Yellow()));
    }

    if(Constants::_cr_type == 1)
    {
        // CRYSTAL
        G4Box* cr_box = new G4Box("crystal",
                                  Constants::_crystal_stf_width/2,
                                  Constants::_crystal_stf_height/2,
                                  Constants::_crystal_stf_thick/2);
        G4LogicalVolume* cr_log = new G4LogicalVolume(cr_box, SILICON, "crystal");
        new G4PVPlacement(0, G4ThreeVector(Constants::_crystal_stf_pos_X,
                                           Constants::_crystal_stf_pos_Y,
                                           Constants::_crystal_stf_pos_Z), cr_log, "crystal", world_log, false, 0);

        cr_log->SetVisAttributes(G4VisAttributes(G4Color::Green()));

        if(Constants::_sw_holder)
        {
            // HOLDER
            G4Box* hld1_box = new G4Box("holder1_stf",5.0*mm/2,31.0*mm/2,14.0*mm/2);
            G4LogicalVolume* hld1_log = new G4LogicalVolume(hld1_box, TITANIUM, "holder1_stf");
            new G4PVPlacement(0, G4ThreeVector(Constants::_crystal_stf_pos_X - 16.5*mm - Constants::_crystal_stf_width/2,Constants::_crystal_stf_pos_Y + 0.0,Constants::_crystal_stf_pos_Z), hld1_log, "holder1", world_log, false, 0);

            G4Box* hld2_box = new G4Box("holder2_stf",33.0*mm/2,12.0*mm/2,40.0*mm/2);
            G4LogicalVolume* hld2_log = new G4LogicalVolume(hld2_box, TITANIUM, "holder2_stf");
            new G4PVPlacement(0, G4ThreeVector(Constants::_crystal_stf_pos_X - 16.5*mm - Constants::_crystal_stf_width/2,Constants::_crystal_stf_pos_Y + 21.5*mm,Constants::_crystal_stf_pos_Z), hld2_log, "holder2", world_log, false, 0);

            G4Box* hld3_box = new G4Box("holder3_stf",33.0*mm/2,12.0*mm/2,40.0*mm/2);
            G4LogicalVolume* hld3_log = new G4LogicalVolume(hld3_box, TITANIUM, "holder3_stf");
            new G4PVPlacement(0, G4ThreeVector(Constants::_crystal_stf_pos_X - 16.5*mm - Constants::_crystal_stf_width/2,Constants::_crystal_stf_pos_Y - 21.5*mm,Constants::_crystal_stf_pos_Z), hld3_log, "holder3", world_log, false, 0);

            hld1_log->SetVisAttributes(G4VisAttributes(G4Color::Gray()));
            hld2_log->SetVisAttributes(G4VisAttributes(G4Color::Gray()));
            hld3_log->SetVisAttributes(G4VisAttributes(G4Color::Gray()));
        }
    }
    else if(Constants::_cr_type == 2)
    {
        // CRYSTAL
        G4Box* cr_box = new G4Box("crystal",
                                  Constants::_crystal_qmp_width/2,
                                  Constants::_crystal_qmp_height/2,
                                  Constants::_crystal_qmp_thick/2);
        G4LogicalVolume* cr_log = new G4LogicalVolume(cr_box, SILICON, "crystal");
        new G4PVPlacement(0, G4ThreeVector(Constants::_crystal_qmp_pos_X,
                                           Constants::_crystal_qmp_pos_Y,
                                           Constants::_crystal_qmp_pos_Z), cr_log, "crystal", world_log, false, 0);

        cr_log->SetVisAttributes(G4VisAttributes(G4Color::Green()));

        if(Constants::_sw_holder)
        {
            // HOLDER
            G4Box* hld_box = new G4Box("holder_qmp",23.0*mm/2,40.0*mm/2,25.0*mm/2);
            G4LogicalVolume* hld_log = new G4LogicalVolume(hld_box, TITANIUM, "holder_qmp");
            new G4PVPlacement(0, G4ThreeVector(Constants::_crystal_qmp_pos_X - 11.5*mm - Constants::_crystal_qmp_width/2,Constants::_crystal_qmp_pos_Y,Constants::_crystal_qmp_pos_Z), hld_log, "holder_qmp", world_log, false, 0);

            hld_log->SetVisAttributes(G4VisAttributes(G4Color::Gray()));
        }
    }

    // FIRST PLASTIC SCINTILLATOR S1
    G4Box* sc1_box = new G4Box("plastic_scintillator1",
                               Constants::_small_plastic_scint_width/2,
                               Constants::_small_plastic_scint_height/2,
                               Constants::_small_plastic_scint_thick/2);
    //G4LogicalVolume* sc1_log = new G4LogicalVolume(sc1_box, PLASTIC, "plastic_scintillator1");
    G4LogicalVolume* sc1_log = new G4LogicalVolume(sc1_box, POLYSTYRENE, "plastic_scintillator1");
    new G4PVPlacement(0, G4ThreeVector((Constants::_small_plastic_scint_gap_X + Constants::_small_plastic_scint_width - 2*Constants::_small_plastic_scint_shift_X)/2,
                                       Constants::_small_plastic_scint_pos_Y,
                                       Constants::_small_plastic_scint_pos_Z), sc1_log, "plastic_scintillator1", world_log, false, 0);

    // SECOND PLASTIC SCINTILLATOR S2
    G4Box* sc2_box = new G4Box("plastic_scintillator2",
                               Constants::_small_plastic_scint_width/2,
                               Constants::_small_plastic_scint_height/2,
                               Constants::_small_plastic_scint_thick/2);
    //G4LogicalVolume* sc2_log = new G4LogicalVolume(sc2_box, PLASTIC, "plastic_scintillator2");
    G4LogicalVolume* sc2_log = new G4LogicalVolume(sc2_box, POLYSTYRENE, "plastic_scintillator2");
    new G4PVPlacement(0, G4ThreeVector(-(Constants::_small_plastic_scint_gap_X + Constants::_small_plastic_scint_width + 2*Constants::_small_plastic_scint_shift_X)/2,
                                       Constants::_small_plastic_scint_pos_Y,
                                       Constants::_small_plastic_scint_pos_Z), sc2_log, "plastic_scintillator2", world_log, false, 0);


    // SENSITIVE DETECTOR S1
    DetectorSD* detectorSD1 = new DetectorSD("plastic_scintillator1");
    G4SDManager* sdMan1 = G4SDManager::GetSDMpointer();
    sdMan1->AddNewDetector(detectorSD1);
    sc1_log->SetSensitiveDetector(detectorSD1);

    // SENSITIVE DETECTOR S2
    DetectorSD* detectorSD2 = new DetectorSD("plastic_scintillator2");
    G4SDManager* sdMan2 = G4SDManager::GetSDMpointer();
    sdMan2->AddNewDetector(detectorSD2);
    sc2_log->SetSensitiveDetector(detectorSD2);

    //world_log->SetVisAttributes(G4VisAttributes::Invisible);
    sc1_log->SetVisAttributes(G4VisAttributes(G4Color::Blue()));
    sc2_log->SetVisAttributes(G4VisAttributes(G4Color::Blue()));

    return world_phys;
}
