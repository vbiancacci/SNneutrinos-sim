#ifndef SNneutrinosDetectorConstruction_h
#define SNneutrinosDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4SystemOfUnits.hh"
#include "TGraph.h"
#include <cmath>


class SNneutrinosDetectorMessenger;

class SNneutrinosDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  SNneutrinosDetectorConstruction();
  ~SNneutrinosDetectorConstruction();

  G4VPhysicalVolume* Construct() override;
  void SetVerbose(G4bool verbose);
  G4bool IsVerbose() const;


  G4double stone       = 500.0;  // Hall wall thickness 5 m
  G4double hallrad     = 600.0;  // Hall cylinder diam 12 m
  G4double hallhheight = 850.0;  // Hall cylinder height 17 m
  G4double offset = 200.0;  // shift cavern floor to keep detector centre at origin
  
  //Water tank with water and air buffer
  G4double water_tank_thickness = 7.0;
  G4double water_tank_height = 8900.0;
  G4double inner_tank_height = water_tank_height - 2 * water_tank_thickness;
  G4double inner_radius = 0.0;
  G4double water_radius = 5000.0;
  G4double water_height = inner_tank_height;
  
  //# Reflective foil
  G4double reflective_foil_thickness = 0.04;

  //cryostat
  G4double cryo_radius = 3976 / 2;
  G4double cryo_wall = 12;
  G4double cryo_tub_height = 3900;
  G4double cryo_top_height = 826;
  G4double cryo_bottom_height = 829;

  G4double cryo_access_radius = 800 / 2;
  G4double cryo_access_wall = 10;
  G4double cryo_access_height = 1720;
  G4double access_overlap = 200;

  //Pillbox
  G4double shielding_foot_or = 2000.0;
  G4double shielding_foot_thickness = 1.2;
  G4double shielding_foot_ir = shielding_foot_or - shielding_foot_thickness - reflective_foil_thickness;
  G4double pillbox_cryo_bottom_height = (inner_tank_height / 2) - (cryo_tub_height / 2) - cryo_bottom_height - reflective_foil_thickness - access_overlap + cryo_wall+ cryo_access_wall;
  G4double pillbox_offset = -water_height / 2 + 0.5 * pillbox_cryo_bottom_height;
  G4double pillbox_tube_foil_offset = pillbox_offset;

  //Manhole
  G4double manhole_outer_radius = 400.0;
  G4double manhole_inner_radius = 0 ; // # No inner radius for the manhole
  G4double manhole_height = 2. * (shielding_foot_or + reflective_foil_thickness);
  G4double manhole_angle =  CLHEP::pi ;// # Half-circle (180 degrees)
  G4double manhole_offset = 0.5 * pillbox_cryo_bottom_height - manhole_outer_radius;

  //# Air buffer
  G4double outer_water_tank_radius = water_radius + water_tank_thickness;
  G4double air_buffer_radius = water_radius - reflective_foil_thickness;
  G4double air_buffer_height = 486.0;
  
  
  //# z-axis Offsets
  G4double air_buffer_offset = 0.5 * (inner_tank_height - air_buffer_height);
  G4double tank_offset = 0.0;
  G4double bottom_foil_offset = -0.5 * water_height + 0.5 * reflective_foil_thickness;
  G4double cryo_z_displacement =-153.0; // # (innertank_height/2-cryo_acess_height-cryo_top_height-access_overlap/2  // inner_tank_height/2.
  
  // PMT
  G4double PMTrad     = 10.2;   // diameter 8 inch
  G4double PMTheight  = 2.0;    //random value
 
  G4double LambdaE = 1239.84193 * nanometer * eV;

  std::vector<G4double> to_E_in_eV(std::vector<G4double> E_in_nm){
    G4double e_in_eV=0;
    std::vector<G4double> E_in_eV;
    for (auto e_in_nm: E_in_nm){
      e_in_eV = LambdaE/ e_in_nm /eV;
      E_in_eV.push_back(e_in_eV);
    }
    return E_in_eV;
  }
  G4double to_e_in_eV(G4double e_in_nm){
      G4double e_in_eV = LambdaE/ e_in_nm /eV;
    return e_in_eV;
  }

  G4double vm2000_calculate_wls_mfp(G4double yield_value){
    //# Set total path length (currently hardcoded to 1 mm)
    G4double total_path = 1.0e-3 * m;

    //# Handle edge cases
    if (yield_value == 0){
        return 10.0 * m; } // # Large mean free path, no absorption
    if (yield_value == 1){
        return 0.01e-3 * m;}//  # 0.01 mm, Very small mean free path, 100% absorption

    //# Calculate mean free path for valid yield values
    G4double help_value = log(1.0 - yield_value);
    return -total_path / help_value;
  }

  std::vector<double> EmissionSpectrum(std::string volume, std::vector<double> energy){
    G4String nameFile;
    if (volume=="WLS_abs")
      nameFile = "TPBAbsorption.dat";
    else if (volume=="WLS_em")
      nameFile = "VM2000_em_spec.dat";
    else if (volume=="WLS_ref")
      nameFile = "Reflectivity_VM200.dat";
    else 
      G4cout << "Wrong ...";

    G4String pathString("../materials/");
    G4String pathFile = pathString + nameFile;

    auto Graph = new TGraph(pathFile.data());

    std::vector <double> table;
 
    
    for (auto e : energy){
      G4cout << e << "energy " <<LambdaE/e/nanometer << G4endl;
      auto em = Graph->Eval(LambdaE/e/nanometer);
      G4cout << volume << " " << em << G4endl;
      table.push_back(em >= 0 ? em : 0);
    }

    return table;
  }

 private:
  void PrintError(G4String);

  SNneutrinosDetectorMessenger* fDetectorMessenger;

  G4bool fVerbose;

  G4Material *waterMat;
  G4MaterialPropertiesTable *waterMPT;


  G4Material *steelMat;
  G4MaterialPropertiesTable *steelMPT;

  G4Material *PMTMat;
  G4MaterialPropertiesTable *CsMPT;

  G4MaterialPropertiesTable *reflectorMPT;
  G4MaterialPropertiesTable *borderMPT;
  
  G4Material *foilMat;
  G4MaterialPropertiesTable *foilMPT;
  
  G4Material *worldMat;
  G4MaterialPropertiesTable *worldMPT;

  G4Material *airMat;
  G4MaterialPropertiesTable *airMPT;

};




#endif
