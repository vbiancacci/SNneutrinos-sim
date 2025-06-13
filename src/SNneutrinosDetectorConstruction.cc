#include "SNneutrinosDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4NistManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include <numeric>


SNneutrinosDetectorConstruction::SNneutrinosDetectorConstruction()
  : G4VUserDetectorConstruction()
{
  fVerbose          = false;
  // create a messenger for this class
  
  waterMPT = new G4MaterialPropertiesTable();
  steelMPT = new G4MaterialPropertiesTable();
  CsMPT = new G4MaterialPropertiesTable();
  reflectorMPT = new G4MaterialPropertiesTable();
  borderMPT = new G4MaterialPropertiesTable();
  foilMPT = new G4MaterialPropertiesTable();
  worldMPT = new G4MaterialPropertiesTable();

}

SNneutrinosDetectorConstruction::~SNneutrinosDetectorConstruction()
{
}


G4VPhysicalVolume* SNneutrinosDetectorConstruction::Construct()
{
  G4bool checkOverlaps = true;
  // ------------- Materials -------------

  G4NistManager *nist = G4NistManager::Instance();
  worldMat = nist->FindOrBuildMaterial("G4_Galactic");
  steelMat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  waterMat = nist->FindOrBuildMaterial("G4_WATER");
  airMat   = nist->FindOrBuildMaterial("G4_AIR");
  PMTMat = nist->FindOrBuildMaterial("G4_Al"); //CathodeMetalAluminium
  foilMat = new G4Material("VM2000", 1.15*g/cm3, 4);
  foilMat->AddElement(nist->FindOrBuildElement("C"),13);
  foilMat->AddElement(nist->FindOrBuildElement("H"),2);
  foilMat->AddElement(nist->FindOrBuildElement("N"),2);
  foilMat->AddElement(nist->FindOrBuildElement("O"),3);



  // https://www.mpi-hd.mpg.de/gerda/public/2008/c08_ndip08_MuonVeto_mk.pdf

  //
  
  // ------------ Generate & Add Material Properties Table ------------
  //
    
  std::vector<G4double> E_in_nm = {600. *nm, 550.*nm, 500.*nm, 450.*nm,  400.*nm, 350.*nm,  300.*nm,  250.*nm,  200.*nm,  150.*nm,  100.* nm};
  std::vector<G4double> E_in_eV= to_E_in_eV(E_in_nm);
  G4int n_energy = size(E_in_eV);
  for (G4double i : E_in_eV){
    G4cout << i << G4endl;
  }
  std::vector<G4double> waterAbsorption = {
    10 * 1000 * mm, // # 10 m
    20 * 1000* mm,  //# 20 m
    50 * 1000* mm,  //# 50 m
    100 * 1000* mm,  //# 100 m
    100 * 1000* mm,  //# 100 m
    100 * 1000* mm,  ///# 100 m
    90 * 1000* mm,  //# 90 m
    20 * 1000* mm,  //# 20 m
    1 * 1000* mm,  //# 1 m
    0.001* mm,  //# 0.001 mm
    0.0001* mm,  //# 0.0001 mm
  };

  std::vector<G4double> waterRIndex (n_energy, 1.33);

  std::vector<G4double>  steelRIndex  (n_energy, 2.86);
  std::vector<G4double>  steelAbsorption (n_energy, 1.e-20*m);
    
  //water
  waterMPT->AddProperty("RINDEX", E_in_eV, waterRIndex); //, false, true);
  waterMPT->AddProperty("ABSLENGTH", E_in_eV, waterAbsorption); //, false, true);
  G4cout << "Water G4MaterialPropertiesTable:" << G4endl;
  waterMPT->DumpTable();
  waterMat->SetMaterialPropertiesTable(waterMPT);
  //GdwaterMat->SetMaterialPropertiesTable(waterMPT);

  //steel
  steelMPT->AddProperty("RINDEX", E_in_eV, steelRIndex);//, false, true);
  steelMPT->AddProperty("ABSLENGTH", E_in_eV, steelAbsorption);//, false, true); //if uncommented there will be some absorption even in the cryostat/water tank (????)
  steelMat->SetMaterialPropertiesTable(steelMPT);
  G4cout << "Steel G4MaterialPropertiesTable:" << G4endl;
  steelMPT->DumpTable();
  


  //foil
  //# VM2000 seem to consist of PMMA and PEN layers https://iopscience.iop.org/article/10.1088/1748-0221/12/06/P06017/pdf
  //legendoptics.pen.pyg4_pen_attach_scintillation(self._vm2000, self.g4_registry)
  //legendoptics.vm2000.pyg4_vm2000_attach_particle_scintillationyields(self._vm2000, self.g4_registry)


  //"""Wavelength-shifting parameters for the reflective foil VM2000."""

  G4double wls_yield = 0.075; //  # 0.6 MaGe, 0.075 XENON paper

  //# Populate VM2000_energy_range array with energy values
  G4double ppsci_high_e = to_e_in_eV(115 * nm); 
  G4double ppsci_low_e = to_e_in_eV(650 * nm);


  int num1 = 251;
  G4double dee = (ppsci_high_e - ppsci_low_e) / (num1 - 2);
  std::vector<G4double> vm2000_energy_range (num1,0.*eV);
  vm2000_energy_range[0] = 1.8 * eV;

  std::vector<int> nn (num1-1); 
  std::iota(nn.begin(), nn.end(), 1);
  for (int i : nn){
      vm2000_energy_range.at(i) = ppsci_low_e + i * dee;
  }   
  //# Create arrays for energy and optical properties
  std::vector<G4double> vm2000_reflectivity(num1, 0.);
  std::vector<double>  vm2000_efficiency (num1, 0.);
  std::vector<double> wls_absorption (num1, 0.*m);
  std::vector<double> wls_emission (num1, 0.*m);

  //# Set reflectivity, absorption, and emission
  std::vector<int> nn_ (num1-1); 
  std::iota(nn_.begin(), nn_.end(), 0);
  for (int i :  nn_){
      if (vm2000_energy_range[i] < to_e_in_eV(370 * nm)){
          vm2000_reflectivity[i] = 0.95;} // # Visible light 0.95, 0.99
      else{
          vm2000_reflectivity[i] = 0.12;} //  # UV light 0.15, 0.3 (paper)

      if (vm2000_energy_range[i] > 3.35 * eV){
          //# depending on path length in foil --> angle
          wls_absorption[i] = vm2000_calculate_wls_mfp(wls_yield);}  //# Absorbs UV
      else{
          wls_absorption[i] = (1.0 * m);} //  # Imperturbed, no absorption of visible light
  }
   //update emission data
  wls_emission=EmissionSpectrum("WLS_em",vm2000_energy_range);

  //# Copy the first element to 0th position
  wls_absorption[0] = wls_absorption[1]; // # depending on path length in foil --> angle
  wls_emission[0] = wls_emission[1];

  std::vector<G4double> FiberRIndex (n_energy, 1.15);
  std::vector<G4double> FiberABS (n_energy, 50*m); //*m  
  foilMPT->AddProperty("RINDEX", E_in_eV, FiberRIndex);//, false, true);
  foilMPT->AddProperty("ABSLENGTH", E_in_eV, FiberABS);//, false, true);
  G4cout << "foil G4MaterialPropertiesTable:" << G4endl;
  foilMPT->DumpTable();

  foilMPT->AddProperty("WLSABSLENGTH",vm2000_energy_range, wls_absorption);//, false, true);
  foilMPT->AddProperty("WLSCOMPONENT",vm2000_energy_range, wls_emission);//, false, true);
  G4double TimeConstant = 0.5 *ns;
  foilMPT->AddConstProperty("WLSTIMECONSTANT", TimeConstant);
  foilMat->SetMaterialPropertiesTable(foilMPT);

  //foilMPT->AddProperty("SCINTILLATIONCOMPONENT1", λ_scint.to("eV"), scint_em);
  //foilMPT->AddProperty("SCINTILLATIONTIMECONSTANT1", pen_scint_timeconstant());
  //foilMPT->AddProperty("RESOLUTIONSCALE", 1);


  reflectorMPT->AddProperty("REFLECTIVITY", vm2000_energy_range, vm2000_reflectivity);//, false, true);
  reflectorMPT->AddProperty("EFFICIENCY", vm2000_energy_range, vm2000_efficiency);//, false, true);

  std::vector<G4double> vm2000_reflectivity_border (num1, 0.);
  std::vector<G4double>  vm2000_transmittance_border (num1, 1.);;
  std::vector<G4double>  vm2000_efficiency_border (num1, 0.);;  

  borderMPT->AddProperty("REFLECTIVITY", vm2000_energy_range, vm2000_reflectivity_border);//, false, true);
  borderMPT->AddProperty("EFFICIENCY", vm2000_energy_range, vm2000_efficiency_border);//, false, true);
  borderMPT->AddProperty("TRANSMITTANCE", vm2000_energy_range, vm2000_transmittance_border);//, false, true);
 

  //vacuum
  worldMPT->AddProperty("ABSLENGTH", E_in_eV, steelAbsorption);//, false, true);
  worldMat->SetMaterialPropertiesTable(worldMPT);

  foilMat->SetMaterialPropertiesTable(foilMPT);



 //------------- Volumes --------------


  // World (hall)
  auto* worldSolid = new G4Tubs("World", 0.0 * cm, (hallrad + stone + 0.1) * cm, (hallhheight + stone + offset + 0.1) * cm, 0.0, CLHEP::twopi);
  auto* fWorldLogical  = new G4LogicalVolume(worldSolid, worldMat, "World_log");
  auto* fWorldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLogical,"World_phys", nullptr, false, 0);

  //Tank
  auto* TankSolid = new G4Tubs("Tank", 0.,outer_water_tank_radius, inner_tank_height/2., 0, CLHEP::twopi);
  auto* fTankLogical = new G4LogicalVolume(TankSolid, steelMat, "Tank_log");
  auto* fTankPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fTankLogical, "Tank_phys", fWorldLogical, false, 0, true);
  
  //Water
  auto* WaterSolid = new G4Tubs("Water", 0.,water_radius, water_height/2., 0, CLHEP::twopi);
  auto* fWaterLogical  = new G4LogicalVolume(WaterSolid, waterMat, "Water_log"); // waterMat 
  auto* fWaterPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWaterLogical, "Water_phys", fTankLogical, false, 0, true);

  //AirBuffer
  auto* AirBufferSolid = new G4Tubs("AirBuffer", cryo_access_radius + cryo_access_wall, air_buffer_radius, air_buffer_height/2., 0, CLHEP::twopi);
  auto* fAirBufferLogical  = new G4LogicalVolume(AirBufferSolid, airMat, "AirBuffer_log"); 
  auto* fAirBufferPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,air_buffer_offset), fAirBufferLogical, "AirBuffer_phys", fWaterLogical, false, 0, true);
  
  //PillBox
    G4RotationMatrix* Rot = new G4RotationMatrix; 
    Rot->rotateZ(-CLHEP::pi / 2.0*rad);
    Rot->rotateY(0);
    Rot->rotateX(-CLHEP::pi /2.0*rad);
    auto* pillbox_tube = new G4Tubs("pillbox_tube", shielding_foot_ir, shielding_foot_or, pillbox_cryo_bottom_height/2., 0, CLHEP::twopi ); //  # outer steel cylinder
    //# Create the half-tube (semi-cylinder) for the manhole
    auto* manhole_pillbox_arc = new G4Tubs("manhole_pillbox", manhole_inner_radius, manhole_outer_radius, manhole_height/2., 0, manhole_angle);
    auto* manholepillbox_box = new G4Box("manholepillbox_box", manhole_outer_radius, manhole_outer_radius/2., manhole_height/2.);
    auto* manhole_pillbox = new G4UnionSolid("manhole_union", manhole_pillbox_arc, manholepillbox_box, 0, G4ThreeVector(0, -0.5 * manhole_outer_radius));
    //# Subtract the first manhole (half-tube) from the pillbox
    auto* PillboxSolid = new G4SubtractionSolid("pillbox_subtraction1", pillbox_tube, manhole_pillbox, Rot, G4ThreeVector(0, 0, 0 - manhole_offset));
    auto* fPillboxLogical = new G4LogicalVolume(PillboxSolid, steelMat, "Pillbox_log");
    auto* fPillboxPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0, pillbox_offset), fPillboxLogical, "Pillbox_phys", fWaterLogical, false, 0, true);

    /*
//Cryostat
  G4double cryostat_r [] = {
    0, 
    cryo_access_radius+cryo_access_wall,
    cryo_access_radius+cryo_access_wall,
    cryo_radius+cryo_wall,
    cryo_radius+cryo_wall,
    0};
  G4double cryostat_z [] = {
    0,
    0,
    -(cryo_access_height+access_overlap),
    -(cryo_access_height+access_overlap),
    -(cryo_access_height+access_overlap+cryo_tub_height+cryo_top_height+cryo_bottom_height),
    -(cryo_access_height+access_overlap+cryo_tub_height+cryo_top_height+cryo_bottom_height)
  } ;

  G4double zeros[6]={0.};
  auto CryostatSolid = new G4Polycone("Cryostat", 0, CLHEP::twopi, 6, cryostat_z, zeros, cryostat_r);
  auto* fCryostatLogical = new G4LogicalVolume(CryostatSolid, steelMat, "Cryostat_log");
  auto* fCryostatPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0, cryo_z_displacement), fCryostatLogical, "Cryostat_phys", fWaterLogical, false, 0, true);
  */

//cryostat from pygeom
auto* cryo_top = new G4Tubs(
  "cryo_top",
  0,
  cryo_radius + cryo_wall,
  cryo_top_height + cryo_wall,
  0, CLHEP::twopi );

auto* cryo_access_tub = new G4Tubs(
  "cryo_access_tub",
  0,
  cryo_access_radius + cryo_access_wall,
  (cryo_access_height + access_overlap)/2.,
  0, CLHEP::twopi );

auto* cryo_bottom = new G4Tubs(
  "cryo_bottom",
  0,
  cryo_radius + cryo_wall,
  cryo_bottom_height + cryo_wall,
  0, CLHEP::twopi );

auto* cryo_tub = new G4Tubs("cryo_tub", 0, cryo_radius + cryo_wall, cryo_tub_height/2.,  0, CLHEP::twopi );

auto* cryo1 = new G4UnionSolid("cryo1", cryo_tub, cryo_top, 0, G4ThreeVector(0, 0, cryo_tub_height / 2));
auto* cryo2 = new G4UnionSolid("cryo2", cryo1, cryo_bottom, 0, G4ThreeVector(0, 0, -cryo_tub_height / 2));
auto* CryostatSolid = new G4UnionSolid(
  "cryostat",
  cryo2,
  cryo_access_tub,
  0,
  G4ThreeVector(0, 0, +cryo_tub_height / 2 + cryo_top_height + cryo_access_height / 2));
  auto* fCryostatLogical = new G4LogicalVolume(CryostatSolid, steelMat, "Cryostat_log");
  auto* fCryostatPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0, cryo_z_displacement), fCryostatLogical, "Cryostat_phys", fWaterLogical, false, 0, true);




//Foil water wall tube
  auto* ReflectionFoilWaterTankTubeSolid = new G4Tubs(
    "ReflectionFoil_WaterTankTube",
    water_radius - reflective_foil_thickness,
    water_radius,
    water_height/2.,
    0, CLHEP::twopi );
  auto* fReflectionFoilWaterTankTubeLogical = new G4LogicalVolume(ReflectionFoilWaterTankTubeSolid, foilMat, "ReflectionFoil_WaterTankTube_log");
  auto* fReflectionFoilWaterTankTubePhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fReflectionFoilWaterTankTubeLogical, "ReflectionFoil_WaterTankTube_phys", fWaterLogical, false, 0, true);

// Foil water bottom
  auto* ReflectionFoilWaterTankBottomSolid = new G4Tubs(
  "ReflectionFoil_WaterTankBottom",
  shielding_foot_or + reflective_foil_thickness,
  water_radius - reflective_foil_thickness,
  reflective_foil_thickness/2.,
  0, CLHEP::twopi );
  auto* fReflectionFoilWaterTankBottomLogical = new G4LogicalVolume(ReflectionFoilWaterTankBottomSolid, foilMat, "ReflectionFoil_WaterTankBottom_log");
  auto* fReflectionFoilWaterTankBottomPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0, bottom_foil_offset), fReflectionFoilWaterTankBottomLogical, "ReflectionFoil_WaterTankBottom_phys", fWaterLogical, false, 0, true);

///////////////////////////////////////////////////////////////////////////////////////////////////

// Foil water bottom
/*auto* ReflectionFoilWaterTankTopSolid = new G4Tubs(
  "ReflectionFoil_WaterTankTop",
  cryo_access_radius + cryo_access_wall,
  air_buffer_radius,
  reflective_foil_thickness/2.,
  0, CLHEP::twopi );
  //auto* fReflectionFoilWaterTankTopLogical = new G4LogicalVolume(ReflectionFoilWaterTankTopSolid, foilMat, "ReflectionFoil_WaterTankTop_log");
  //auto* fReflectionFoilWaterTankTopPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0, water_height/2. - air_buffer_height-reflective_foil_thickness/2.), fReflectionFoilWaterTankTopLogical, "ReflectionFoil_WaterTankTop_phys", fWaterLogical, false, 0, true);
*/
///////////////////////////////////////////////////////////////////////////////////////////////////

//Foil pillbox outer
  auto* pillbox_outer_reflection_foil_tube_subtraction1 = new G4Tubs(
  "pillbox_outer_reflection_foil_tube_subtraction1",
  shielding_foot_or,
  shielding_foot_or + reflective_foil_thickness,
  pillbox_cryo_bottom_height/2.,
  0, CLHEP::twopi );

  auto* ReflectionFoilPillboxOuterSolid = new G4SubtractionSolid(
  "pillbox_outer_reflection_foil_tube ",
  pillbox_outer_reflection_foil_tube_subtraction1,
  manhole_pillbox,Rot, G4ThreeVector(0, 0, 0 - manhole_offset));  
  auto* fReflectionFoilPillboxOuterLogical = new G4LogicalVolume(ReflectionFoilPillboxOuterSolid, foilMat, "ReflectionFoil_PillboxOuter_log");
  auto* fReflectionFoilPillboxOuterPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,pillbox_offset), fReflectionFoilPillboxOuterLogical, "ReflectionFoil_PillboxOuter_log_phys",  fWaterLogical, false, 0, true);

//Foil pillbox inner
  auto* pillbox_inner_reflection_foil_tube_subtraction1 = new G4Tubs(
  "pillbox_inner_reflection_foil_tube_subtraction1",
  shielding_foot_ir - reflective_foil_thickness,
  shielding_foot_ir,
  (pillbox_cryo_bottom_height - shielding_foot_thickness)/2.,
  0, CLHEP::twopi );

  auto* ReflectionFoilPillboxInnerSolid = new G4SubtractionSolid(
  "pillbox_inner_reflection_foil_tube ",
  pillbox_inner_reflection_foil_tube_subtraction1,
  manhole_pillbox, Rot, G4ThreeVector(0, 0, 0 - manhole_offset));
  auto* fReflectionFoilPillboxInnerLogical = new G4LogicalVolume(ReflectionFoilPillboxInnerSolid, foilMat, "ReflectionFoil_PillboxInner_log");
  auto* fReflectionFoilPillboxInnerPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,pillbox_offset),  fReflectionFoilPillboxInnerLogical, "ReflectionFoil_PillboxInner_phys", fWaterLogical, false, 0, true);

//Foil pillbox top
  auto* ReflectionFoilPillboxTopSolid = new G4Tubs(
  "pillbox_reflection_foil_top",
  inner_radius,
  shielding_foot_ir,
  reflective_foil_thickness/2.,
  0, CLHEP::twopi );
  auto* fReflectionFoilPillboxTopLogical = new G4LogicalVolume(ReflectionFoilPillboxTopSolid, foilMat, "ReflectionFoil_PillboxTop_log");
  auto* fReflectionFoilPillboxTopPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0, bottom_foil_offset + pillbox_cryo_bottom_height - reflective_foil_thickness),  fReflectionFoilPillboxTopLogical, "ReflectionFoil_PillboxTop_phys", fWaterLogical, false, 0, true);

//Foil pillbox bottom
  auto* ReflectionFoilPillboxBottomSolid = new G4Tubs(
  "pillbox_reflection_foil_bottom",
  inner_radius,
  shielding_foot_ir,
  reflective_foil_thickness/2.,
  0, CLHEP::twopi );
  auto* fReflectionFoilPillboxBottomLogical = new G4LogicalVolume(ReflectionFoilPillboxBottomSolid, foilMat, "ReflectionFoil_PillboxBottom_log");
  auto* fReflectionFoilPillboxBottomPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,bottom_foil_offset),  fReflectionFoilPillboxBottomLogical, "ReflectionFoil_PillboxBottom_log", fWaterLogical, false, 0, true);

//Foil cryostat
auto* ReflectionFoilCryostatSolid = new G4Tubs(
  "cryo_reflection_foil",
  cryo_radius + cryo_wall,
  cryo_radius + cryo_wall + reflective_foil_thickness,
  (cryo_tub_height + cryo_top_height + cryo_bottom_height + 2 * cryo_wall)/2.,
  0, CLHEP::twopi );
  auto* fReflectionFoilCryostatLogical = new G4LogicalVolume(ReflectionFoilCryostatSolid, foilMat, "ReflectionFoil_Cryostat_log");
  auto* fReflectionFoilCryostatPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,cryo_z_displacement),  fReflectionFoilCryostatLogical, "ReflectionFoil_Cryostat_phys", fWaterLogical, false, 0, true);



// ------------- Surfaces --------------

  //Reflective surface for VM2000
  G4OpticalSurface* VM2000_opSurface = new G4OpticalSurface("VM2000_opSurface",
    unified,  
    polished,          
    dielectric_metal,
    0.3
);

VM2000_opSurface->SetMaterialPropertiesTable(reflectorMPT);

G4OpticalSurface* WaterToVM2000_opSurface = new G4OpticalSurface("WaterToVM2000_opSurface",
  unified,  
  polished,          
  dielectric_metal,
  0.3
);
WaterToVM2000_opSurface->SetMaterialPropertiesTable(borderMPT);


//Border Surface

  auto* ReflectionFoilPillboxOuterSurface = new G4LogicalBorderSurface("ReflectionFoilPillboxOuterSurface", fReflectionFoilPillboxOuterPhysical,fWaterPhysical, WaterToVM2000_opSurface);

  auto* ReflectionFoilPillboxInnerSurface = new G4LogicalBorderSurface("ReflectionFoilPillboxOuterSurface",fReflectionFoilPillboxInnerPhysical,fWaterPhysical, WaterToVM2000_opSurface);

  auto* ReflectionFoilPillboxBottomSurface = new G4LogicalBorderSurface("ReflectionFoilPillboxBottomSurface",fReflectionFoilPillboxBottomPhysical,fWaterPhysical, WaterToVM2000_opSurface);

  auto* ReflectionFoilPillboxTopSurface = new G4LogicalBorderSurface("ReflectionFoilPillboxTopSurface",fReflectionFoilPillboxTopPhysical,fWaterPhysical, WaterToVM2000_opSurface);

  auto* ReflectionFoilWaterTankTubeSurface = new G4LogicalBorderSurface("ReflectionFoilWaterTankTubeSurface",fReflectionFoilWaterTankTubePhysical,fWaterPhysical, WaterToVM2000_opSurface);

  auto* ReflectionFoilWaterTankBottomSurface = new G4LogicalBorderSurface("ReflectionFoilWaterTankBottomSurface",fReflectionFoilWaterTankBottomPhysical,fWaterPhysical, WaterToVM2000_opSurface);
  //auto* ReflectionFoilWaterTankTopSurface = new G4LogicalBorderSurface("ReflectionFoilWaterTankTopSurface",fReflectionFoilWaterTankTopPhysical,fWaterPhysical, WaterToVM2000_opSurface);

  auto* ReflectionFoilCryostatSurface = new G4LogicalBorderSurface("ReflectionFoilCryostatSurface",fReflectionFoilCryostatPhysical,fWaterPhysical, WaterToVM2000_opSurface);


//Skin Surface

  auto* ReflectionFoilPillboxOuterSkin = new  G4LogicalSkinSurface("ReflectionFoilPillboxOuterSkin", fReflectionFoilPillboxOuterLogical, VM2000_opSurface);

  auto* ReflectionFoilPillboxInnerSkin = new  G4LogicalSkinSurface("ReflectionFoilPillboxInnerSkin", fReflectionFoilPillboxInnerLogical, VM2000_opSurface);

  auto* ReflectionFoilPillboxBottomSkin = new  G4LogicalSkinSurface("ReflectionFoilPillboxBottomSkin", fReflectionFoilPillboxBottomLogical, VM2000_opSurface);

  auto* ReflectionFoilPillboxTopSkin = new  G4LogicalSkinSurface("ReflectionFoilPillboxTopSkin", fReflectionFoilPillboxTopLogical, VM2000_opSurface);

  auto* ReflectionFoilWaterTankTubeSkin = new  G4LogicalSkinSurface("ReflectionFoilWaterTankTubeSkin", fReflectionFoilWaterTankTubeLogical, VM2000_opSurface);

  auto* ReflectionFoilWaterTankBottomSkin = new  G4LogicalSkinSurface("ReflectionFoilWaterTankBottomSkin", fReflectionFoilWaterTankBottomLogical, VM2000_opSurface);
  //auto* ReflectionFoilWaterTankTopSkin = new  G4LogicalSkinSurface("ReflectionFoilWaterTankTopSkin", fReflectionFoilWaterTankTopLogical, VM2000_opSurface);

  auto* ReflectionFoilCryostatSkin = new  G4LogicalSkinSurface("ReflectionFoilCryostatSkin", fReflectionFoilCryostatLogical, VM2000_opSurface);




  //
  // PMT
  //
  
  G4double pos_x = 0.;
  G4double pos_y = 0.;
  G4double pos_z = 0.;
  G4double pos_r = 0.;
  G4double pos_theta = 0.;

  G4int PMT_ID = 0.;

  auto* PMTSolid =
    new G4Tubs("PMT", 0.0 * cm, PMTrad * cm, PMTheight *  cm, 0.0, CLHEP::twopi);
  auto* fPMTLogical  = new G4LogicalVolume(PMTSolid, PMTMat, "PMT_log");
 
  //Bottom PMT
 G4int bottom_circles = 3; 
 G4int n_bottom_PMT=0;
 G4int n_PMT_ring [3] = {12,8,14}; //in the first ring, there are only 10 PMTs, two spot are left empty  //{12,5,8};only currently not working PMT (6 in first ring)
 for (int j=1; j<=bottom_circles; j++){ 
    G4int PMT_bottom_circle = G4int(n_PMT_ring[j-1]/(bottom_circles - j +1));
    for (int i=1;i<=n_PMT_ring[j-1]; i++) {
      if (j==1 && (i==6 || i==12 )){
        continue;
      }
      else{
        pos_r = G4double((water_radius+500.)  * j /(bottom_circles+1.)) ;
        pos_theta = G4double(CLHEP::twopi*i/n_PMT_ring[j-1]);
        pos_x = pos_r * cos(pos_theta);
        pos_y = pos_r * sin(pos_theta);
        pos_z = -water_height/2.+10.*PMTheight+reflective_foil_thickness;
        PMT_ID = 8000+j*100+i;  //ID starting with 8 is set to the bottom and the second number indicates the ring - number 1 is the most internal ring
        new G4PVPlacement(nullptr, G4ThreeVector(pos_x, pos_y, pos_z), fPMTLogical, "PMT_phys", fWaterLogical, false, PMT_ID);
        n_bottom_PMT++;
        G4cout << PMT_ID <<  pos_r  <<" r " << " pos_x " << pos_x << " pos_y " << pos_y << " pos_z " << pos_z << G4endl;
      }
      }
  } 
  G4cout <<"bottom PMT  " << n_bottom_PMT << G4endl;

  G4int n_PMT=0;
  G4int rings = 4;       
  G4int PMT_lateral [4] ={10, 10, 10, 1}; //{3,4,4,1}; only currently working PMTs
  G4double ring_lateral_pos [4] = {2000, 3500, 5000, 6500};
  for (int j=1; j<=rings; j++) { 
    for (int i=0;i<PMT_lateral[j-1]; i++) {
      if (j==2){
        pos_theta = double(CLHEP::twopi)*(i+1/2.)/PMT_lateral[j-1];
      }
      else{
        pos_theta = double(CLHEP::twopi)*i/PMT_lateral[j-1];
      }
      G4RotateY3D rotTheta; //(0);
      G4int quad=i/25;
      G4double pos_theta_mod =pos_theta-quad*CLHEP::pi/2.;
      if (quad==0 || quad==2){
        rotTheta= G4RotateY3D(CLHEP::pi/2.+pos_theta_mod);}
      else{
        rotTheta= G4RotateY3D(pos_theta_mod);}
      G4RotateX3D rotPhi(CLHEP::pi/2.);
      pos_r = water_radius-10.*PMTheight-reflective_foil_thickness;
      pos_x = pos_r * cos(pos_theta);
      pos_y = pos_r * sin(pos_theta);
      G4Translate3D shift(pos_x, pos_y, -water_height/2+ ring_lateral_pos[j-1]);
      auto transform = shift*rotPhi*rotTheta; 
      PMT_ID = (rings-j+1)*100+i;   //ID starting with 1 belogs to the first ring from the top
      new G4PVPlacement(transform, fPMTLogical, "PMT_phys", fWaterLogical, false, PMT_ID);
      n_PMT+=1;
    } 
  }
  G4cout <<"lateral PMT  " << n_PMT << G4endl;


  auto* water_only_solid1 = new G4SubtractionSolid("water_only1", WaterSolid, CryostatSolid, 0, G4ThreeVector(0,0, cryo_z_displacement));
  auto* water_only_solid2 = new G4SubtractionSolid("water_only2", water_only_solid1, AirBufferSolid, 0, G4ThreeVector(0,0, air_buffer_offset));
  auto* water_only_solid3 = new G4SubtractionSolid("water_only3", water_only_solid2, PillboxSolid, 0, G4ThreeVector(0,0, pillbox_offset));
  auto* water_only_solid4 = new G4SubtractionSolid("water_only4", water_only_solid3, ReflectionFoilWaterTankTubeSolid, 0, G4ThreeVector(0,0, 0));
  auto* water_only_solid5 = new G4SubtractionSolid("water_only5", water_only_solid4, ReflectionFoilWaterTankBottomSolid, 0, G4ThreeVector(0,0, bottom_foil_offset));
  auto* water_only_solid6 = new G4SubtractionSolid("water_only6", water_only_solid5, ReflectionFoilPillboxOuterSolid, 0, G4ThreeVector(0,0, pillbox_offset));
  auto* water_only_solid7 = new G4SubtractionSolid("water_only7", water_only_solid6, ReflectionFoilPillboxInnerSolid, 0, G4ThreeVector(0,0, pillbox_offset));
  auto* water_only_solid8 = new G4SubtractionSolid("water_only8", water_only_solid7, ReflectionFoilPillboxTopSolid, 0, G4ThreeVector(0,0, bottom_foil_offset + pillbox_cryo_bottom_height - reflective_foil_thickness));
  auto* water_only_solid9 = new G4SubtractionSolid("water_only9", water_only_solid8, ReflectionFoilPillboxBottomSolid, 0, G4ThreeVector(0,0, bottom_foil_offset));
  auto* water_only_solid10 = new G4SubtractionSolid("water_only10", water_only_solid9, ReflectionFoilCryostatSolid, 0, G4ThreeVector(0,0, cryo_z_displacement));
  
  auto* fwater_only_logical  = new G4LogicalVolume(water_only_solid10, waterMat, "Water_only_log");
  //auto* fWaterOnlyPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), fwater_only_logical, "Water_only_phys", fTankLogical, false, 0, true);

  

  /* 
  //PMT photocathode aluminium
  std::vector <double> PMT_emission;
  for (int i==0; i< size(photonEnergy); i++){
    PMT_emission.push_back(EmissionSpectrum("PMT",photonEnergy[i]));
  }
  
  const G4int CathMetal_NUMENTRIES = 5;
  G4double CathodeMetalAluminium_PP[CathMetal_NUMENTRIES]   = {1.8*eV, 2.0*eV , 3.06*eV, 6.2*eV, 10.5*eV };
  G4double CathodeMetalAluminium_RIND[CathMetal_NUMENTRIES] = { 1.74, 1.3, 0.523, 0.119, 0.040 };     // refraction index 
  G4double CathodeMetalAluminium_ABSL[CathMetal_NUMENTRIES] = { 1.e-30*m, 1.e-30*m, 1.e-30*m,  1.e-30*m,  1.e-30*m };// attenuation length
  G4MaterialPropertiesTable *CathodeMetalAluminium_mt = new G4MaterialPropertiesTable();
  CathodeMetalAluminium_mt->AddProperty("RINDEX", CathodeMetalAluminium_PP, CathodeMetalAluminium_RIND,CathMetal_NUMENTRIES);
  CathodeMetalAluminium_mt->AddProperty("ABSLENGTH", CathodeMetalAluminium_PP, CathodeMetalAluminium_ABSL, CathMetal_NUMENTRIES);
  CathodeMetalAluminium_mt->AddProperty("EFFICIENCY", CathodeMetalAluminium_PP, PMT_emission, CathMetal_NUMENTRIES);
  //PMTMat->SetMaterialPropertiesTable(CathodeMetalAluminium_mt);
  */


 

  //
  // Visualization
  //
  G4Color testColor_water(170 / 255., 191 / 225., 219 / 225.); //170, 191, 219
  auto*   testVisAtt_water = new G4VisAttributes(testColor_water);
  testVisAtt_water->SetVisibility(true);
  auto* greyVisAtt = new G4VisAttributes(G4Colour::Grey());
  auto* cyanVisAtt = new G4VisAttributes(G4Colour::Cyan());
  auto* redVisAtt = new G4VisAttributes(G4Colour::Red());
  auto* greenVisAtt = new G4VisAttributes(G4Colour::Green());
  greyVisAtt->SetVisibility(true);

  fTankLogical->SetVisAttributes(greyVisAtt);
  //fWaterLogical->SetVisAttributes(redVisAtt);
  //fCryostatLogical->SetVisAttributes(redVisAtt);
  fReflectionFoilPillboxOuterLogical->SetVisAttributes(cyanVisAtt);
  //fwater_only_logical->SetVisAttributes(redVisAtt);



  return fWorldPhysical;
}


void SNneutrinosDetectorConstruction::SetVerbose(G4bool val) { fVerbose = val; }

G4bool SNneutrinosDetectorConstruction::IsVerbose() const { return fVerbose; }

void SNneutrinosDetectorConstruction::PrintError(G4String ed)
{
  G4Exception("SNneutrinosDetectorConstruction:MaterialProperty test", "op001",
              FatalException, ed);
}


