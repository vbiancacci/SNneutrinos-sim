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
  PMTMat = nist->FindOrBuildMaterial("G4_Al"); //CathodeMetalAluminium
  foilMat = new G4Material("VM200", 2.201*g/cm3, 2);
  foilMat->AddElement(nist->FindOrBuildElement("C"),2);
  foilMat->AddElement(nist->FindOrBuildElement("H"),4);

  /*
  G4Element* elGd = new G4Element("Gadolinium", "Gd", 64, 157.25 * CLHEP::g / CLHEP::mole);
  G4Element* elS = new G4Element("Sulfur", "S", 16., 32.066 * CLHEP::g / CLHEP::mole);
  G4Element* O = new G4Element("Oxygen", "O", 8., 16.00 * CLHEP::g / CLHEP::mole);

  G4Material* gadoliniumSulfate = new G4Material("GadoliniumSulfate",  3.01 * g /cm3, 3); // Gd2(SO4)3
  gadoliniumSulfate->AddElement(elGd, 2);
  gadoliniumSulfate->AddElement(elS, 3);
  gadoliniumSulfate->AddElement(O, 12);

  G4Material* GdwaterMat = new G4Material("GdLoadedWater", 1.000000 * g /cm3, 2);
  GdwaterMat->AddMaterial(waterMat, 1. - 0.002);
  GdwaterMat->AddMaterial(gadoliniumSulfate, 0.002);
  */
  


  // https://www.mpi-hd.mpg.de/gerda/public/2008/c08_ndip08_MuonVeto_mk.pdf

  //
  
  // ------------ Generate & Add Material Properties Table ------------
  //
  
  std::vector<G4double> E_in_nm = {600. *nm, 550.*nm, 500.*nm, 450.*nm,  400.*nm, 350.*nm,  300.*nm,  250.*nm,  200.*nm,  150.*nm,  100.* nm};
  std::vector<G4double> E_in_eV= to_E_in_eV(E_in_nm);
  G4int n_energy = size(E_in_eV);
  
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
  

  //reflective foil
  std::vector<G4double> tyvek_E_in_nm = {800. *nm, 500.*nm, 440.*nm, 400.*nm,  355.*nm, 300.*nm,  250.*nm};
  std::vector<G4double> tyvek_E_in_eV= to_E_in_eV(tyvek_E_in_nm);
  G4int tyvek_n_energy = size(tyvek_E_in_eV);
  std::vector<G4double> tyvek_reflectivity = { 0.96, 0.96, 0.96, 0.96, 0.945, 0.90, 0.80}; //from pygeom.optics
  std::vector<double>  tyvek_efficiency (tyvek_n_energy, 0.);

  reflectorMPT->AddProperty("REFLECTIVITY", tyvek_E_in_eV, tyvek_reflectivity);//, false, true);
  reflectorMPT->AddProperty("EFFICIENCY", tyvek_E_in_eV, tyvek_efficiency);//, false, true);

  /*std::vector<G4double> tyvek_reflectivity_border (tyvek_n_energy, 0.);
  std::vector<G4double>  tyvek_transmittance_border (tyvek_n_energy, 1.);;
  std::vector<G4double>  tyvek_efficiency_border (tyvek_n_energy, 0.);;  

  borderMPT->AddProperty("REFLECTIVITY", tyvek_E_in_eV, tyvek_reflectivity_border);//, false, true);
  borderMPT->AddProperty("EFFICIENCY", tyvek_E_in_eV, tyvek_efficiency_border);//, false, true);
  borderMPT->AddProperty("TRANSMITTANCE", tyvek_E_in_eV, tyvek_transmittance_border);//, false, true); 

  //foil
  std::vector<G4double> FiberRIndex (tyvek_n_energy, 1.60);  
  foilMPT->AddProperty("RINDEX", tyvek_E_in_eV, FiberRIndex, false, true);
  foilMat->SetMaterialPropertiesTable(foilMPT);
*/

  //vacuum
  worldMPT->AddProperty("ABSLENGTH", E_in_eV, steelAbsorption, false, true);
  worldMat->SetMaterialPropertiesTable(worldMPT);

/*
  std::vector<G4double> absorption = { //from GEANT4
  
    3.448 * m,  4.082 * m,  6.329 * m,  9.174 * m,  12.346 * m, 13.889 * m,
    15.152 * m, 17.241 * m, 18.868 * m, 20.000 * m, 26.316 * m, 35.714 * m,
    45.455 * m, 47.619 * m, 52.632 * m, 52.632 * m, 55.556 * m, 52.632 * m,
    52.632 * m, 47.619 * m, 45.455 * m, 41.667 * m, 37.037 * m, 33.333 * m,
    30.000 * m, 28.500 * m, 27.000 * m, 24.500 * m, 22.000 * m, 19.500 * m,
    17.500 * m, 14.500 * m
  };

   std::vector<G4double> absorption = { //from MAGE
  
    12*m, 12*m, 12*m, 12*m, 12*m, 12*m,
    12*m, 12*m, 12*m, 12*m, 12*m, 12*m,
    12*m, 12*m, 12*m, 12*m, 12*m, 12*m,
    12*m, 12*m, 12*m, 12*m, 12*m, 12*m,
    12*m, 12*m, 12*m, 12*m, 12*m, 12*m,
    12*m, 12*m
  };
  
     std::vector<G4double> absorption = {  // from article
  
    3.78*m, 4.5*m, 7.4*m, 11.16*m, 14.39*m, 16.15*m,
    17.7*m, 21*m, 23*m, 24.5*m, 30.7*m, 49*m,
    66.7*m, 78.7*m, 94.3*m, 102.1*m, 108.45*m, 157.4*m,
    202*m, 220*m, 211*m, 150*m, 117*m, 87*m//,
    //12*m, 12*m, 12*m, 12*m, 12*m, 12*m,
    //12*m, 12*m
  };

   std::vector<G4double> absorption = { //from website
  
    4.01 * m,  4.28 * m,  5.88 * m,  9.12 * m,  13.86 * m, 17.03 * m,
    22.33 * m, 25.21 * m, 29.12 * m, 32.94 * m, 35.97 * m, 39.78 * m,
    40.03 * m, 40.29 * m, 39.28 * m, 37.12 * m, 35.10 * m, 30.93 * m,
    27.5* m, 23.67 * m, 19.94 * m, 17.11 * m, 12.13* m, 9.53 * m,
    7.18* m, 5.4 * m, 4.28 * m, 3.29 * m, 2.64 * m, 12.15 * m,
    1.77* m, 1.49 * m//,
    //1.21*m, 1.01*m, 0.84*m, 0.70*m, 0.59*m, 0.48*m, 0.40*m, 0.29*m, 0.20*m, 0.14*m
  };
*/




 //------------- Volumes --------------


  //
  // Foil
  //   
 
      G4double foil_in_cryo[] = {1214.9, 1214.9, 1252.9, 1455.7, 1794.5, 2469.8,
        2815.9, 3092.2, 3217.2, 3280.0, 3289.9, 3289.9, 3280.4, 3232.4, 3129.7,
        2969.6, 2815.9, 2507.9, 1983.9, 1489.3, 923.3, 407.5, 208.0, 149.9, 99.9, 0.00, 0.00, 0.00};
          
      G4double foil_out_cryo[] = {1215.0, 1215.0, 1253.0, 1455.8, 1794.6, 2469.9,
        2816.0, 3092.3, 3217.3, 3280.1, 3290.0, 3290.0, 3280.5, 3232.5, 3129.8,
        2969.7, 2816.0, 2508.0, 1984.0, 1489.4, 923.4, 407.6, 208.1, 150., 100., 100., 0.00, 0.00};
      
      G4double foil_z_cryo[] = {4357.8, 3640.1, 3628.3, 3565.3, 3455.4,
        3155.3, 2953.9, 2695.7, 2470.4, 2234.1, 2092.0, -2092.0,
      -2230.7, -2429.6, -2641.0, -2834.5, -2953.9, -3134.8, -3382.6,
      -3555.5, -3688.5, -3753.6, -3768.7, -3783.7, -3787.6,-3787.6, -3787.7, -4357.7}; //-3783.71, -4357.7};
      
      G4double foil_in_tank[] = {1214.9, 1214.9, 6206.0, 6206.0, 0., 0.};
      G4double foil_out_tank[] = {6206.1,6206.1, 6206.1, 6206.1, 6206.1,6206.1}; 
      G4double foil_z_tank[] = {4357.8, 4357.7,4357.7, -4357.7, -4357.7, -4357.8};
      


G4double outerCryo_r_in[] = {0,0,669,909,1127,1293,1491,1661,1824,1999,2190,2460,
  2742,2912,2993,3065,3126,3178,3224,3265,3301,3332,3360,3385,3404,
  3416,3429,3440,3454,3464,3472,3478,3481,3482,3488,3488,3482,3481,
  3478,3472,3464,3454,3440,3429,3416,3404,3385,3360,3332,3301,3265,
  3224,3178,3126,3065,2993,2912,2742,2460,2190,1999,1824,1661,1491,
  1293,1127,909,669,0,0};

G4double outerCryo_r_out[] = {0,1216,1216,1216,1216,1367,1559,1723,1879,2046,2234,
  2501,2778,2946,3024,3094,3152,3204,3247,3287,3322,3353,3380,3404,
  3423,3435,3448,3459,3472,3482,3490,3496,3500,3500,3500,3500,3500,
  3500,3496,3490,3482,3472,3459,3448,3435,3423,3404,3380,3353,3322,
  3287,3247,3204,3152,3094,3024,2946,2778,2501,2234,2046,1879,1723,
  1559,1367,1216,1015,808,461,0};

G4double outerCryo_z[] = {water_h_base/2.,water_h_base/2.,3833,3793,3759,3718,3681,3630,3580,3527,3464,3386,
  3262,3113,3013,2964,2913,2863,2812,2763,2713,2663,2613,2563,2513,
  2470,2438,2400,2363,2313,2263,2213,2163,2111,2000,2000,-2000,-2000,
  -2111,-2163,-2213,-2263,-2313,-2363,-2400,-2438,-2470,-2513,-2563,
  -2613,-2663,-2713,-2763,-2812,-2863,-2913,-2964,-3013,-3113,-3262,
  -3386,-3464,-3527,-3580,-3630,-3681,-3718,-3759,-3793,-3833,-3848};

  int n = sizeof(outerCryo_r_out) / sizeof(outerCryo_r_out[0]); 
  
  //
  // World (hall)
  //
  auto* worldSolid =
    new G4Tubs("World", 0.0 * cm, (hallrad + stone + 0.1) * cm,
               (hallhheight + stone + offset + 0.1) * cm, 0.0, CLHEP::twopi);
  auto* fWorldLogical  = new G4LogicalVolume(worldSolid, worldMat, "World_log");
  auto* fWorldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLogical,
                                           "World_phys", nullptr, false, 0);

  //Tank
  G4double tank_h_base= tank_base_height-tank_pit_height+ tank_top_height-tank_base_height;
  G4double tank_h_base2= tank_pit_height+ tank_top_height;
  G4cout << " tank_h_base " << tank_h_base << G4endl;
  G4cout << " tank_h_base2 " << tank_h_base2 << G4endl;

  auto* TankBase = new G4Tubs("TankBase", 0, tank_base_radius, tank_h_base/2., 0., CLHEP::twopi);
  auto* TankBottom = new G4Tubs("TankBottom", 0, tank_pit_radius, (tank_pit_height)/2., 0., CLHEP::twopi);
  auto* TankSolid = new G4UnionSolid("TankUnion", TankBase, TankBottom, 0, G4ThreeVector(0, 0,-(tank_h_base+tank_pit_height)/2.));
  auto* fTankLogical = new G4LogicalVolume(TankSolid, steelMat, "Tank_log");
  auto* fTankPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), fTankLogical, "Tank_phys", fWorldLogical, false, 0, true);

  //water
  
  auto* WaterBase = new G4Tubs("TankBase", 0,r_water_base, water_h_base/2., 0., CLHEP::twopi);
  auto* WaterBottom = new G4Tubs("TankBottom", 0, r_water_bottom, (tank_pit_height)/2., 0., CLHEP::twopi);
  auto* WaterSolid = new G4UnionSolid("WaterSolid", WaterBase, WaterBottom, 0, G4ThreeVector(0, 0,-(water_h_base+tank_pit_height)/2.));
  auto* fWaterLogical = new G4LogicalVolume(WaterSolid, waterMat, "Water_log");
  auto* fWaterPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), fWaterLogical, "Water_phys", fTankLogical, false, 0, true);

 
  const int n_out = 70;
  G4double zeros[n_out]={0.};
  G4double FoilCryo_r_out[] = {0,1216,1216,1216,1216,1367,1559,1723,1879,2046,2234,
    2501,2778,2946,3024,3094,3152,3204,3247,3287,3322,3353,3380,3404,
    3423,3435,3448,3459,3472,3482,3490,3496,3500,3500,3500,3500,3500,
    3500,3496,3490,3482,3472,3459,3448,3435,3423,3404,3380,3353,3322,
    3287,3247,3204,3152,3094,3024,2946,2778,2501,2234,2046,1879,1723,
    1559,1367,1216,1015,808,461,0};
  G4double FoilCryo_z[] = {water_h_base/2,water_h_base/2.,3833,3793,3759,3718,3681,3630,3580,3527,3464,3386,
    3262,3113,3013,2964,2913,2863,2812,2763,2713,2663,2613,2563,2513,
    2470,2438,2400,2363,2313,2263,2213,2163,2111,2000,2000,-2000,-2000,
    -2111,-2163,-2213,-2263,-2313,-2363,-2400,-2438,-2470,-2513,-2563,
    -2613,-2663,-2713,-2763,-2812,-2863,-2913,-2964,-3013,-3113,-3262,
    -3386,-3464,-3527,-3580,-3630,-3681,-3718,-3759,-3793,-3836};
  int nn = sizeof(FoilCryo_r_out) / sizeof(FoilCryo_r_out[0]);

int numF = 500+tyvek_thickness; 
  std::transform(FoilCryo_r_out, FoilCryo_r_out + nn, FoilCryo_r_out, [numF](int x) { return x + numF; });
  auto* FoilCryostatSolid = new G4Polycone("FoilCryostat", 0, CLHEP::twopi, nn, FoilCryo_z, zeros, FoilCryo_r_out);
  auto* fFoilCryostatLogical  = new G4LogicalVolume(FoilCryostatSolid, foilMat, "FoilCryostat_log");
  auto* fFoilCryostatPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0.), fFoilCryostatLogical,"FoilCryostat_phys", fWaterLogical, false, 0, true);

  auto* FoilBottomSolid = new G4Tubs("FoilBottom", 0, r_water_bottom, tyvek_thickness/2., 0., CLHEP::twopi);
  auto* fFoilBottomLogical  = new G4LogicalVolume(FoilBottomSolid, foilMat, "FoilBottom_log");
  auto* fFoilBottomPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,-water_h_base/2.-tank_pit_height+tyvek_thickness/2.), fFoilBottomLogical,"FoilBottom_phys", fWaterLogical, false, 0, true);
  
  //reflector foil at the tank wall
  //auto* FoilWallSolid = new G4Tubs("FoilWall", r_water_base-tyvek_thickness, r_water_base, (water_h_base)/2., 0.0, CLHEP::twopi);
  //auto* fFoilWallLogical  = new G4LogicalVolume(FoilWallSolid, foilMat, "FoilWall_log");
  //auto* fFoilWallPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), fFoilWallLogical,"FoilWall_phys", fWaterLogical, false, 0, true);
  
  //reflector foil at the PMT wall
  auto* FoilWallSolid = new G4Tubs("FoilWall", tyvek_effective_radius, tyvek_effective_radius +tyvek_thickness, (water_h_base+tank_pit_height)/2., 0.0, CLHEP::twopi);
  auto* fFoilWallLogical  = new G4LogicalVolume(FoilWallSolid, foilMat, "FoilWall_log");
  auto* fFoilWallPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,-tank_pit_height/2.), fFoilWallLogical,"FoilWall_phys", fWaterLogical, false, 0, true);
  


  int num = 500; 
  std::transform(outerCryo_r_out, outerCryo_r_out + n, outerCryo_r_out, [num](int x) { return x + num; });
  std::transform(outerCryo_r_in, outerCryo_r_in + n, outerCryo_r_in, [num](int x) { return x + num; });
  auto CryostatSolid = new G4Polycone("Cryostat", 0, CLHEP::twopi, n, outerCryo_z, zeros, outerCryo_r_out);
  auto* fCryostatLogical  = new G4LogicalVolume(CryostatSolid, worldMat, "Cryostat_log");
  auto* fCryostatPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), fCryostatLogical,"Cryostat_phys", fFoilCryostatLogical, false, 0, true);
  
  //PillBox
  G4RotationMatrix* Rot = new G4RotationMatrix; 
  Rot->rotateZ(-CLHEP::pi / 2.0*rad);
  Rot->rotateY(0);
  Rot->rotateX(-CLHEP::pi /2.0*rad);
  auto* pillbox_tube = new G4Tubs("pillbox_tube", shielding_foot_ir, shielding_foot_or, pillbox_cryo_bottom_height/2., 0, CLHEP::twopi ); //  # outer steel cylinder
  auto* manhole_pillbox_arc = new G4Tubs("manhole_pillbox", manhole_inner_radius, manhole_outer_radius, manhole_height/2., 0, manhole_angle);
  auto* manholepillbox_box = new G4Box("manholepillbox_box", manhole_outer_radius, manhole_outer_radius/2., manhole_height/2.);
  auto* manhole_pillbox = new G4UnionSolid("manhole_union", manhole_pillbox_arc, manholepillbox_box, 0, G4ThreeVector(0, -0.5 * manhole_outer_radius));
  auto* PillboxSolid = new G4SubtractionSolid("pillbox_subtraction1", pillbox_tube, manhole_pillbox, Rot, G4ThreeVector(0, 0, 0 - manhole_offset));
  auto* fPillboxLogical = new G4LogicalVolume(PillboxSolid, steelMat, "Pillbox_log");
  auto* fPillboxPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,pillbox_offset2), fPillboxLogical, "Pillbox_phys", fWaterLogical, false, 0, true);

  auto* FoilPillboxInnerSolid1 = new G4Tubs("FoilPillboxInner", shielding_foot_ir - tyvek_thickness, shielding_foot_ir, (pillbox_cryo_bottom_height - shielding_foot_thickness)/2., 0., CLHEP::twopi);
  auto* FoilPillboxInnerSolid = new G4SubtractionSolid(
    "pillbox_inner_reflection_foil_tube ",
    FoilPillboxInnerSolid1,
    manhole_pillbox, Rot, G4ThreeVector(0, 0, 0 - manhole_offset));
  auto* fFoilPillboxInnerLogical = new G4LogicalVolume(FoilPillboxInnerSolid, foilMat, "FoilPillboxInner_log");
  auto* fFoilPillboxInnerPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,pillbox_offset2),  fFoilPillboxInnerLogical, "FoilPillboxInner_phys", fWaterLogical, false, 0, true);

  auto* FoilPillboxOuterSolid1 = new G4Tubs(
    "FoilPillboxIOuter",
    shielding_foot_or,
    shielding_foot_or + tyvek_thickness,
    pillbox_cryo_bottom_height/2.,
    0, CLHEP::twopi );
  
    auto* FoilPillboxOuterSolid = new G4SubtractionSolid(
    "pillbox_outer_reflection_foil_tube ",
    FoilPillboxOuterSolid1,
    manhole_pillbox,Rot, G4ThreeVector(0, 0, 0 - manhole_offset));  
    auto* fFoilPillboxOuterLogical = new G4LogicalVolume(FoilPillboxOuterSolid, foilMat, "FoilPillboxOuter_log");
    auto* fFoilPillboxOuterPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,pillbox_offset2), fFoilPillboxOuterLogical, "FoilPillboxOuter_phys",  fWaterLogical, false, 0, true);


  // ------------- Surfaces --------------


  // Water-Tyvek surface
  G4OpticalSurface* Tyvek_opSurface = new G4OpticalSurface("Tyvek_opSurface",
      unified,  
      groundfrontpainted,          
      dielectric_dielectric 
  );
Tyvek_opSurface->SetMaterialPropertiesTable(reflectorMPT);
  
  G4OpticalSurface* WaterToTyvek_opSurface = new G4OpticalSurface("WaterToTyvek_opSurface",
    unified,  
    groundfrontpainted,          
    dielectric_dielectric 
);
WaterToTyvek_opSurface->SetMaterialPropertiesTable(borderMPT);

/*
auto* ReflectionFoilCryostatSurface = new G4LogicalBorderSurface("ReflectionFoilCryostatSurface", fFoilCryostatPhysical,fWaterPhysical, WaterToTyvek_opSurface);
auto* ReflectionFoilWallSurface = new G4LogicalBorderSurface("ReflectionFoilWallSurface", fFoilWallPhysical,fWaterPhysical, WaterToTyvek_opSurface);
auto* ReflectionFoilBottomSurface = new G4LogicalBorderSurface("ReflectionFoilBottomSurface", fFoilBottomPhysical,fWaterPhysical, WaterToTyvek_opSurface);
auto* ReflectionFoilPillobxOuterSurface = new G4LogicalBorderSurface("ReflectionFoilPillboxOuterSurface", fFoilPillboxOuterPhysical,fWaterPhysical, WaterToTyvek_opSurface);
auto* ReflectionFoilPillobxInnerSurface = new G4LogicalBorderSurface("ReflectionFoilPillboxInnerSurface", fFoilPillboxInnerPhysical,fWaterPhysical, WaterToTyvek_opSurface);
*/

auto* ReflectionFoilCryostatSkin = new  G4LogicalSkinSurface("ReflectionFoilCryostatSkin", fFoilCryostatLogical, Tyvek_opSurface);
auto* ReflectionFoilWallSkin = new  G4LogicalSkinSurface("ReflectionFoilWallSkin", fFoilWallLogical, Tyvek_opSurface);
auto* ReflectionFoilBottomSkin = new  G4LogicalSkinSurface("ReflectionFoilBottomSkin", fFoilBottomLogical, Tyvek_opSurface);
auto* ReflectionFoilPillboxOuterSkin = new  G4LogicalSkinSurface("ReflectionFoilPillboxOuterkin", fFoilPillboxOuterLogical, Tyvek_opSurface);
auto* ReflectionFoilPillboxInnerSkin = new  G4LogicalSkinSurface("ReflectionFoilPillboxInnerkin", fFoilPillboxInnerLogical, Tyvek_opSurface);



 


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
 
 G4int bottom_circles = 4;  //in the previous versione =5
 G4int n_bottom_PMT=0;
 for (int j=1; j<=bottom_circles; j++){ //150 PMT at the bottom // in the previous versione they were 100, updated from the latest Josef's presentation given in TC call in 10-12-24
    G4int PMT_bottom_circle = G4int(45/(bottom_circles - j +1));
    for (int i=1;i<=PMT_bottom_circle; i++) {
      pos_r =  G4double((tank_pit_radius +400.)  * j /(bottom_circles+1.)) ;//G4double(water_r_out[0]  * j /(bottom_circles+1.)) ;
      pos_theta = G4double(CLHEP::twopi*i/PMT_bottom_circle);
      pos_x = pos_r * cos(pos_theta);
      pos_y = pos_r * sin(pos_theta);
      pos_z =  10.*PMTheight+tyvek_thickness-water_h_base/2.-tank_pit_height;
      PMT_ID = 8000+j*100+i;  //ID starting with 8 is set to the bottom and the second number indicates the ring - number 1 is the most internal ring
      new G4PVPlacement(nullptr, G4ThreeVector(pos_x, pos_y, pos_z), fPMTLogical, "PMT_phys", fWaterLogical, false, PMT_ID);
      n_bottom_PMT++;
      }
    if (j==1){
      PMT_ID=8000;
      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, pos_z), fPMTLogical, "PMT_phys", fWaterLogical, false, PMT_ID);
    }  
  } 
  G4cout <<"bottom PMT  " << n_bottom_PMT << G4endl;

  G4int n_PMT=0;
  int rings = 11;       //rings=7, PMT_later=71 in the previous version, updated to 11 from the latest Josef's presentation given in TC call in 10-12-24
  int PMT_lateral=45;
  for (int j=1; j<=rings; j++) { //500 PMT at the later surface 80cm from the cryostat.
    for (int i=0;i<PMT_lateral; i++) {
      pos_theta = double(CLHEP::twopi)*i/PMT_lateral;
      G4RotateY3D rotTheta; //(0);
      G4int quad=i/25;
      G4double pos_theta_mod =pos_theta-quad*CLHEP::pi/2.;
      if (quad==0 || quad==2){
        rotTheta= G4RotateY3D(CLHEP::pi/2.+pos_theta_mod);}
      else{
        rotTheta= G4RotateY3D(pos_theta_mod);}
      G4RotateX3D rotPhi(CLHEP::pi/2.);
      pos_r = 800. + outerCryo_r_out[34]- PMTheight*10;
      pos_x = pos_r * cos(pos_theta);
      pos_y = pos_r * sin(pos_theta);
      G4Translate3D shift(pos_x, pos_y, (outerCryo_z[2]+4500)*j/(rings+1)-5000);// (2*water_z[0] * j/(rings+1)-water_z[0]));
      auto transform = shift*rotPhi*rotTheta; 
      PMT_ID = (rings-j+1)*100+i;   //ID starting with 1 belogs to the first ring from the top
      new G4PVPlacement(transform, fPMTLogical, "PMT_phys", fWaterLogical, false, PMT_ID);
      n_PMT+=1;
    } 
  }
  G4cout <<"lateral PMT  " << n_PMT << G4endl;

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

  auto* water_only_solid1 = new G4SubtractionSolid("water_only1", WaterSolid, FoilCryostatSolid, 0, G4ThreeVector(0,0, 0));
  auto* water_only_solid2 = new G4SubtractionSolid("water_only2", water_only_solid1, FoilBottomSolid, 0, G4ThreeVector(0,0, -water_h_base/2.-tank_pit_height+tyvek_thickness/2.));
  auto* water_only_solid3 = new G4SubtractionSolid("water_only3", water_only_solid2, FoilWallSolid, 0, G4ThreeVector(0,0, -tank_pit_height/2.));
  auto* water_only_solid4 = new G4SubtractionSolid("water_only4", water_only_solid3, FoilPillboxOuterSolid, 0, G4ThreeVector(0,0, pillbox_offset2));
  auto* water_only_solid5 = new G4SubtractionSolid("water_only5", water_only_solid4, FoilPillboxInnerSolid, 0, G4ThreeVector(0,0, pillbox_offset2));

  auto* fwater_only_logical  = new G4LogicalVolume(water_only_solid5, waterMat, "Water_only_log");
  //auto* fWaterOnlyPhysical = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), fwater_only_logical, "Water_only_phys", fTankLogical, false, 0, true);

  

  //
  // Visualization
  //
  G4Color testColor_water(170 / 255., 191 / 225., 219 / 225.); //170, 191, 219
  auto*   testVisAtt_water = new G4VisAttributes(testColor_water);
  testVisAtt_water->SetVisibility(true);
  auto* greyVisAtt = new G4VisAttributes(G4Colour::Grey());
  auto* greenVisAtt = new G4VisAttributes(G4Colour::Green());
  auto* cyanVisAtt = new G4VisAttributes(G4Colour::Cyan());
  auto* redVisAtt = new G4VisAttributes(G4Colour::Red());
  greyVisAtt->SetVisibility(true);

  //fTankLogical->SetVisAttributes(greyVisAtt);
  fWaterLogical->SetVisAttributes(greyVisAtt);
  //fwater_only_logical->SetVisAttributes(cyanVisAtt);
  fFoilCryostatLogical->SetVisAttributes(redVisAtt);
  //fPillboxLogical->SetVisAttributes(greenVisAtt);
  fPMTLogical->SetVisAttributes(cyanVisAtt);
  


  return fWorldPhysical;
}


void SNneutrinosDetectorConstruction::SetVerbose(G4bool val) { fVerbose = val; }

G4bool SNneutrinosDetectorConstruction::IsVerbose() const { return fVerbose; }

void SNneutrinosDetectorConstruction::PrintError(G4String ed)
{
  G4Exception("SNneutrinosDetectorConstruction:MaterialProperty test", "op001",
              FatalException, ed);
}


