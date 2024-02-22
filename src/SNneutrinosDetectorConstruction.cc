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




  // https://www.mpi-hd.mpg.de/gerda/public/2008/c08_ndip08_MuonVeto_mk.pdf

  //
  
  // ------------ Generate & Add Material Properties Table ------------
  //
  
  std::vector<G4double> photonEnergy = {
    2.034 * eV, 2.068 * eV, 2.103 * eV, 2.139 * eV, 2.177 * eV, 2.216 * eV,
    2.256 * eV, 2.298 * eV, 2.341 * eV, 2.386 * eV, 2.433 * eV, 2.481 * eV,
    2.532 * eV, 2.585 * eV, 2.640 * eV, 2.697 * eV, 2.757 * eV, 2.820 * eV,
    2.885 * eV, 2.954 * eV, 3.026 * eV, 3.102 * eV, 3.181 * eV, 3.265 * eV,
    3.353 * eV, 3.446 * eV, 3.545 * eV, 3.649 * eV, 3.760 * eV, 3.877 * eV,
    4.002 * eV, 4.136 * eV//, 
    //4.275 * eV, 4.428* eV, 4.592* eV, 4.769 * eV, 4.959* eV, 5.166* eV, 5.391* eV, 5.636* eV, 5.904 *eV, 6.199*eV
  };

  // Water
  std::vector<G4double> waterRIndex = {
    1.3435, 1.344,  1.3445, 1.345,  1.3455, 1.346,  1.3465, 1.347,
    1.3475, 1.348,  1.3485, 1.3492, 1.35,   1.3505, 1.351,  1.3518,
    1.3522, 1.3530, 1.3535, 1.354,  1.3545, 1.355,  1.3555, 1.356,
    1.3568, 1.3572, 1.358,  1.3585, 1.359,  1.3595, 1.36,   1.3608//,
    //1.3610, 1.3612, 1.3614, 1.3617, 1.3620, 1.3664, 1.3708, 1.3776,  1.3868, 1.3960
    ////1.351, 1.3530, 1.3556, 1.3588, 1.3620, 1.3664, 1.3708, 1.3776,  1.3868, 1.3960
  };

  std::vector<G4double> waterAbsorption = { //from GEANT4
  
    3.448 * m,  4.082 * m,  6.329 * m,  9.174 * m,  12.346 * m, 13.889 * m,
    15.152 * m, 17.241 * m, 18.868 * m, 20.000 * m, 26.316 * m, 35.714 * m,
    45.455 * m, 47.619 * m, 52.632 * m, 52.632 * m, 55.556 * m, 52.632 * m,
    52.632 * m, 47.619 * m, 45.455 * m, 41.667 * m, 37.037 * m, 33.333 * m,
    30.000 * m, 28.500 * m, 27.000 * m, 24.500 * m, 22.000 * m, 19.500 * m,
    17.500 * m, 14.500 * m
  };


  std::vector<G4double>  steelRIndex  (size(photonEnergy), 2.86);
  std::vector<G4double>  steelAbsorption (size(photonEnergy), 1.e-20*m);

    
  //water
  waterMPT->AddProperty("RINDEX", photonEnergy, waterRIndex, false, true);
  waterMPT->AddProperty("ABSLENGTH", photonEnergy, waterAbsorption, false, true);
  G4cout << "Water G4MaterialPropertiesTable:" << G4endl;
  waterMPT->DumpTable();
  waterMat->SetMaterialPropertiesTable(waterMPT);

  //steel
  steelMPT->AddProperty("RINDEX", photonEnergy, steelRIndex, false, true);
  steelMPT->AddProperty("ABSLENGTH", photonEnergy, steelAbsorption, false, true); //if uncommented there will be some absorption even in the cryostat/water tank (????)
  steelMat->SetMaterialPropertiesTable(steelMPT);
  



  std::vector <double> WLS_emission;
  WLS_emission=EmissionSpectrum("WLS_em",photonEnergy);
  
  std::vector <double> WLS_absorption;
  WLS_absorption=EmissionSpectrum("WLS_abs",photonEnergy);

  std::vector<G4double>  refReflectivity;
  refReflectivity=EmissionSpectrum("WLS_ref",photonEnergy);
  reflectorMPT->AddProperty("REFLECTIVITY", photonEnergy, refReflectivity, false, true);


  //foil
  std::vector<G4double> FiberRIndex (size(photonEnergy), 1.60);  
  foilMPT->AddProperty("RINDEX", photonEnergy, FiberRIndex, false, true);
  foilMPT->AddProperty("WLSABSLENGTH",photonEnergy, WLS_absorption, false, true);
  foilMPT->AddProperty("WLSCOMPONENT",photonEnergy, WLS_emission, false, true);
  G4double TimeConstant = 0.01 *ns;
  foilMPT->AddConstProperty("WLSTIMECONSTANT", TimeConstant);
  
  foilMat->SetMaterialPropertiesTable(foilMPT);


  //vacuum
  worldMPT->AddProperty("ABSLENGTH", photonEnergy, steelAbsorption, false, true);
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
  // World (hall)
  //
  auto* worldSolid =
    new G4Tubs("World", 0.0 * cm, (hallrad + stone + 0.1) * cm,
               (hallhheight + stone + offset + 0.1) * cm, 0.0, CLHEP::twopi);
  auto* fWorldLogical  = new G4LogicalVolume(worldSolid, worldMat, "World_log");
  auto* fWorldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLogical,
                                           "World_phys", nullptr, false, 0);



  //
  // Water
  //
        const int n_in = 27;

      G4double water_r_in[] = {1215.0, 1215.0, 1253.0, 1455.8, 1794.6, 2469.9,
        2816.0, 3092.3, 3217.3, 3280.1, 3290.0, 3290.0, 3280.5, 3232.5, 3129.8,
        2969.7, 2816.0, 2508.0, 1984.0, 1489.4, 923.4, 407.6, 208.1, 150., 100., 0.00, 0.00};

      G4double water_r_out[] = {6206, 6206, 6206, 6206, 6206, 6206,
        6206, 6206, 6206, 6206, 6206, 6206, 6206, 6206, 6206, 6206,
        6206, 6206, 6206, 6206, 6206, 6206, 6206, 6206, 6206, 6206, 6206};

      G4double water_z[] = {4357.7, 3640.1, 3628.3, 3565.3, 3455.4,
         3155.3, 2953.9, 2695.7, 2470.4, 2234.1, 2092.0, -2092.0,
        -2230.7, -2429.6, -2641.0, -2834.5, -2953.9, -3134.8, -3382.6,
        -3555.5, -3688.5, -3753.6, -3768.7, -3783.7, -3785.7, -3787.7, -4357.7}; //-3783.71, -4357.7};
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
      

  //
  // Cryostat
  //
      G4double outerCryo_r_in[] = {1190.0, 1190.0, 1190.0, 1390.4, 1740.1,
        2427.8, 2777.9, 3063.6, 3191.1, 3254.9, 3265.0, 3265.0, 3255.4, 3206.6,
        3101.9, 2937.6, 2777.9, 2466.4, 1933.8, 1425.4, 823.1, 33.6, 0.0, 0.0};
  
      G4double outerCryo_r_out[] = {1205.0, 1205.0, 1243.0, 1445.8, 1784.6, 2459.9,
        2806.0, 3082.3, 3207.3, 3270.1, 3280.0, 3280.0, 3270.5, 3222.5, 3119.8,
        2959.7, 2806.0, 2498.0, 1974.0, 1479.4, 913.4, 397.6, 198.1, 0.0};

      G4double outerCryo_z[] = {4807.0, 3640.1, 3628.3, 3565.3, 3455.4, 3155.3,
        2953.9, 2695.7, 2470.4, 2234.1, 2092.0, -2092.0, -2230.7, -2429.6, -2641.0,
        -2834.5, -2953.9, -3134.8, -3382.6, -3555.5, -3688.5, -3753.6, -3768.7, -3783.7};


  auto* EnvelopeSolid =   new G4Tubs("Envelope", 0.0*cm , water_r_out[0]+20.,outerCryo_z[0]+20., 0.0*cm, CLHEP::twopi);
  auto* fEnvelopeLogical  = new G4LogicalVolume(EnvelopeSolid, steelMat, "Envelope_log");
  auto* fEnvelopePhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fEnvelopeLogical,
                                           "Envelope_phys", fWorldLogical, false, 0);

  auto* WaterSolid = new G4Polycone("Water", 0, CLHEP::twopi, n_in, water_z, water_r_in, water_r_out);
  auto* fWaterLogical  = new G4LogicalVolume(WaterSolid, waterMat, "Water_log");
  auto* fWaterPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWaterLogical,
                                           "Water_phys", fEnvelopeLogical, false, 0, true);

  const int n_out = 24;
  auto* FoilSolidCryostat = new G4Polycone("FoilCryostat", 0, CLHEP::twopi, std::size(foil_in_cryo), foil_z_cryo, foil_in_cryo, foil_out_cryo);
  auto* FoilSolidTank = new G4Polycone("FoilTank",  0, CLHEP::twopi,6,foil_z_tank,foil_in_tank,foil_out_tank );
  auto* FoilSolid = new G4UnionSolid("FoilCryostat+FoilTank", FoilSolidCryostat, FoilSolidTank);
  auto* fFoilLogical  = new G4LogicalVolume(FoilSolid, foilMat, "Foil_log");
  auto* fFoilPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fFoilLogical,"Foil_phys", fEnvelopeLogical, false, 0, true);
  


  G4double zeros[n_out]={0.};
  auto CryostatSolid = new G4Polycone("Cryostat", 0, CLHEP::twopi, n_out, outerCryo_z, zeros, outerCryo_r_in);
  auto* fCryostatLogical  = new G4LogicalVolume(CryostatSolid, worldMat, "Cryostat_log");
  auto* fCryostatPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fCryostatLogical,"Cryostat_phys", fEnvelopeLogical, false, 0, true);
    
 

  // ------------- Surfaces --------------

  // Foil-Steel surface
  G4OpticalSurface* opSurface_FoilSteel = new G4OpticalSurface("opSurface_FoilSteel",
      unified,  
      polished,          
      dielectric_metal  
  );

  opSurface_FoilSteel->SetMaterialPropertiesTable(reflectorMPT);
   G4LogicalBorderSurface *FoilSteelSurface = new G4LogicalBorderSurface("FoilSteelSurface", fFoilPhysical,fEnvelopePhysical, opSurface_FoilSteel);
   //G4LogicalBorderSurface *outSurface = new G4LogicalBorderSurface("CryoSurface", fEnvelopePhysical,fFoilPhysical, opSurface);
  G4OpticalSurface* opticalSurfaceFoilSteel = dynamic_cast<G4OpticalSurface*>(
    FoilSteelSurface->GetSurface(fFoilPhysical,fEnvelopePhysical)->GetSurfaceProperty());
    //outSurface->GetSurface(fFoilPhysical, fWaterPhysical)->GetSurfaceProperty());

if(opSurface_FoilSteel)
    opSurface_FoilSteel->DumpInfo();

/*
  //Foil-Water surface
    G4OpticalSurface* opSurface_FoilWater = new G4OpticalSurface("opSurface_WaterFoil",
      unified,  
      ground,          //polishedfrontpainted
      dielectric_dielectric  // dielectric_metal
  );
   G4LogicalBorderSurface *FoilWaterSurface = new G4LogicalBorderSurface("FoilWaterSurface",fFoilPhysical, fWaterPhysical, opSurface_FoilWater);
   G4OpticalSurface* opticalSurfaceFoilWater = dynamic_cast<G4OpticalSurface*>(
    FoilWaterSurface->GetSurface(fFoilPhysical,fWaterPhysical)->GetSurfaceProperty());
  
  opticalSurfaceFoilWater->SetSigmaAlpha(1.);
  if(opSurface_FoilWater)
    opSurface_FoilWater->DumpInfo();
  
*/

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
 
 G4int bottom_circles = 5;
 G4int n_bottom_PMT=0;
 for (int j=1; j<=bottom_circles; j++){ //100 PMT at the bottom
    G4int PMT_bottom_circle = G4int(45/(bottom_circles - j +1));
    for (int i=1;i<=PMT_bottom_circle; i++) {
      pos_r = G4double(water_r_out[0]  * j /(bottom_circles+1.)) ;
      pos_theta = G4double(CLHEP::twopi*i/PMT_bottom_circle);
      pos_x = pos_r * cos(pos_theta);
      pos_y = pos_r * sin(pos_theta);
      pos_z = -water_z[0]+10.*PMTheight;
      PMT_ID = 8000+j*100+i;  //ID starting with 8 is set to the bottom and and second number indicates the ring - number 1 is the most internal ring
      new G4PVPlacement(nullptr, G4ThreeVector(pos_x, pos_y, pos_z), fPMTLogical, "PMT_phys", fWaterLogical, false, PMT_ID);
      n_bottom_PMT++;
      }
  } 
  G4cout <<"bottom PMT  " << n_bottom_PMT << G4endl;

  G4int n_PMT=0;
  for (int j=1; j<=7; j++) { //500 PMT at the later surface 80cm from the cryostat
    for (int i=0;i<71; i++) {
      pos_theta = double(CLHEP::twopi)*i/71;
      G4RotateY3D rotTheta; //(0);
      G4int quad=i/25;
      G4double pos_theta_mod =pos_theta-quad*CLHEP::pi/2.;
      if (quad==0 || quad==2){
        rotTheta= G4RotateY3D(CLHEP::pi/2.+pos_theta_mod);}
      else{
        rotTheta= G4RotateY3D(pos_theta_mod);}
      G4RotateX3D rotPhi(CLHEP::pi/2.);
      pos_r = 800. + outerCryo_r_out[10];
      pos_x = pos_r * cos(pos_theta);
      pos_y = pos_r * sin(pos_theta);
      G4Translate3D shift(pos_x, pos_y,  (2*water_z[0] * j/8.-water_z[0]));
      auto transform = shift*rotPhi*rotTheta; 
      PMT_ID = (7-j+1)*100+i;   //ID starting with 1 belogs to the first ring from the top
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


  

  //
  // Visualization
  //
  G4Color testColor_water(170 / 255., 191 / 225., 219 / 225.); //170, 191, 219
  auto*   testVisAtt_water = new G4VisAttributes(testColor_water);
  testVisAtt_water->SetVisibility(true);
  auto* greyVisAtt = new G4VisAttributes(G4Colour::Grey());
  auto* greenVisAtt = new G4VisAttributes(G4Colour::Cyan());
  auto* redVisAtt = new G4VisAttributes(G4Colour::Red());
  greyVisAtt->SetVisibility(true);

  //fTankLogical->SetVisAttributes(greyVisAtt);
  fWaterLogical->SetVisAttributes(redVisAtt);
  fCryostatLogical->SetVisAttributes(greenVisAtt);
  fEnvelopeLogical->SetVisAttributes(redVisAtt);
  fPMTLogical->SetVisAttributes(greenVisAtt);
  


  return fWorldPhysical;
}


void SNneutrinosDetectorConstruction::SetVerbose(G4bool val) { fVerbose = val; }

G4bool SNneutrinosDetectorConstruction::IsVerbose() const { return fVerbose; }

void SNneutrinosDetectorConstruction::PrintError(G4String ed)
{
  G4Exception("SNneutrinosDetectorConstruction:MaterialProperty test", "op001",
              FatalException, ed);
}


