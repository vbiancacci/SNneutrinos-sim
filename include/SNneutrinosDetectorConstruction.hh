#ifndef SNneutrinosDetectorConstruction_h
#define SNneutrinosDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericPolycone.hh"
#include "TGraph.h"


class SNneutrinosDetectorMessenger;

class SNneutrinosDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  SNneutrinosDetectorConstruction();
  ~SNneutrinosDetectorConstruction();

  G4VPhysicalVolume* Construct() override;
  void SetVerbose(G4bool verbose);
  G4bool IsVerbose() const;

/*

    if (file.is_open())
	  {
      G4double gg,hh;      
      for (G4int i=0;i<npoints;i++)
        {
          file >> gg >> hh;
          frequencyV[i] = gg*nanometer;
          efficiencyV[i] = hh;
        }
      file.close();
      G4bool successfulInitialization = true;
      G4cout << " data successfully read from file" << endlog;
    }
    else
	  {
      G4err << "Data file not found!" << endlog;
      successfulInitialization = false;  
    }
    if (successfulInitialization)
      {
        //Here the data is loaded
        static const G4double LambdaE = twopi *1.973269602e-16 * m * GeV;
        G4double targetf = LambdaE/energy;

        if (targetf < frequencyV[0] || targetf > frequencyV[npoints-1])
          return 0.0;

        for(j=0;j<npoints-1;j++)
          {
            if (frequencyV[j]<targetf && targetf <= frequencyV[j+1]) break;
          }
        G4double eff = efficiencyV[j] + (targetf-frequencyV[j])
          *(efficiencyV[j+1]-efficiencyV[j])/(frequencyV[j+1]-frequencyV[j]);
        return eff;
      }
    else
      {
        return 0.2;
      }
  }
  */

  G4double stone       = 500.0;  // Hall wall thickness 5 m
  G4double hallrad     = 600.0;  // Hall cylinder diam 12 m
  G4double hallhheight = 850.0;  // Hall cylinder height 17 m
  G4double offset = 200.0;  // shift cavern floor to keep detector centre at origin
 
  // water tank
  /*G4double tankwalltop = 0.6;  // water tank thickness at top 6 mm
  G4double tankwallbot = 0.8;  // water tank thickness at bottom 8 mm
  G4double tankrad     = 550;  // water tank diam 11 m
  G4double tankheight = 650;  // water tank height 13 m
  */
  G4double tank_pit_radius = 9950.0 / 2;  // Radius of the outer tank wall inside icarus pit
  G4double tank_vertical_wall = 10.0;
  G4double tank_horizontal_wall = 20.0; //  # If i read the drawing correctly the horizontal wall is thicker
  G4double tank_base_radius = 12000.0 / 2; //  # Radius of the base of the tank
  G4double tank_pit_height = 800.0; //  # Height of the icarus pit
  G4double tank_base_height =  8877.6 + tank_pit_height; //# height being the z position with regards to the tank bottom at z = 0

  G4double tank_top_height = 9409.8 + tank_pit_height; //  # This value is therefore equal to the entire tank height.
  G4double tank_top_bulge_width = 3330.0; //  # width of the bulged rectangle on top of the tank
  G4double tank_top_bulge_depth = 169.2; //
  G4double tank_top_bulge_radius = 3025.0; //  # radius of the bulged sections

  G4double tank_top_bulge_hwidth = tank_top_bulge_width / 2;
  G4double tank_top_bulge_height = tank_top_height - tank_top_bulge_depth;
  
  
  
  G4double teflon_outer_radius = 5000; //  # rough estimation, the real radius should be smaller than this value
  G4double  offset2 = tank_horizontal_wall;
  G4double out = tank_base_radius - tank_vertical_wall - teflon_outer_radius;
  G4double h_diff = tank_top_height - tank_base_height;
  G4double inner = tank_base_radius - offset2 - tank_top_bulge_width / 2;
  G4double teflon_effective_height = tank_base_height - 4 * offset2 + out * h_diff / inner;
  //  # Accurate would be 2*offset, to be safe we take 4*offset
  
  
  // PMT
  G4double PMTrad     = 10.2;   // diameter 8 inch
  G4double PMTheight  = 2.0;    //random value

  G4GenericPolycone* create_base(std::string name, G4double v_wall=0.0, G4double h_wall=0.0){
    G4double r_base[] = {
      0,
      tank_pit_radius - v_wall,
      tank_pit_radius - v_wall,
      tank_base_radius - v_wall,
      tank_base_radius - v_wall,
      tank_top_bulge_hwidth + v_wall,
      tank_top_bulge_hwidth + v_wall,
      0,
    };
    G4double z_base[] = {
      h_wall,
      h_wall,
      tank_pit_height + h_wall,
      tank_pit_height + h_wall,
      tank_base_height - h_wall,
      tank_top_height - h_wall,
      tank_top_bulge_height - h_wall,
      tank_top_bulge_height - h_wall,
    };
    return new G4GenericPolycone(name, 0, CLHEP::twopi, 8, r_base, z_base);
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
    static const G4double LambdaE = 1239.84193 * nanometer * eV;
    
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
  
  G4Material *foilMat;
  G4MaterialPropertiesTable *foilMPT;
  
  G4Material *worldMat;
  G4MaterialPropertiesTable *worldMPT;

  G4Material *teflonMat;
  G4MaterialPropertiesTable *teflonMPT;
  

};




#endif
