/// \file DetectorConstruction.hh
/// \brief Definition of the sim::DetectorConstruction class

#ifndef SIM_DetectorConstruction_h
#define SIM_DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;

namespace sim
{

/// HPGe detector geometry.
///
/// The active germanium crystal is a cylinder, radius 15 mm and height 10 mm.
/// It is placed so that its coordinate frame matches the SolidStateDetectors
/// (SSD) configuration used in ~/ssd_projects: crystal axis = z, bottom face
/// at z = 0, top face at z = 10 mm. A hit recorded at world position
/// (x, y, z) in mm therefore maps directly onto an SSD CartesianPoint
/// (x, y, z) * 1e-3 in metres -- no coordinate transform is needed.
class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    G4VPhysicalVolume* Construct() override;

    G4LogicalVolume* GetCrystalVolume() const { return fCrystalVolume; }

  private:
    G4LogicalVolume* fCrystalVolume = nullptr;
};

}  // namespace sim

#endif
