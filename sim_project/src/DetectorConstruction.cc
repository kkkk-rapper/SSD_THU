/// \file DetectorConstruction.cc
/// \brief Implementation of the sim::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

namespace sim
{

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();
  G4bool checkOverlaps = true;

  // --- Crystal dimensions (must stay in sync with the SSD config.yaml) ---
  const G4double crystalRadius = 15. * mm;
  const G4double crystalHeight = 10. * mm;

  //
  // World
  //
  const G4double worldSize = 20. * cm;
  G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");

  auto solidWorld =
    new G4Box("World", 0.5 * worldSize, 0.5 * worldSize, 0.5 * worldSize);
  auto logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
  auto physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World",
                                     nullptr, false, 0, checkOverlaps);

  //
  // HPGe crystal
  //
  // A G4Tubs is centred on its own origin, so we place it at z = +height/2.
  // The crystal then spans world-z in [0, 10] mm, matching the SSD frame.
  //
  G4Material* geMat = nist->FindOrBuildMaterial("G4_Ge");

  auto solidCrystal = new G4Tubs("Crystal", 0., crystalRadius, 0.5 * crystalHeight,
                                 0. * deg, 360. * deg);
  fCrystalVolume = new G4LogicalVolume(solidCrystal, geMat, "Crystal");
  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.5 * crystalHeight),
                    fCrystalVolume, "Crystal", logicWorld, false, 0, checkOverlaps);

  // Visualization attributes.
  auto crystalVis = new G4VisAttributes(G4Colour(0.2, 0.6, 1.0, 0.4));
  crystalVis->SetForceSolid(true);
  fCrystalVolume->SetVisAttributes(crystalVis);
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

  // NOTE: a real HPGe detector also has a cryostat, an aluminium end-cap
  // window and an n+ dead layer. Add them here as extra volumes when needed;
  // they do not change the crystal coordinate frame used for the SSD hand-off.

  return physWorld;
}

}  // namespace sim
