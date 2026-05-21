/// \file SteppingAction.cc
/// \brief Implementation of the sim::SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"

namespace sim
{

SteppingAction::SteppingAction(RunAction* runAction) : fRunAction(runAction) {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.) return;

  // Resolve the crystal logical volume once, lazily.
  if (!fCrystalVolume) {
    const auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fCrystalVolume = detConstruction->GetCrystalVolume();
  }

  // Only score steps inside the germanium crystal.
  G4LogicalVolume* volume =
    step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  if (volume != fCrystalVolume) return;

  // Interaction position = step mid-point (world coordinates == SSD frame).
  const G4ThreeVector pos =
    0.5 * (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition());

  const G4int eventID =
    G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  fRunAction->RecordHit(eventID, pos.x() / mm, pos.y() / mm, pos.z() / mm,
                        edep / keV);
}

}  // namespace sim
