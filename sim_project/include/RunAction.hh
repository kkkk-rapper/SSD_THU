/// \file RunAction.hh
/// \brief Definition of the sim::RunAction class

#ifndef SIM_RunAction_h
#define SIM_RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <fstream>

class G4Run;

namespace sim
{

/// Owns the hits output file (hits.csv).
///
/// Writes one CSV row per energy-depositing step inside the crystal:
///
///     event,x_mm,y_mm,z_mm,edep_keV
///
///   - event   : Geant4 event id = one primary gamma
///   - x,y,z   : interaction position in mm, in the SSD crystal frame
///   - edep    : energy deposited in that step, in keV
///
/// All hits sharing an event id are coincident (one detector pulse), so the
/// downstream signal simulation builds one SSD Event per Geant4 event.
///
/// The run is serial, so a single RunAction instance owns the stream and no
/// locking is required.
class RunAction : public G4UserRunAction
{
  public:
    RunAction() = default;
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    void RecordHit(G4int eventID, G4double x_mm, G4double y_mm, G4double z_mm,
                   G4double edep_keV);

  private:
    std::ofstream fOut;
    G4long fHitCount = 0;
};

}  // namespace sim

#endif
