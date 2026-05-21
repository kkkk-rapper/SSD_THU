/// \file SteppingAction.hh
/// \brief Definition of the sim::SteppingAction class

#ifndef SIM_SteppingAction_h
#define SIM_SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4LogicalVolume;

namespace sim
{

class RunAction;

/// Records every energy-depositing step inside the germanium crystal,
/// forwarding (position, energy, time) to the RunAction output stream.
class SteppingAction : public G4UserSteppingAction
{
  public:
    explicit SteppingAction(RunAction* runAction);
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step* step) override;

  private:
    RunAction* fRunAction = nullptr;
    G4LogicalVolume* fCrystalVolume = nullptr;
};

}  // namespace sim

#endif
