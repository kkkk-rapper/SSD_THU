/// \file ActionInitialization.hh
/// \brief Definition of the sim::ActionInitialization class

#ifndef SIM_ActionInitialization_h
#define SIM_ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

namespace sim
{

/// Registers the user action classes.
///
/// Only Build() is defined: the project runs serially (see sim_project.cc),
/// so there is no master/worker split and BuildForMaster() is not needed.
class ActionInitialization : public G4VUserActionInitialization
{
  public:
    void Build() const override;
};

}  // namespace sim

#endif
