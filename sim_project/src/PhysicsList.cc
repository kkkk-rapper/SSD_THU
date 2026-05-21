/// \file PhysicsList.cc
/// \brief Implementation of the sim::PhysicsList class

#include "PhysicsList.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4SystemOfUnits.hh"

namespace sim
{

PhysicsList::PhysicsList()
{
  SetVerboseLevel(1);

  // A 0.1 mm production cut keeps energy deposits well localised, so the hit
  // positions handed to the signal simulation stay sharp. Override at run
  // time with /run/setCut if needed.
  SetDefaultCutValue(0.1 * mm);

  // EM physics only: a mono-energetic gamma source needs no decay processes.
  // To simulate a radioactive isotope source later, also register
  // G4DecayPhysics + G4RadioactiveDecayPhysics here.
  RegisterPhysics(new G4EmStandardPhysics_option4());
}

}  // namespace sim
