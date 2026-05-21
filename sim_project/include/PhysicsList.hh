/// \file PhysicsList.hh
/// \brief Definition of the sim::PhysicsList class

#ifndef SIM_PhysicsList_h
#define SIM_PhysicsList_h 1

#include "G4VModularPhysicsList.hh"

namespace sim
{

/// Physics list for HPGe gamma spectroscopy.
///
/// Registers G4EmStandardPhysics_option4 -- the most accurate EM option, with
/// low-energy models for the photo-electric effect, Compton scattering and
/// Rayleigh scattering, essential for realistic HPGe spectra.
///
/// Only EM physics is needed for a mono-energetic gamma source. To simulate a
/// radioactive isotope instead, also register G4DecayPhysics and
/// G4RadioactiveDecayPhysics in the constructor.
class PhysicsList : public G4VModularPhysicsList
{
  public:
    PhysicsList();
    ~PhysicsList() override = default;
};

}  // namespace sim

#endif
