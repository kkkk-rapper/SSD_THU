/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the sim::PrimaryGeneratorAction class

#ifndef SIM_PrimaryGeneratorAction_h
#define SIM_PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4GeneralParticleSource;
class G4Event;

namespace sim
{

/// Gamma point-source primary generator.
///
/// Uses G4GeneralParticleSource so the particle, energy, source position and
/// emission geometry are all configured from /gps/... macro commands.
/// See run.mac for the default 122 keV setup.
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

  private:
    G4GeneralParticleSource* fGPS = nullptr;
};

}  // namespace sim

#endif
