/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the sim::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4GeneralParticleSource.hh"

namespace sim
{

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  fGPS = new G4GeneralParticleSource();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGPS;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // Particle type, energy, position and emission direction all come from the
  // /gps/... commands in the macro (see run.mac).
  fGPS->GeneratePrimaryVertex(event);
}

}  // namespace sim
