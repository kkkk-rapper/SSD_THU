/// \file sim_project.cc
/// \brief Main program of the HPGe gamma-spectroscopy simulation.
///
/// Simulates a mono-energetic gamma point source illuminating a high-purity
/// germanium crystal and writes every energy-depositing step inside the
/// crystal to hits.csv. Those hits (interaction position + deposited energy)
/// are the input for the SolidStateDetectors.jl signal simulations in
/// ~/ssd_projects (ppc / ppc_2).

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

using namespace sim;

int main(int argc, char** argv)
{
  // Interactive mode when launched without a macro argument.
  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  // Serial run manager: a single thread, so one RunAction owns one output
  // stream and hits.csv cannot be interleaved/corrupted. See README.md for
  // notes on switching to multithreaded mode.
  auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new PhysicsList());
  runManager->SetUserInitialization(new ActionInitialization());

  auto visManager = new G4VisExecutive(argc, argv);
  visManager->Initialize();

  auto UImanager = G4UImanager::GetUIpointer();

  if (!ui) {
    // Batch mode: execute the supplied macro.
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }
  else {
    // Interactive mode.
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;
}
