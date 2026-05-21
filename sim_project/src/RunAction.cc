/// \file RunAction.cc
/// \brief Implementation of the sim::RunAction class

#include "RunAction.hh"

#include "G4Run.hh"

namespace sim
{

namespace
{
// Output file, created in the current working directory at /run/beamOn.
// Overwritten by every run -- rename it between runs to keep results.
const char* kOutputFile = "hits.csv";
}  // namespace

void RunAction::BeginOfRunAction(const G4Run*)
{
  fHitCount = 0;
  fOut.open(kOutputFile, std::ios::out | std::ios::trunc);
  fOut << "event,x_mm,y_mm,z_mm,edep_keV\n";
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  fOut.close();
  G4cout << "\n--------------------------------------------------------------\n"
         << " Run finished: " << run->GetNumberOfEvent() << " primary gammas\n"
         << " Recorded " << fHitCount << " hits -> " << kOutputFile << "\n"
         << "--------------------------------------------------------------\n"
         << G4endl;
}

void RunAction::RecordHit(G4int eventID, G4double x_mm, G4double y_mm,
                          G4double z_mm, G4double edep_keV)
{
  fOut << eventID << ',' << x_mm << ',' << y_mm << ',' << z_mm << ','
       << edep_keV << '\n';
  ++fHitCount;
}

}  // namespace sim
