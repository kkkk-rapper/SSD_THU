/// \file ActionInitialization.cc
/// \brief Implementation of the sim::ActionInitialization class

#include "ActionInitialization.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

namespace sim
{

void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction());

  auto runAction = new RunAction();
  SetUserAction(runAction);

  // The stepping action writes hits through the run action's output stream.
  SetUserAction(new SteppingAction(runAction));
}

}  // namespace sim
