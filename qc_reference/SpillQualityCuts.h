// B. Howard - 2023
//   howard <at> fnal.gov
// The cuts borrow heavily from NOvA, thanks NOvA!

#pragma once

#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/SBNAna/Vars/BeamExposureVars.h"

namespace ana
{
  extern const SpillCut kNuMISpillHornCurrentCut;
  extern const SpillCut kNuMISpillPOTCut;
  extern const SpillCut kNuMISpillBeamPosCut;
  extern const SpillCut kNuMISpillBeamWidthCut;

  //extern const SpillCut kNuMIDeltaTimeCut;

  // Combined cut for analyzer to use
  const SpillCut kNuMISpillQualityCut = kNuMISpillHornCurrentCut && kNuMISpillPOTCut && kNuMISpillBeamPosCut && kNuMISpillBeamWidthCut;
}
