#include "sbnana/SBNAna/Cuts/SpillQualityCuts.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRNuMIInfo.h"

namespace ana
{

  const SpillCut kNuMISpillHornCurrentCut ( [](const caf::SRSpillProxy *sr)
  {
    return HornCurrentCut( sr->hdr.spillnumiinfo, true );
  });

  const SpillCut kNuMISpillPOTCut ( [](const caf::SRSpillProxy *sr)
  {
    return POTInSpillCut( sr->hdr.spillnumiinfo, true );
  });

  // CHANGED TO STRINGENT MODE ON 4 APRIL 2025 (BH)
  const SpillCut kNuMISpillBeamPosCut ( [](const caf::SRSpillProxy *sr)
  {
    return BeamPositionAtTargetCut( sr->hdr.spillnumiinfo, sr->hdr.run, true );
  });

  const SpillCut kNuMISpillBeamWidthCut ( [](const caf::SRSpillProxy *sr)
  {
    return BeamWidthCut( sr->hdr.spillnumiinfo );
  });

}
