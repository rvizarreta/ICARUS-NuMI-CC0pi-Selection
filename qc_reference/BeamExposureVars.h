// B. Howard - 2023
//   howard <at> fnal.gov
// The cuts borrow heavily from NOvA, thanks NOvA!

#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/MultiVar.h"

#include "sbnanaobj/StandardRecord/SRNuMIInfo.h"

namespace ana
{
  double HornCurrentVal ( const caf::SRNuMIInfoProxy& info );
  bool HornCurrentCut ( const caf::SRNuMIInfoProxy& info, const bool stringentMode );

  double POTInSpillVal ( const caf::SRNuMIInfoProxy& info );
  bool POTInSpillCut ( const caf::SRNuMIInfoProxy& info, const bool stringentMode );

  std::pair<double,double> BeamPositionMonitors( const caf::SRNuMIInfoProxy& info );

  double extrapToLoc ( const double var1, const double loc1, const double var2, const double loc2, const double loc3 );
  std::pair<double,double> BeamPositionAtTargetVal( const caf::SRNuMIInfoProxy& info, const unsigned int runNumber );
  bool BeamPositionAtTargetCut ( const caf::SRNuMIInfoProxy& info, const unsigned int runNumber, const bool stringentMode );

  std::pair<double,double> BeamWidthVal ( const caf::SRNuMIInfoProxy& info );
  bool BeamWidthCut ( const caf::SRNuMIInfoProxy& info );

  // For all spills -- only should return info for first event in subrun
  extern const SpillMultiVar kTRTGTDAll;
  extern const SpillMultiVar kTR101DAll;
  extern const SpillMultiVar kHornCurrentAll;
  extern const SpillMultiVar kPOTInSpillAll;
  extern const SpillMultiVar kBeamHPTGT;
  extern const SpillMultiVar kBeamVPTGT;
  extern const SpillMultiVar kBeamPosHAll;
  extern const SpillMultiVar kBeamPosVAll;
  extern const SpillMultiVar kBeamWidthHAll;
  extern const SpillMultiVar kBeamWidthVAll;
  extern const SpillMultiVar kEventAll;

  // For our "triggering" spills
  extern const SpillVar kTRTGTD;
  extern const SpillVar kTR101D;
  extern const SpillVar kHornCurrent;
  extern const SpillVar kPOTInSpill;
  extern const SpillVar kBeamPosH;
  extern const SpillVar kBeamPosV;
  extern const SpillVar kBeamWidthH;
  extern const SpillVar kBeamWidthV;
  extern const SpillVar kDeltaBeamTimeDAQTime;

  extern const SpillVar kRunNumber;
  extern const SpillVar kSubrunNumber;
  extern const SpillVar kEventNumber;

  // Exposure accounting with and without cuts...
  extern const SpillVar kDummyVarForPOTCounting;
  extern const SpillVar kSummedPOT_NuMI_All;
  extern const SpillVar kSummedPOT_NuMI_TRTGTD_All;
  extern const SpillVar kSummedPOT_NuMI_Cuts;
  // Few extras
  extern const SpillVar	kSummedPOT_NuMI_Cuts_TRTGTDgt0p02;
  extern const SpillVar kSummedPOT_NuMI_Cuts_TRTGTDgt0;
  extern const SpillVar kSummedPOT_NuMI_Cuts_TRTGTDall;
  //////////////
  extern const SpillVar kSummedCuts_NuMI_Cuts;

  // NSpills Info
  extern const SpillMultiVar kNuMI_NSpills_Run;
  extern const SpillMultiVar kNuMI_NSpills_Subrun;
  extern const SpillMultiVar kNuMI_NSpills_Event;
  extern const SpillMultiVar kNuMI_NSpills_NSpills;
  // -- and add a few more >> but not sure how to recover these without a lot of work... maybe just using the RUN/EVENT numbers we can figure out where to remove POT?
  //extern const SpillMultiVar kNuMI_NSpills_Time;
  //extern const SpillMultiVar kNuMI_NSpills_SpillTimeS;
  //extern const SpillMultiVar kNuMI_NSpills_SpillTimeNS;
}
