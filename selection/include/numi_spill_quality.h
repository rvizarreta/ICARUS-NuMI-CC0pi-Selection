/**
 * @file numi_spill_quality.h
 * @brief NuMI spill quality cuts (spill POT, horn current, beam position at
 * the target, beam width) for ICARUS.
 * @details Verbatim copy of the cut logic of B. Howard's BeamExposureVars
 * (SBNAna), which is not part of sbnana v10_01_02_01. Only the numerical
 * cut code is kept. The only changes are that the functions are inline
 * (header-only), live in the numi_qc namespace, and the ROOT objects used
 * for the beam width fit are owned by unique_ptr instead of being leaked.
 * Reference copy of the original: qc_reference/BeamExposureVars.cxx.
 */
#pragma once

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TGraph.h"
#include "TF1.h"
#include "TMath.h"

#include <cmath>
#include <limits>
#include <memory>
#include <utility>

namespace numi_qc
{
  // Horn Current:
  // NB: current corrections from old NOvA code
  //////////////////////////////////////////////////////////////
  inline double HornCurrentVal ( const caf::SRNuMIInfoProxy& info )
  {
    return ( ((info.NSLINA - 0.01)/0.9951) + ((info.NSLINB - (-0.14))/0.9957) + ((info.NSLINC - (-0.05))/0.9965) + ((info.NSLIND - (-0.07))/0.9945) );
  }

  inline bool HornCurrentCut ( const caf::SRNuMIInfoProxy& info, const bool stringentMode )
  {
    double iHorn = HornCurrentVal(info);

    if ( stringentMode ) return (iHorn > -202. && iHorn <= -196.4);
    return (iHorn > -205. && iHorn <= -195.);
  }

  // POT in spill: NOvA style
  // Try TRTGTD first, switch to TR101D if it's too small, then pin to 0 if necessary
  //////////////////////////////////////////////////////////////
  inline double POTInSpillVal ( const caf::SRNuMIInfoProxy& info )
  {
    double pot = info.TRTGTD;
    // TODO: 0.02 or 0.02e13?
    if ( info.TRTGTD < 0.02 ) {
      pot = info.TR101D;
      if ( pot < 0. ) pot = 0.;
    }

    return pot;
  }

  inline bool POTInSpillCut ( const caf::SRNuMIInfoProxy& info, const bool stringentMode )
  {
    double pot = POTInSpillVal(info);

    if ( stringentMode ) return pot > 2.0e12;
    return pot > 0.5e12;
  }

  // Some weighted averages of beam position monitors
  //////////////////////////////////////////////////////////////
  // NOTE: this function does averages THEN projects, I believe NOvA projects then averages. Maybe that's better...
  inline std::pair<double,double> BeamPositionMonitors ( const caf::SRNuMIInfoProxy& info )
  {
    // If any monitor doesn't have 7 then return -999, -999
    if ( info.HPTGT.size() != 7 || info.HITGT.size() != 7 ||
         info.VPTGT.size() != 7 || info.VITGT.size() != 7 ) {
      std::pair<double,double> ret = std::make_pair<double,double>(-999.,-999.);
      return ret;
    }

    double ave_HPTGT(0.), ave_VPTGT(0.), wgt_HPTGT(0.), wgt_VPTGT(0.);
    for ( unsigned int idxSpillPart=0; idxSpillPart<7; ++idxSpillPart ) {
      ave_HPTGT+=info.HPTGT[idxSpillPart]*info.HITGT[idxSpillPart];
      wgt_HPTGT+=info.HITGT[idxSpillPart];
      ave_VPTGT+=info.VPTGT[idxSpillPart]*info.VITGT[idxSpillPart];
      wgt_VPTGT+=info.VITGT[idxSpillPart];
    }

    if ( wgt_HPTGT < std::numeric_limits<double>::epsilon() ||
         wgt_VPTGT < std::numeric_limits<double>::epsilon() ) {
      std::pair<double,double> ret = std::make_pair<double,double>(-999.,-999.);
      return ret;
    }

    std::pair<double,double> ret = std::make_pair<double,double>( ave_HPTGT/wgt_HPTGT, ave_VPTGT/wgt_VPTGT );
    return ret;
  }

  // Beam Pos at target
  //////////////////////////////////////////////////////////////
  inline double extrapToLoc ( const double var1, const double loc1, const double var2, const double loc2, const double loc3 )
  {
    double slope = (var2-var1)/(loc2-loc1);
    double proj = var1 + ((loc3-loc1)*slope);
    return proj;
  }

  // NOTE: this function does averages THEN projects, I believe NOvA projects then averages. Maybe that's better...
  inline std::pair<double,double> BeamPositionAtTargetVal ( const caf::SRNuMIInfoProxy& info, const unsigned int runNumber )
  {
    double mmPerFoot = 12.*2.54*10.; // inch per foot * cm per inch * mm per cm

    // From NOvA, converted to mm
    double z_hp121 = -68.04458 * mmPerFoot;
    double z_vp121 = -66.99283 * mmPerFoot;
    double z_hptgt = -31.25508 * mmPerFoot;
    double z_vptgt = -30.16533 * mmPerFoot;
    double z_targ  = 0.; 

    // If any monitor doesn't have 7 then return -999, -999
    if ( info.HP121.size() != 7 || info.VP121.size() != 7 || info.HPTGT.size() != 7 ||
         info.HITGT.size() != 7 || info.VPTGT.size() != 7 || info.VITGT.size() != 7 ) {
      std::pair<double,double> ret = std::make_pair<double,double>(-999.,-999.);
      return ret;
    }

    double ave_HP121(0.), ave_VP121(0.), ave_HPTGT(0.), ave_VPTGT(0.), wgt_HPTGT(0.), wgt_VPTGT(0.);
    for ( unsigned int idxSpillPart=0; idxSpillPart<7; ++idxSpillPart ) {
      ave_HP121+=info.HP121[idxSpillPart];
      ave_VP121+=info.VP121[idxSpillPart];
      ave_HPTGT+=info.HPTGT[idxSpillPart]*info.HITGT[idxSpillPart];
      wgt_HPTGT+=info.HITGT[idxSpillPart];
      ave_VPTGT+=info.VPTGT[idxSpillPart]*info.VITGT[idxSpillPart];
      wgt_VPTGT+=info.VITGT[idxSpillPart];
    }

    ave_HP121/=7.;
    ave_VP121/=7.;
    if ( wgt_HPTGT < std::numeric_limits<double>::epsilon() ||
         wgt_VPTGT < std::numeric_limits<double>::epsilon() ) {
      std::pair<double,double> ret = std::make_pair<double,double>(-999.,-999.);
      return ret;
    }
    ave_HPTGT/=wgt_HPTGT;
    ave_VPTGT/=wgt_VPTGT;

    // TODO: double check that these offsets make sense...
    // Based on values provided by NuMI Target group. Thanks!
    double xCorr(-9999.);
    double yCorr(-9999.);
    // Run 1 offsets -- for now apply same cuts to everything before Run1 as well
    if ( runNumber <= 8553 ) {
      double x_nom121(0.), x_nomtgt(0.), y_nom121(-1.5), y_nomtgt(-1.5), x_nomTarg(0.398), y_nomTarg(-0.39);
      
      xCorr = extrapToLoc(x_nom121,z_hp121,x_nomtgt,z_hptgt,z_targ) - x_nomTarg;
      yCorr = extrapToLoc(y_nom121,z_hp121,y_nomtgt,z_hptgt,z_targ) - y_nomTarg;
    }
    // Run 2 offsets -- will have to edit this to give the max run at the end of Run 2, or when conditions change.
    else if ( runNumber > 8553 ) {
      double x_nom121(1.2), x_nomtgt(0.7), y_nom121(-0.4), y_nomtgt(-0.94), x_nomTarg(0.03), y_nomTarg(-0.59);
      
      xCorr = extrapToLoc(x_nom121,z_hp121,x_nomtgt,z_hptgt,z_targ) - x_nomTarg;
      yCorr = extrapToLoc(y_nom121,z_hp121,y_nomtgt,z_hptgt,z_targ) - y_nomTarg;
    }

    std::pair<double,double> ret = std::make_pair<double,double>( extrapToLoc(ave_HP121,z_hp121,ave_HPTGT,z_hptgt,z_targ) - xCorr,
                                                                  extrapToLoc(ave_VP121,z_vp121,ave_VPTGT,z_vptgt,z_targ) - yCorr );

    return ret;
  }

  inline bool BeamPositionAtTargetCut ( const caf::SRNuMIInfoProxy& info, const unsigned int runNumber, const bool stringentMode )
  {
    std::pair<double,double> beamPos = BeamPositionAtTargetVal(info,runNumber);

    if ( stringentMode ) {
      return fabs(beamPos.first) <= 1. && fabs(beamPos.second) <= 1.;
    }

    return fabs(beamPos.first) <= 2. && fabs(beamPos.second) <= 2.;
  }

  // Beam Widths
  // This also motivated by/borrows what is done in the NOvA modules
  // is it [103,103+48) or (103,103+48] ? -- I think the first?
  //////////////////////////////////////////////////////////////
  inline std::pair<double,double> BeamWidthVal ( const caf::SRNuMIInfoProxy& info )
  {
    if ( info.MTGTDS.size() < 199 ) {
      std::pair<double,double> ret = std::make_pair<double,double>(-999.,-999.);
      return ret;
    }

    double x_pos[48];
    double h_chan[48];
    double v_chan[48];

    for ( unsigned int idxX=0; idxX<48; ++idxX ) {
      x_pos[ idxX ] = 0.5 * (idxX+1-24.5);
    }

    double h_sum = 0.;
    double h_min = 99999;
    for ( unsigned int idxH=0; idxH<48; ++idxH ) {
      h_chan[ idxH ] = -1.*info.MTGTDS.at(103+idxH);
      h_sum += h_chan[ idxH ];
      if(h_chan[ idxH ]<h_min){
        h_min = h_chan[ idxH ];
      }
    }

    double v_sum = 0.;
    double v_min = 99999;
    for ( unsigned int idxV=0; idxV<48; ++idxV ) {
      v_chan[ idxV ] = -1.*info.MTGTDS.at(151+idxV);
      v_sum += v_chan[ idxV ];
      if(v_chan[ idxV ]<v_min){
        v_min = v_chan[ idxV ];
      }
    }

    // H width
    std::unique_ptr<TGraph> hGraph(new TGraph(48, x_pos, h_chan));
    std::unique_ptr<TF1> hGauss(new TF1("hGauss","[0] + [1]*TMath::Exp(-((x-[2])*(x-[2]))/(2*[3]*[3]))"));
    hGauss->SetParameter(0,-1.*h_sum/48.);
    hGauss->SetParameter(1,-1.*h_min);
    hGauss->SetParameter(2,0.);
    hGauss->SetParameter(3,2.);
    hGauss->SetParLimits(3, 0., 10.);
/*
    // Remove 0s
    std::vector<int> to_remove_h;
    for ( int idx=0; idx<hGraph->GetN(); ++idx ) {
      double x, y;
      hGraph->GetPoint( idx, x, y );
      if ( fabs(y) < std::numeric_limits<double>::epsilon() ) to_remove_h.push_back( idx );
    }
    for ( unsigned int idx=0; idx<to_remove_h.size(); ++idx ) hGraph->RemovePoint( to_remove_h.at( to_remove_h.size()-1-idx ) );
*/
    hGraph->Fit(hGauss.get(), "Q");

    // V width
    std::unique_ptr<TGraph> vGraph(new TGraph(48, x_pos, v_chan));
    std::unique_ptr<TF1> vGauss(new TF1("vGauss","[0] + [1]*TMath::Exp(-((x-[2])*(x-[2]))/(2*[3]*[3]))"));
    vGauss->SetParameter(0,-1.*v_sum/48.);
    vGauss->SetParameter(1,-1.*v_min);
    vGauss->SetParameter(2,0.);
    vGauss->SetParameter(3,2.);
    vGauss->SetParLimits(3, 0., 10.);
/*
    // Remove 0s
    std::vector<int> to_remove_v;
    for ( int idx=0; idx<vGraph->GetN(); ++idx ) {
      double x, y;
      vGraph->GetPoint( idx, x, y );
      if ( fabs(y) < std::numeric_limits<double>::epsilon() ) to_remove_v.push_back( idx );
    }
    for ( unsigned int idx=0; idx<to_remove_v.size(); ++idx ) vGraph->RemovePoint( to_remove_v.at( to_remove_v.size()-1-idx ) );
*/

    vGraph->Fit(vGauss.get(),"Q");

    std::pair<double,double> ret = std::make_pair<double,double>(hGauss->GetParameter(3),vGauss->GetParameter(3));


    return ret;
  }

  inline bool BeamWidthCut ( const caf::SRNuMIInfoProxy& info )
  {
    std::pair<double,double> beamWidths = BeamWidthVal(info);

    return beamWidths.first > 0.57 && beamWidths.first <= 1.88 &&
           beamWidths.second > 0.57 && beamWidths.second <= 1.88;
  }

}
