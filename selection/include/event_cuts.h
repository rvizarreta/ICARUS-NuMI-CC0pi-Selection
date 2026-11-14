/**
 * @file event_cuts.h
 * @brief Definitions of analysis cuts which can be applied to the
 * StandardRecord object.
 * @details This file contains definitions of analysis cuts which can be
 * applied to the StandardRecord object. Each cut is implemented as a
 * function which takes a StandardRecord object as an argument and returns
 * a boolean indicating whether the event passes the cut.
 * @author mueller@fnal.gov
 */
#ifndef EVENT_CUTS_H
#define EVENT_CUTS_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "numi_spill_quality.h"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "framework.h"
#include "utilities.h"
#include "event_variables.h"

/**
 * @namespace ecut
 * @brief Namespace for organizing cuts which act on events.
 * @details This namespace is intended to be used for organizing cuts
 * which act on events. Each cut is implemented as a function which takes
 * a StandardRecord object as an argument and returns a boolean.
 */
namespace ecut
{
    /**
     * @brief Apply no cut to the event.
     * @details This cut always returns true, indicating that the event passes
     * the cut. It is useful for testing purposes or when no cut is desired.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @return true, indicating that the event passes the cut.
     */
    template<typename T>
    bool no_cut(const T & sr) { return true; }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, no_cut, no_cut);

    /**
     * @brief Apply a cut enforcing that there are nonzero reco interactions.
     * @details This cut checks if the number of reco interactions in the event
     * is greater than zero. It is useful for ensuring that the file isn't
     * corrupted or empty.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @return true if there are nonzero reco interactions, false otherwise.
     */
    template<typename T>
    bool nonzero_reco_interactions(const T & sr)
    {
        // Check if the number of reco interactions is greater than zero.
        return (sr.ndlp > 0);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, nonzero_reco_interactions, nonzero_reco_interactions);

    /**
     * @brief Apply a cut enforcing a CRT-PMT veto.
     * @details This cut rejects events that do not have at least one flash in-
     * time with the beam that is not associated with a CRT hit.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @param params a vector of parameters for the cut, which can be used to
     * specify the time window for the flash association.
     * @return true if the event has at least one flash in-time with the beam
     * that is not associated with a CRT hit, false otherwise.
     */
    template<typename T>
    bool crtpmt_veto(const T & sr, std::vector<double> params={})
    {
        if(params.size() < 2)
        {
            throw std::invalid_argument("crtpmt_veto requires at least two parameters: time window start and end.");
        }
        bool crtpmt_matched(sr.ncrtpmt_matches == 0);
        for(auto const & c : sr.crtpmt_matches)
        {
            if(c.flashGateTime > params[0] && c.flashGateTime < params[1] && c.flashClassification == 0)
                crtpmt_matched = true;
        }
        return crtpmt_matched;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, crtpmt_veto, crtpmt_veto);

    /**
     * @brief Apply a cut on the global trigger time.
     * @details This cut checks if the global trigger time is within some
     * interval of time. This is useful for partitioning the dataset into
     * subsets based on the time of the global trigger. If no time interval is
     * specified, the cut will always return true, allowing all events to pass.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @return true if the global trigger time is within the specified interval,
     */
    template<typename T>
    bool global_trigger_time_cut(const T & sr, std::vector<double> params={})
    {
        if(params.empty())
        {
            // If no parameters are provided, we assume the cut passes for all events.
            return true;
        }
        else if(params.size() == 2)
        {
            // Check if the global trigger time is within the specified interval.
            double global_trigger_time = sr.hdr.triggerinfo.global_trigger_time;
            return (global_trigger_time >= params[0] && global_trigger_time <= params[1]);
        }
        else
        {
            throw std::invalid_argument("global_trigger_time_cut requires either no parameters or exactly two parameters.");
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, global_trigger_time_cut, global_trigger_time_cut);

    /**
     * @brief Apply a data quality cut on the event metadata.
     * @details This cut rejects events that belong to runs that are marked as
     * bad. This can happen due to various reasons, such as hardware issues,
     * non-physics configuration, or other data quality concerns.
     * @param sr the StandardRecord to apply the cut on.
     * @return true if the event is from a good run, false otherwise.
     */
    template<typename T>
    bool data_quality_cut(const T & sr)
    {
        if(sr.hdr.det == caf::kICARUS && !sr.hdr.ismc)
        {
            return utilities::is_icarus_good_run(sr.hdr.run);
        }
        else
        {
            // For other detectors, we assume all runs are good.
            // This can be extended to include other detectors' run quality checks.
            return true;
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, data_quality_cut, data_quality_cut);

    /**
     * @brief Apply the NuMI beam quality cuts to the spill that triggered the
     * event.
     * @details This cut rejects events whose own spill (hdr.spillnumiinfo)
     * fails any of the four NuMI spill quality cuts of the mesonless analysis
     * (Howard, SBNAna BeamExposureVars): spill POT, horn current, beam
     * position at the target (stringent mode, run dependent) and beam width.
     * The cuts are evaluated in this order and the beam width is only
     * evaluated if the first three pass. It is applied on ICARUS data only;
     * off-beam events (no NuMI spill, noffbeamnumi > 0) pass.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @return true if the spill passes all four cuts or the event is MC,
     * false otherwise.
     */
    template<typename T>
    bool numi_beam_quality_cut(const T & sr)
    {
        if(sr.hdr.det == caf::kICARUS && !sr.hdr.ismc)
        {
            // Off-beam events have no NuMI spill (empty numiinfo) and are
            // counted through noffbeamnumi (see the livetime in framework.cc),
            // so there is no beam to check: they pass. Every other data
            // event is cut, including an on-beam event with no spill record.
            if(sr.hdr.numiinfo.size() == 0 && sr.hdr.noffbeamnumi > 0) return true;

            // The beam width needs two ROOT fits and this cut is evaluated
            // once per tree (and per exposure tree) for every event, so the
            // result is cached per event.
            static std::unordered_map<uint64_t, bool> cache;
            const uint64_t key = ((uint64_t)sr.hdr.run << 44) ^
                                 ((uint64_t)sr.hdr.subrun << 24) ^
                                 (uint64_t)sr.hdr.evt;
            const auto hit = cache.find(key);
            if(hit != cache.end()) return hit->second;

            const auto & info = sr.hdr.spillnumiinfo;
            const bool pass = numi_qc::POTInSpillCut(info, true) &&
                              numi_qc::HornCurrentCut(info, true) &&
                              numi_qc::BeamPositionAtTargetCut(info, sr.hdr.run, true) &&
                              numi_qc::BeamWidthCut(info);
            cache.emplace(key, pass);
            return pass;
        }
        else
        {
            return true;
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, numi_beam_quality_cut, numi_beam_quality_cut);

    /**
     * @brief Apply the trigger quality cut using an external list of bad
     * events.
     * @details The trigger quality check compares, for consecutive majority
     * gates, the change of the NuMI gate counter with the number of IFBeam
     * spills between the gates. This needs the trigger database and the
     * IFBeam spill times, which are not in the flat CAFs, so the check is done
     * outside of medulla and written to a text file with one "run event" pair
     * per line for every event that fails it. This cut rejects the events in
     * that file. The file path is read from the environment variable
     * MEDULLA_TRIGGER_BAD_LIST. It is applied on ICARUS data only, and the cut
     * throws if the variable is not set or the file cannot be read, so it can
     * never pass silently without the list. The join is on (run, event),
     * never on subrun.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @return true if the event is not in the bad-event list or is MC, false
     * otherwise.
     */
    template<typename T>
    bool trigger_quality_cut(const T & sr)
    {
        if(sr.hdr.det == caf::kICARUS && !sr.hdr.ismc)
        {
            static const std::unordered_set<uint64_t> bad_events = []() {
                const char * path = std::getenv("MEDULLA_TRIGGER_BAD_LIST");
                if(path == nullptr)
                    throw std::runtime_error("trigger_quality_cut: MEDULLA_TRIGGER_BAD_LIST is not set.");
                std::ifstream in(path);
                if(!in)
                    throw std::runtime_error(std::string("trigger_quality_cut: cannot read ") + path);
                std::unordered_set<uint64_t> events;
                uint64_t run, evt;
                while(in >> run >> evt)
                    events.insert((run << 32) | evt);
                if(events.empty())
                    throw std::runtime_error(std::string("trigger_quality_cut: no events read from ") + path);
                return events;
            }();
            const uint64_t key = ((uint64_t)sr.hdr.run << 32) | (uint64_t)sr.hdr.evt;
            return bad_events.find(key) == bad_events.end();
        }
        else
        {
            return true;
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, trigger_quality_cut, trigger_quality_cut);

    /**
     * @brief Cut that implements trigger emulation on simulated data.
     * @details This cut emulates the trigger logic on simulated data by
     * checking if the event has a valid trigger time.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @return true if the event has a valid trigger time, false otherwise.
     */
    template<typename T>
    bool trigger_emulation_cut(const T & sr)
    {
        if(sr.hdr.ismc && std::isnan(sr.hdr.triggerinfo.global_trigger_det_time))
        {
            // If the event is simulated and the global trigger time is NaN,
            // we assume the trigger emulation cut fails.
            return false;
        }
        else
        {
            // For real data or if the trigger time is valid, we pass the cut.
            return true;
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, trigger_emulation_cut, trigger_emulation_cut);

    /**
     * @brief A cut that places a threshold on the BNB Figure of Merit 2 (FoM2).
     * @details This cut checks if the FoM2 value is above a certain threshold,
     * which can be configured by the user. This is intended to be a handle on
     * the quality of the beam-target overlap of the spill that is associated
     * with the event
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @param params a vector of parameters for the cut, which can be used to
     * specify the threshold for the FoM2 value.
     * @return true if the FoM2 value is above the threshold, false otherwise.
     */
    /*template<typename T>
    bool bnb_fom2_cut(const T & sr, std::vector<double> params={})
    {
        if(params.empty())
        {
            throw std::invalid_argument("bnb_fom2_cut requires at least one parameter for the threshold (recommended 0.98).");
        }
        if(sr.hdr.ismc)
        {
            // If the event is simulated, do not apply the cut.
            return true;
        }
        else
        {
            double threshold = params[0];
            return (evar::bnb_fom2(sr) >= threshold);
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, bnb_fom2_cut, bnb_fom2_cut);*/

    /**
     * @brief A cut that checks if the event is the first in the subrun.
     * @details This cut checks if the event is the first event in the subrun.
     * This can be relevant for certain attributes (e.g. spill info) that is
     * stored only once per subrun in the first event.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @return true if the event is the first in the subrun, false otherwise.
     */
    template<typename T>
    bool is_first_in_subrun_cut(const T & sr)
    {
        // This cut checks if the event is the first in the subrun.
        // It returns true if it is the first event in the subrun, false otherwise.
        return sr.hdr.first_in_subrun;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, is_first_in_subrun_cut, is_first_in_subrun_cut);
}

#endif // EVENT_CUTS_H