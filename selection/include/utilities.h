/**
 * @file utilities.h
 * @brief Header file for definitions of utility functions acting on
 * interactions.
 * @details This file contains definitions of utility functions which are used
 * to support the implementation of analysis variables and cuts. These functions
 * are intended to be used to simplify the implementation of variables and cuts
 * by providing common functionality which can be reused across multiple
 * variables and cuts.
 * @author mueller@fnal.gov
 */
#ifndef UTILITIES_H
#define UTILITIES_H

#define PI0_MASS 134.9768

#include <vector>
#include <array>

#include "framework.h"
#include "include/particle_variables.h"
#include "include/particle_cuts.h"

/**
 * @namespace utilities
 * @brief Namespace for organizing utility functions for supporting analysis
 * variables and cuts.
 * @details This namespace is intended to be used for organizing utility
 * functions which are used to support the implementation of analysis variables
 * and cuts. These functions are intended to be used to simplify the
 * implementation of variables and cuts by providing common functionality which
 * can be reused across multiple variables and cuts.
 * @note The namespace is intended to be used in conjunction with the
 * vars and cuts namespaces, which are used for organizing variables and cuts
 * which act on interactions.
 */
namespace utilities
{
    /**
     * @brief Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to find the topology of.
     * @return the count of primaries of each particle type within the
     * interaction.
     */
    template<class T>
        std::vector<uint32_t> count_primaries(const T & obj)
        {
            std::vector<uint32_t> counts(5, 0);
            for(auto &p : obj.particles)
            {
                if(pcuts::final_state_signal(p))
                    ++counts[pvars::pid(p)];
            }
            return counts;
        }

    /**
     * @brief List of good runs for ICARUS Run 2 (NuMI on-beam).
     * @details This list contains the 92 good runs used by the
     * collaboration for the NuMI Run 2 on-beam data (the same list as the
     * prescaled sample). Data events from runs not in this list are removed
     * by the data quality cut.
     */
    constexpr std::array<unsigned int, 92>
    icarus_good_runs_run2{  9593,  9594,  9595,  9597,  9599,  9602,  9610,  9642,  9646,  9648,
                             9649,  9688,  9690,  9691,  9692,  9693,  9694,  9695,  9696,  9699,
                             9700,  9704,  9715,  9716,  9717,  9721,  9723,  9725,  9726,  9728,
                             9729,  9730,  9731,  9732,  9733,  9735,  9743,  9744,  9745,  9746,
                             9747,  9750,  9752,  9753,  9755,  9762,  9763,  9764,  9765,  9781,
                             9791,  9792,  9794,  9796,  9807,  9834,  9835,  9837,  9838,  9840,
                             9844,  9847,  9849,  9851,  9854,  9855,  9860,  9862,  9870,  9892,
                             9894,  9944,  9945,  9951,  9953,  9954,  9956,  9959,  9960,  9972,
                             9974,  9977,  9979,  9981, 10054, 10059, 10060, 10061, 10064, 10065,
                             10066, 10084};
    
    /**
     * @brief Boolean function to check if a run is in the list of good runs
     * for ICARUS Run 2.
     * @param run the run number to check.
     * @return true if the run is in the list of good runs, false otherwise.
     */
    bool is_icarus_good_run(unsigned int run)
    {
        return std::find(icarus_good_runs_run2.begin(), icarus_good_runs_run2.end(), run) != icarus_good_runs_run2.end();
    }

    /**
     * @brief Utility function to find the index of the first optical flash
     * @details This function returns the index of the first optical flash in
     * the event, shifted by a given amount. It is intended to be used to find
     * the flash nearest to the trigger of the event. The optional shift is
     * used to allow for tuning to correct for the natural offset of the
     * reconstructed flash time from the trigger time, which is not zero
     * despite all systems being referenced to the trigger. This function uses
     * the `firsttime` field of the optical flash.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @param shift the amount to shift the time of the flash by, default is 0.
     * @return size_t the index of the first optical flash in the event, or
     * kNoMatch if no flash is found.
     */
    template<typename T>
    size_t first_opflash_firsttime(const T & sr, double shift=0.0)
    {
        // This variable returns the index of the first optical flash in the
        // event, shifted by a given amount.
        if(sr.opflashes.empty()) return kNoMatch;
        size_t first_flash_index = 0;
        double min_time_diff = std::abs(sr.opflashes[0].firsttime - shift);
        for(size_t i = 1; i < sr.opflashes.size(); ++i)
        {
            double time_diff = std::abs(sr.opflashes[i].firsttime - shift);
            if(time_diff < min_time_diff)
            {
                min_time_diff = time_diff;
                first_flash_index = i;
            }
        }
        return first_flash_index;
    }

    /**
     * @brief Utility function to find the index of the first optical flash
     * @details This function returns the index of the first optical flash in
     * the event, shifted by a given amount. It is intended to be used to find
     * the flash nearest to the trigger of the event. The optional shift is
     * used to allow for tuning to correct for the natural offset of the
     * reconstructed flash time from the trigger time, which is not zero
     * despite all systems being referenced to the trigger. This function uses
     * the `time` field of the optical flash.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @param shift the amount to shift the time of the flash by, default is 0.
     * @return size_t the index of the first optical flash in the event, or
     * kNoMatch if no flash is found.
     */
    template<typename T>
    size_t first_opflash_rawtime(const T & sr, double shift=0.0)
    {
        // This variable returns the index of the first optical flash in the
        // event, shifted by a given amount.
        if(sr.opflashes.empty()) return kNoMatch;
        size_t first_flash_index = 0;
        double min_time_diff = std::abs(sr.opflashes[0].time - shift);
        for(size_t i = 1; i < sr.opflashes.size(); ++i)
        {
            double time_diff = std::abs(sr.opflashes[i].time - shift);
            if(time_diff < min_time_diff)
            {
                min_time_diff = time_diff;
                first_flash_index = i;
            }
        }
        return first_flash_index;
    }

    /**
     * @brief Group true pi0 daughters by parent ID
     * @param obj the interaction from which to count
     * @return an unordered map with parent ID as key and daughter particle indices as value
     */
    std::unordered_map< int, std::vector<size_t> > get_true_pi0s(const caf::SRInteractionTruthDLPProxy & obj, bool primaries = true, std::vector<double> params = {0.0,})
    {
        std::unordered_map< int, std::vector<size_t> > true_pi0s;

	/**
	 * Loop over particles, grouping photons
	 * and electrons that share a parent pi0.
	 */
	for(size_t i(0); i < obj.particles.size(); ++i)
	{
	    const auto & p = obj.particles[i];


	    if(primaries == true and !p.is_primary) continue;
	    if(primaries == false and p.is_primary) continue;

	    if(p.parent_pdg_code == 111 && (p.pdg_code == 22 || p.pdg_code == 11 || p.pdg_code == -11))
	    {
	        true_pi0s[p.parent_track_id].push_back(i);
	    }
	}

	/**
	 * Remove pi0 candidates that are either
	 * subthreshold or have fewer than two daughters.
	 */
	std::vector<size_t> bad_pi0_ids;
        for(auto const & _pi0 : true_pi0s)
        {

	  int num_pi0_daughters(0);
	  TVector3 pi0_momentum(0,0,0);
	  double pi0_ke(0);
	  for(auto & _pidx : _pi0.second)
	    {
	      const auto & _p = obj.particles[_pidx];
	      TVector3 _p_momentum(_p.momentum[0], _p.momentum[1], _p.momentum[2]);
	      pi0_momentum += _p_momentum;
	      num_pi0_daughters++;
	    }
	    pi0_ke = std::sqrt(std::pow(PI0_MASS, 2) + std::pow(pi0_momentum.Mag(), 2)) - PI0_MASS;

	    // TEST: RELAX SIGNAL DEFINITION
	    //if (pi0_ke < params[0]) bad_pi0_ids.push_back(_pi0.first);
            if(num_pi0_daughters < 2 || pi0_ke < params[0]) bad_pi0_ids.push_back(_pi0.first);
	}
	for(size_t i=0; i<bad_pi0_ids.size(); i++)
	{
            true_pi0s.erase(bad_pi0_ids[i]);
	}
	return true_pi0s;
    }

    /**
     * @brief Count the number of true primary neutral pions in interaction.
     * @param obj the interaction from which to count.
     * @param params the neutral pion kinetic energy threshold to consider.
     * @return the number of true primary neutral pions in interaction.
     */
    double true_primary_pi0_multiplicity(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={0.0,})
    {
        std::unordered_map< int, std::vector<size_t> > true_primary_pi0s = get_true_pi0s(obj, true, params);
	return true_primary_pi0s.size();
    }
}
#endif // UTILITIES_H
