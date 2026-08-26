 /**
 * @file cuts.h
 * @brief Header file for definitions of analysis cuts.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions. Each cut is implemented as a function which takes an
 * interaction object as an argument and returns a boolean. These are the
 * building blocks for defining more complex selections.
 * @author mueller@fnal.gov
*/
#ifndef CUTS_H
#define CUTS_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "utilities.h"
#include "framework.h"

/**
 * @namespace cuts
 * @brief Namespace for organizing generic cuts which act on interactions.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions. Each cut is implemented as a function which takes an
 * interaction object as an argument and returns a boolean. The function should
 * be templated on the type of interaction object if the cut is intended to be
 * used on both true and reconstructed interactions.
 */
namespace cuts
{   
    /**
     * @brief Apply a cut on the validity of the flash match.
     * @details A "valid" flash match is defined as a flash-interaction
     * association with a flash time that is not NaN and a flash match
     * status of 1. The upstream flash matching algorithm (OpT0Finder) has a
     * flash filter that restricts candidate flashes to near the beam window,
     * which means that the majority of cosmogenic interactions are not
     * flash matched. If no flash match is found, the flash time is NaN. This
     * cut is intended to be applied as a preselection cut to reduce comparisons
     * to NaN values, which tend to be noisy on stderr.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction on which to place the flash validity cut.
     * @return true if the interaction is flash matched and the time is valid.
     */
    template<class T>
    bool valid_flashmatch(const T & obj)
    {
        return obj.flash_times.size() > 0 && obj.is_flash_matched == 1 && !std::isnan(obj.flash_times[0]);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, valid_flashmatch, valid_flashmatch);

    /**
     * @brief Apply no cut; all interactions passed.
     * @details This is a placeholder function for a cut which does not apply
     * any selection criteria. It is intended to be used in cases where a cut
     * function is required, but no selection is desired.
     * @tparam T the type of object (true or reco).
     * @param obj the interaction to select on.
     * @return true (always).
     */
    template<class T>
    bool no_cut(const T & obj) { return true; }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_cut, no_cut);

    /**
     * @brief Apply a cut to select neutrinos.
     * @details This function applies a cut to select neutrinos. This cut
     * makes use of the is_neutrino flag in the true interaction object and is
     * intended to be used to identify signal neutrinos.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction is a neutrino.
     * @note This cut is intended to be used for identifying neutrinos in
     * truth, which is useful for making signal definitions.
     */
    template<class T>
    bool neutrino(const T & obj) { return obj.nu_id >= 0; }
    REGISTER_CUT_SCOPE(RegistrationScope::True, neutrino, neutrino);

    /**
     * @brief Apply a cut to select cosmogenic interactions.
     * @details This function applies a cut to select cosmogenic interactions.
     * This cut makes use of the is_neutrino flag in the true interaction
     * object and is intended to be used to identify cosmogenic interactions.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction is a cosmogenic interaction.
     * @note This cut is intended to be used for identifying cosmogenic
     * interactions in truth, which is useful for making background definitions.
     */
    template<class T>
    bool cosmic(const T & obj) { return !neutrino(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, cosmic, cosmic);

    /**
     * @brief Apply a cut to select charged current interactions.
     * @details This function applies a cut to select charged current
     * interactions. This cut makes use of the `current_type` attribute in the
     * true interaction object and is intended to be used to identify charged
     * current interactions.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction is a charged current interaction.
     */
    template<class T>
    bool iscc(const T & obj) { return obj.current_type == 0; }
    REGISTER_CUT_SCOPE(RegistrationScope::True, iscc, iscc);

    /**
     * @brief Apply a cut on the interaction mode.
     * @details This function applies a cut to select interactions based on
     * the interaction mode. The interaction mode is stored by Genie as an
     * enumerated category.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this is a vector
     * of interaction modes to select on.
     * @return true if the interaction mode is one of the specified modes.
     */
    template<class T>
    bool is_interaction_mode(const T & obj, std::vector<double> params={})
    {
        if(params.empty())
            return true; // No cut applied if no parameters are given.
        return std::find(params.begin(), params.end(), obj.interaction_mode) != params.end();
    }
    REGISTER_CUT_SCOPE(RegistrationScope::True, is_interaction_mode, is_interaction_mode);

    /**
     * @brief Apply a fiducial volume cut; the interaction vertex must be
     * reconstructed within the fiducial volume.
     * @details The fiducial volume cut is applied on the reconstructed
     * interaction vertex upstream in SPINE. The fiducial volume is defined
     * (in a SPINE post-processor) as a 25 cm border around the x and y
     * detector faces, a 50 cm border around the downstream (+) z face, and a
     * 30 cm border around the upstream (-) z face. The fiducial volume is
     * intended to reduce the impact of detector edge effects on the analysis.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the vertex is in the fiducial volume.
     */
    template<class T>
    bool fiducial_cut(const T & obj)
    {
        return obj.is_fiducial && !(obj.vertex[0] > 210.215 && obj.vertex[1] > 60 && (obj.vertex[2] > 290 && obj.vertex[2] < 390));
//        return obj.is_fiducial
//            && !(obj.vertex[0] > 210.215 && obj.vertex[1] > 60 && (obj.vertex[2] > 290 && obj.vertex[2] < 390))
//            && !(obj.vertex[2] > -100 && obj.vertex[2] < 100);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, fiducial_cut, fiducial_cut);

    /**
     * @brief Veto interactions with true primary neutral pions.
     * @details Counts the number of true primary neutral pions in the interaction
     * using the pi0 utility function, which groups photon and electron daughters
     * by parent track ID and reconstructs the pi0 kinetic energy from its daughters.
     * Pi0 candidates with fewer than two daughters or with kinetic energy below
     * the threshold defined by params[0] are discarded. This cut is applied at
     * truth level only and is intended to align the signal definition with the
     * NUISANCE ICARUS_1muNp0pi_IsSignal flag, which vetoes events containing
     * neutral pions in the final state.
     * @param obj the interaction to select on.
     * @param params pi0 kinetic energy threshold in MeV (default 0.0, no threshold).
     * @return true if no true primary neutral pions are present in the interaction.
     */
    template<class T>
    bool no_pi0s(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params = {0.0,})
    {
        double num_primary_pi0s = utilities::true_primary_pi0_multiplicity(obj, params);
        return num_primary_pi0s == 0;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::True, no_pi0s, no_pi0s);
    
    /**
     * @brief Apply a containment cut on the entire interaction.
     * @details The containment cut is applied on the entire interaction. The
     * interaction is considered contained if all particles and all spacepoints
     * are contained within 5cm of the detector edges (configured in a SPINE 
     * post-processor). Additionally, no spacepoints are allowed to be
     * reconstructed in a TPC that did not create it. This is an unphysical
     * condition that can occur when a cosmic muon is moved according to an
     * assumed t0 that is very out-of-time.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the vertex is contained.
     */
    template<class T>
    bool containment_cut(const T & obj) { return obj.is_contained; }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, containment_cut, containment_cut);

    /**
     * @brief Apply a flash time cut on the interaction.
     * @details The flash time cut is applied on the interaction. The flash time
     * is required to be within the beam window, which is expected to be
     * [0 us, 1.6 us] for BNB and [0 us, 9.6 us] for NuMI. This cut is intended
     * to reduce the impact of cosmogenic interactions on analyses.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut.
     * @return true if the interaction has been matched to an in-time flash.
     * @note The switch to the NuMI beam window is applied by the definition of
     * a preprocessor macro (BEAM_IS_NUMI).
     * @note The cut window has been widened to reconcile the beam window as
     * observed in data and simulation.
     */
    template<class T>
    bool flash_cut(const T & obj, std::vector<double> params={})
    {
        if(!valid_flashmatch(obj))
            return false;
        else if(params.size() == 2 && obj.flash_times[0] >= params[0] && obj.flash_times[0] <= params[1])
            return true;
        else if(params.size() !=2)
            return true;
        else
            return false;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, flash_cut, flash_cut);

    /**
     * @brief Base particle multiplicity cut for a single particle.
     * @details This function applies a cut to select interactions with a
     * multiplicity of 1 for a specific particle type. The particle type is
     * specified by the `particle_species` parameter, which corresponds to the
     * index in the @ref utilities::count_primaries function.
     * @param obj the interaction to select on.
     * @param particle_species the index of the particle species to count.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for the particle to count towards the
     * multiplicity.
     * @return true if the interaction has a multiplicity of 1 for the specified
     * particle species.
     */
    template<class T>
    bool single_particle_multiplicity(const T & obj, size_t particle_species, std::vector<double> params={})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == particle_species && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ++count;
            if(count > 1)
                break; // No need to count further, we only care about multiplicity of 1.
        }
        return count == 1;
    }

    /**
     * @brief Binding for a single particle photon multiplicity cut.
     * @details This function binds the single particle multiplicity cut for
     * photons, which corresponds to the index 0 in the
     * @ref utilities::count_primaries function.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a photon to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return true if the interaction has a single primary photon.
     */
    template<class T>
    bool single_photon(const T & obj, std::vector<double> params={25.0,})
    {
        return single_particle_multiplicity(obj, 0, params);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_photon, single_photon);

    /**
     * @brief Binding for a single particle electron multiplicity cut.
     * @details This function binds the single particle multiplicity cut for
     * electrons, which corresponds to the index 1 in the
     * @ref utilities::count_primaries function.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for an electron to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return true if the interaction has a single primary electron.
     */
    template<class T>
    bool single_electron(const T & obj, std::vector<double> params={25.0,})
    {
        return single_particle_multiplicity(obj, 1, params);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_electron, single_electron);

    /**
     * @brief Binding for a single particle muon multiplicity cut.
     * @details This function binds the single particle multiplicity cut for
     * muons, which corresponds to the index 2 in the
     * @ref utilities::count_primaries function.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a muon to count towards the multiplicity.
     * Defaults to 143.425 MeV, which corresponds to a muon of length 50 cm
     * (assuming the muon stops).
     * @return true if the interaction has a single primary muon.
     */
    template<class T>
    bool single_muon(const T & obj, std::vector<double> params={143.425,})
    {
        return single_particle_multiplicity(obj, 2, params);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_muon, single_muon);

    /**
     * @brief Binding for a single particle pion multiplicity cut.
     * @details This function binds the single particle multiplicity cut for
     * charged pions, which corresponds to the index 3 in the
     * @ref utilities::count_primaries function.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a pion to count towards the multiplicity.
     * Defaults to 25 MeV.
     * @return true if the interaction has a single primary charged pion.
     */
    template<class T>
    bool single_pion(const T & obj, std::vector<double> params={25.0,})
    {
        return single_particle_multiplicity(obj, 3, params);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_pion, single_pion);

    /**
     * @brief Binding for a single particle proton multiplicity cut.
     * @details This function binds the single particle multiplicity cut for
     * protons, which corresponds to the index 4 in the
     * @ref utilities::count_primaries function.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a proton to count towards the multiplicity.
     * Defaults to 50 MeV.
     * @return true if the interaction has a single primary proton.
     */
    template<class T>
    bool single_proton(const T & obj, std::vector<double> params={50.0,})
    {
        return single_particle_multiplicity(obj, 4, params);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_proton, single_proton);

    /**
     * @brief Base particle multiplicity cut for nonzero particle multiplicity.
     * @details This function applies a cut to select interactions with a
     * nonzero multiplicity for a specific particle type. The particle type is
     * specified by the `particle_species` parameter, which corresponds to the
     * index in the @ref utilities::count_primaries function.
     * @param obj the interaction to select on.
     * @param particle_species the index of the particle species to count.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for the particle to count towards the
     * multiplicity.
     * @return true if the interaction has a nonzero multiplicity for the
     * specified particle species.
     */
    template<class T>
    bool nonzero_particle_multiplicity(const T & obj, size_t particle_species, std::vector<double> params={})
    {
        size_t count(0);
        double ke_min = params.size() > 0 ? params[0] : 0.0;
        double ke_max = params.size() > 1 ? params[1] : std::numeric_limits<double>::max();

        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == particle_species && pvars::primary_classification(p)
               && pvars::ke(p) >= ke_min && pvars::ke(p) < ke_max)
                ++count;
            if(count > 0)
                break;
        }
        return count > 0;
    }

    /**
     * @brief Binding for zero particle photon multiplicity cut (negation of
     * nonzero_particle_multiplicity).
     * @details This function binds the nonzero particle multiplicity cut for
     * photons, which corresponds to the index 0 in the
     * @ref utilities::count_primaries function. The negation of this
     * function is used to select interactions with no primary photons.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a photon to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return true if the interaction has a nonzero primary photon.
     */
    template<class T>
    bool no_photons(const T & obj, std::vector<double> params={25.0,})
    {
        return !nonzero_particle_multiplicity(obj, 0, params);
    }

    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_photons, no_photons);

    template<class T>
    bool no_photons_fsi(const T & obj, std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == 0 && pvars::ke(p) >= params[0])
                return false;
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_photons_fsi, no_photons_fsi);

    /**
     * @brief Binding for zero particle electron multiplicity cut (negation of
     * nonzero_particle_multiplicity).
     * @details This function binds the nonzero particle multiplicity cut for
     * electrons, which corresponds to the index 1 in the
     * @ref utilities::count_primaries function. The negation of this
     * function is used to select interactions with no primary electrons.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for an electron to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return true if the interaction has a nonzero primary electron.
     */
    template<class T>
    bool no_electrons(const T & obj, std::vector<double> params={25.0,})
    {
        return !nonzero_particle_multiplicity(obj, 1, params);
    }

    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_electrons, no_electrons);

    template<class T>
    bool no_electrons_fsi(const T & obj, std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == 1 && pvars::ke(p) >= params[0])
                return false;
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_electrons_fsi, no_electrons_fsi);

    /**
     * @brief Binding for zero particle muon multiplicity cut (negation of
     * nonzero_particle_multiplicity).
     * @details This function binds the nonzero particle multiplicity cut for
     * muons, which corresponds to the index 2 in the
     * @ref utilities::count_primaries function. The negation of this
     * function is used to select interactions with no primary muons.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a muon to count towards the
     * multiplicity. Defaults to 143.425 MeV, which corresponds to a muon of
     * length 50 cm (assuming the muon stops).
     * @return true if the interaction has a nonzero primary muon.
     */

    template<class T>
    bool no_muons(const T & obj, std::vector<double> params={143.425,})
    {
        return !nonzero_particle_multiplicity(obj, 2, params);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_muons, no_muons);

    /**
     * @brief Binding for zero particle pion multiplicity cut (negation of
     * nonzero_particle_multiplicity).
     * @details This function binds the nonzero particle multiplicity cut for
     * charged pions, which corresponds to the index 3 in the
     * @ref utilities::count_primaries function. The negation of this
     * function is used to select interactions with no primary charged pions.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a pion to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return true if the interaction has a nonzero primary charged pion.
     */
    template<class T>
    bool no_charged_pions(const T & obj, std::vector<double> params={25.0,})
    {
        return !nonzero_particle_multiplicity(obj, 3, params);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_charged_pions, no_charged_pions);

    /**
     * @brief Binding for zero particle pion multiplicity cut including FSI pions.
     * @details Same as no_charged_pions but vetoes ALL final state charged pions
     * regardless of whether they are primary or secondary (FSI products).
     * @param obj the interaction to select on.
     * @param params kinetic energy threshold, defaults to 25 MeV.
     * @return true if the interaction has no final state charged pions above threshold.
     */
    template<class T>
    bool no_charged_pions_fsi(const T & obj, std::vector<double> params={25.0,})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == 3 && pvars::ke(p) >= params[0])
                ++count;
            if(count > 0)
                break;
        }
        return count == 0;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_charged_pions_fsi, no_charged_pions_fsi);

    /**
     * @brief Cut to select interactions with at least one final state proton
     * above threshold, including FSI secondary protons.
     * @param obj the interaction to select on.
     * @param params kinetic energy threshold, defaults to 50 MeV.
     * @return true if the interaction has at least one final state proton above threshold.
     */
    template<class T>
    bool has_proton_fsi(const T & obj, std::vector<double> params={50.0,})
    {
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == 4 && pvars::ke(p) >= params[0])
                return true;
        }
        return false;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, has_proton_fsi, has_proton_fsi);

    /**
     * @brief Binding for zero particle proton multiplicity cut (negation of
     * nonzero_particle_multiplicity).
     * @details This function binds the nonzero particle multiplicity cut for
     * protons, which corresponds to the index 4 in the
     * @ref utilities::count_primaries function. The negation of this
     * function is used to select interactions with no primary protons.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a proton to count towards the
     * multiplicity. Defaults to 50 MeV.
     * @return true if the interaction has a nonzero primary proton.
     */
    template<class T>
    bool no_protons(const T & obj, std::vector<double> params={50.0,})
    {
        return !nonzero_particle_multiplicity(obj, 4, params);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_protons, no_protons);

    /**
     * @brief Cut to select interactions with a single Michel electron.
     * @details This function applies a cut to select interactions with a
     * single Michel electron.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a single Michel electron.
    */
    template<class T>
    bool single_michel(const T & obj)
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::semantic_type(p) == 2)
                ++count;
            if(count > 1)
                break; // No need to count further, we only care about multiplicity of 1.
        }
        return count == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_michel, single_michel);

/**
     * @brief Cut to select interactions with contained muons only.
     * @details This function applies a cut to select interactions where all
     * muons above the energy threshold are contained within the detector
     * volume. A muon is considered contained if its containment flag is set
     * and it meets the energy and primary particle requirements. This cut
     * is useful for analyses that require accurate momentum reconstruction
     * of muons, which is only possible when the entire track is contained
     * within the detector.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a muon to be considered in the selection.
     * Defaults to 143.425 MeV, which corresponds to a muon of length 50 cm
     * (assuming the muon stops).
     * @return true if all muons above threshold in the interaction are contained.
     */
    template<class T>
    bool contained_muons_only(const T & obj, std::vector<double> params={143.425,})
    {
        for(const auto & p : obj.particles)
        {
            // Check if particle is a muon above threshold and primary
            if(pvars::pid(p) == 2 &&
               pvars::primary_classification(p) &&
               pvars::ke(p) >= params[0])
            {
                // If we find a muon above threshold, check if it's contained
                if(p.is_contained != 1)
                {
                    return false; // Found an uncontained muon, fail the cut
                }
            }
        }
        return true; // All muons above threshold are contained (or no muons found)
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, contained_muons_only, contained_muons_only);

    /**
     * @brief Cut to select interactions with contained protons only.
     * @details This function applies a cut to select interactions where all
     * primary protons above the kinetic energy threshold are contained
     * within the detector volume. A proton is considered contained if its
     * containment flag is set. This cut mirrors @ref contained_muons_only,
     * but for protons: since proton momentum in this analysis is
     * reconstructed from range (CSDA), the reconstruction is only reliable
     * when the full trajectory, from the interaction vertex to the stopping
     * point, is contained within the active volume. Requiring containment
     * therefore removes events where the leading proton candidate exits the
     * detector, which would otherwise bias its reconstructed momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a proton to be considered in the
     * selection. Defaults to 50 MeV, matching the proton threshold used
     * elsewhere in the selection.
     * @return true if all primary protons above threshold in the
     * interaction are contained.
     */
    template<class T>
    bool contained_protons_only(const T & obj, std::vector<double> params={50.0,})
    {
        for(const auto & p : obj.particles)
        {
            // Check if particle is a proton above threshold and primary
            if(pvars::pid(p) == 4 &&
               pvars::primary_classification(p) &&
               pvars::ke(p) >= params[0])
            {
                // If we find a proton above threshold, check if it's contained
                if(p.is_contained != 1)
                {
                    return false; // Found an uncontained proton, fail the cut
                }
            }
        }
        return true; // All protons above threshold are contained (or no protons found)
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, contained_protons_only, contained_protons_only);

    /**
     * @brief Cut to reject events whose leading muon is exiting and has an
     * unreliable MCS energy estimate.
     * @details When the leading muon is not contained, its kinetic energy is
     * taken from the MCS (multiple Coulomb scattering) fit, @ref
     * pvars::mcs_ke, rather than from range (CSDA). The MCS fit is bounded
     * to a maximum kinetic energy of 100 GeV in SPINE, and a sizable
     * fraction of exiting muons in this selection pile up at that bound,
     * indicating the fit failed to converge on a physical value. Even below
     * that hard bound, the MCS bias and resolution both degrade rapidly as
     * a function of the reported kinetic energy (checked against truth in a
     * dedicated study), which biases kinematic imbalance variables built
     * from the muon momentum (e.g. dpT). This cut removes events in which
     * an exiting primary muon's MCS kinetic energy exceeds a configurable
     * ceiling, above which the estimate is no longer considered reliable.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. params[0] sets the maximum
     * allowed MCS kinetic energy (MeV) for an exiting primary muon, and
     * defaults to 4000 MeV (4 GeV). params[1] sets the kinetic energy
     * threshold (MeV) for a muon to be considered in the selection, and
     * defaults to 143.425 MeV, matching the muon threshold used elsewhere
     * in the selection.
     * @return true unless a primary muon above threshold is exiting and its
     * MCS kinetic energy exceeds the configured ceiling.
     */
    template<class T>
    bool exiting_muon_energy_cut(const T & obj, std::vector<double> params={4000.0, 143.425,})
    {
        for(const auto & p : obj.particles)
        {
            // Check if particle is a muon above threshold and primary
            if(pvars::pid(p) == 2 &&
               pvars::primary_classification(p) &&
               pvars::ke(p) >= params[1])
            {
                // Only exiting muons rely on the MCS energy estimate.
                if(p.is_contained != 1 && pvars::mcs_ke(p) > params[0])
                {
                    return false; // Exiting muon with unreliable MCS energy, fail the cut
                }
            }
        }
        return true; // No exiting muon exceeded the MCS energy ceiling (or no muons found)
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, exiting_muon_energy_cut, exiting_muon_energy_cut);

    /**
     * @brief Apply a cut on the pion group code.
     * @details This function applies a cut to select interactions based on
     * the fate of the leading pion. The pion group code classifies what happened
     * to the leading pion in the interaction: DECAY=1 (π→μν), CAPTURE=2 (π⁻ absorbed),
     * INELASTIC=3 (π scattering with hadron production), ELASTIC=4 (π elastic scatter),
     * UNKNOWN=0 (no leading pion or unclear).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this is a vector
     * of pion group codes to select on (e.g., {1.0} for decay only, {1.0, 2.0} for
     * decay or capture).
     * @return true if the pion group code is one of the specified codes.
     * @note This function will be properly linked when variables.h is included elsewhere.
     */
    template<class T>
    bool is_pion_group_code(const T & obj, std::vector<double> params={});

    REGISTER_CUT_SCOPE(RegistrationScope::True, is_pion_group_code, is_pion_group_code);

}
#endif
