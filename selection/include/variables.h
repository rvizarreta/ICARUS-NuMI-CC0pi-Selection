/**
 * @file variables.h
 * @brief Header file for definitions of analysis variables.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from interactions. Each variable is implemented
 * as a function which takes an interaction object as an argument and returns a
 * double. These are the building blocks for producing high-level plots of the
 * selected interactions.
 * @author mueller@fnal.gov
*/
#ifndef VARIABLES_H
#define VARIABLES_H
#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813
#define NUCLEON_MASS 938.9187473

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"
#include "sbnanaobj/StandardRecord/Proxy/EpilogFwd.h"

#include "include/particle_variables.h"
#include "include/particle_cuts.h"
#include "include/cuts.h"
#include "include/utilities.h"
#include "include/particle_utilities.h"
#include "include/selectors.h"
#include "framework.h"

/**
 * @namespace vars
 * @brief Namespace for organizing generic variables which act on interactions.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions. Each variable is implemented as a function which takes
 * an interaction object as an argument and returns a double. The function
 * should be templated on the type of interaction object if the variable is
 * intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the
 * pvars namespace, which is used for organizing variables which act on single
 * particles.
 */
namespace vars
{
    /**
     * @brief Variable for the neutrino ID of the interaction.
     * @details This variable is intended to provide a unique identifier for
     * each parent neutrino within the event record. This number is assigned
     * starting at 0 for the first neutrino in the event and is incremented
     * for each subsequent neutrino. Non-neutrino interactions are assigned
     * a value of -1.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the neutrino ID.
     */
    template<class T>
    double neutrino_id(const T & obj) { return obj.nu_id; }
    REGISTER_VAR_SCOPE(RegistrationScope::True, neutrino_id, neutrino_id);

    /**
     * @brief Variable for the interaction ID.
     * @details This variable is intended to provide a unique identifier for
     * each interaction within the event record. This number is assigned
     * starting at 0 for the first interaction in the event and is incremented
     * for each subsequent interaction. This assignment is done upstream in the
     * SPINE reconstruction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the interaction ID.
     */
    template<class T>
    double interaction_id(const T & obj) { return obj.id; }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, interaction_id, interaction_id);

    /**
     * @brief Variable for the best-match IoU of the interaction.
     * @details The best-match IoU is the intersection over union of the
     * points belonging to a pair of reconstructed and true interactions. The
     * best-match IoU is calculated upstream in the SPINE reconstruction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the best-match IoU of the interaction.
     */
    template<class T>
    double iou(const T & obj)
    {
        if(obj.match_ids.size() > 0)
            return obj.match_overlaps[0];
        else 
            return PLACEHOLDERVALUE;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, iou, iou);

    /**
     * @brief Variable for the containment status of the interaction.
     * @details The containment status is determined upstream in the SPINE
     * reconstruction and is based on the set of all points in the interaction,
     * which must be contained within the volume of the TPC that created them.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the containment status of the interaction.
     */
    template<class T>
    double containment(const T & obj) { return cuts::containment_cut(obj); }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, containment, containment);

    /**
     * @brief Variable for the fiducial volume status of the interaction.
     * @details The fiducial volume status is determined upstream in the SPINE
     * reconstruction and is a requirement that the interaction vertex is within
     * the fiducial volume of the TPC.
     */
    template<class T>
    double fiducial(const T & obj) { return cuts::fiducial_cut(obj); }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, fiducial, fiducial);

    /**
     * @brief Variable for total visible energy of interaction.
     * @details This function calculates the total visible energy of the
     * interaction by summing the energy of all particles that are identified
     * as counting towards the final state of the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj interaction to apply the variable on.
     * @return the total visible energy of the interaction.
     */
    template<class T>
    double visible_energy(const T & obj)
    {
        double energy(0);
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                energy += pvars::energy(p);
                if(pvars::pid(p) == pvars::kProton) energy -= pvars::mass(p) - PROTON_BINDING_ENERGY;
            }
        }
        return energy/1000.0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, visible_energy, visible_energy);

    /**
     * @brief Variable for energy reconstruction assuming CCQE kinematics using
     * the lepton.
     * @details This function calculates the neutrino energy assuming CCQE
     * kinematics using the leading lepton in the interaction. The leading
     * lepton is defined as the highest kinetic energy electron or muon in the
     * interaction. If no electron or muon is found, the function returns
     * PLACEHOLDERVALUE. This does not check that the interaction is actually
     * QE-like.
     * @tparam T the type of interaction (true or reco).
     * @param obj interaction to apply the variable on.
     * @return the reconstructed neutrino energy assuming CCQE kinematics.
     */
    template<class T>
    double energy_qel(const T & obj)
    {
        size_t ei = selectors::leading_electron(obj);
        size_t mi = selectors::leading_muon(obj);
        size_t li = kNoMatch;

        if(ei != kNoMatch && mi != kNoMatch)
            li = (pvars::ke(obj.particles[ei]) > pvars::ke(obj.particles[mi])) ? ei : mi;
        else if(ei != kNoMatch)
            li = ei;
        else if(mi != kNoMatch)
            li = mi;
        else
            return PLACEHOLDERVALUE;

        double Mn = 939.565;
        double Mp = 938.272;
        double Ml = (li == ei) ? ELECTRON_MASS : MUON_MASS;
        double EB = PROTON_BINDING_ENERGY;

        double El = pvars::energy(obj.particles[li]);
        double pz = pvars::pz(obj.particles[li]);

        double numerator   = 2*(Mn - EB)*El - ((Mn - EB)*(Mn - EB) + Ml*Ml - Mp*Mp);
        double denominator = 2*((Mn - EB) - El + pz);

        return (numerator / denominator) / 1000.0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, energy_qel, energy_qel);

    /**
     * @brief Variable for energy reconstruction assuming CCQE kinematics using
     * the proton.
     * @details This function calculates the neutrino energy assuming CCQE
     * kinematics using the leading proton in the interaction. The leading
     * proton is defined as the highest kinetic energy proton in the
     * interaction. If no proton is found, the function returns 
     * PLACEHOLDERVALUE. This does not check that the interaction is actually
     * QE-like.
     * @tparam T the type of interaction (true or reco).
     * @param obj interaction to apply the variable on.
     * @return the reconstructed neutrino energy assuming CCQE kinematics.
     */
    template<class T>
    double energy_qep(const T & obj)
    {
        size_t ei = selectors::leading_electron(obj);
        size_t mi = selectors::leading_muon(obj);
        size_t pi = selectors::leading_proton(obj);
        size_t li = kNoMatch;

        if(ei != kNoMatch && mi != kNoMatch)
            li = (pvars::ke(obj.particles[ei]) > pvars::ke(obj.particles[mi])) ? ei : mi;
        else if(ei != kNoMatch)
            li = ei;
        else if(mi != kNoMatch)
            li = mi;
        else
            return PLACEHOLDERVALUE;

        if(pi == kNoMatch)
            return PLACEHOLDERVALUE;

        double Mn = 939.565;
        double Mp = 938.272;
        double Ml = (li == ei) ? ELECTRON_MASS : MUON_MASS;
        double EB = PROTON_BINDING_ENERGY;

        double Ep = pvars::energy(obj.particles[pi]);
        double pz = pvars::pz(obj.particles[pi]);

        double numerator   = 2*(Mn - EB)*Ep - ((Mn - EB)*(Mn - EB) + Mp*Mp - Ml*Ml);
        double denominator = 2*((Mn - EB) - Ep + pz);

        return (numerator / denominator) / 1000.0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, energy_qep, energy_qep);

    /**
     * @brief Variable for total visible energy of interaction, including
     * sub-threshold particles.
     * @details This function calculates the total visible energy of the
     * interaction by summing the energy of all particles that are identified
     * as counting towards the final state of the interaction. Sub-threshold
     * particles are included calorimetrically.
     * @tparam T the type of interaction (true or reco).
     * @param obj interaction to apply the variable on.
     * @return the total visible energy of the interaction.
     */
    template<class T>
    double visible_energy_calosub(const T & obj)
    {
        double energy(0);
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                energy += pvars::energy(p);
                if(pvars::pid(p) == pvars::kProton) energy -= PROTON_MASS - PROTON_BINDING_ENERGY;
            }
            else if(pcuts::is_primary(p))
                energy += p.calo_ke;
        }
        return energy/1000.0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, visible_energy_calosub, visible_energy_calosub);

    /**
     * @brief Variable for the flash time of the interaction.
     * @details The flash time is the time of the flash observed in the PMTs
     * and associated with the charge deposition in the interaction using the
     * OpT0Finder likelihood method.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash time of the interaction.
     */
    template<class T>
    double flash_time(const T & obj)
    {
        if(obj.flash_times.size() > 0)
            return obj.flash_times[0];
        return PLACEHOLDERVALUE;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, flash_time, flash_time);

    /**
     * @brief Variable for the flash score of the interaction.
     * @details The flash score is the likelihood score of the flash observed
     * in the PMTs and associated with the charge deposition in the interaction
     * by the OpT0Finder likelihood method.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash score of the interaction.
     */
    template<class T>
    double flash_score(const T & obj)
    {
        if(obj.flash_scores.size() > 0)
            return obj.flash_scores[0];
        return PLACEHOLDERVALUE;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, flash_score, flash_score);

    /**
     * @brief Variable for the flash total photoelectron count of the
     * interaction.
     * @details The flash total photoelectron count is the total number of
     * photoelectrons observed in the PMTs and associated with the charge
     * deposition in the interaction using the OpT0Finder likelihood method.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash total photoelectron count of the interaction.
     */
    template<class T>
    double flash_total_pe(const T & obj) { return obj.flash_total_pe; }
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, flash_total_pe, flash_total_pe);

    /**
     * @brief Variable for the flash hypothesis total photoelectron count of
     * the interaction.
     * @details The flash hypothesis total photoelectron count is the total
     * number of photoelectrons predicted by OpT0Finder for the interaction in
     * the flash associated with the charge deposition in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash hypothesis total photoelectron count of the interaction.
     */
    template<class T>
    double flash_hypothesis(const T & obj) { return obj.flash_hypo_pe; }
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, flash_hypothesis, flash_hypothesis);

    /**
     * @brief Variable for the x-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the x-coordinate of the interaction vertex.
     */
    template<class T>
    double vertex_x(const T & obj) { return obj.vertex[0]; }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, vertex_x, vertex_x);

    /**
     * @brief Variable for the y-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the y-coordinate of the interaction vertex.
     */
    template<class T>
    double vertex_y(const T & obj) { return obj.vertex[1]; }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, vertex_y, vertex_y);

    /**
     * @brief Variable for the z-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the z-coordinate of the interaction vertex.
     */
    template<class T>
    double vertex_z(const T & obj) { return obj.vertex[2]; }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, vertex_z, vertex_z);

    /**
     * @brief Variable for the transverse momentum of the interaction counting
     * only particles identified as contributing to the final state.
     * @details This function calculates the transverse momentum of the
     * interaction by summing the transverse momentum of all particles that are
     * identified as counting towards the final state of the interaction. The
     * neutrino direction is assumed to either be the BNB axis direction
     * (z-axis) or the unit vector pointing from the NuMI target to the
     * interaction vertex. See @ref utilities::transverse_momentum for details
     * on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the transverse momentum of the primary particles.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double dpT(const T & obj)
    {
        utilities::three_vector pt = {0, 0, 0};
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                // Sum up the transverse momentum of all final state particles
                utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
                pt = utilities::add(pt, this_pt);
            }
        }
        return utilities::magnitude(pt);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dpT, dpT);

    /**
     * @brief Variable for the transverse momentum of the interaction counting
     * only the leading charged lepton and proton.
     * @details This function calculates the transverse momentum of the
     * interaction by summing the transverse momentum of the leading charged
     * lepton and proton. The neutrino direction is assumed to either be the
     * BNB axis direction (z-axis) or the unit vector pointing from the NuMI
     * target to the interaction vertex. See @ref utilities::transverse_momentum
     * for details on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the transverse momentum of the leading charged lepton and proton.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double dpT_lp(const T & obj)
    {

        utilities::three_vector l_pt = {0, 0, 0};
        utilities::three_vector p_pt = {0, 0, 0};
        double l_ke(0), p_ke(0);
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                // Find the leading charged lepton and proton
                if((pvars::pid(p) == pvars::kElectron || pvars::pid(p) == pvars::kMuon) && pvars::ke(p) > l_ke)
                {
                    l_ke = pvars::ke(p);
                    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    l_pt = utilities::transverse_momentum(momentum, vtx);
                }
                else if(pvars::pid(p) == pvars::kProton && pvars::ke(p) > p_ke)
                {
                    p_ke = pvars::ke(p);
                    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    p_pt = utilities::transverse_momentum(momentum, vtx);
                }
            }
        }
        if(l_ke == 0 || p_ke == 0)
            return PLACEHOLDERVALUE;
        else
            return utilities::magnitude(utilities::add(l_pt, p_pt));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dpT_lp, dpT_lp);

    /**
     * @brief Variable for dphi_T of the interaction.
     * @details dphi_T is a transverse kinematic imbalance variable defined
     * using the transverse momentum of the leading muon and the total hadronic
     * system. This variable is sensitive to the presence of F.S.I. The
     * neutrino direction is assumed to either be the BNB axis direction
     * (z-axis) or the unit vector pointing from the NuMI target to the
     * interaction vertex. See @ref utilities::transverse_momentum for details
     * on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the phi_T of the interaction.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double dphiT(const T & obj)
    {
        utilities::three_vector lepton_pt = {0, 0, 0};
        utilities::three_vector hadronic_pt = {0, 0, 0};
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                // There should only be one lepton, so replace the lepton
                // transverse momentum if the particle is a lepton.
                utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
                if(pvars::pid(p) == pvars::kElectron || pvars::pid(p) == pvars::kMuon)
                    lepton_pt = this_pt;
                // The total hadronic system is treated as a single object.
                else if(pvars::pid(p) > 2)
                    hadronic_pt = utilities::add(hadronic_pt, this_pt);
            }
        }
        return std::acos(-1 * utilities::dot_product(lepton_pt, hadronic_pt) / (utilities::magnitude(lepton_pt) * utilities::magnitude(hadronic_pt)));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dphiT, dphiT);

    /**
     * @brief Variable for dphi_T of the interaction in degrees.
     * @details Same as dphiT but converted to degrees for easier interpretation.
     */
    template<class T>
    double dphiT_deg(const T & obj)
    {
        return dphiT(obj) * 180.0 / M_PI;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dphiT_deg, dphiT_deg);

    /**
     * @brief Variable for dalpha_T of the interaction.
     * @details dalpha_T is a transverse kinematic imbalance variable defined
     * using the transverse momentum of the total hadronic system and the
     * outgoing lepton. The neutrino direction is assumed to either be the BNB
     * axis direction (z-axis) or the unit vector pointing from the NuMI target
     * to the interaction vertex. See @ref utilities::transverse_momentum for
     * details on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the alpha_T of the interaction.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double dalphaT(const T & obj)
    {
        utilities::three_vector lepton_pt = {0, 0, 0};
        utilities::three_vector total_pt = {0, 0, 0};
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                // There should only be one lepton, so replace the lepton
                // transverse momentum if the particle is a lepton.
                utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
                if(pvars::pid(p) == pvars::kElectron || pvars::pid(p) == pvars::kMuon)
                    lepton_pt = this_pt;
                total_pt = utilities::add(total_pt, this_pt);
            }
        }
        return std::acos(-1 * utilities::dot_product(total_pt, lepton_pt) / (utilities::magnitude(total_pt) * utilities::magnitude(lepton_pt)));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dalphaT, dalphaT);

    /**
     * @brief Variable for dalpha_T of the interaction in degrees.
     * @details Same as dalphaT but converted to degrees for easier interpretation.
     */
    template<class T>
    double dalphaT_deg(const T & obj)
    {
        return dalphaT(obj) * 180.0 / M_PI;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dalphaT_deg, dalphaT_deg);

    /**
     * @brief Variable for the missing longitudinal momentum of the
     * interaction.
     * @details The missing longitudinal momentum is calculated as the
     * difference between the total longitudinal momentum of the final state
     * particles and the best estimate of the neutrino energy. The neutrino
     * energy is calculated using @ref vars::visible_energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the missing longitudinal momentum of the interaction.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double dpL(const T & obj)
    {
        utilities::three_vector lepton_pl = {0, 0, 0};
        utilities::three_vector hadronic_pl = {0, 0, 0};
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                // There should only be one lepton, so replace the lepton
                // transverse momentum if the particle is a lepton.
                utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                utilities::three_vector this_pl = utilities::longitudinal_momentum(momentum, vtx);
                if(pvars::pid(p) == pvars::kElectron || pvars::pid(p) == pvars::kMuon)
                    lepton_pl = this_pl;
                // The total hadronic system is treated as a single object.
                else if(pvars::pid(p) > 2)
                    hadronic_pl = utilities::add(hadronic_pl, this_pl);
            }
        }
        return utilities::magnitude(utilities::add(hadronic_pl, lepton_pl)) - 1000*vars::visible_energy(obj);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dpL, dpL);

    /**
     * @brief Variable for the missing longitudinal momentum of the interaction
     * counting only the leading charged lepton and proton.
     * @details The missing longitudinal momentum is calculated as the
     * difference between the total longitudinal momentum of the leading charged
     * lepton and proton and the best estimate of the neutrino energy. The
     * neutrino energy is calculated using @ref vars::visible_energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the missing longitudinal momentum of the leading charged lepton
     * and proton.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double dpL_lp(const T & obj)
    {
        utilities::three_vector l_pl = {0, 0, 0};
        utilities::three_vector p_pl = {0, 0, 0};
        double l_ke(0), p_ke(0);
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                // Find the leading charged lepton and proton
                if((pvars::pid(p) == pvars::kElectron || pvars::pid(p) == pvars::kMuon) && pvars::ke(p) > l_ke)
                {
                    l_ke = pvars::ke(p);
                    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    l_pl = utilities::longitudinal_momentum(momentum, vtx);
                }
                else if(pvars::pid(p) == pvars::kProton && pvars::ke(p) > p_ke)
                {
                    p_ke = pvars::ke(p);
                    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    p_pl = utilities::longitudinal_momentum(momentum, vtx);
                }
            }
        }
        if(l_ke == 0 || p_ke == 0)
            return PLACEHOLDERVALUE;
        else
            return utilities::magnitude(utilities::add(l_pl, p_pl)) - 1000*vars::visible_energy(obj);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dpL_lp, dpL_lp);

    /**
     * @brief Variable for the estimate of the momentum of the struck nucleon.
     * @details The estimate of the momentum of the struck nucleon is calculated
     * as the quadrature sum of the transverse momentum (see @ref vars::dpT) and
     * the missing longitudinal momentum (see @ref vars::dpL).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the estimate of the momentum of the struck nucleon.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double pn(const T & obj) { return std::sqrt(std::pow(vars::dpT(obj), 2) + std::pow(vars::dpL(obj), 2)); }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pn, pn);

    /**
     * @brief Variable for the estimate of the momentum of the struck nucleon
     * counting only the leading charged lepton and proton.
     * @details The estimate of the momentum of the struck nucleon is calculated
     * as the quadrature sum of the transverse momentum (see @ref vars::dpT_lp)
     * and the missing longitudinal momentum (see @ref vars::dpL_lp).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the estimate of the momentum of the struck nucleon.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double pn_lp(const T & obj) { return std::sqrt(std::pow(vars::dpT_lp(obj), 2) + std::pow(vars::dpL_lp(obj), 2)); }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pn_lp, pn_lp);

    /**
     * @brief Variable for the opening angle between leading muon and proton.
     * @details The leading muon and proton are defined as the particles with the
     * highest kinetic energy. The opening angle is defined as the arccosine of
     * the dot product of the momentum vectors of the leading muon and proton.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the opening angle between the leading muon and
     * proton.
     */
    template<class T>
    double opening_angle(const T & obj)
    {
        size_t mi = selectors::leading_muon(obj);
        size_t pi = selectors::leading_proton(obj);
        if(mi == kNoMatch || pi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & m(obj.particles[mi]);
            auto & p(obj.particles[pi]);
            return std::acos(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, opening_angle, opening_angle);

    /**
     * @brief Variable for the opening angle between leading muon and proton in degrees.
     * @details Same as opening_angle but converted to degrees for easier interpretation.
     * The leading muon and proton are defined as the particles with the
     * highest kinetic energy. The opening angle is defined as the arccosine of
     * the dot product of the momentum vectors of the leading muon and proton.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the opening angle between the leading muon and proton in degrees.
     */
    template<class T>
    double opening_angle_deg(const T & obj)
    {
        double angle_rad = opening_angle(obj);
        if(std::isnan(angle_rad))
            return kNoMatchValue; // No leading muon or proton found.
        else
            return angle_rad * 180.0 / M_PI;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, opening_angle_deg, opening_angle_deg);

    /**
     * @brief Variable for the cosine of the opening angle between leading muon and proton.
     * @details The leading muon and proton are defined as the particles with the
     * highest kinetic energy. This variable returns the cosine of the opening angle,
     * which is simply the dot product of the direction vectors of the leading muon
     * and proton (assuming normalized direction vectors).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the cosine of the opening angle between the leading muon and proton.
     */
    template<class T>
    double opening_angle_cos(const T & obj)
    {
        size_t mi = selectors::leading_muon(obj);
        size_t pi = selectors::leading_proton(obj);
        if(mi == kNoMatch || pi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & m(obj.particles[mi]);
            auto & p(obj.particles[pi]);
            return m.start_dir[0] * p.start_dir[0] +
                   m.start_dir[1] * p.start_dir[1] +
                   m.start_dir[2] * p.start_dir[2];
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, opening_angle_cos, opening_angle_cos);

    /**
     * @brief Variable for the (primary) photon multiplicity of the
     * interaction.
     * @details This function calculates the multiplicity of primary
     * photons in the interaction by counting the number of primary particles
     * that are identified as photons and have a kinetic energy above a
     * threshold. The threshold is set by the `params` vector, which defaults
     * to 25 MeV. The function returns the number of primary photons in the
     * interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a photon to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return the multiplicity of primary photons in the interaction.
     */
    template<class T>
    double photon_multiplicity(const T & obj, std::vector<double> params={25.0,})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == pvars::kPhoton && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ++count;
        }
        return count;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, photon_multiplicity, photon_multiplicity);

    /**
     * @brief Variable for the (primary) electron multiplicity of the
     * interaction.
     * @details This function calculates the multiplicity of primary electrons
     * in the interaction by counting the number of primary particles that are
     * identified as electrons and have a kinetic energy above a threshold. The
     * threshold is set by the `params` vector, which defaults to 25 MeV. The
     * function returns the number of primary electrons in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for an electron to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return the multiplicity of primary electrons in the interaction.
     */
    template<class T>
    double electron_multiplicity(const T & obj, std::vector<double> params={25.0,})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == pvars::kElectron && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ++count;
        }
        return count;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, electron_multiplicity, electron_multiplicity);

    /**
     * @brief Variable for the (primary) muon multiplicity of the
     * interaction.
     * @details This function calculates the multiplicity of primary muons in
     * the interaction by counting the number of primary particles that are
     * identified as muons and have a kinetic energy above a threshold. The
     * threshold is set by the `params` vector, which defaults to 25 MeV. The
     * function returns the number of primary muons in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a muon to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return the multiplicity of primary muons in the interaction.
     */
    template<class T>
    double muon_multiplicity(const T & obj, std::vector<double> params={25.0,})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == pvars::kMuon && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ++count;
        }
        return count;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, muon_multiplicity, muon_multiplicity);

    /**
     * @brief Variable for the (primary) pion multiplicity of the
     * interaction.
     * @details This function calculates the multiplicity of primary pions in
     * the interaction by counting the number of primary particles that are
     * identified as pions and have a kinetic energy above a threshold. The
     * threshold is set by the `params` vector, which defaults to 25 MeV. The
     * function returns the number of primary pions in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a pion to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return the multiplicity of primary pions in the interaction.
     */
    template<class T>
    double pion_multiplicity(const T & obj, std::vector<double> params={25.0,})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == pvars::kPion && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ++count;
        }
        return count;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pion_multiplicity, pion_multiplicity);

    /**
     * @brief Variable for the (primary) proton multiplicity of the
     * interaction.
     * @details This function calculates the multiplicity of primary protons in
     * the interaction by counting the number of primary particles that are
     * identified as protons and have a kinetic energy above a threshold. The
     * threshold is set by the `params` vector, which defaults to 25 MeV. The
     * function returns the number of primary protons in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a proton to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return the multiplicity of primary protons in the interaction.
     */
    template<class T>
    double proton_multiplicity(const T & obj, std::vector<double> params={25.0,})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == pvars::kProton && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ++count;
        }
        return count;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, proton_multiplicity, proton_multiplicity);

    /**
     * @brief Variable for the distance between the interaction vertex and the
     * leading muon start point.
     * @details This function calculates the distance from the leading muon
     * start point to the interaction vertex. The leading muon is defined as
     * the particle with the highest kinetic energy that is identified as a
     * muon. If no leading muon is found, the function returns the usual 
     * PLACEHOLDERVALUE.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the distance from the leading muon start point to the
     * interaction vertex.
     */
    template<class T>
    double leading_muon_vertex_gap(const T & obj)
    {
        // Find the leading muon in the interaction.
        size_t mi = selectors::leading_muon(obj);
        if(mi == kNoMatch) return PLACEHOLDERVALUE;
        auto & m(obj.particles[mi]);
        
        // Calculate the distance from the leading muon start point to the
        // interaction vertex.
        utilities::three_vector vtx = {obj.vertex[0], obj.vertex[1], obj.vertex[2]};
        utilities::three_vector muon_start = {pvars::start_x(m), pvars::start_y(m), pvars::start_z(m)};
        return utilities::magnitude(utilities::subtract(muon_start, vtx));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, leading_muon_vertex_gap, leading_muon_vertex_gap);

    /**
     * @brief Variable for four-momentum transfer Q**2.
     * @details Variable for four-momentum transfer from incoming
     * neutrino to incident nucleon.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the four-momentum Q**2.
     */
    template<class T>
    double Q2(const T & obj)
    {
      // Find the leading muon in the interaction
      size_t mi = selectors::leading_muon(obj);
      if(mi == kNoMatch) return PLACEHOLDERVALUE;
      auto & m(obj.particles[mi]);
      return 2*visible_energy(obj)*((pvars::energy(m)/1000.0) - (pvars::momentum(m)/1000.0)*pvars::beam_costheta(m)) - std::pow(MUON_MASS/1000.0, 2);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, Q2, Q2);

    /**
     * @brief Variable for hadronic invariant mass.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the hadronic invariant mass W.
     */
    template<class T>
    double W(const T & obj)
    {
        // Find the leading muon in the interaction
        size_t mi = selectors::leading_muon(obj);
	if(mi == kNoMatch) return PLACEHOLDERVALUE;
	auto & m(obj.particles[mi]);
	return std::sqrt( std::pow(NUCLEON_MASS/1000.0, 2) + 2*(NUCLEON_MASS/1000.0)*(visible_energy(obj) - (pvars::energy(m)/1000.0)) - Q2(obj) );
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, W, W);







    // --- Keep this helper: resolve leading-pion children PDGs ---
    template <class T>
    std::vector<int64_t> get_pion_children_pdgs(const T& obj){
        std::vector<int64_t> v;
        const size_t pi = selectors::leading_pion(obj);
        if (pi == kNoMatch) return v;

        const auto& pion = obj.particles[pi];

        std::unordered_map<int64_t,size_t> id2i;
        id2i.reserve(obj.particles.size());
        for (size_t i=0;i<obj.particles.size();++i)
            id2i.emplace((int64_t)obj.particles[i].id, i);

        for (auto cid_raw: pion.children_id){
            auto it = id2i.find((int64_t)cid_raw);
            if (it==id2i.end()) continue;
            v.push_back((int64_t)obj.particles[it->second].pdg_code);
        }
        return v;
    }

    // Build: parent_id -> list of child indices
    template <class T>
    static inline std::unordered_map<int64_t, std::vector<size_t>>
    build_parent_children_map(const T& obj)
    {
        std::unordered_map<int64_t, std::vector<size_t>> pc;
        pc.reserve(obj.particles.size());
        for (size_t i = 0; i < obj.particles.size(); ++i)
        {
            if (obj.particles[i].parent_id == obj.particles[i].id) continue; // <-- skip self-parented entries
            pc[obj.particles[i].parent_id].push_back(i);
        }
        return pc;
    }

    // Descendant scan flags must be templated consistently
    template <class T>
        struct DescendantFlags {
            bool hasMuon        = false;
            bool hasMichelStrict= false; // (shape==2 && |pdg|==11)
            bool hasChargedPi   = false; // |pdg|==211
            bool hasPi0         = false; // 111
            bool hasNuclear     = false; // p/n or PDG >= 1e9
        };

    template <class T>
        static inline DescendantFlags<T>
        scan_direct_children(const T& obj, size_t pion_idx)
        {
            DescendantFlags<T> f;  // same fields: hasMuon, hasMichelStrict, hasChargedPi, hasPi0, hasNuclear
            if (pion_idx == kNoMatch || pion_idx >= obj.particles.size()) return f;

            const auto pc = build_parent_children_map(obj);
            const int64_t root = obj.particles[pion_idx].id;

            auto it = pc.find(root);
            if (it == pc.end()) return f; // no direct children

            for (size_t idx : it->second) {
                const auto& ch = obj.particles[idx];
                if (ch.id == root) continue; // <-- ignore self as a "child"
                const int pdg = ch.pdg_code;
                const int ap  = std::abs(pdg);

                if (ap == 13)                     f.hasMuon = true;
                if (ch.shape == 2 && ap == 11)    f.hasMichelStrict = true; // Michel-like electron as a *direct* daughter only
                if (ap == 211)                    f.hasChargedPi = true;
                if (pdg == 111)                   f.hasPi0 = true;
                if (ap == 2212 || ap == 2112 || ap >= 1000000000) f.hasNuclear = true; // p/n or ion
            }
            return f;
        }

    enum class PionGroupCode:int{ UNKNOWN=0, DECAY=1, CAPTURE=2, INELASTIC=3, ELASTIC=4 };

    template <class T>
        int classify_pion_group_code_int(const T& obj)
        {
            const size_t pi = selectors::leading_pion(obj);
            if (pi == kNoMatch) return (int)PionGroupCode::UNKNOWN;

            // Direct children PDGs (needed for gamma-only and counting nucleons)
            const auto ch = get_pion_children_pdgs(obj);

            auto cntN = [&](){
                int c = 0;
                for (auto x: ch) if (x==2212 || x==2112 || std::abs(x)>=1000000000) ++c;
                return c;
            };

                auto cntGammas = [&](){
                int c = 0;
                for (auto x: ch) if (x == 22) ++c;
                return c;
            };

            auto onlyGammasAndNucleons = [&](){
                if (ch.empty()) return false;
                for (auto x: ch) {
                    if (x != 22 && x != 2212 && x != 2112 && std::abs(x) < 1000000000)
                        return false;
                }
                return true;
            };

            const int  Nnuc_dir = cntN();
            const bool anyN_dir = (Nnuc_dir > 0);
            const int  Ngamma = cntGammas();

            // *** Direct-only scan ***
            const auto df = scan_direct_children(obj, pi);

            // -------- DECAY (direct-only) ----------
            // Require a direct muon to call DECAY.
            // (Optional) Allow a direct Michel-shaped e± ONLY if there are no direct pions/pi0/nuclear products.
            if (df.hasMuon)
                return (int)PionGroupCode::DECAY;

            if (df.hasMichelStrict && !df.hasChargedPi && !df.hasPi0 && !anyN_dir)
                return (int)PionGroupCode::DECAY;

            // -------- INELASTIC (check BEFORE capture) ----------
            // π0 among direct daughters → inelastic
            if (df.hasPi0)
                return (int)PionGroupCode::INELASTIC;

            // Charged π with nucleons → inelastic
            if (df.hasChargedPi && anyN_dir)
                return (int)PionGroupCode::INELASTIC;

            // **NEW**: Exactly 2 gammas (likely π⁰ → γγ from charge-exchange)
            // This catches: π± + nucleus → π⁰ + X → γγ + X
            if (Ngamma == 2 && onlyGammasAndNucleons())
                return (int)PionGroupCode::INELASTIC;

            // -------- CAPTURE (now more restrictive) --------
            // Multiple gammas (NOT exactly 2) suggests nuclear de-excitation
            if (Ngamma > 0 && Ngamma != 2 && onlyGammasAndNucleons() && !df.hasChargedPi)
                return (int)PionGroupCode::CAPTURE;

            // Nucleons with no pions and no gammas → capture
            if (anyN_dir && !df.hasChargedPi && !df.hasPi0 && Ngamma == 0)
                return (int)PionGroupCode::CAPTURE;

            // -------- ELASTIC ----------
            // Exactly one direct charged π and nothing else → elastic
            {
                int nCpi = 0, nOther = 0;
                for (auto x: ch) {
                    if (std::abs(x) == 211) ++nCpi;
                    else if (x != 0) ++nOther;
                }
                if (nCpi == 1 && nOther == 0)
                    return (int)PionGroupCode::ELASTIC;
            }

            // Any charged π at all (and no stronger signature above) → inelastic bucket
            if (df.hasChargedPi)
                return (int)PionGroupCode::INELASTIC;

            // -------- Fallback ----------
            return (int)PionGroupCode::UNKNOWN;
        }

    // --- Scalar wrappers (double) for your binder ---
    template <class T> double pion_group_code(const T& obj){ return (double)classify_pion_group_code_int(obj); }
    template <class T> double pion_is_decay(const T& obj){ return pion_group_code(obj)==(double)PionGroupCode::DECAY ? 1.0:0.0; }
    template <class T> double pion_is_capture(const T& obj){ return pion_group_code(obj)==(double)PionGroupCode::CAPTURE ? 1.0:0.0; }
    template <class T> double pion_is_inelastic(const T& obj){ return pion_group_code(obj)==(double)PionGroupCode::INELASTIC ? 1.0:0.0; }

    // (Optional) tiny debug helpers
    template <class T> double pion_children_count(const T& obj){ return (double)get_pion_children_pdgs(obj).size(); }
    template <class T> double pion_child0_pdg(const T& obj){

        auto v=get_pion_children_pdgs(obj);

        return v.empty()? -9999.0 : (double)v[0]; }

    // --- Register once (in a single TU) ---
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_group_code,     vars::pion_group_code);
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_is_decay,       vars::pion_is_decay);
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_is_capture,     vars::pion_is_capture);
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_is_inelastic,   vars::pion_is_inelastic);
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_children_count, vars::pion_children_count);
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_child0_pdg,     vars::pion_child0_pdg);

    template <class T>
        std::vector<int64_t> get_muon_children_pdgs(const T& obj){
            std::vector<int64_t> v;
            const size_t mu = selectors::leading_muon(obj);
            if (mu == kNoMatch) return v;

            const auto& muon = obj.particles[mu];

            std::unordered_map<int64_t,size_t> id2i;
            id2i.reserve(obj.particles.size());
            for (size_t i=0;i<obj.particles.size();++i)
                id2i.emplace((int64_t)obj.particles[i].id, i);

            for (auto cid_raw: muon.children_id){
                auto it = id2i.find((int64_t)cid_raw);
                if (it==id2i.end()) continue;
                v.push_back((int64_t)obj.particles[it->second].pdg_code);
            }
            return v;
        }

    // --- Muon group codes ---
    enum class MuonGroupCode:int{ UNKNOWN=0, DECAY=1, CAPTURE=2, INELASTIC=3 };

    template <class T>
        int classify_muon_group_code_int(const T& obj){
            auto ch = get_muon_children_pdgs(obj);
            if (ch.empty()) return (int)MuonGroupCode::UNKNOWN;

            auto has = [&](std::initializer_list<int64_t> s){
                for (auto x: ch) for (auto y: s) if (x==y) return true; return false;
            };
            auto cntN = [&](){
                int c=0; for (auto x: ch) if (x==2212||x==2112) ++c; return c;
            };
            auto onlyGammas = [&](){
                if (ch.empty()) return false;
                for (auto x: ch) if (x!=22) return false; return true;
            };

            const bool hasE  = has({11,-11});
            const bool anyN  = (cntN()>0) || has({1000010020,1000010030,1000020040});

            // DECAY (μ± → e± ν ν)
            if (hasE) return (int)MuonGroupCode::DECAY;

            // CAPTURE (μ− absorption in nucleus)
            if (onlyGammas())       return (int)MuonGroupCode::CAPTURE;
            if (anyN && !hasE)      return (int)MuonGroupCode::CAPTURE;

            // INELASTIC (μ scatters and survives with extra hadrons)
            if (has({13,-13}))      return (int)MuonGroupCode::INELASTIC;

            return (int)MuonGroupCode::UNKNOWN;
        }

    // --- Scalar wrappers (double) ---
    template <class T> double muon_group_code(const T& obj){ return (double)classify_muon_group_code_int(obj); }
    template <class T> double muon_is_decay(const T& obj){ return muon_group_code(obj)==(double)MuonGroupCode::DECAY ? 1.0:0.0; }
    template <class T> double muon_is_capture(const T& obj){ return muon_group_code(obj)==(double)MuonGroupCode::CAPTURE ? 1.0:0.0; }
    template <class T> double muon_is_inelastic(const T& obj){ return muon_group_code(obj)==(double)MuonGroupCode::INELASTIC ? 1.0:0.0; }

    // (Optional debug helpers)
    template <class T> double muon_children_count(const T& obj){ return (double)get_muon_children_pdgs(obj).size(); }
    template <class T> double muon_child0_pdg(const T& obj){ auto v=get_muon_children_pdgs(obj); return v.empty()? -9999.0 : (double)v[0]; }

    // --- Register once (in a single TU) ---
    REGISTER_VAR_SCOPE(RegistrationScope::True, muon_group_code,     vars::muon_group_code);
    REGISTER_VAR_SCOPE(RegistrationScope::True, muon_is_decay,       vars::muon_is_decay);
    REGISTER_VAR_SCOPE(RegistrationScope::True, muon_is_capture,     vars::muon_is_capture);
    REGISTER_VAR_SCOPE(RegistrationScope::True, muon_is_inelastic,   vars::muon_is_inelastic);


}

namespace cuts {
    template<class T>
    bool is_pion_group_code(const T & obj, std::vector<double> params)
    {
        if(params.empty())
            return true;

        double code = vars::pion_group_code(obj);
        return std::find(params.begin(), params.end(), code) != params.end();
    }
}


#endif // VARIABLES_H
