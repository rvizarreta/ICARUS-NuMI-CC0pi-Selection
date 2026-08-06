/**
 * @file mctruth.h
 * @brief Definitions of analysis variables which can extract information from
 * the SRTrueInteraction object.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from the SRTrueInteraction object. Each variable
 * is implemented as a function which takes an SRTrueInteraction object as an
 * argument and returns a double. The association of an SRInteractionTruthDLP
 * object to an SRTrueInteraction object is handled upstream in the SpineVar
 * functions.
 * @author mueller@fnal.gov
 */
#ifndef MCTRUTH_H
#define MCTRUTH_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

#include "framework.h"

/**
 * @namespace mctruth
 * @brief Namespace for organizing variables which act on true interactions.
 * @details This namespace is intended to be used for organizing variables
 * which act on true interactions. Each variable is implemented as a function
 * which takes an SRTrueInteraction object as an argument and returns a double.
 */
namespace mctruth
{
    /**
     * @brief Variable for the true neutrino energy.
     * @details This variable is intended to provide the true energy of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino energy.
     */
    template<typename T>
        double neutrino_energy(const T & obj) { return obj.E; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, neutrino_energy, neutrino_energy);

    /**
     * @brief Variable for the true neutrino baseline.
     * @details This variable is intended to provide the true baseline of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino baseline.
     */
    template<typename T>
        double baseline(const T & obj) { return obj.baseline; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, baseline, baseline);

    /**
     * @brief Variable for the true neutrino PDG code.
     * @details This variable is intended to provide the true PDG code of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino PDG code.
     */
    template<typename T>
        double pdg(const T & obj) { return obj.pdg; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pdg, pdg);

    /**
     * @brief Variable for the PDG code of the parent of the neutrino.
     * @details This variable is intended to provide the PDG code of the
     * parent of the neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the PDG code of the parent of the neutrino.
     */
    template<typename T>
        double parent_pdg(const T & obj) { return obj.parent_pdg; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, parent_pdg, parent_pdg);

    /**
     * @brief Variable for the true neutrino current value.
     * @details This variable is intended to provide the true current value of
     * the parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino current value.
     */
    template<typename T>
        double cc(const T & obj) { return obj.iscc; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, cc, cc);

    /**
     * @brief Variable for the interaction mode of the interaction.
     * @details This variable is intended to provide the interaction mode of the
     * interaction. This is based on the GENIE interaction mode enumeration 
     * defined in the LArSoft MCNeutrino class.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the interaction mode.
     * Here for numeration: https://code-doc.larsoft.org/docs/latest/html/MCNeutrino_8h_source.html
     */
    template<typename T>
        double interaction_mode(const T & obj) { return obj.genie_mode; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, interaction_mode, interaction_mode);

    /**
     * @brief Variable for the interaction type of the interaction.
     * @details This variable is intended to provide the interaction type of the
     * interaction. This is based on the GENIE interaction type enumeration 
     * defined in the LArSoft MCNeutrino class.
     * @param T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the interaction type.
     */
    template<typename T>
        double interaction_type(const T & obj) { return obj.genie_inttype; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, interaction_type, interaction_type);

    /**
     * @brief Cut for exactly one true final state muon above threshold.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 143.425 MeV.
     * @return true if exactly one muon above threshold.
     */
    template<typename T>
    bool single_muon_srtruth(const T & obj, std::vector<double> params={143.425,})
    {
        int num_muons(0);
        for(const auto & p : obj.prim)
        {
            if(abs(p.pdg) == 13)
            {
                double ke = 1000. * (p.genE - (MUON_MASS/1000.));
                if(ke >= params[0])
                    num_muons++;
            }
        }
        return num_muons == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, single_muon_srtruth, single_muon_srtruth);

    /**
     * @brief Cut for at least one true final state proton above threshold.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 50 MeV.
     * @return true if at least one proton above threshold.
     */
    template<typename T>
    bool has_proton_srtruth(const T & obj, std::vector<double> params={50.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 2212)
            {
                double ke = 1000. * (p.genE - (PROTON_MASS/1000.));
                if(ke >= params[0])
                    return true;
            }
        }
        return false;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, has_proton_srtruth, has_proton_srtruth);

    /**
     * @brief Cut for exactly one true final state proton above threshold.
     */
    template<typename T>
    bool single_proton_srtruth(const T & obj, std::vector<double> params={50.0,})
    {
        int count = 0;
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 2212)
            {
                double ke = 1000. * (p.genE - (PROTON_MASS/1000.));
                if(ke >= params[0])
                    ++count;
            }
        }
        return count == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, single_proton_srtruth, single_proton_srtruth);

    /**
     * @brief Cut for no true final state protons above threshold.
     */
    template<typename T>
    bool no_protons_srtruth(const T & obj, std::vector<double> params={50.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 2212)
            {
                double ke = 1000. * (p.genE - (PROTON_MASS/1000.));
                if(ke >= params[0])
                    return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_protons_srtruth, no_protons_srtruth);

    /**
     * @brief Cut for zero true final state charged pions above threshold.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 25 MeV.
     * @return true if no charged pions above threshold.
     */
    template<typename T>
    bool no_charged_pions_srtruth(const T & obj, std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(abs(p.pdg) == 211)
            {
                double ke = 1000. * (p.genE - (PION_MASS/1000.));
                if(ke >= params[0])
                    return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_charged_pions_srtruth, no_charged_pions_srtruth);

    /**
     * @brief Cut for zero true final state neutral pions.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 0 MeV.
     * @return true if no neutral pions above threshold.
     */
    template<typename T>
    bool no_neutral_pions(const T & obj)
    {
        for(const auto & p : obj.prim)
        {
            if(p.start_process != 0) continue;
            if(p.pdg == 111)
            {
                return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_neutral_pions, no_neutral_pions);

    /**
     * @brief Returns 1 if no extra mesons are present in the event.
     * @details Checks for the absence of:
     *   - Kaons (charged and neutral)
     *   - Eta mesons
     *   - K* resonances
     * Uses GENIE truth variables (obj.prim) to match the NUISANCE FlatTree
     * signal definition for ICARUS_1muNp0pi.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return 1 if no extra mesons are found, 0 otherwise.
     */
    template<typename T>
    double no_extra_mesons(const T & obj)
    {
        unsigned int nMesons(0);

        for(const auto & p : obj.prim)
        {
            if(p.start_process != 0) continue;

            int pdg = p.pdg;

            // ── Mesons (kaons, eta, neutral kaons, K*) ────────────────────────
            if(std::abs(pdg) == 321 || std::abs(pdg) == 323 ||
               pdg == 130            || pdg == 310            ||
               pdg == 311            || pdg == 313            ||
               std::abs(pdg) == 221  || std::abs(pdg) == 331)
                nMesons++;
        }

        if(nMesons > 0) return false;

        return true;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, no_extra_mesons, no_extra_mesons);

    /**
     * @brief Returns 1 if no heavy baryons are present in the event.
     * @details Checks for the absence of:
     *   - Strange baryons (Sigma+/-, Lambda, Sigma0)
     *   - Charmed baryons and D mesons
     * Uses GENIE truth variables (obj.prim) to match the NUISANCE FlatTree
     * signal definition for ICARUS_1muNp0pi.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return 1 if no heavy baryons are found, 0 otherwise.
     */
    template<typename T>
    double no_extra_baryons(const T & obj)
    {
        unsigned int nBaryons(0);

        for(const auto & p : obj.prim)
        {
            if(p.start_process != 0) continue;

            int pdg = p.pdg;

            // ── Heavy baryons (strange, charmed) ──────────────────────────────
            if(pdg == 3112 || pdg == 3122 || pdg == 3212 || pdg == 3222 ||
               pdg == 4112 || pdg == 4122 || pdg == 4212 || pdg == 4222 ||
               pdg == 411  || pdg == 421)
                nBaryons++;
        }

        if(nBaryons > 0) return false;

        return true;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, no_extra_baryons, no_extra_baryons);

    /**
     * @brief Returns 1 if no photons, extra mesons, or heavy baryons/pi0 are
     *        present in the event, replicating the MINERvA-style signal definition.
     * @details Checks (exclusively, in order) for the absence of:
     *   - Photons with E > 10 MeV (PDG 22)
     *   - Mesons: charged pions, kaons (charged/neutral), eta, pi0, K* (PDG 211,
     *             321, 323, 111, 130, 310, 311, 313, 221, 331)
     *   - Heavy baryons + pi0: strange/charmed baryons, D mesons, pi0 (PDG 3112,
     *             3122, 3212, 3222, 4112, 4122, 4212, 4222, 411, 421, 111)
     * Each particle increments at most one counter (else-if chain), matching the
     * NUISANCE FlatTree signal definition for ICARUS_1muNp0pi.
     * Inspired by MINERvA: https://github.com/NUISANCEMC/nuisance/blob/main/src/MINERvA/MINERvA_SignalDef.cxx#L469-L481
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return 1 if none of the above are found, 0 otherwise.
     */
    template<typename T>
    double no_extra_particles_minerva(const T & obj)
    {
        unsigned int nPhotons(0);
        unsigned int nMesons(0);
        unsigned int nBaryonsAndPi0(0);

        for(const auto & p : obj.prim)
        {
            if(p.start_process != 0) continue;

            int pdg = p.pdg;

            // ── Photons with E > 10 MeV ───────────────────────────────────────
            if(std::abs(pdg) == 22 && p.startE > 0.01)
                nPhotons++;

            // ── Mesons (charged pions, kaons, eta, pi0, K*) ───────────────────
            else if(std::abs(pdg) == 211 || std::abs(pdg) == 321 || std::abs(pdg) == 323 ||
                    pdg == 111           || pdg == 130            || pdg == 310            ||
                    pdg == 311           || pdg == 313            ||
                    std::abs(pdg) == 221 || std::abs(pdg) == 331)
                nMesons++;

            // ── Heavy baryons + pi0 (strange, charmed, D mesons) ─────────────
            else if(pdg == 3112 || pdg == 3122 || pdg == 3212 || pdg == 3222 ||
                    pdg == 4112 || pdg == 4122 || pdg == 4212 || pdg == 4222 ||
                    pdg == 411  || pdg == 421  || pdg == 111)
                nBaryonsAndPi0++;
        }

        if(nPhotons       > 0) return false;
        if(nMesons        > 0) return false;
        if(nBaryonsAndPi0 > 0) return false;

        return true;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, no_extra_particles_minerva, no_extra_particles_minerva);

    /**
     * @brief Cut for zero true final state photons above threshold.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params energy threshold in MeV, defaults to 25 MeV.
     * @return true if no photons above threshold.
     */
    template<typename T>
    bool no_photons_srtruth(const T & obj, std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 22)
            {
                double energy = 1000. * p.genE;
                if(energy >= params[0])
                    return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_photons_srtruth, no_photons_srtruth);

    /**
     * @brief Cut for zero true final state electrons above threshold.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 25 MeV.
     * @return true if no electrons above threshold.
     */
    template<typename T>
    bool no_electrons_srtruth(const T & obj, std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(abs(p.pdg) == 11)
            {
                double ke = 1000. * (p.genE - (ELECTRON_MASS/1000.));
                if(ke >= params[0])
                    return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_electrons_srtruth, no_electrons_srtruth);

    /**
     * @brief Variable for the generator-level four-momentum transfer Q**2.
     * @details Returns the four-momentum transfer squared as computed by the
     * event generator (GENIE) and stored directly on the truth interaction
     * object, rather than reconstructing it from final-state kinematics as
     * Q2() does. This is the quantity the Stowell et al. (arXiv:1903.01558)
     * low-Q2 suppression parameterization was derived against, and should be
     * used in place of Q2() whenever the suppression weight is computed at
     * truth level. Returns a placeholder if the stored value is non-positive
     * or non-finite.
     * @tparam T the type of interaction (true only).
     * @param obj the interaction to apply the variable on.
     * @return the generator-level four-momentum transfer Q**2.
     */
    template<class T>
    double generator_q2(const T & obj)
    {
        if(obj.Q2 <= 0.0 || !std::isfinite(obj.Q2))
            return PLACEHOLDERVALUE;
        return obj.Q2;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, generator_q2, generator_q2);

    /**
     * @brief Replicates NUISANCE ICARUS_1muNp0pi_IsSignal definition exactly.
     * @details This variable is implemented using GENIE truth variables (obj.prim)
     * to match the NUISANCE FlatTree signal definition for ICARUS_1muNp0pi.
     * The selection requires:
     *   - Exactly one muon with momentum > 226 MeV/c
     *   - At least one proton with leading momentum > 310 MeV/c
     *   - No charged pions with momentum > 87 MeV/c
     *   - No neutral pions (no threshold)
     *   - No photons with energy > 10 MeV
     *   - No electrons with energy > 10 MeV
     *   - No mesons (kaons, eta, neutral kaons, K*)
     *   - No heavy baryons (strange, charmed)
     * This is intended to be saved as a branch for direct comparison with
     * NUISANCE cross-section predictions.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return 1 if the event passes the NUISANCE signal definition, 0 otherwise.
     */
    template<typename T>
    double nuisance_is_signal(const T & obj)
    {
        unsigned int nMu(0), nP(0), nPi(0);
        unsigned int nPhoton(0), nElectron(0), nMesons(0), nBaryons(0);
        double maxMomentumP = -999.;
        bool passProtonPCut = false;

        for(const auto & p : obj.prim)
        {
            if(p.start_process != 0) continue;

            // Compute momentum magnitude in GeV/c from GENIE genp
            double px = p.genp.x, py = p.genp.y, pz = p.genp.z;
            double momentum = std::sqrt(px*px + py*py + pz*pz);  // GeV/c

            int pdg = p.pdg;

            // ── Muon ──────────────────────────────────────────────────────────
            if(std::abs(pdg) == 13 && momentum > 0.226)
                nMu++;

            // ── Proton ────────────────────────────────────────────────────────
            if(std::abs(pdg) == 2212)
            {
                nP++;
                if(momentum > maxMomentumP)
                {
                    maxMomentumP   = momentum;
                    passProtonPCut = (momentum > 0.310);
                }
            }

            // ── Charged pion (with threshold) ─────────────────────────────────
            if(std::abs(pdg) == 211 && momentum > 0.087)
                nPi++;

            // ── Neutral pion ─────────────────────────────────
            if(std::abs(pdg) == 111)
                nPi++;

            // Photon — 25 MeV threshold
            if(pdg == 22 && p.genE * 1000. > 25.0)
                nPhoton++;

            // Electron — 25.5 MeV/c momentum threshold
            else if(std::abs(pdg) == 11 && momentum > 0.0255)
                nElectron++;

            // Mesons — only if NOT photon/electron (else if chain)
            else if(std::abs(pdg) == 211 || std::abs(pdg) == 321 || std::abs(pdg) == 323 ||
                    pdg == 111            || pdg == 130            || pdg == 310            ||
                    pdg == 311            || pdg == 313            ||
                    std::abs(pdg) == 221  || std::abs(pdg) == 331)
                nMesons++;

            else if(pdg == 3112 || pdg == 3122 || pdg == 3212 || pdg == 3222 ||
                    pdg == 4112 || pdg == 4122 || pdg == 4212 || pdg == 4222 ||
                    pdg == 411  || pdg == 421  || pdg == 111)
                nBaryons++;
        }

        // ── Signal condition (exact NUISANCE logic) ───────────────────────────
        if(nMu != 1)         return 0.;
        if(nP  == 0)         return 0.;
        if(!passProtonPCut)  return 0.;
        if(nPi > 0)          return 0.;
        if(nPhoton > 0)      return 0.;
        if(nElectron > 0)    return 0.;
        if(nMesons > 0)      return 0.;
        if(nBaryons > 0)     return 0.;

        return 1.;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nuisance_is_signal, nuisance_is_signal);

    // ── Diagnostic counters for nuisance_is_signal debugging ─────────────────

    template<typename T>
    double nuisance_nMu(const T & obj) {
        unsigned int n = 0;
        for(const auto & p : obj.prim) {
            double px = p.genp.x, py = p.genp.y, pz = p.genp.z;
            double momentum = std::sqrt(px*px + py*py + pz*pz);
            if(std::abs(p.pdg) == 13 && momentum > 0.226) n++;
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nuisance_nMu, nuisance_nMu);

    template<typename T>
    double nuisance_nP(const T & obj) {
        unsigned int n = 0;
        bool pass = false;
        double maxP = -999.;
        for(const auto & p : obj.prim) {
            double px = p.genp.x, py = p.genp.y, pz = p.genp.z;
            double momentum = std::sqrt(px*px + py*py + pz*pz);
            if(std::abs(p.pdg) == 2212) {
                n++;
                if(momentum > maxP) { maxP = momentum; pass = (momentum > 0.310); }
            }
        }
        return pass ? n : -n;  // negative if leading proton fails threshold
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nuisance_nP, nuisance_nP);

    template<typename T>
    double nuisance_nPi(const T & obj) {
        unsigned int n = 0;
        for(const auto & p : obj.prim) {
            double px = p.genp.x, py = p.genp.y, pz = p.genp.z;
            double momentum = std::sqrt(px*px + py*py + pz*pz);
            if(std::abs(p.pdg) == 211 && momentum > 0.087) n++;
            if(p.pdg == 111) n++;
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nuisance_nPi, nuisance_nPi);

    template<typename T>
    double nuisance_nPiZero(const T & obj) {
        unsigned int n = 0;
        for(const auto & p : obj.prim) {
            double px = p.genp.x, py = p.genp.y, pz = p.genp.z;
            double momentum = std::sqrt(px*px + py*py + pz*pz);
            if(p.pdg == 111) n++;
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nuisance_nPiZero, nuisance_nPiZero);

    template<typename T>
    double nuisance_nPhoton(const T & obj) {
        unsigned int n = 0;
        for(const auto & p : obj.prim)
            if(p.pdg == 22 && p.genE * 1000. > 25.0) n++;
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nuisance_nPhoton, nuisance_nPhoton);

    template<typename T>
    double nuisance_nElectron(const T & obj) {
        unsigned int n = 0;
        for(const auto & p : obj.prim) {
            double px = p.genp.x, py = p.genp.y, pz = p.genp.z;
            double momentum = std::sqrt(px*px + py*py + pz*pz);
            if(std::abs(p.pdg) == 11 && momentum > 0.0255) n++;
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nuisance_nElectron, nuisance_nElectron);

    template<typename T>
    double nuisance_nMesons(const T & obj) {
        unsigned int n = 0;
        for(const auto & p : obj.prim) {
            double px = p.genp.x, py = p.genp.y, pz = p.genp.z;
            double momentum = std::sqrt(px*px + py*py + pz*pz);
            if(p.pdg == 22 && p.genE * 1000. > 25.0) continue;
            else if(std::abs(p.pdg) == 11 && momentum > 0.0255) continue;
            else if(std::abs(p.pdg) == 211 || std::abs(p.pdg) == 321 || std::abs(p.pdg) == 323 ||
                    p.pdg == 111            || p.pdg == 130            || p.pdg == 310            ||
                    p.pdg == 311            || p.pdg == 313            ||
                    std::abs(p.pdg) == 221  || std::abs(p.pdg) == 331) n++;
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nuisance_nMesons, nuisance_nMesons);

    template<typename T>
    double nuisance_nBaryons(const T & obj) {
        unsigned int n = 0;
        for(const auto & p : obj.prim) {
            double px = p.genp.x, py = p.genp.y, pz = p.genp.z;
            double momentum = std::sqrt(px*px + py*py + pz*pz);
            if(p.pdg == 22 && p.genE * 1000. > 25.0) continue;
            else if(std::abs(p.pdg) == 11 && momentum > 0.0255) continue;
            else if(p.pdg == 3112 || p.pdg == 3122 || p.pdg == 3212 || p.pdg == 3222 ||
                    p.pdg == 4112 || p.pdg == 4122 || p.pdg == 4212 || p.pdg == 4222 ||
                    p.pdg == 411  || p.pdg == 421  || p.pdg == 111) n++;
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nuisance_nBaryons, nuisance_nBaryons);

    template<class T>
    double dpT_lp_genie(const T & obj)
    {
        double mu_px = 0, mu_py = 0, mu_pz = 0;
        double lp_px = 0, lp_py = 0, lp_pz = 0;
        double max_mu_p = -1, max_p_p = -1;

        for(const auto & p : obj.prim)
        {
            double px = p.genp.x, py = p.genp.y, pz = p.genp.z;
            double mom = std::sqrt(px*px + py*py + pz*pz);

            if(std::abs(p.pdg) == 13 && mom > max_mu_p)
            {
                max_mu_p = mom;
                mu_px = px; mu_py = py; mu_pz = pz;
            }
            else if(p.pdg == 2212 && mom > max_p_p)
            {
                max_p_p = mom;
                lp_px = px; lp_py = py; lp_pz = pz;
            }
        }

        if(max_mu_p < 0 || max_p_p < 0) return PLACEHOLDERVALUE;

        // Use true neutrino momentum direction (consistent with NUISANCE)
        double nu_mag = std::sqrt(obj.momentum.x*obj.momentum.x +
                                  obj.momentum.y*obj.momentum.y +
                                  obj.momentum.z*obj.momentum.z);
        if(nu_mag == 0) return PLACEHOLDERVALUE;

        utilities::three_vector nu_dir = std::make_tuple(obj.momentum.x/nu_mag,
                                                          obj.momentum.y/nu_mag,
                                                          obj.momentum.z/nu_mag);

        auto transverse = [&](utilities::three_vector v) {
            double dot = std::get<0>(v)*std::get<0>(nu_dir) +
                         std::get<1>(v)*std::get<1>(nu_dir) +
                         std::get<2>(v)*std::get<2>(nu_dir);
            return std::make_tuple(std::get<0>(v) - dot*std::get<0>(nu_dir),
                                   std::get<1>(v) - dot*std::get<1>(nu_dir),
                                   std::get<2>(v) - dot*std::get<2>(nu_dir));
        };

        utilities::three_vector mu_pt = transverse(std::make_tuple(mu_px, mu_py, mu_pz));
        utilities::three_vector lp_pt = transverse(std::make_tuple(lp_px, lp_py, lp_pz));

        return 1000.0 * utilities::magnitude(utilities::add(mu_pt, lp_pt));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, dpT_lp_genie, dpT_lp_genie);

} // namespace mctruth
#endif