// physics_constants.h
#ifndef PHYSICS_CONSTANTS_H
#define PHYSICS_CONSTANTS_H

namespace PhysicsConstants
{
    // ------------------------------
    // Sensitivity Values (Previous Limits)
    // ------------------------------
    constexpr double sigma_prev_germanium = 2.2e-16; // e·cm
    constexpr double sigma_prev_silicium = 5.6e-16; // e·cm

    // ------------------------------
    // Cross Sections (13 TeV, in barn)
    // ------------------------------
    constexpr double sigma_bb = 495e-6;   // b b̄
    constexpr double sigma_cc = 2369e-6;  // c c̄

    // ------------------------------
    // Fragmentation Fractions
    // ------------------------------
    constexpr double f_c_to_Lambda_c = 0.18;
    constexpr double f_q_to_Xib = 0.054;

    // ------------------------------
    // Branching Ratios: Lambda_c channels
    // ------------------------------
    constexpr double BR_Lambda_1 = 0.3067;
    constexpr double BR_Lambda_2 = 0.1772;
    constexpr double BR_Lambda_3 = 0.1772;
    constexpr double BR_Lambda_4 = 0.0632;
    constexpr double BR_Lambda_5 = 0.2182;
    constexpr double BR_Lambda_6 = 0.1772;

    // ------------------------------
    // Branching Ratios: Xi_b channels
    // ------------------------------
    constexpr double BR_Xib_Jpsi_Xi = 1.02e-5;
    constexpr double BR_Xib_Jpsi_Lambda_K = 2.5e-6;
    constexpr double BR_Xib_pbar_K_K = 2.3e-6;
    constexpr double BR_Xib_pbar_pi_pi = 1.3e-6;
    constexpr double BR_Xib_pbar_K_pi = 2.3e-6;
    constexpr double BR_Xib_Lambdab_pi = 7.0e-4;
    constexpr double BR_Xib_Xic_Ds = 1.9e-3;
    constexpr double BR_Xib_Xi_gamma = 1.3e-4;
}

#endif // PHYSICS_CONSTANTS_H
