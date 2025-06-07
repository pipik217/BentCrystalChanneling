#include "properties.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <cmath>
#include <iostream>

// Global theta_c variable (placeholder)
double g_theta_c = 0.01; // Crystal bending angle (radians)

// Silicium sensitivity value σδ (×10⁻¹⁶ e·cm)
const double sigma_prev_germanium = 2.2e-16;
const double sigma_prev_silicium = 5.6e-16;

namespace eric
{
    double getPV(const TLorentzVector& p)
    {
        return p.P() * p.Beta();
    }

    double getThetaLindhard(const TLorentzVector& p)
    {
        return std::sqrt(2 * Config::Uxc / getPV(p));
    }

    double getThetaY(const TLorentzVector& p)
    {
        return p.Py() / p.Pz();
    }

    double getAStraight(double Xc, double Dp, double thetaY, double thetaLindhard)
    {
        return ((2 * Xc / Dp) * std::sqrt(1 - pow(thetaY / thetaLindhard, 2)));
    }
}

namespace original {
    Bool_t inLindhard(const TLorentzVector& p)
    {
        Double_t thetaY = p.Py() / p.Pz();
        return (TMath::Abs(thetaY) < eric::getThetaLindhard(p));
    }

    Bool_t decayAfterCrystal(Double_t l)
    {
        return (l > Config::lengthCrystal);
    }
}

Bool_t inCriticalRadius(const TLorentzVector& p)
{
    Double_t pv = eric::getPV(p);
    Double_t Rc = pv / Config::Uprime;
    Double_t R = (Config::lengthCrystal / Config::thetaCrystal) * 100.0; // cm
    return Rc <= R;
}

Bool_t passesTrappingProbability(const TLorentzVector& p)
{
    Double_t pv = eric::getPV(p);
    Double_t Rc = pv / Config::Uprime;
    Double_t R = (Config::lengthCrystal / Config::thetaCrystal) * 100.0; // cm
    Double_t A_bent = eric::getAStraight(Config::xc, Config::d, eric::getThetaY(p), eric::getThetaLindhard(p)) * (1.0 - Rc / R);
    if (A_bent < 0.0 || A_bent > 1.0) return kFALSE;
    return (gRandom->Uniform() <= A_bent);
}

Bool_t passesDechanneling(const TLorentzVector& p, Double_t pv)
{
    Double_t aTF_cm = Config::aTF * Constants::angstromcmRatio;
    Double_t d_cm = Config::d * Constants::angstromcmRatio;

    Double_t gamma = p.Gamma();
    Double_t N = Constants::NA * Config::rho / Config::A;

    Double_t logArg = (2.0 * Constants::me_c2 * gamma) / Config::I;
    if (logArg <= 1.0) {
        printf("[Dech] logTerm invalid: logArg=%.3e\n", logArg);
        return kFALSE;
    }

    Double_t logTerm = log(logArg);
    Double_t thetaD = Constants::prefactor * N * Config::Z * aTF_cm * pow(d_cm, 2) / (logTerm - 1);

    Double_t Rc = pv / Config::Uprime;
    Double_t R = (Config::lengthCrystal / Config::thetaCrystal) * 100.0; // cm

    if (Rc >= R) return kFALSE;

    Double_t frac = Rc / R;
    Double_t exponent = Config::thetaCrystal / (thetaD * frac * pow(1.0 - frac, 2));
    Double_t w = exp(-exponent);

    return (gRandom->Uniform() <= w);
}

// Compute expected number of detected events, including channeling efficiency
double computeExpectedEvents(
    double sigma_qq,
    double f_q_to_H,
    double B_decay,
    double channeling_efficiency
) {
    double N_ideal = sigma_qq
        * f_q_to_H
        * B_decay;

    double N_detected = N_ideal * channeling_efficiency;
    return N_detected;
}

double computeSensitivity(
    double gamma_prev,
    double theta_c_prev,
    double N_prev,
    double gamma_target,
    double theta_c_target,
    double N_target,
    double sigma_prev  // old sensitivity value (from paper)
) {
    if (N_prev <= 0 || N_target <= 0) {
        return -1.0; // invalid input
    }

    double param_ratio = (gamma_prev * theta_c_prev) / (gamma_target * theta_c_target);
    double stat_ratio = std::sqrt(N_target / N_prev);

    double sig_target = param_ratio * stat_ratio * sigma_prev;
    return sig_target;
}

// Compute channeling efficiency and accumulate gamma for a given file
struct Result {
    double channelingEfficiency;
    double gammaSum;
    Long64_t entries;
};

Result processFile(const char* filename)
{
    TFile f(filename, "READ");
    if (f.IsZombie()) {
        std::cerr << "Error opening file " << filename << std::endl;
        return { 0.0, 0.0, 0 };
    }

    TTree* t = (TTree*)f.Get("tree");
    if (!t) {
        std::cerr << "Tree not found in file " << filename << std::endl;
        return { 0.0, 0.0, 0 };
    }

    Double_t px, py, pz, E, l;
    t->SetBranchAddress("px", &px);
    t->SetBranchAddress("py", &py);
    t->SetBranchAddress("pz", &pz);
    t->SetBranchAddress("E", &E);
    t->SetBranchAddress("l", &l);

    Int_t count = 0;
    Double_t gammaSum = 0.0;
    Long64_t n = t->GetEntries();

    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        TLorentzVector p(px, py, pz, E);
        Double_t pv = eric::getPV(p);

        bool failed = false;

        if (!original::inLindhard(p)) failed = true;
        if (!original::decayAfterCrystal(l)) failed = true;
        if (!inCriticalRadius(p)) failed = true;
        if (!passesTrappingProbability(p)) failed = true;
        if (!passesDechanneling(p, pv)) failed = true;

        if (!failed) {
            count++;
            gammaSum += p.Gamma();
        }
    }
    f.Close();

    if (n == 0) return { 0.0, 0.0, 0 };
    double efficiency = static_cast<double>(count) / static_cast<double>(n);
    return { efficiency, gammaSum, n };
}

int EdizEndFirstTake()
{
    const char* fileLambda = "events.root";
    const char* fileXib = "eventsXib.root";

    // Process both files to get channeling efficiencies, sum of gamma, and total entries
    Result resLambda = processFile(fileLambda);
    Result resXib = processFile(fileXib);

    std::cout << "Lambda channeling efficiency: " << resLambda.channelingEfficiency << std::endl;
    std::cout << "Xi_b channeling efficiency: " << resXib.channelingEfficiency << std::endl;

    // Calculate average gamma for Lambda and Xi_b (sum gamma / total entries)
    double gamma_prev = (resLambda.entries > 0) ? (resLambda.gammaSum / resLambda.entries) : 0.0;
    double gamma_target = (resXib.entries > 0) ? (resXib.gammaSum / resXib.entries) : 0.0;

    std::cout << "Average gamma Lambda: " << gamma_prev << std::endl;
    std::cout << "Average gamma Xi_b: " << gamma_target << std::endl;

    // Cross sections at 13 TeV
    double sigma_bb = 495e-6;   // σ(pp→b b̄ X) = 495 μb
    double sigma_cc = 2369e-6;  // σ(pp→c c̄ X) in 0<pT<8, 2<y<4.5 = 2369 μb

    std::cout << "σ(pp→b b̄ X) at 13 TeV = " << sigma_bb << " b" << std::endl;
    std::cout << "σ(pp→c c̄ X) (charm)    = " << sigma_cc << " b" << std::endl;

    // c→Λ_c branching and fragmentation fractions these were the closest ones I could find from the PDG not sure if correct
    double f_c_to_Lambda_c = 0.18;
    double BR_Lambda_1 = 0.3067;  // Λ₁ → Λ K⁻ π⁺
    double BR_Lambda_2 = 0.1772;  // Λ₁ → Ξ⁻ π⁺ π⁺ π⁻
    double BR_Lambda_3 = 0.1772;  // Λ₁ → Λ π⁺ π⁺ π⁻
    double BR_Lambda_4 = 0.0632;  // Λ₁ → Ξ⁻ π⁺
    double BR_Lambda_5 = 0.2182;  // Λ₁ → Λ K⁻ π⁺ π⁺
    double BR_Lambda_6 = 0.1772;  // Λ₁ → Ξ⁻ π⁺ π⁺


    double N_Lambda_1 = computeExpectedEvents(sigma_cc, f_c_to_Lambda_c, BR_Lambda_1, resLambda.channelingEfficiency);
    double N_Lambda_2 = computeExpectedEvents(sigma_cc, f_c_to_Lambda_c, BR_Lambda_2, resLambda.channelingEfficiency);
    double N_Lambda_3 = computeExpectedEvents(sigma_cc, f_c_to_Lambda_c, BR_Lambda_3, resLambda.channelingEfficiency);
    double N_Lambda_4 = computeExpectedEvents(sigma_cc, f_c_to_Lambda_c, BR_Lambda_4, resLambda.channelingEfficiency);
    double N_Lambda_5 = computeExpectedEvents(sigma_cc, f_c_to_Lambda_c, BR_Lambda_5, resLambda.channelingEfficiency);
    double N_Lambda_6 = computeExpectedEvents(sigma_cc, f_c_to_Lambda_c, BR_Lambda_6, resLambda.channelingEfficiency);

    double N_Lambda_c = N_Lambda_1 + N_Lambda_2 + N_Lambda_3 + N_Lambda_4 + N_Lambda_5 + N_Lambda_6;

    std::cout << "Expected detected events Λ_c (sum of 6 channels): " << N_Lambda_c << std::endl;

    // Xi_b branching ratios and fragmentation ratio
    double f_q_to_H = 0.054;
    double B_Xib_Jpsi_Xi = 1.02e-5;
    double B_Xib_Jpsi_Lambda_K = 2.5e-6;
    double B_Xib_pbar_K_K = 2.3e-6;
    double B_Xib_pbar_pi_pi = 1.3e-6;
    double B_Xib_pbar_K_pi = 2.3e-6;
    double B_Xib_Lambdab_pi = 7.0e-4;
    double B_Xib_Xic_Ds = 1.9e-3;
    double B_Xib_Xi_gamma = 1.3e-4;

    // Calculate detected events per Xi_b channel
    double N_Xib_1 = computeExpectedEvents(sigma_bb, f_q_to_H, B_Xib_Jpsi_Xi, resXib.channelingEfficiency);
    double N_Xib_2 = computeExpectedEvents(sigma_bb, f_q_to_H, B_Xib_Jpsi_Lambda_K, resXib.channelingEfficiency);
    double N_Xib_3 = computeExpectedEvents(sigma_bb, f_q_to_H, B_Xib_pbar_K_K, resXib.channelingEfficiency);
    double N_Xib_4 = computeExpectedEvents(sigma_bb, f_q_to_H, B_Xib_pbar_pi_pi, resXib.channelingEfficiency);
    double N_Xib_5 = computeExpectedEvents(sigma_bb, f_q_to_H, B_Xib_pbar_K_pi, resXib.channelingEfficiency);
    double N_Xib_6 = computeExpectedEvents(sigma_bb, f_q_to_H, B_Xib_Lambdab_pi, resXib.channelingEfficiency);
    double N_Xib_7 = computeExpectedEvents(sigma_bb, f_q_to_H, B_Xib_Xic_Ds, resXib.channelingEfficiency);
    double N_Xib_8 = computeExpectedEvents(sigma_bb, f_q_to_H, B_Xib_Xi_gamma, resXib.channelingEfficiency);

    double N_Xib = N_Xib_1 + N_Xib_2 + N_Xib_3 + N_Xib_4 + N_Xib_5 + N_Xib_6 + N_Xib_7 + N_Xib_8;

    std::cout << "Expected detected events Xi_b (sum of 8 channels): " << N_Xib << std::endl;

    double sigma_target = computeSensitivity(
        gamma_prev, g_theta_c, N_Lambda_c,
        gamma_target, g_theta_c, N_Xib,
        sigma_prev_silicium
    );

    std::cout << "Estimated sensitivity for target particle: " << sigma_target << " e·cm" << std::endl;

    return 0;
}
