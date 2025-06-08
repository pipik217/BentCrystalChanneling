#include "properties.h"
#include "physics_constants.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <cmath>
#include <iostream>

// Use global theta_c from constants
const double g_theta_c = 0.1;

namespace eric {
    double getPV(const TLorentzVector& p) {
        return p.P() * p.Beta();
    }

    double getThetaLindhard(const TLorentzVector& p) {
        return std::sqrt(2 * Config::Uxc / getPV(p));
    }

    double getThetaY(const TLorentzVector& p) {
        return p.Py() / p.Pz();
    }

    double getAStraight(double Xc, double Dp, double thetaY, double thetaLindhard) {
        return ((2 * Xc / Dp) * std::sqrt(1 - pow(thetaY / thetaLindhard, 2)));
    }
}

namespace original {
    Bool_t inLindhard(const TLorentzVector& p) {
        Double_t thetaY = p.Py() / p.Pz();
        return (TMath::Abs(thetaY) < eric::getThetaLindhard(p));
    }

    Bool_t decayAfterCrystal(Double_t l) {
        return (l > Config::lengthCrystal);
    }
}

Bool_t inCriticalRadius(const TLorentzVector& p) {
    Double_t pv = eric::getPV(p);
    Double_t Rc = pv / Config::Uprime;
    Double_t R = (Config::lengthCrystal / Config::thetaCrystal) * 100.0; // cm
    return Rc <= R;
}

Bool_t passesTrappingProbability(const TLorentzVector& p) {
    Double_t pv = eric::getPV(p);
    Double_t Rc = pv / Config::Uprime;
    Double_t R = (Config::lengthCrystal / Config::thetaCrystal) * 100.0; // cm
    Double_t A_bent = eric::getAStraight(Config::xc, Config::d, eric::getThetaY(p), eric::getThetaLindhard(p)) * (1.0 - Rc / R);
    if (A_bent < 0.0 || A_bent > 1.0) return kFALSE;
    return (gRandom->Uniform() <= A_bent);
}

Bool_t passesDechanneling(const TLorentzVector& p, Double_t pv) {
    Double_t aTF_cm = Config::aTF * Constants::angstromcmRatio;
    Double_t d_cm = Config::d * Constants::angstromcmRatio;
    Double_t gamma = p.Gamma();
    Double_t N = Constants::NA * Config::rho / Config::A;
    Double_t logArg = (2.0 * Constants::me_c2 * gamma) / Config::I;
    if (logArg <= 1.0) return kFALSE;
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

double computeExpectedEvents(double sigma_qq, double f_q_to_H, double B_decay, double channeling_efficiency) {
    double N_ideal = sigma_qq * f_q_to_H * B_decay;
    return N_ideal * channeling_efficiency;
}

double computeSensitivity(double gamma_prev, double theta_c_prev, double N_prev,
    double gamma_target, double theta_c_target, double N_target,
    double sigma_prev) {
    if (N_prev <= 0 || N_target <= 0) return -1.0;
    double param_ratio = (gamma_prev * theta_c_prev) / (gamma_target * theta_c_target);
    double stat_ratio = std::sqrt(N_target / N_prev);
    return param_ratio * stat_ratio * sigma_prev;
}

struct Result {
    double channelingEfficiency;
    double gammaSum;
    Long64_t entries;
};

Result processFile(const char* filename) {
    TFile f(filename, "READ");
    if (f.IsZombie()) return { 0.0, 0.0, 0 };
    TTree* t = (TTree*)f.Get("tree");
    if (!t) return { 0.0, 0.0, 0 };

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
        if (!original::inLindhard(p) ||
            !original::decayAfterCrystal(l) ||
            !inCriticalRadius(p) ||
            !passesTrappingProbability(p) ||
            !passesDechanneling(p, pv)) continue;
        count++;
        gammaSum += p.Gamma();
    }
    f.Close();
    if (n == 0) return { 0.0, 0.0, 0 };
    return { static_cast<double>(count) / static_cast<double>(n), gammaSum, n };
}

int Sensitivity() {
    const char* fileLambda = "events.root";
    const char* fileXib = "eventsXib.root";
    Result resLambda = processFile(fileLambda);
    Result resXib = processFile(fileXib);

    double gamma_prev = (resLambda.entries > 0) ? (resLambda.gammaSum / resLambda.entries) : 0.0;
    double gamma_target = (resXib.entries > 0) ? (resXib.gammaSum / resXib.entries) : 0.0;

    double N_Lambda_c =
        computeExpectedEvents(PhysicsConstants::sigma_cc, PhysicsConstants::f_c_to_Lambda_c, PhysicsConstants::BR_Lambda_1, resLambda.channelingEfficiency) +
        computeExpectedEvents(PhysicsConstants::sigma_cc, PhysicsConstants::f_c_to_Lambda_c, PhysicsConstants::BR_Lambda_2, resLambda.channelingEfficiency) +
        computeExpectedEvents(PhysicsConstants::sigma_cc, PhysicsConstants::f_c_to_Lambda_c, PhysicsConstants::BR_Lambda_3, resLambda.channelingEfficiency) +
        computeExpectedEvents(PhysicsConstants::sigma_cc, PhysicsConstants::f_c_to_Lambda_c, PhysicsConstants::BR_Lambda_4, resLambda.channelingEfficiency) +
        computeExpectedEvents(PhysicsConstants::sigma_cc, PhysicsConstants::f_c_to_Lambda_c, PhysicsConstants::BR_Lambda_5, resLambda.channelingEfficiency) +
        computeExpectedEvents(PhysicsConstants::sigma_cc, PhysicsConstants::f_c_to_Lambda_c, PhysicsConstants::BR_Lambda_6, resLambda.channelingEfficiency);

    double N_Xib =
        computeExpectedEvents(PhysicsConstants::sigma_bb, PhysicsConstants::f_q_to_Xib, PhysicsConstants::BR_Xib_Jpsi_Xi, resXib.channelingEfficiency) +
        computeExpectedEvents(PhysicsConstants::sigma_bb, PhysicsConstants::f_q_to_Xib, PhysicsConstants::BR_Xib_Jpsi_Lambda_K, resXib.channelingEfficiency) +
        computeExpectedEvents(PhysicsConstants::sigma_bb, PhysicsConstants::f_q_to_Xib, PhysicsConstants::BR_Xib_pbar_K_K, resXib.channelingEfficiency) +
        computeExpectedEvents(PhysicsConstants::sigma_bb, PhysicsConstants::f_q_to_Xib, PhysicsConstants::BR_Xib_pbar_pi_pi, resXib.channelingEfficiency) +
        computeExpectedEvents(PhysicsConstants::sigma_bb, PhysicsConstants::f_q_to_Xib, PhysicsConstants::BR_Xib_pbar_K_pi, resXib.channelingEfficiency) +
        computeExpectedEvents(PhysicsConstants::sigma_bb, PhysicsConstants::f_q_to_Xib, PhysicsConstants::BR_Xib_Lambdab_pi, resXib.channelingEfficiency) +
        computeExpectedEvents(PhysicsConstants::sigma_bb, PhysicsConstants::f_q_to_Xib, PhysicsConstants::BR_Xib_Xic_Ds, resXib.channelingEfficiency) +
        computeExpectedEvents(PhysicsConstants::sigma_bb, PhysicsConstants::f_q_to_Xib, PhysicsConstants::BR_Xib_Xi_gamma, resXib.channelingEfficiency);

    double sigma_target = computeSensitivity(
        gamma_prev, g_theta_c, N_Lambda_c,
        gamma_target, g_theta_c, N_Xib,
        PhysicsConstants::sigma_prev_silicium
    );

    std::cout << "Estimated sensitivity for target particle: " << sigma_target << " e·cm" << std::endl;
    return 0;
}