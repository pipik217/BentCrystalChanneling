#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <cmath>
#include <iostream>

namespace Config {
    constexpr double thetaCrystal = 0.010;   // rad
    constexpr double lengthCrystal = 0.10;   // m
    constexpr double Uprime = 5.7;          // GeV/cm
    constexpr double A_straight = 0.75;      // this value is not dependant on anything
    constexpr double thetaLindhard = 7e-6;   // rad

    // material properties Si(110)
    constexpr double Uxc = 16e-9;            // GeV
    constexpr double aTF = 0.194;             // Å
    constexpr double d = 1.92;               // Å
    constexpr double rho = 2.329;            // g/cm³
    constexpr double Z = 14.0;
    constexpr double A = 28.0855;            // g/mol
    constexpr double I = 173e-9;             // GeV
}

Bool_t inLindhard(const TLorentzVector& p) {
    Double_t thetaY = p.Py() / p.Pz();
    return (TMath::Abs(thetaY) < Config::thetaLindhard);
}

Bool_t decayAfterCrystal(Double_t l) {
    return (l > Config::lengthCrystal);
}

Bool_t inCriticalRadius(const TLorentzVector& p) {
    Double_t pv = p.P() * p.Beta(); // GeV
    Double_t Rc = pv / Config::Uprime;      // cm
    Double_t R = (Config::lengthCrystal / Config::thetaCrystal) * 100.0; // m to cm
    return Rc <= R;
}

Bool_t passesTrappingProbability(const TLorentzVector& p) {
    Double_t pv = p.P() * p.Beta(); // GeV
    Double_t Rc = pv / Config::Uprime;
    Double_t R = (Config::lengthCrystal / Config::thetaCrystal) * 100.0; // cm
    Double_t A_bent = Config::A_straight * (1.0 - Rc / R);
    if (A_bent < 0.0 || A_bent > 1.0) return kFALSE;
    return (gRandom->Uniform() <= A_bent);
}

Bool_t passesDechanneling(Double_t E, Double_t pv) {
    const Double_t me_c2 = 0.000511;     // GeV
    const Double_t NA = 6.022e23;        // mol⁻¹
    const Double_t pi = TMath::Pi();
    const Double_t prefactor = 256.0 / (9.0 * pi);

    Double_t aTF_cm = Config::aTF * 1e-8;
    Double_t d_cm = Config::d * 1e-8;

    Double_t gamma = E / me_c2;
    Double_t N = NA * Config::rho / Config::A;

    Double_t logArg = (2.0 * me_c2 * gamma) / Config::I;
    if (logArg <= 1.0) {
        printf("[Dech] logTerm invalid: logArg=%.3e\n", logArg);
        return kFALSE;
    }

    Double_t logTerm = log(logArg);
    Double_t thetaD = prefactor * N * Config::Z * aTF_cm * pow(d_cm, 2) / logTerm;
    Double_t thetaC = sqrt(2.0 * Config::Uxc / E);     // rad

    Double_t Rc = pv / Config::Uprime;                 // cm
    Double_t R = (Config::lengthCrystal / Config::thetaCrystal) * 100.0; // cm

    if (Rc >= R) return kFALSE;

    Double_t frac = Rc / R;
    Double_t exponent = (thetaC / thetaD) * frac * pow(1.0 - frac, 2);
    Double_t w = exp(-exponent);

    printf("[Dech] thetaC=%.2e, thetaD=%.2e, Rc=%.2e, R=%.2e, w=%.4f\n",
        thetaC, thetaD, Rc, R, w);

    return (gRandom->Uniform() <= w);
}

// false (more data) true (faster)
void EdizChanneling(bool skipOnFail = true) {
    TFile f("events.root", "READ");
    TTree* t = (TTree*)f.Get("tree");

    Double_t px, py, pz, E, l;
    t->SetBranchAddress("px", &px);
    t->SetBranchAddress("py", &py);
    t->SetBranchAddress("pz", &pz);
    t->SetBranchAddress("E", &E);
    t->SetBranchAddress("l", &l);

    Int_t count = 0;
    Int_t fail_Lindhard = 0;
    Int_t fail_Decay = 0;
    Int_t fail_CriticalRad = 0;
    Int_t fail_TrapProb = 0;
    Int_t fail_Dechannel = 0;

    Long64_t n = t->GetEntries();

    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        TLorentzVector p(px, py, pz, E);
        Double_t pv = p.P() * p.Beta();

        bool failed = false;

        if (!inLindhard(p)) {
            fail_Lindhard++;
            if (skipOnFail) continue;
            failed = true;
        }
        if (!decayAfterCrystal(l)) {
            fail_Decay++;
            if (skipOnFail) continue;
            failed = true;
        }
        if (!inCriticalRadius(p)) {
            fail_CriticalRad++;
            if (skipOnFail) continue;
            failed = true;
        }
        if (!passesTrappingProbability(p)) {
            fail_TrapProb++;
            if (skipOnFail) continue;
            failed = true;
        }
        if (!passesDechanneling(E, pv)) {
            fail_Dechannel++;
            if (skipOnFail) continue;
            failed = true;
        }

        if (!failed) {
            count++;
        }
    }

    printf("\n--- Channeling Summary ---\n");
    printf("Total particles: %lld\n", n);
    printf("Channeled:       %d (%.6f)\n", count, (double)count / (double)n);
    printf("Failed Lindhard: %d\n", fail_Lindhard);
    printf("Failed Decay:    %d\n", fail_Decay);
    printf("Failed Rc > R:   %d\n", fail_CriticalRad);
    printf("Failed TrapProb: %d\n", fail_TrapProb);
    printf("Failed Dech:     %d\n", fail_Dechannel);
}