#include "properties.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <cmath>
#include <iostream>

namespace original{
//function returns whether a particle is trapped or not given its angle and the Lindhardt angle
    Bool_t inLindhard(const TLorentzVector& p) 
    {
        Double_t thetaY = p.Py() / p.Pz();
        return (TMath::Abs(thetaY) < Config::thetaLindhard);
    }

//function returns whether a particle doesn't decay before exciting the crystal
    Bool_t decayAfterCrystal(Double_t l) {
        return (l > Config::lengthCrystal);
    }
}

namespace eric
{
    double getThetaY(const TLorentzVector &p)
    {
        return p.Py() / p.Pz();
    }

    double getAStraight(double Xc, double Dp, double thetaY, double thetaLindhard)
    {
        return ((2 * Xc / Dp) * std::sqrt(1 - pow(thetaY / thetaLindhard,2)));
    }
}

//function returns whether a particle isn't too fast to be able to be bend by the crystal
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
    Double_t A_bent = eric::getAStraight(Config::xc, Config::d, eric::getThetaY(p), Config::thetaLindhard) * (1.0 - Rc / R);
    if (A_bent < 0.0 || A_bent > 1.0) return kFALSE;
    return (gRandom->Uniform() <= A_bent);
}

Bool_t passesDechanneling(Double_t E, Double_t pv) {

    Double_t aTF_cm = Config::aTF * Constants::angstromcmRatio;
    Double_t d_cm = Config::d * Constants::angstromcmRatio;

    Double_t gamma = E / Constants::me_c2;
    Double_t N = Constants::NA * Config::rho / Config::A;

    Double_t logArg = (2.0 * Constants::me_c2 * gamma) / Config::I;
    if (logArg <= 1.0) {
        printf("[Dech] logTerm invalid: logArg=%.3e\n", logArg);
        return kFALSE;
    }

    Double_t logTerm = log(logArg);
    Double_t thetaD = Constants::prefactor * N * Config::Z * aTF_cm * pow(d_cm, 2) / logTerm;
    Double_t thetaC = sqrt(2.0 * Config::Uxc / E);     // rad

    Double_t Rc = pv / Config::Uprime;                 // cm
    Double_t R = (Config::lengthCrystal / Config::thetaCrystal) * 100.0; // cm

    if (Rc >= R) return kFALSE;

    Double_t frac = Rc / R;
    Double_t exponent = (thetaC / thetaD) * frac * pow(1.0 - frac, 2);
    Double_t w = exp(-exponent);

    //printf("[Dech] thetaC=%.2e, thetaD=%.2e, Rc=%.2e, R=%.2e, w=%.4f\n",thetaC, thetaD, Rc, R, w);

    return (gRandom->Uniform() <= w);
}



// false (more data) true (faster)
void EdizChannelingImproved(bool skipOnFail = false) {
    TFile f("eventsXib.root", "READ");
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

        if (!original::inLindhard(p)) {
            fail_Lindhard++;
            if (skipOnFail) continue;
            failed = true;
        }
        if (!original::decayAfterCrystal(l)) {
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
