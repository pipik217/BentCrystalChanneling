#include "propertiesNeutral.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <cmath>
#include <fstream>
#include <iostream>

namespace eric
{
    double getPV(const TLorentzVector &p)
    {
        return p.P() * p.Beta();
    }    

    double getThetaLindhard(const TLorentzVector& p, bool silicium)
    {
    double Uxc{};    
    if (silicium)
    {
        Uxc = silicium::Uxc;
    }
    else
    {
        Uxc = germanium::Uxc;
    }
    return std::sqrt(2 * Uxc / getPV(p));
    }
}

namespace original{
//function returns whether a particle is trapped or not given its angle and the Lindhardt angle
    Bool_t inLindhard(const TLorentzVector& p, bool silicium) //1 
    {
        Double_t thetaY = p.Py() / p.Pz();
        return (TMath::Abs(thetaY) < eric::getThetaLindhard(p, silicium));
    }

//function returns whether a particle doesn't decay before exciting the crystal
    Bool_t decayAfterCrystal(Double_t l, Double_t lengthCrystal) //2 
    {
        return (l > lengthCrystal);
    }
}

namespace eric
{
    double getThetaY(const TLorentzVector &p)
    {
        return p.Py() / p.Pz();
    }

    double getAStraight(const TLorentzVector& p, double Xc, double Dp, double thetaY, double thetaLindhard, bool silicium)
    {
        return ((2 * Xc / Dp) * std::sqrt(1 - pow(thetaY / eric::getThetaLindhard(p, silicium),2)));
    }
}

//function returns whether a particle isn't too fast to be able to be bend by the crystal
Bool_t inCriticalRadius(const TLorentzVector& p, double lengthCrystal, double thetaCrystal, bool silicium) //3
    {
    Double_t pv = eric::getPV(p); // GeV = momentum(GeV/c) * beta(fraction times c)

    double Uprime{};
    if (silicium)
    {
        Uprime = silicium::Uprime;
    }
    else
    {
        Uprime = germanium::Uprime;
    }

    Double_t Rc = pv / Uprime;      // cm
    Double_t R = (lengthCrystal / thetaCrystal) * 100.0; // m to cm
    return Rc <= R;
}

Bool_t passesTrappingProbability(const TLorentzVector& p, double lengthCrystal, double thetaCrystal, bool silicium) //4
{
    Double_t pv = eric::getPV(p); // GeV

    double Uprime{};
    double xc{};
    double d{};
    if (silicium)
    {
        Uprime = silicium::Uprime;
        xc = silicium::xc;
        d = silicium::d;
    }
    else
    {
        Uprime = germanium::Uprime;
        xc = germanium::xc;
        d = germanium::d;
    }

    Double_t Rc = pv / Uprime;
    Double_t R = (lengthCrystal / thetaCrystal) * 100.0; // cm
    Double_t A_bent = eric::getAStraight(p, xc, d, eric::getThetaY(p), eric::getThetaLindhard(p, silicium) * (1.0 - Rc / R), silicium);
    if (A_bent < 0.0 || A_bent > 1.0) return kFALSE;
    return (gRandom->Uniform() <= A_bent);
}

Bool_t passesDechanneling(const TLorentzVector& p, Double_t pv, double lengthCrystal, double thetaCrystal, bool silicium) //5
{
    double aTF{};
    double d{};
    double rho{};
    double A{};
    double I{};
    double Z{};
    double Uprime{};
    if (silicium)
    {
        aTF = silicium::aTF;
        d = silicium::d;
        rho = silicium::rho;
        A = silicium::A;
        I = silicium::I;
        Z = silicium::Z;
        Uprime = silicium::Uprime;
    }
    else
    {
        aTF = germanium::aTF;
        d = germanium::d;
        rho = germanium::rho;
        A = germanium::A;
        I = germanium::I;
        Z = germanium::Z;
        Uprime = germanium::Uprime;
    }
    Double_t aTF_cm = aTF * Constants::angstromcmRatio;
    Double_t d_cm = d * Constants::angstromcmRatio;

    Double_t gamma = p.Gamma();
    Double_t N = Constants::NA * rho / A;

    Double_t logArg = (2.0 * Constants::me_c2 * gamma) / I;
    if (logArg <= 1.0) {
        printf("[Dech] logTerm invalid: logArg=%.3e\n", logArg);
        return kFALSE;
    }

    Double_t logTerm = log(logArg);
    Double_t thetaD = Constants::prefactor * N * Z * aTF_cm * pow(d_cm, 2) / (logTerm - 1);
    //Double_t thetaC = sqrt(2.0 * Config::Uxc / E);     // rad

    Double_t Rc = pv / Uprime;                 // cm
    Double_t R = (lengthCrystal / thetaCrystal) * 100.0; // cm

    if (Rc >= R) return kFALSE;

    Double_t frac = Rc / R;
    Double_t exponent = thetaCrystal / (thetaD * frac * pow(1.0 - frac, 2));
    Double_t w = exp(-exponent);

    //printf("[Dech] thetaC=%.2e, thetaD=%.2e, Rc=%.2e, R=%.2e, w=%.4f\n",thetaC, thetaD, Rc, R, w);

    return (gRandom->Uniform() <= w);
}



// false (more data) true (faster)
void matrixMaker(bool skipOnFail = false) 
{
    TFile f("eventsXib.root", "READ");
    TTree* t = (TTree*)f.Get("tree");

    Double_t px, py, pz, E, l;
    t->SetBranchAddress("px", &px);
    t->SetBranchAddress("py", &py);
    t->SetBranchAddress("pz", &pz);
    t->SetBranchAddress("E", &E);
    t->SetBranchAddress("l", &l);

    Int_t count{};
    Int_t fail_Lindhard{};
    Int_t fail_Decay{};
    Int_t fail_CriticalRad{};
    Int_t fail_TrapProb{};
    Int_t fail_Dechannel{};

    bool silicium{};
    

    constexpr double thetaStep{0.0005};
    constexpr int nThetaStep{21};
    constexpr double lengthStep{0.01};
    constexpr int nLengthStep{21};

    double thetaCrystal = 0.01 - thetaStep;   // rad
    double lengthCrystal{};   // m

    int siliciumMatrix[nThetaStep][nLengthStep];
    int germaniumMatrix[nThetaStep][nLengthStep];

    for(int iii{0}; iii < nThetaStep; ++iii)
    {
        thetaCrystal += thetaStep;
        lengthCrystal = 0.0;
        for(int jjj{0}; jjj < nLengthStep; ++jjj)
        {
            lengthCrystal += lengthStep;
            for(int binaryoption{0}; binaryoption <= 1; ++binaryoption)
            {
                silicium = static_cast<bool>(binaryoption);
                count = 0;
                fail_Lindhard = 0;
                fail_Decay = 0;
                fail_CriticalRad = 0;
                fail_TrapProb = 0;
                fail_Dechannel = 0;
                Long64_t n = t->GetEntries();

                for (Long64_t i = 0; i < n; ++i) 
                {
                    t->GetEntry(i);
                    TLorentzVector p(px, py, pz, E);
                    Double_t pv{eric::getPV(p)};
       
                    bool failed = false;

                    if (!original::inLindhard(p, silicium)) 
                    {
                        fail_Lindhard++;
                        if (skipOnFail) continue;
                        failed = true;
                    }
                    if (!original::decayAfterCrystal(l, lengthCrystal)) 
                    {
                        fail_Decay++;
                        if (skipOnFail) continue;
                        failed = true;
                    }
                    if (!inCriticalRadius(p, lengthCrystal, thetaCrystal, silicium)) 
                    {
                        fail_CriticalRad++;
                        if (skipOnFail) continue;
                        failed = true;
                    }

                    if (!passesTrappingProbability(p, lengthCrystal, thetaCrystal, silicium)) 
                    {
                        fail_TrapProb++;
                        if (skipOnFail) continue;
                        failed = true;
                    }
                    if (!passesDechanneling(p, pv, lengthCrystal, thetaCrystal, silicium)) 
                    {
                        fail_Dechannel++;
                        if (skipOnFail) continue;
                        failed = true;
                    }

                    if (!failed) 
                    {
                        count++;
                    }
                }

                printf("\n--- Channeling Summary ---\n");
                std::cout << "Material was " << (silicium ? "silicium" : "germanium") << ".\n";
                std::cout << "Crystallength was " << lengthCrystal << " and angle was " << thetaCrystal << ".\n";
                printf("Total particles: %lld\n", n);
                printf("Channeled:       %d (%.6f)\n", count, (double)count / (double)n);
                printf("Failed Lindhard: %d\n", fail_Lindhard);
                printf("Failed Decay:    %d\n", fail_Decay);
                printf("Failed Rc > R:   %d\n", fail_CriticalRad);
                printf("Failed TrapProb: %d\n", fail_TrapProb);
                printf("Failed Dech:     %d\n", fail_Dechannel);

                if (silicium)
                {
                    siliciumMatrix[iii][jjj] = count;
                }
                else
                {
                    germaniumMatrix[iii][jjj] = count;
                }
            }
        }
    }

    ofstream myfile;
    myfile.open("silicium.dat");
    std::cout << "silicium.\n";
    for(int iii{0}; iii < nLengthStep; ++iii)
    {
        for(int jjj{0}; jjj < nLengthStep; ++jjj)
        {
            std::cout << siliciumMatrix[iii][jjj] << "\t";
            myfile << siliciumMatrix[iii][jjj] << "\t";
        }
        std::cout << "\n";
        myfile << "\n";
    }
    myfile.close();
    myfile.open("germanium.dat");
    std::cout << "\ngermanium\n";
    for(int iii{0}; iii < nLengthStep; ++iii)
    {
        for(int jjj{0}; jjj < nLengthStep; ++jjj)
        {
            std::cout << germaniumMatrix[iii][jjj] << "\t";
            myfile << germaniumMatrix[iii][jjj] << "\t";
        }
        std::cout << "\n";
        myfile << "\n";
    }
    myfile.close();
}
