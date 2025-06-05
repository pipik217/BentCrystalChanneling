#include "propertiesImproved.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <cmath>
#include <fstream>
#include <iostream>

namespace Eric
{
    Double_t getThetaY(const TLorentzVector& p);
    Double_t getPV(const TLorentzVector& p);
    Double_t getThetaLindhard(Double_t pv, Bool_t silicium);
    Double_t getAStraight(Double_t thetaY, Double_t thetaLindhard, Bool_t isSilicium);

    Bool_t inLindhard(Double_t thetaY, Double_t thetaLindhard); //condition 1
    Bool_t decayAfterCrystal(Double_t l, Double_t lengthCrystal); //condition 2
    Bool_t inCriticalRadius(Double_t Rc, Double_t R); // condition 3
    Bool_t passesTrappingProbability(Double_t thetaY, Double_t thetaLindhard, Double_t RRatio, Bool_t isSilicium); //condition 4
    Bool_t passesDechanneling(Double_t gamma, Double_t RRatio, Double_t thetaCrystal, Bool_t isSilicium); //condition 5
}

void matrixMakerImproved(bool skipOnFail = true)
{   
    //this block reads in the data for analysis.

    //change name here to decide which file to read the events from
    TFile f("eventsMultiplied.root", "READ");
    TTree* t = (TTree*)f.Get("tree");
    auto outFile = std::unique_ptr<TFile>(TFile::Open("outfile.root", "RECREATE"));

    constexpr Int_t lengthMin{1};           //cm
    constexpr Double_t lengthMin_m{lengthMin / 100.0}; //m
    constexpr Int_t lengthMax{20};          //cm
    constexpr Double_t lengthMax_m{(lengthMax) / 100.0};          //m
    constexpr Double_t dLengthStep{1.0};    //cm
    constexpr Int_t nLengthStep{Eric::returnNumberOfSteps(lengthMin, lengthMax, dLengthStep)}; //cm   

    constexpr Int_t thetaMin{5};        //mrad
    constexpr Double_t thetaMin_rad{thetaMin / 1000.0}; //rad
    constexpr Int_t thetaMax{20};       //mrad
    constexpr Double_t thetaMax_rad{(thetaMax) / 1000.0}; //rad
    constexpr Double_t dThetaStep{1.0}; //mrad
    constexpr Int_t nThetaStep{Eric::returnNumberOfSteps(thetaMin, thetaMax, dThetaStep)};


    
    auto hSilicium = new TH2I("histogram Silicium", "XiB Silicium", nLengthStep, lengthMin_m, lengthMax_m, nThetaStep, thetaMin_rad, thetaMax_rad);
    auto hGermanium = new TH2I("histogram Germanium", "XiB Germanium", nLengthStep, lengthMin_m, lengthMax_m, nThetaStep, thetaMin_rad, thetaMax_rad);
    
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(33);
    gStyle->SetFrameFillColor(18);  
    hSilicium->SetFillColor(46);
    hGermanium->SetFillColor(46);
    TPaveLabel pl;
    
    Double_t px, py, pz, E, l;
    t->SetBranchAddress("px", &px);
    t->SetBranchAddress("py", &py);
    t->SetBranchAddress("pz", &pz);
    t->SetBranchAddress("E", &E);
    t->SetBranchAddress("l", &l);

    //line stores the amount of events we have
    Long64_t n = t->GetEntries();

    //first loop goes over all the events
    for (Long64_t i = 0; i < n; ++i) 
    {
        if(!(i%(n/10))) std::cout << i / static_cast<Double_t>(n) << "\n";        
        t->GetEntry(i);
        TLorentzVector p(px, py, pz, E);

        Bool_t fail{false};

        //these values don't depend on the crystal and only the particle so they are calculated before the other loops
        Double_t thetaY{Eric::getThetaY(p)};
        Double_t pv{Eric::getPV(p)};

        for(Int_t binaryOption{0}; binaryOption <= 1; ++binaryOption)
        {
            //
            Bool_t isSilicium{static_cast<Bool_t>(binaryOption)};
            Double_t thetaLindhard{Eric::getThetaLindhard(pv, isSilicium)};
            
            Bool_t lindhardFail{!Eric::inLindhard(thetaY, thetaLindhard)};
            if(lindhardFail)
            {        
                if (skipOnFail) continue;
            }

                 
            
            for(Int_t iii{0}; iii < nLengthStep; ++iii)
            {
                Double_t lengthCrystal{(lengthMin + iii * dLengthStep) / 100.0}; //cm to m
                Bool_t decayFail{!Eric::decayAfterCrystal(l, lengthCrystal)};
                if (decayFail)
                {                                      
                    if (skipOnFail) continue;
                }



                for(Int_t jjj{0}; jjj <= nThetaStep; ++jjj)
                {
                    Double_t thetaCrystal{(thetaMin + jjj * dThetaStep) / 1000.0}; //mrad to rad
                    Double_t Rc{pv / (isSilicium ? Silicium::Uprime : Germanium::Uprime)};
                    Double_t R{(lengthCrystal / thetaCrystal) * 100.0}; // m to cm
                    Double_t RRatio{Rc / R};

                    Bool_t criticalRadiusFail{!Eric::inCriticalRadius(Rc, R)};
                    if (criticalRadiusFail)
                    {
                        if (skipOnFail) continue;
                    }
                    
                    Bool_t trappingFail{!Eric::passesTrappingProbability(thetaY, thetaLindhard, RRatio, isSilicium)};
                    if (trappingFail)
                    {
                        if (skipOnFail) continue;
                    }

                    Bool_t channelingFail{!Eric::passesDechanneling(p.Gamma(), RRatio, thetaCrystal, isSilicium)};
                    if (channelingFail)
                    {
                        if (skipOnFail) continue;
                    }
                    if (isSilicium)
                    {
                        hSilicium->Fill(lengthCrystal, thetaCrystal);
                    }
                    else
                    {
                        hGermanium->Fill(lengthCrystal, thetaCrystal);
                    }
                    //std::cout << lengthCrystal << " " << thetaCrystal << " " << gRandom->Uniform() << "\n";
                }
            }
        }
    }
    auto c2h = new TCanvas("c2h", "Xib", 1200, 600);
    
    c2h->Divide(2,1);
    c2h->cd(1);
    hSilicium->DrawClone();
    c2h->cd(2);
    hGermanium->DrawClone();
    c2h->Update();
    /*
    outFile->WriteObject(hSilicium, hSilicium->GetName());
    outFile->WriteObject(hGermanium, hGermanium->GetName());
    */ 
}

namespace Eric
{
    Double_t getThetaY(const TLorentzVector& p)
    {
        return TMath::Abs(p.Py() / p.Pz());
    }

    Double_t getPV(const TLorentzVector& p)
    {
        return p.P() * p.Beta();
    }

    Double_t getThetaLindhard(Double_t pv, Bool_t isSilicium)
    {
        return std::sqrt(2 * (isSilicium ? Silicium::Uxc : Germanium::Uxc) / pv);
    }

    Double_t getAStraight(Double_t thetaY, Double_t thetaLindhard, Bool_t isSilicium)
    {
        return 2 * (isSilicium ? (Silicium::xc / Silicium::d) : (Germanium::xc / Germanium::d)) * std::sqrt(1 - std::pow(thetaY / thetaLindhard, 2));   
    }

    Bool_t inLindhard(Double_t thetaY, Double_t thetaLindhard) // condition 1
    {
        return thetaY < thetaLindhard;
    }

    Bool_t decayAfterCrystal(Double_t l, Double_t lengthCrystal) // condition 2 
    {
        return (l > lengthCrystal);
    }

    Bool_t inCriticalRadius(Double_t Rc, Double_t R) //condition 3
    {
        return Rc <= R;
    }

    Bool_t passesTrappingProbability(Double_t thetaY, Double_t thetaLindhard, Double_t RRatio, Bool_t isSilicium) //condition 4
    {
        Double_t ABent = Eric::getAStraight(thetaY, thetaLindhard, isSilicium) * (1.0 - RRatio);
        if(ABent < 0.0 || ABent > 1.0) return kFALSE;
        return (gRandom->Uniform() <= ABent);
    }

    Bool_t passesDechanneling(Double_t gamma, Double_t RRatio, Double_t thetaCrystal, Bool_t isSilicium) //condition 5
    {
        Double_t I{};
        Double_t N{};
        Double_t Z{};
        Double_t aTF_cm{};
        Double_t d_cm{};        
        if (isSilicium)
        {
            I = Silicium::I;
            N = Silicium::N;
            Z = Silicium::Z;
            aTF_cm = Silicium::aTF_cm;
            d_cm = Silicium::d_cm;
        }
        else
        {
            I = Germanium::I;
            N = Germanium::N;
            Z = Germanium::Z;
            aTF_cm = Germanium::aTF_cm;
            d_cm = Germanium::d_cm;
        }

        Double_t logArg{(2.0 * Constants::me_c2 * gamma) / I};
        if (logArg <= 1.0) 
        {
            printf("[Dech] logTerm invalid: logArg=%.3e\n", logArg);
            return kFALSE;
        }

        Double_t thetaD{Constants::prefactor * N * Z * aTF_cm * std::pow(d_cm, 2) / (std::log(logArg) - 1.0)};
        Double_t exponent{thetaCrystal / (thetaD * RRatio * std::pow(1.0 - RRatio, 2))};
        Double_t w{std::exp(-exponent)};

        return (gRandom->Uniform() <= w);
    }
}
