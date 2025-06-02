#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <random>
#include <cmath>

//setup to use some random behaviour I need later
std::random_device rd;
std::mt19937 gen(rd());

Bool_t inLindhard(const TLorentzVector &p, Double_t lengthCrystal, Double_t thetaCrystal){
  Double_t thetaY = p.Py()/p.Pz();
  // thetaLindhard must be calculated for every event!
  Double_t thetaLindhard = 0.000007; // 7 microrad typical 
  return (TMath::Abs(thetaY) < thetaLindhard);
}

Bool_t decayAfterCrystal(Double_t l, Double_t lengthCrystal){
  return (l > lengthCrystal);
}

namespace eric
{
    //more setup for the random behavior    
    std::uniform_real_distribution<> dis(0.0, 1.0);    

    //function to calculate the theta
    double getThetaY(const TLorentzVector &p)
    {
        return p.Py() / p.Pz();
    }    

    function to calculate the Xc according to theory in paper
    constexpr double calculateXc(double Dp, double Atf)
    {
        return (Dp / 2 - Atf);
    }

    //defining some values we need
    constexpr double siliciumDp{1.92}; //angstrom
    constexpr double germaniumDp{2.00};

    constexpr double siliciumAtf{0.194};//angstom
    constexpr double germaniumAtf{0.148};

    constexpr double siliciumXc{eric::calculateXc(eric::siliciumDp, eric::siliciumAtf)};
    constexpr double germaniumXc{calculateXc(eric::germaniumDp, eric::germaniumAtf)};
    //function calculates the trapping efficiency or chance for parallel beam according to paper
    Double_t getTrappingEfficiency(double Xc, double Dp, double thetaY, double thetaLinhard)
    {
        return static_cast<Double_t>((2 * Xc / Dp) * std::sqrt(1 - pow(thetaY / thetaLinhard,2)));
    }
}

void channeling(){

  // Crystal parameters
  constexpr double thetaCrystal  = 0.010;   // Bending angle (rad)
  constexpr double lengthCrystal = 0.10;    // Crystal length (m)

  TFile  f("events.root","READ");
  TTree *t = (TTree*)f.Get("tree");

  Double_t px, py, pz, E, l;
  t->SetBranchAddress("px",&px);
  t->SetBranchAddress("py",&py);
  t->SetBranchAddress("pz",&pz);
  t->SetBranchAddress("E" ,&E);
  t->SetBranchAddress("l" ,&l);

  Int_t count = 0;
  Long64_t n = t->GetEntries();
  for(Long64_t i=0; i<n; ++i){
    t->GetEntry(i);
    TLorentzVector p(px,py,pz,E);
    

    // test channeling conditions
    Bool_t trapped = inLindhard(p, lengthCrystal, thetaCrystal);
    Double_t trappingEfficiency = eric::getTrappingEfficiency(eric::siliciumXc, eric::siliciumDp, eric::getThetaY(p), 0.000007);
    Bool_t expandedTrapped{eric::dis(gen) >= trappingEfficiency};
    Bool_t decaysAfter = decayAfterCrystal(l, lengthCrystal);
    //Bool_t channeled = trapped * decaysAfter * expandedTrapped; 
    Bool_t channeled = decaysAfter * expandedTrapped;
    
    if(channeled){ count++; }
  }
  
  printf("\nEfficiency of channeling: %f \n", (double) count / (double) n);
}

