#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

Bool_t inLindhard(const TLorentzVector &p, Double_t lengthCrystal, Double_t thetaCrystal){
  Double_t thetaY = p.Py()/p.Pz();
  // thetaLindhard must be calculated for every event!
  Double_t thetaLindhard = 0.000007; // 7 microrad typical 
  return (TMath::Abs(thetaY) < thetaLindhard);
}

Bool_t decayAfterCrystal(Double_t l, Double_t lengthCrystal){
  return (l > lengthCrystal);
}

void channeling(){

  // Crystal parameters
  const double thetaCrystal  = 0.010;   // Bending angle (rad)
  const double lengthCrystal = 0.10;    // Crystal length (m)

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
    Bool_t decaysAfter = decayAfterCrystal(l, lengthCrystal);
    // ...
    Bool_t channeled = trapped * decaysAfter; 
    
    if(channeled){ count++; }
  }
  
  printf("\nEfficiency of channeling: %f \n", (double) count / (double) n);
}

