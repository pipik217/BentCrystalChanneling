/////////////////////////////////
/////////   CONCLUSION //////////
/////////////////////////////////
//  Using this KDE for a 5D data
// doesn't work. It is too many dim
// I think it is bc the integral.
// But plotting the pdf looks really strange
// with perfect gaussians. It does not look like one smeared kernel per event. its weid
/////////////////////////////////









// multiplyEvents.C
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <vector>
#include <string>

#include <RooRealVar.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooNDKeysPdf.h>
#include <RooFit.h>

#include <chrono>
#include <iostream>

void multiplyevents() {
  // ───── OPTIONS ────────────────────────────────────────────────
  const char* inputFile   = "events.root";
  const char* treeName    = "tree";
  const char* outputFile  = "eventsMultiplied.root";
  const Int_t  factor     = 7;     // how many× more events
  const Double_t rho      = 2;    // kernel bandwidth (0<rho<1)
  // ─────────────────────────────────────────────────────────────

  // Open input TFile & TTree
  TFile inF(inputFile, "READ");
  TTree* t = (TTree*)inF.Get(treeName);

  // 1) Discover double branches
  TObjArray* blist = t->GetListOfBranches();
  std::vector<std::string> names;
  std::vector<Double_t>     vals;
  for (int i = 0; i < blist->GetEntries(); ++i) {
    TBranch* b = (TBranch*)blist->At(i);
    TLeaf*   L = b->GetLeaf(b->GetName());
    if (L && std::string(L->GetTypeName()) == "Double_t") {
      names.push_back(b->GetName());
      vals.push_back(0.);
    }
  }
  const int nvar = names.size();

  // 2) Create RooRealVars (auto-ranged from tree min/max)
  RooArgList varList;
  for (int i = 0; i < nvar; ++i) {
    // attach temporarily to get min/max
    t->SetBranchAddress(names[i].c_str(), &vals[i]);
    Double_t mn = t->GetMinimum(names[i].c_str());
    Double_t mx = t->GetMaximum(names[i].c_str());
    t->ResetBranchAddresses();

    RooRealVar* rv = new RooRealVar(names[i].c_str(),
                                    names[i].c_str(),
                                    mn, mx);
    varList.add(*rv);
  }

  // 3) Import TTree into RooDataSet
  //    Re-attach branches for reading
  for (int i = 0; i < nvar; ++i)
    t->SetBranchAddress(names[i].c_str(), &vals[i]);
  RooDataSet data("data", "input dataset", varList,
                  RooFit::Import(*t));


  // -------------------------------
  // Print out the first few rows with numeric values:
  Int_t toPrint = 2;
  for (Int_t i = 0; i < toPrint; ++i) {
    const RooArgSet* row = data.get(i);
    std::cout << "  row[" << i << "]  ";
    for (int j = 0; j < varList.getSize(); ++j) {
      RooRealVar* rv = dynamic_cast<RooRealVar*>(varList.at(j));
      if (!rv) continue;
      const char* name = rv->GetName();
      Double_t val = row->getRealValue(name);
      std::cout << name << " = " << val;
      if (j < varList.getSize()-1) std::cout << ",  ";
    }
    std::cout << "\n";
  }
  std::cout << "────────────────────────────────────────────\n";





  // 4) Build the ND-kernel pdf and generate new events

  RooNDKeysPdf pdf("pdfND", "smoothed PDF",
                   varList, data, 
                   "ma", 
                   rho, 3);

pdf.getParameters(RooArgSet())->Print();
pdf.Print();



//////////////////////////////////////////
//////////   TEST  ///////////////////////

RooArgSet allVars(varList);
Double_t integral = pdf.createIntegral(allVars)->getVal();
std::cout << "PDF integral over full space = " << integral << "\n";

//Int_t toPrint = std::min<Int_t>(10, data.numEntries());
for (Int_t i = 0; i < toPrint; ++i) {
  const RooArgSet* row = data.get(i);
  std::cout << "  row[" << i << "] = ";
  row->Print(); 
}

// 1) Pick one variable to check, e.g. the first RooRealVar in varList:
RooRealVar* xvar = dynamic_cast<RooRealVar*>(varList.at(4));
if (!xvar) {
  std::cerr << "[!] Error: varList.at(0) is not a RooRealVar.\n";
  return;
}

// 2) Build the 1D projection of your current N‐dimensional KDE
//    (here “pdf” is whatever RooNDKeysPdf object you already constructed)
RooPlot* frame1 = xvar->frame();
//data.plotOn(frame1, RooFit::Binning(50), RooFit::LineColor(kBlack), RooFit::Normalization(1.0, RooAbsReal::Relative));
pdf.plotOn(frame1, RooFit::LineColor(kRed), RooFit::Normalization(1.0, RooAbsReal::Relative));

double nEvents = data.numEntries();  
// pdf.plotOn(frame1,
//            RooFit::LineColor(kRed),
//            RooFit::Normalization(nEvents, RooAbsReal::NumEvent));

// 3) Draw it and save it to a PNG so you can visually inspect edge/tail behavior:
TCanvas c1("c1", "KDE 1D projection", 800, 600);
frame1->Draw();
c1.SaveAs("kde_projection_rho.png");

// 4) Now repeat for an “extreme” tiny‐kernel version (ρ→1e−6, nSigma=0):
const Double_t tinyRho = 10;
RooNDKeysPdf pdf_tiny(
  "pdf_tiny", "tiny‐kernel KDE",
  varList, data,
  "ma",            // no mirroring
  tinyRho,        // ρ
  3
//   3,            // nSigma = 0
//   false,          // rotate = false
//   false           // sortInput = false
);

RooPlot* frame2 = xvar->frame();
//data.plotOn(frame2, RooFit::Binning(50), RooFit::LineColor(kBlack));
pdf_tiny.plotOn(frame2, RooFit::LineColor(kBlue), RooFit::Normalization(nEvents, RooAbsReal::NumEvent));

TCanvas c2("c2", "Tiny‐kernel KDE projection", 800, 600);
frame2->Draw();
c2.SaveAs("kde_projection_tiny.png");


///////////////////////////////////////////7


auto t0 = std::chrono::high_resolution_clock::now();

RooAbsData* gen = pdf.generate(varList,
                    factor * t->GetEntries(), true, false);

  std::cout 
    << "Build+Generate took " 
    << std::chrono::duration_cast<std::chrono::milliseconds>(
         std::chrono::high_resolution_clock::now() - t0
       ).count() 
    << " ms\n";

  // 5) Write out to new TTree
  TFile outF(outputFile, "RECREATE");
  TTree* tout = new TTree(treeName, "multiplied events");

  std::vector<Double_t> newvals(nvar);
  for (int i = 0; i < nvar; ++i) {
    tout->Branch(names[i].c_str(), &newvals[i]);
  }

  // Loop over generated dataset
  for (int i = 0; i < gen->numEntries(); ++i) {
    const RooArgSet* row = gen->get(i);
    for (int j = 0; j < nvar; ++j) {
      newvals[j] = row->getRealValue(names[j].c_str());
    }
    tout->Fill();
  }

  tout->Write();
  outF.Close();
  inF.Close();
}
