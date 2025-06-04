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

void Edizmultiplyevents() {
    const char* inputFile = "events.root";
    const char* treeName = "tree";
    const char* outputFile = "eventsMultiplied.root";
    const Int_t factor = 7;
    const Double_t rho = 0.5;
    const Int_t maxEvents = 10;

    TFile inF(inputFile, "READ");
    TTree* t = (TTree*)inF.Get(treeName);

    TObjArray* blist = t->GetListOfBranches();
    std::vector<std::string> names;
    std::vector<Double_t> vals;
    for (int i = 0; i < blist->GetEntries(); ++i) {
        TBranch* b = (TBranch*)blist->At(i);
        TLeaf* L = b->GetLeaf(b->GetName());
        if (L && std::string(L->GetTypeName()) == "Double_t") {
            names.push_back(b->GetName());
            vals.push_back(0.);
        }
    }
    const int nvar = names.size();

    RooArgList varList;
    for (int i = 0; i < nvar; ++i) {
        t->SetBranchAddress(names[i].c_str(), &vals[i]);
        Double_t mn = t->GetMinimum(names[i].c_str());
        Double_t mx = t->GetMaximum(names[i].c_str());
        t->ResetBranchAddresses();
        RooRealVar* rv = new RooRealVar(names[i].c_str(), names[i].c_str(), mn, mx);
        varList.add(*rv);
    }

    for (int i = 0; i < nvar; ++i)
        t->SetBranchAddress(names[i].c_str(), &vals[i]);

    RooDataSet data("data", "input dataset", varList);

    Long64_t nEntries = std::min((Long64_t)maxEvents, t->GetEntries());
    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);
        RooArgSet row(varList);
        for (int j = 0; j < nvar; ++j) {
            RooRealVar* rv = dynamic_cast<RooRealVar*>(varList.at(j));
            if (!rv) continue;
            rv->setVal(vals[j]);
        }
        data.add(row);
    }

    std::cout << "Loaded " << data.numEntries() << " events for KDE\n";

    auto t0 = std::chrono::high_resolution_clock::now();

    RooNDKeysPdf pdf("pdfND", "smoothed PDF", varList, data, "ma", rho, 3);
    RooAbsData* gen = pdf.generate(varList, factor * data.numEntries(), true, false);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Build+Generate took "
        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
        << " ms\n";

    if (!gen || gen->numEntries() == 0) {
        std::cerr << "Error: Generated dataset is null or empty.\n";
        inF.Close();
        return;
    }

    TFile outF(outputFile, "RECREATE");
    TTree* tout = new TTree(treeName, "multiplied events");
    std::vector<Double_t> newvals(nvar);
    for (int i = 0; i < nvar; ++i) {
        tout->Branch(names[i].c_str(), &newvals[i]);
    }

    for (int i = 0; i < gen->numEntries(); ++i) {
        const RooArgSet* row = gen->get(i);
        if (!row) continue;
        for (int j = 0; j < nvar; ++j) {
            if (row->find(names[j].c_str())) {
                newvals[j] = row->getRealValue(names[j].c_str());
            }
        }
        tout->Fill();
    }

    std::cout << "Tree has " << tout->GetEntries() << " entries before writing.\n";

    tout->Write();
    outF.Close();
    inF.Close();

    std::cout << "Wrote " << tout->GetEntries() << " multiplied events to " << outputFile << "\n";
}
