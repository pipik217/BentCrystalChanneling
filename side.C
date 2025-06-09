void side() {

    // Create canvas and divide into 2 pads
    TCanvas *c1 = new TCanvas("c1", "Energy Distribution Comparison of Xi_{b} Particles", 1600,900);
    c1->SetTopMargin(0.05);
    c1->Divide(2,1);    
    // Two plots side by side    c1->Divide(2,1);

// Create top pad for the title
  // ===============================
    // === Plot 1: All Particles ====
    // ===============================
    TFile *file = TFile::Open("eventsMultiplied.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open eventsMultiplied.root!" << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: TTree 'tree' not found in eventsMultiplied.root!" << std::endl;
        return;
    }

    c1->cd(1);  // First pad
    gPad->SetGrid();
    gPad->SetLeftMargin(0.16);
    gPad->SetTopMargin(0.1);
    TH1F *hAll = new TH1F("hAll", "", 100, 0, 4500);
    hAll->SetXTitle("Energy [GeV]");
    hAll->SetYTitle("Relative Frequency");

    tree->Draw("E>>hAll");
      // Normalize
    if (hAll->Integral() > 0)
        hAll->Scale(1.0 / hAll->Integral());
    hAll->SetMinimum(0);
    hAll->SetMaximum(0.07);
    hAll->SetLineColor(kBlue + 1);
    hAll->SetLineWidth(2);
    hAll->SetFillColor(kBlue - 10);
    hAll->SetFillStyle(3004);  // Transparent fill
    hAll->Draw("HIST");
// Create right pad for Channeled particles


    // ===============================
    // === Plot 2: Channeled Particles ====
    // ===============================
    TFile *channeledFile = TFile::Open("channeled_particles_Germanium.root");
    if (!channeledFile || channeledFile->IsZombie()) {
        std::cerr << "Error: Cannot open channeled_particles.root!" << std::endl;
        return;
    }

    TTree *channeledTree = (TTree*)channeledFile->Get("channeled");
    if (!channeledTree) {
        std::cerr << "Error: TTree 'channeled' not found in channeled_particles.root!" << std::endl;
        return;
    }
    c1->cd(2);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.16);
    gPad->SetTopMargin(0.1);
    TH1F *hChan = new TH1F("hChan", "", 100, 0, 4500);
    hChan->SetXTitle("Energy [GeV]");
    hChan->SetYTitle("Relative Frequency");

    channeledTree->Draw("E>>hChan");
    // Normalize
    if (hChan->Integral() > 0)
        hChan->Scale(1.0 / hChan->Integral());
    hChan->SetMinimum(0);
    hChan->SetMaximum(0.07);
    hChan->SetLineColor(kRed + 1);
    hChan->SetLineWidth(2);
    hChan->SetFillColor(kRed - 10);
    hChan->SetFillStyle(3005);
    hChan->Draw("HIST");
// Ensure both pads have been drawn first
    c1->cd(1);
    gPad->Update();  // Force update of pad 1
    c1->cd(2);
    gPad->Update();  // Force update of pad 2
    c1->cd();
// Create a transparent pad just for the title
    TPad *titlePad = new TPad("titlePad", "titlePad", 0, 0.93, 1, 1.0);  // top strip
    titlePad->SetFillStyle(4000);  // Transparent
    titlePad->SetFrameFillStyle(4000);
    titlePad->Draw();
    titlePad->cd();
    TLatex *title = new TLatex();
    title->SetNDC();               // Use normalized device coordinates (0–1)
    title->SetTextSize(0.8);      // Adjust the size (try 0.05–0.07)
    title->SetTextAlign(22);       // Centered (horizontally and vertically)
    title->DrawLatex(0.5, 0.5, "Energy Distribution Comparison of All vs. Channeled #Xi_{b} Baryons");
   
// Create a custom title
    // Save canvas as PNG
    c1->SaveAs("energy_distributions_side_by_side.png");

    // Cleanup
    file->Close();
    channeledFile->Close();
}
