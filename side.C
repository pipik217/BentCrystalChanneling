void log() {
    // Create a single canvas
    TCanvas *c1 = new TCanvas("c1", "Energy Distribution of Xi_{b} Particles", 1600, 900);
    c1->SetTopMargin(0.1);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    c1->SetGrid();
    c1->SetLogy();  // Set logarithmic scale for y-axis
    
    // Open first file for all Xib baryons
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

    // Find maximum energy value to set x-axis range
    Double_t maxE_all = tree->GetMaximum("E");
    Double_t maxE_chan = 0;
    
    // Open second file for channeled particles to compare max energy
    TFile *channeledFile = TFile::Open("channeled_particles_Germanium.root");
    if (channeledFile && !channeledFile->IsZombie()) {
        TTree *channeledTree = (TTree*)channeledFile->Get("channeled");
        if (channeledTree) {
            maxE_chan = channeledTree->GetMaximum("E");
        }
        channeledFile->Close();
    }
    
    // Use the maximum of both maximums for x-axis range
    Double_t xmax = TMath::Max(maxE_all, maxE_chan) * 1.05; // Add 5% margin

    // Create and draw first histogram (all particles)
    TH1F *hAll = new TH1F("hAll", "Energy Distribution of #Xi_{b} Baryons;Energy [GeV];Counts", 100, 0, xmax);
    tree->Draw("E>>hAll");
    hAll->SetLineColor(kBlue + 1);
    hAll->SetLineWidth(2);
    hAll->SetFillColor(kBlue - 10);
    hAll->SetFillStyle(3004); 
    hAll->Draw("HIST");

    // Reopen channeled file for actual plotting
    channeledFile = TFile::Open("channeled_particles_Germanium.root");
    if (!channeledFile || channeledFile->IsZombie()) {
        std::cerr << "Error: Cannot open channeled_particles.root!" << std::endl;
        return;
    }

    TTree *channeledTree = (TTree*)channeledFile->Get("channeled");
    if (!channeledTree) {
        std::cerr << "Error: TTree 'channeled' not found in channeled_particles.root!" << std::endl;
        return;
    }

    // Create and draw second histogram (channeled particles)
    TH1F *hChan = new TH1F("hChan", "", 100, 0, xmax);
    channeledTree->Draw("E>>hChan");
    hChan->SetLineColor(kRed + 1);
    hChan->SetLineWidth(2);
    hChan->SetFillColor(kRed - 10);
    hChan->SetFillStyle(3005);
    hChan->Draw("HIST SAME");

    // Add legend
    TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hAll, "All #Xi_{b}", "f");
    leg->AddEntry(hChan, "Channeled #Xi_{b}", "f");
    leg->Draw();

    // Adjust y-axis range to look good in log scale
    Double_t ymin = 0.1;  // Small non-zero value for log scale
    Double_t ymax = TMath::Max(hAll->GetMaximum(), hChan->GetMaximum()) * 10;
    hAll->SetMinimum(ymin);
    hAll->SetMaximum(ymax);

    // Update canvas
    c1->Update();
    
    // Save canvas as PNG
    c1->SaveAs("log.png");

    // Cleanup
    file->Close();
    channeledFile->Close();
}
