void jer(TH2D* h){

    //TODO: figure out how to take a file as an input and extract histograms from it.
    Int_t num_of_bins = h->GetXaxis()->GetNbins();
    TH1D *histos[num_of_bins];
    Double_t* xbar = new Double_t [num_of_bins];
    Double_t* xbarerr = new Double_t [num_of_bins];
    Double_t* truthpt = new Double_t [num_of_bins]; //in empty x-bins, graph a point (0,0) with no error.
    Double_t* pterr = new Double_t [num_of_bins];
    Double_t* sigma = new Double_t [num_of_bins];
    Double_t* sigmaerr = new Double_t [num_of_bins];
    TF1::InitStandardFunctions();
    //auto gaus = gROOT->GetFunction("gaus");
    TF1 *f1 = new TF1("f1", "gaus", 0);
    TFitResultPtr r;
    TFile *f = new TFile("JERandJES.root","RECREATE");
    TCanvas *canvas = new TCanvas("canvas");

    //Store the projection of each xbin from *h and fit to a gaussian.
    for (Int_t xbin = 1; xbin <= num_of_bins; xbin++){
        Double_t bin_energy = h->GetXaxis()->GetBinLowEdge(xbin);
        Int_t i = bin_energy; //There may be a more general way to do this.
        TH1D* h_proj = h->ProjectionY(Form("h_Jet_Reconstruction_TruthPt_%i",i),xbin,xbin);
        h_proj->SetTitle(Form("%i_Truth_Pt",i));
        h_proj->GetXaxis()->SetTitle("Reconstructed Pt / Truth Pt");
        h_proj->GetYaxis()->SetTitle("counts");
        h_proj->GetYaxis()->SetTitleSize(0.045);
        h_proj->GetXaxis()->SetTitleSize(0.045);
        histos[xbin-1] = h_proj;
        if (h_proj->GetEntries() != 0){ 
            try {
                Double_t std = h_proj->GetStdDev();
                Double_t mean = h_proj->GetMean();
                r = h_proj->Fit(f1,"SQR", "", mean-std, mean+std);
                r = h_proj->Fit(f1,"SQR", "", r->Parameter(1) - r->Parameter(2), r->Parameter(1) + r->Parameter(2));
                xbar[xbin] = r->Parameter(1);
                xbarerr[xbin] = r->ParError(1);
                sigma[xbin] = r->Parameter(2);
                sigmaerr[xbin] = r->ParError(2);
                truthpt[xbin] = bin_energy;
                pterr[xbin] = 0;
                h_proj->Write();
            } catch (...) {
                xbar[xbin] = 0;
                xbarerr[xbin] = 0;
                sigma[xbin] = 0;
                sigmaerr[xbin] = 0;
                truthpt[xbin] = bin_energy;
                pterr[xbin] = 0;
                h_proj->Write();
            }
            h_proj->Draw("col");
            if (xbin == 1){
                canvas->Print("Projected Histograms.pdf(");
            } else if (xbin == num_of_bins){
                canvas->Print("Projected Histograms.pdf)");
            } else {
                canvas->Print("Projected Histograms.pdf");
            }
        } else {
            xbar[xbin] = 0;
            xbarerr[xbin] = 0;
            sigma[xbin] = 0;
            sigmaerr[xbin] = 0;
            truthpt[xbin] = 0;
            pterr[xbin] = 0;
        }
    }
    canvas->Close();

    //remove 0's from the data for a clearer graph
    /* Removed -- causes segmantation violation
    Int_t num_zeros = 0;
    for (Int_t i; i <= num_of_bins; i++){
        if (xbar[i] == 0){
            num_zeros++;
        }
    }
    Int_t j = 0;
    Double_t *xbar2 = new Double_t [num_zeros];
    Double_t *truthpt2 = new Double_t [num_zeros];
    Double_t *pterr2 = new Double_t [num_zeros];
    Double_t *xbarerr2 = new Double_t [num_zeros];
    for (Int_t i; i <= num_of_bins; i++){
        if (xbar[i] != 0){
            xbar2[j] = xbar[i];
            truthpt2[j] = truthpt[i];
            pterr2[j] = pterr[i];
            xbarerr2[j] = xbarerr[i];
            j++;
        }
    } */

    //Graph mean (JES) and sigma (JER) vs truth pt
    TGraphErrors *JES = new TGraphErrors(num_of_bins, truthpt, xbar, pterr, xbarerr);
    JES->SetTitle("Jet Energy Scale vs Pt;Truth Pt (GeV);JES");
    JES->SetName("JES");
    JES->Write();
    TGraphErrors *JER = new TGraphErrors(num_of_bins, truthpt, sigma, pterr, sigmaerr);
    JER->SetTitle("Jet Energy Resolution vs Pt;Truth Pt (GeV);JER");
    JER->SetName("JER");
    JER->Write();
    f->Close();
}
