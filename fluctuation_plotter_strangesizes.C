void fluctuation_plotter_strangesizes(int nsizes, string sizes[], const char * datafile = "testsizes_unsubtracted.root", string destination = "fluctuation_hists_strangesizes"){
    //Read file
    TFile *f = TFile::Open(datafile);

    string size_no[nsizes];
    for (int i=0; i<nsizes; i++){
        size_no[i] = to_string(i);
    }
    string layers[3] = {"EMCal","IHCal","OHCal"};
    string layer_no[3] = {"0","1","2"};
    string avg_type[2] = {"AvgE","STD"};
    string collision_parameter[3] = {"Cent", "IP","TE"};

    TCanvas *canvas = new TCanvas("canvas","canvas",800,600);
    for (int window = 0; window < nsizes; window++){
        for (int average = 0; average < 2; average++){
            for (int param = 0; param < 3; param++){
                string name = Form("h_%stot%s_%s", avg_type[average].c_str(), size_no[window].c_str(), collision_parameter[param].c_str());
                TH2D *h = (TH2D*) f->Get(name.c_str());
                h->Draw("col");
                h->GetXaxis()->SetTitle(Form("%s (GeV)",collision_parameter[param].c_str()));
                h->GetXaxis()->SetTitleSize(0.045);
                h->GetYaxis()->SetTitle(Form("%s (GeV)",avg_type[average].c_str()));
                h->GetYaxis()->SetTitleSize(0.045);
                canvas->Print(Form("%s/h_%stot%s_%s.pdf",destination.c_str(), avg_type[average].c_str(), sizes[window].c_str(), collision_parameter[param].c_str()));
                for (int calo = 0; calo < 3; calo++){
                    string name = Form("h_%s%s_layer%s_%s", avg_type[average].c_str(), size_no[window].c_str(), layer_no[calo].c_str(), collision_parameter[param].c_str());
                    TH2D *h = (TH2D*) f->Get(name.c_str());
                    h->Draw("col");
                    h->GetXaxis()->SetTitle(Form("%s (GeV)",collision_parameter[param].c_str()));
                    h->GetXaxis()->SetTitleSize(0.045);
                    h->GetYaxis()->SetTitle(Form("%s (GeV)",avg_type[average].c_str()));
                    h->GetYaxis()->SetTitleSize(0.045);
                    canvas->Print(Form("%s/h_%s%s_layer%s_%s.pdf",destination.c_str(), avg_type[average].c_str(), sizes[window].c_str(), layers[calo].c_str(), collision_parameter[param].c_str()));
                }
            }
        }
    }
    
    canvas->Close();
    f->Close();
}
