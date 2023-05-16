float subtractTower(float towerE, float towerUE, float towerPhi, float v2, float psi2, int doFlow = 1){
  
  if(doFlow){
    towerUE += towerUE * (2*v2*TMath::Cos(2*(towerPhi - psi2)));
  }
  
  return towerE-towerUE;

}

int getEtaBin(float eta){
  for(int i = 0; i < 24; i++){
    if(eta < -1.1+(i+1)*2.2/24) return i;
  }
  return -1;
}

int getPhiBin(float phi){
  for(int i = 0; i < 64; i++){
    if(phi < (i+1)*2*TMath::Pi()/64) return i;
  }
  return -1;
}

vector<vector<int>> chimney_sectors = {{0,13},{0,14},{0,15},{0,16},{0,17},{0,18},{1,13},{1,14},{1,15},{1,16},{1,17},{1,18},{2,13},{2,14},{2,15},{2,16},{2,17},{2,18},{3,13},{3,14},{3,15},{3,16},{3,17},{3,18}, {0,8},{1,8},{22,8},{23,8},{0,22},{1,22},{22,22},{23,22}}; //also includes towers affected by TCP supports.

void getFluctuations(const int nsizes, const string sizes[], string in = "outputformax.root", string out ="fluctuations", int subtraction = 0, bool allow_partial_exclusions = false, int nEvents = 0, bool Et = true, vector<vector<int>> exclude = chimney_sectors){
    //"nsizes" is the number of window sizes to be analysed.
    //"sizes" is an array of window sizes, expressed as strings of
    //the form "<width in eta>x<width in phi>".
    //"in" is the input file name, which contains event data.
    //"out" is the output file name prefix.
    //"subtraction" determines subtraction. Subtraction: 0= no subtraction, 1=subtraction without flow, 2= subtraction with flow.
    //"exclude" is a list of towers (as {etabin, phibin} pairs) to exclude from the analysis. Follows the
    //same bin indexing as the getEtaBin and getPhiBin functions.
    //"allow_partial_exclusions" determines whether windows with some excluded
    //towers and some non-excluded towers have their energy scaled to match the
    //average energy of the window's non-excluded towers, or whether these
    //windows are discarded entirely.
    //"nEvents" gives the number of events to be included in the analysis. A
    //value of 0 allows all events in the input file to be analyzed.
    //"Et" determines whether the transverse energy or energy is analysed.
  int doFlow = (subtraction == 2) ? 1 : 0;

  //Define geometric parameters of towers
  float etaWidth = 2.2 / 24;
  float phiWidth = 2*TMath::Pi() / 64;
  
  //read in file

  TChain* chain = new TChain("tree");
  chain->Add(in.c_str());
  if (chain == 0) return;
  Long64_t nentries = (nEvents > 0) ? nEvents+1 : chain->GetEntries();
  if(nentries == 0){
    std::cout<<"No entries"<<std::endl;
    return;
  }

  int m_event;
  std::vector<float> *m_towerCEMC_e = 0;
  std::vector<float> *m_towerCEMC_eta = 0;
  std::vector<float> *m_towerCEMC_phi = 0;
  std::vector<int> *m_towerCEMC_etaBin = 0;
  std::vector<int> *m_towerCEMC_phiBin = 0;

  std::vector<float> *m_towerIHCAL_e = 0;
  std::vector<float> *m_towerIHCAL_eta = 0;
  std::vector<float> *m_towerIHCAL_phi = 0;
  std::vector<int> *m_towerIHCAL_etaBin = 0;
  std::vector<int> *m_towerIHCAL_phiBin = 0;

  std::vector<float> *m_towerOHCAL_e = 0;
  std::vector<float> *m_towerOHCAL_eta = 0;
  std::vector<float> *m_towerOHCAL_phi = 0;
  std::vector<int> *m_towerOHCAL_etaBin = 0;
  std::vector<int> *m_towerOHCAL_phiBin = 0;

  std::vector<float> *m_UE0 = 0;
  std::vector<float> *m_UE1 = 0;
  std::vector<float> *m_UE2 = 0;

  float totalCEMC;
  float totalIHCAL;
  float totalOHCAL;

  float m_v2;
  double m_b;
  double m_cent;
  float m_psi2;

  vector<float> *eta = 0;
  vector<float> *phi = 0;
  vector<float> *pt = 0;
  vector<float> *e = 0;
  vector<float> *truthEta = 0;
  vector<float> *truthPhi = 0;
  vector<float> *truthPt = 0;
  vector<float> *truthE = 0;

  chain->SetBranchAddress("m_event", &m_event);
  chain->SetBranchAddress("m_towerCEMC_e", &m_towerCEMC_e);
  chain->SetBranchAddress("m_towerCEMC_eta", &m_towerCEMC_eta);
  chain->SetBranchAddress("m_towerCEMC_phi", &m_towerCEMC_phi);
  chain->SetBranchAddress("m_towerCEMC_etaBin", &m_towerCEMC_etaBin);
  chain->SetBranchAddress("m_towerCEMC_phiBin", &m_towerCEMC_phiBin);

  chain->SetBranchAddress("m_towerIHCAL_e", &m_towerIHCAL_e);
  chain->SetBranchAddress("m_towerIHCAL_eta", &m_towerIHCAL_eta);
  chain->SetBranchAddress("m_towerIHCAL_phi", &m_towerIHCAL_phi);
  chain->SetBranchAddress("m_towerIHCAL_etaBin", &m_towerIHCAL_etaBin);
  chain->SetBranchAddress("m_towerIHCAL_phiBin", &m_towerIHCAL_phiBin);

  chain->SetBranchAddress("m_towerOHCAL_e", &m_towerOHCAL_e);
  chain->SetBranchAddress("m_towerOHCAL_eta", &m_towerOHCAL_eta);
  chain->SetBranchAddress("m_towerOHCAL_phi", &m_towerOHCAL_phi);
  chain->SetBranchAddress("m_towerOHCAL_etaBin", &m_towerOHCAL_etaBin);
  chain->SetBranchAddress("m_towerOHCAL_phiBin", &m_towerOHCAL_phiBin);

  chain->SetBranchAddress("m_UE0", &m_UE0);
  chain->SetBranchAddress("m_UE1", &m_UE1);
  chain->SetBranchAddress("m_UE2", &m_UE2);

  chain->SetBranchAddress("m_totalCEMC",&totalCEMC);
  chain->SetBranchAddress("m_totalIHCAL",&totalIHCAL);
  chain->SetBranchAddress("m_totalOHCAL",&totalOHCAL);

  chain->SetBranchAddress("m_v2", &m_v2);
  chain->SetBranchAddress("m_b", &m_b);
  chain->SetBranchAddress("m_cent", &m_cent);
  chain->SetBranchAddress("m_psi2", &m_psi2);

  chain->SetBranchAddress("eta",&eta);
  chain->SetBranchAddress("phi",&phi);
  chain->SetBranchAddress("pt",&pt);
  chain->SetBranchAddress("e",&e);
  chain->SetBranchAddress("truthEta",&truthEta);
  chain->SetBranchAddress("truthPhi",&truthPhi);
  chain->SetBranchAddress("truthPt",&truthPt);
  chain->SetBranchAddress("truthE",&truthE);


  //create histograms
  TH2D *h_AvgE[nsizes][3];//average energy for a given window size: n window sizes, 3 subdetectors
  TH2D *h_STD[nsizes][3];//STD of energy for a given window size: n window sizes, 3 subdetectors  
  TH2D *h_AvgEtot[nsizes];
  TH2D *h_STDtot[nsizes];
  TH3D *h_etaPhiEt[nsizes][3];

  TH2D *h_AvgE_TE[nsizes][3];
  TH2D *h_AvgE_IP[nsizes][3];
  TH2D *h_AvgE_Cent[nsizes][3];
  TH2D *h_STD_TE[nsizes][3];
  TH2D *h_STD_IP[nsizes][3];
  TH2D *h_STD_Cent[nsizes][3];

  TH2D *h_AvgEtot_TE[nsizes];
  TH2D *h_STDtot_TE[nsizes];
  TH2D *h_AvgEtot_IP[nsizes];
  TH2D *h_STDtot_IP[nsizes];
  TH2D *h_AvgEtot_Cent[nsizes];
  TH2D *h_STDtot_Cent[nsizes];

  //TH2D *h_windowenergy[nsizes]; //to examine where the bloody energy goes
  
  TH2D *h_windowEtot[nsizes][3];
  TH2I *h_windowExclusions[nsizes][3]; //Tracks # of bin exclusions in each window in each event.
  double meantot[nsizes][3];

  //Count # of excluded windows for each window size (in each event)
  int ExcludedWindows[nsizes];

  double meanEtsub[nsizes][3];
  double meanEt[nsizes];

  string energy_str = (Et) ? "Et" : "Energy";
  int Deta[nsizes];
  int Dphi[nsizes];
  int numineta[nsizes];
  int numinphi[nsizes];
  int numbins[nsizes];
  int *windowsize = new int[nsizes];

  float *AvgE_Per_Tower = new float[4] {1.5, 0.15, 0.2, 2};
  float *STD_Per_Tower = new float[4] {0.7, 0.12, 0.5, 0.8};

  string detectors[] = {"EMCal","IHCal","OHCal"};
  for(int i = 0; i < nsizes; i++){
      //preliminary work for assigning window boundaries.
      int size = 0; string Detastr = ""; string Dphistr = "";
      do{
        Detastr = Form("%s%c",Detastr.c_str(),sizes[i][size]);
        size++;
      } while (sizes[i][size] != 'x');
      size++;
      do{
        Dphistr = Form("%s%c",Dphistr.c_str(),sizes[i][size]);
        size++;
      } while (sizes[i][size] != '\0');
      Deta[i] = stoi(Detastr);
      Dphi[i] = stoi(Dphistr);
      numineta[i] = 24 / Deta[i];
      numinphi[i] = 64 / Dphi[i];
      numbins[i] = numineta[i]*numinphi[i];
      windowsize[i] = Deta[i]*Dphi[i];

      //Histogram Bounds: TODO make ure this system captures all data adequately
      //TODO: IHCal, OHCal boundaries much too large on 7x7 windows (avg & std).
        //Avg:
          //EMCal: 1.5GeV per tower
          //IHCal: 0.15GeV pwer tower
          //OHCal: 0.2GeV per tower
          //Total: 2GeV per tower
        //STD: (Using windowsize^0.7 to scale up to larger windows)
          //EMCal: 0.7GeV per tower
          //IHCal: 0.12GeV per tower
          //OHCal: 0.5 GeV per tower
          //Total: 0.8GeV per tower
    for(int j = 0; j < 3; j++){
      //TODO: label bins  of centrality histograms with the centrality ranges.
      h_AvgE_Cent[i][j]  = new TH2D(Form("h_AvgE%i_layer%i_Cent",i,j),Form("Mean Window %s vs Centrality. Size %s %s",energy_str.c_str(), sizes[i].c_str(),detectors[j].c_str()),5,0,100,8000,0,AvgE_Per_Tower[j]*windowsize[i]);
      h_STD_Cent[i][j]  = new TH2D(Form("h_STD%i_layer%i_Cent",i,j),Form("Standard Deviation Window %s vs Centrality. Size %s %s",energy_str.c_str(), sizes[i].c_str(),detectors[j].c_str()),5,0,100,1400,0,STD_Per_Tower[j]*TMath::Power(windowsize[i],0.7));
      h_AvgE_IP[i][j]  = new TH2D(Form("h_AvgE%i_layer%i_IP",i,j),Form("Mean Window %s vs Impact Parameter. Size %s %s",energy_str.c_str(), sizes[i].c_str(),detectors[j].c_str()),200,0,12,8000,0,AvgE_Per_Tower[j]*windowsize[i]);
      h_STD_IP[i][j]  = new TH2D(Form("h_STD%i_layer%i_IP",i,j),Form("Standard Deviation Window %s vs Impact Parameter. Size %s %s",energy_str.c_str(), sizes[i].c_str(),detectors[j].c_str()),200,0,12,1400,0,STD_Per_Tower[j]*TMath::Power(windowsize[i],0.7));
      h_AvgE_TE[i][j]  = new TH2D(Form("h_AvgE%i_layer%i_TE",i,j),Form("Mean Window %s vs Total %s. Size %s %s",energy_str.c_str(), energy_str.c_str(), sizes[i].c_str(),detectors[j].c_str()),200,0,2500,8000,0,AvgE_Per_Tower[j]*windowsize[i]);
      h_STD_TE[i][j]  = new TH2D(Form("h_STD%i_layer%i_TE",i,j),Form("Standard Deviation Window %s vs Total %s. Size %s %s",energy_str.c_str(), energy_str.c_str(), sizes[i].c_str(),detectors[j].c_str()),200,0,2500,1400,0,STD_Per_Tower[j]*TMath::Power(windowsize[i],0.7));
      h_windowEtot[i][j] = new TH2D(Form("h_window%sEtot_sub%i",sizes[i].c_str(),j),Form("%s Window %s Subdetector %i",energy_str.c_str(), sizes[i].c_str(),j),numineta[i], -1.1+(1.1*(24%Deta[i])/24), 1.1-(1.1*(24%Deta[i])/24), numinphi[i], 0, (2-1.0*(64%Dphi[i])/32)*TMath::Pi()); 
      h_windowExclusions[i][j] = new TH2I(Form("h_window%sExclusions%i",sizes[i].c_str(),j),Form("%s number of Excluded Towers in Windows Subdetector %i",sizes[i].c_str(),j),numineta[i], -1.1+(1.1*(24%Deta[i])/24), 1.1-(1.1*(24%Deta[i])/24), numinphi[i], 0, (2-1.0*(64%Dphi[i])/32)*TMath::Pi()); 
    }

        //h_windowenergy[i] = new TH2D(Form("h_window%stotalenergy (sum over %i events)",sizes[i].c_str(),nEvents),Form("Size %s Window Energy",sizes[i].c_str()),numineta[i], -1.1+(1.1*(24%Deta[i])/24), 1.1-(1.1*(24%Deta[i])/24), numinphi[i], 0, (2-1.0*(64%Dphi[i])/32)*TMath::Pi()); //to examine where the bloody energy goes.

    h_AvgEtot_Cent[i] = new TH2D(Form("h_AvgEtot%i_Cent",i),Form("Mean Winodw %s vs Centrality. Window Size %s",energy_str.c_str(), sizes[i].c_str()),5,0,100,200,0,AvgE_Per_Tower[3]*windowsize[i]);
    h_STDtot_Cent[i] = new TH2D(Form("h_STDtot%i_Cent",i),Form("Standard Deviation Winodw %s vs Centrality. Window Size %s",energy_str.c_str(), sizes[i].c_str()),5,0,100,1400,0,STD_Per_Tower[3]*TMath::Power(windowsize[i],0.7));
    h_AvgEtot_IP[i] = new TH2D(Form("h_AvgEtot%i_IP",i),Form("Mean Winodw %s vs Impact Parameter. Window Size %s",energy_str.c_str(), sizes[i].c_str()),200,0,12,200,0,AvgE_Per_Tower[3]*windowsize[i]);
    h_STDtot_IP[i] = new TH2D(Form("h_STDtot%i_IP",i),Form("Standard Deviation Winodw %s vs Impact Parameter. Window Size %s",energy_str.c_str(), sizes[i].c_str()),200,0,12,1400,0,STD_Per_Tower[3]*TMath::Power(windowsize[i],0.7));
    h_AvgEtot_TE[i] = new TH2D(Form("h_AvgEtot%i_TE",i),Form("Mean Winodw %s vs Total %s. Window Size %s",energy_str.c_str(), energy_str.c_str(), sizes[i].c_str()),200,0,2500,200,0,AvgE_Per_Tower[3]*windowsize[i]);
    h_STDtot_TE[i] = new TH2D(Form("h_STDtot%i_TE",i),Form("Standard Deviation Winodw %s vs Total %s. Window Size %s",energy_str.c_str(), energy_str.c_str(), sizes[i].c_str()),200,0,2500,1400,0,STD_Per_Tower[3]*TMath::Power(windowsize[i],0.7));
  }
  delete[] windowsize;
  delete[] AvgE_Per_Tower;
  delete[] STD_Per_Tower;
  //loop through all the entries in our file and fill the histograms
  int nancount = 0;

  std::cout<<"running "<<nentries<<" events"<<std::endl;
  for(Long64_t jentry=0; jentry<nentries; jentry++){ //TODO: See if we can handle exclusions before the jentry loop
    if(jentry%1000 == 0) std::cout<<"Entry: "<<jentry<<std::endl;
    
    chain->GetEntry(jentry);

    //get the total energy in all three subdetectors
    float totalEnergy = totalCEMC+totalIHCAL+totalOHCAL;
    float towerEt[64][24] = {}; //there are 64 towers in phi x 24 towers in eta
    
    //initialize a few things.
    for (int i = 0; i < nsizes; i++){
        ExcludedWindows[i] = 0;
        for (int j=0; j<3; j++){
            meantot[i][j] = 0;
        }
    }
    //loop through each subdetector: 0: EMCal, 1: inner HCal, 2: outer HCal
    for(int isys = 0; isys < 3; isys++){
        
      //get the number of towers in an event for the given subdetector 
      int nTower = 0;
      if(isys == 0) nTower = m_towerCEMC_e->size();
      else if(isys == 1) nTower = m_towerIHCAL_e->size();
      else if(isys == 2) nTower = m_towerOHCAL_e->size();

      //Mark all the excluded windows with an energy of -1 in towerEt
      for (vector<int> excluded : exclude){
        towerEt[excluded.at(1)][excluded.at(0)] = -1;
        for (int n = 0; n<nsizes; n++){
          h_windowExclusions[n][isys]->Fill( -1.1 + etaWidth*(0.5+excluded.at(0)), phiWidth*(0.5+excluded.at(1)) ); //Mark the excluded tower at its center
        }
      }

      //get the energy for each tower (and subtract the underlying event if using subtraction)
      for(int i = 0; i < nTower; i++){
	float eta;
	float phi;
	int etabin;
	int phibin;
    float energy;
    if(isys == 0){
      eta = m_towerCEMC_eta->at(i);
      phi = m_towerCEMC_phi->at(i);
      energy = (subtraction == 0) ? m_towerCEMC_e->at(i) : subtractTower(m_towerCEMC_e->at(i), m_UE0->at(getEtaBin(eta)), phi, m_v2, m_psi2, doFlow);
    }
    else if(isys== 1){
          eta = m_towerIHCAL_eta->at(i);
          phi = m_towerIHCAL_phi->at(i);
      energy = (subtraction == 0) ? m_towerIHCAL_e->at(i) : subtractTower(m_towerIHCAL_e->at(i), m_UE1->at(getEtaBin(eta)), phi, m_v2, m_psi2, doFlow);
    }
    else if(isys== 2){
          eta = m_towerOHCAL_eta->at(i);
          phi = m_towerOHCAL_phi->at(i);
      energy = (subtraction == 0) ? m_towerOHCAL_e->at(i) : subtractTower(m_towerOHCAL_e->at(i), m_UE2->at(getEtaBin(eta)), phi, m_v2, m_psi2, doFlow);
    }

    etabin = getEtaBin(eta);
    phibin = getPhiBin(phi);
    bool BinNotExcluded = true;

    if (towerEt[phibin][etabin] != -1) towerEt[phibin][etabin] = (Et) ? energy/cosh(eta) : energy;
    else BinNotExcluded = false;

    for (int n = 0; n<nsizes; n++){
      if (BinNotExcluded) {
        h_windowEtot[n][isys]->Fill(eta, phi, towerEt[phibin][etabin]);
        //h_windowenergy[n]->Fill(eta, phi, towerEt[phibin][etabin]); //to examine where the bloody energy goes.
      }
    }
      }

    //calculate mean and STD of window Et in each subdetector
    double win[64][24];
    for (int n=0; n<nsizes; n++){
        //Compute STD and mean
        double STDEt = 0;
        for (int i=0; i<numineta[n]; i++){
            for (int j=0; j<numinphi[n]; j++){
                //First, handle exclusions.
                int TowerExclusions = h_windowExclusions[n][isys]->GetBinContent(i+1, j+1);
                if (TowerExclusions > 0){
                    //If window is to be excluded, set its energy to be -1.
                    if (allow_partial_exclusions){
                        if (TowerExclusions == windowsize[n]) h_windowEtot[n][isys]->SetBinContent(i+1,j+1,-1); //Mark entire window for exclusion
                        else h_windowEtot[n][isys]->SetBinContent(i+1, j+1, ((double) windowsize[n]) / (windowsize[n] - TowerExclusions) * h_windowEtot[n][isys]->GetBinContent(i+1,j+1)); //Scale window energy based on # of exclusions
                    }
                    else{
                        h_windowEtot[n][isys]->SetBinContent(i+1, j+1, -1); //Mark entire window for exclusions
                    }
                }

                //Start the calculations
                win[i][j] = h_windowEtot[n][isys]->GetBinContent(i+1,j+1);
                if (win[i][j] == -1){
                  if (isys==0) ExcludedWindows[n]++; //Count exlcuded windows
                }
                else meantot[n][isys] += win[i][j];
            }
        }
        meanEtsub[n][isys] = meantot[n][isys] / (numbins[n] - ExcludedWindows[n]); //Average Et in each window
        for (int i=0; i<numineta[n]; i++){
            for (int j=0; j<numinphi[n]; j++){
                if (win[i][j] != -1) STDEt += TMath::Sq(meanEtsub[n][isys] - win[i][j]);
            }
        }
        STDEt = TMath::Sqrt(STDEt / (numbins[n] - ExcludedWindows[n]));

        h_AvgE_Cent[n][isys]->Fill(m_cent, meanEtsub[n][isys]);
        h_STD_Cent[n][isys]->Fill(m_cent, STDEt);
        h_AvgE_IP[n][isys]->Fill(m_b, meanEtsub[n][isys]);
        h_STD_IP[n][isys]->Fill(m_b, STDEt);
        h_AvgE_TE[n][isys]->Fill(totalEnergy, meanEtsub[n][isys]);
        h_STD_TE[n][isys]->Fill(totalEnergy, STDEt);

        meanEt[n] += meanEtsub[n][isys];
    }
    }
    //Calculate mean and STD of total Et windows
    //TODO (later) Try not to calculate this all again! We did it above already.
    for (int n=0; n<nsizes; n++){
        double rrmss=0;
        double STDEt = 0;
        double wintot[64][24];
        for (int i=0; i<numineta[n]; i++){
            for (int j=0; j<numinphi[n]; j++){
                wintot[i][j] = 0;
                for (int isys=0; isys<3; isys++){
                    double win2 = h_windowEtot[n][isys]->GetBinContent(i+1,j+1);
                    if (win2 != -1) wintot[i][j] += win2;
                    else{
                        wintot[i][j] = -1;
                        break;
                    }
                }
                if (wintot[i][j] != -1) rrmss += wintot[i][j];
            }
        }
        rrmss /= (numbins[n] - ExcludedWindows[n]);
        for (int i=0; i<numineta[n]; i++){
            for (int j=0; j<numinphi[n]; j++){
                if (wintot[i][j] != -1) STDEt += TMath::Sq(rrmss - wintot[i][j]);
            }
        }
        
        STDEt = TMath::Sqrt(STDEt/(numbins[n] - ExcludedWindows[n]));

        h_AvgEtot_Cent[n]->Fill(m_cent, meanEt[n]);
        h_STDtot_Cent[n]->Fill(m_cent, STDEt);
        h_AvgEtot_IP[n]->Fill(m_b, meanEt[n]);
        h_STDtot_IP[n]->Fill(m_b, STDEt);
        h_AvgEtot_TE[n]->Fill(totalEnergy, meanEt[n]);
        h_STDtot_TE[n]->Fill(totalEnergy, STDEt );
    
        meanEt[n] = 0;
        for(int k = 0; k < 3; k++){
              if (jentry != nentries-1){ //Allow writing of the final histogram of window energy and exclusions to check that the code is doing what it should be doing.
              h_windowExclusions[n][k]->Reset();
              h_windowEtot[n][k]->Reset();
          }
        }
    }
  }

  string s = "";
  if(subtraction == 0) s = "_unsubtracted";
  else if(subtraction == 1) s = "_subtracted_noflow";
  else if(subtraction == 2) s = "_subtracted_withflow";
  TFile *f = new TFile(Form("%s%s.root",out.c_str(),s.c_str()),"RECREATE");
 
  for(int i = 0; i < nsizes; i++){
    h_AvgEtot_Cent[i]->Write();
    h_STDtot_Cent[i]->Write();
    h_AvgEtot_IP[i]->Write();
    h_STDtot_IP[i]->Write();
    h_AvgEtot_TE[i]->Write();
    h_STDtot_TE[i]->Write();
    delete h_AvgEtot_Cent[i];
    delete h_STDtot_Cent[i];
    delete h_AvgEtot_IP[i];
    delete h_STDtot_IP[i];
    delete h_AvgEtot_TE[i];
    delete h_STDtot_TE[i];
        //h_windowenergy[i]->Scale(1.0 / nentries);
        //h_windowenergy[i]->Write(); //to examine where the bloody energy goes.
    for(int j = 0; j < 3; j++){
      h_STD_Cent[i][j]->Write();
      h_AvgE_Cent[i][j]->Write();
      h_STD_IP[i][j]->Write();
      h_AvgE_IP[i][j]->Write();
      h_STD_TE[i][j]->Write();
      h_AvgE_TE[i][j]->Write();
      h_windowExclusions[i][j]->Write(); //Debugging
      h_windowEtot[i][j]->Write(); //Debugging
      delete h_STD_Cent[i][j];
      delete h_AvgE_Cent[i][j];
      delete h_STD_IP[i][j];
      delete h_AvgE_IP[i][j];
      delete h_STD_TE[i][j];
      delete h_AvgE_TE[i][j];
      delete h_windowEtot[i][j];
      delete h_windowExclusions[i][j];
      /*for(int k=0; k<nentries; k++){
      h_windowEtot[i][j][k]->Write();
      }*/
    }
  }

  f->Close();
}
