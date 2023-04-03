

////// This is the main function ///////
void timing_hit_maps(int run = 2587){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TString data_dir = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/";
  TString plots_dir = "../plots/";

  TChain *T = new TChain("T");
  T->Add(Form(data_dir + "gen_replayed_%i_1000k_events.root",run)); 
  

  int maxhits = 5000;

  int nhits;
  double module[maxhits];
  double isamp[maxhits];
  double xhit[maxhits];
  double yhit[maxhits];


  T->SetBranchAddress("Ndata.bb.gem.hit.isampmaxUclust",&nhits);
  T->SetBranchAddress("bb.gem.hit.module",module);
  T->SetBranchAddress("bb.gem.hit.isampmaxUclust",isamp);
  T->SetBranchAddress("bb.gem.hit.xlocal",xhit);
  T->SetBranchAddress("bb.gem.hit.ylocal",yhit);
  
  const int nbinsx = 50;
  const int nbinsy = 50;
  const int nmodules = 8;

  TH2F *h[nmodules];

  for(int imod=0; imod<nmodules;imod++){
    if(imod < 4) h[imod] = new TH2F(Form("h%i",imod),Form("Module %i Average Peak Time Sample;y pos (m);x pos (m)",imod),nbinsx,-0.25,0.25,nbinsy,-0.8,0.8);
    else
      h[imod] = new TH2F(Form("h%i",imod),Form("Module %i Average Peak Time Sample;y pos (m);x pos (m)",imod),nbinsx,-0.35,0.35,nbinsy,-0.3,0.3);

  }

  int nentries[nmodules][nbinsx*nbinsy] = {0};
  int isamp_tot[nmodules][nbinsx*nbinsy] = {0};
  int ievent = 0;


  while(T->GetEntry(ievent++)){

    for(int ihit = 0; ihit < nhits; ihit++){
      	isamp_tot[(int)module[ihit]][h[(int)module[ihit]]->FindBin(yhit[ihit],xhit[ihit])] += isamp[ihit];
	nentries[(int)module[ihit]][h[(int)module[ihit]]->FindBin(yhit[ihit],xhit[ihit])] += 1;
    }

  }

  for(int imod=0; imod<nmodules;imod++)
    for(int ibin=0; ibin < nbinsx*nbinsy; ibin++)
      if(nentries[imod][ibin] > 0) 
	h[imod]->SetBinContent(ibin,isamp_tot[imod][ibin]*1.0/nentries[imod][ibin]);

  TCanvas *c = new TCanvas();

  TString outputfile = Form("../plots/timing_hit_map_%i.pdf",run);

  for(int imod=0; imod<nmodules;imod++){
    h[imod]->GetZaxis()->SetRangeUser(0,5);
    h[imod]->Draw("colz");

    if(imod == 0) c->Print(outputfile + "(");
    else if(imod == nmodules - 1) c->Print(outputfile + ")");
    else c->Print(outputfile);
  }
  

}
