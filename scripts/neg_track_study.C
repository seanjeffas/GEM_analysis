////////////////////////////////////////////////////////////////
//  Written by Sean Jeffas
//  sj9ry@virginia.edu
//  Last updated March 30, 2023
//
//  The purpose of this script is to take full readout GEM data
//  generated from SBS-offline and do a tracking analysis of the
//  negative signals. It will compare the negative signal occupancy,
//  efficiency, and distribution to positive signals to look for
//  any correlations.
//
//  Usage: root neg_track_study.C
//
////////////////////////////////////////////////////////////////
const int nfiber = 24;
const int nadc = 16;
const int nchan = 128;
const int nlayer = 5;
const int nruns = 4;
int runs[nruns] = {2813,2815,2817,2820};
int currents[nruns] = {5,15,30,45};


void neg_track_study(){

  gStyle->SetOptStat(0);

  TString titles[nruns] = {"5 uA on LD2", "10 uA on LD2", "12 uA on LD2", "20 uA on LD2"};

  TString output = "../plots/neg_track_study.pdf";

  TH1F *hADC_track_neg[nlayer][nruns];
  TH1F *hADC_all_neg[nlayer][nruns];

  TH1F *hADC_track_pos[nlayer][nruns];
  TH1F *hADC_all_pos[nlayer][nruns];
  
  TH1F *hstrip_track_neg[nlayer][nruns];
  TH1F *hstrip_all_neg[nlayer][nruns];

  TH1F *hstrip_SampMax_neg_good[nlayer][nruns];
  TH1F *hstrip_SampMax_neg_all[nlayer][nruns];
  TH1F *hstrip_SampMax_pos_good[nlayer][nruns];
  TH1F *hstrip_SampMax_pos_all[nlayer][nruns];
  
  TH1F *hstrip_track_pos[nlayer][nruns];
  TH1F *hstrip_all_pos[nlayer][nruns];

  TH2F *htrack_nstrip_miss[nlayer][nruns];
  TH2F *htrack_nstrip_hit[nlayer][nruns];

  double y_neg_frac[nlayer][nruns],y_pos_frac[nlayer][nruns],y_pos_eff[nlayer][nruns], y_neg_diff[nlayer][nruns];  
  double y_track_miss[nlayer][nruns],y_track_hit[nlayer][nruns];
  double x_neg_occu[nlayer][nruns], x_pos_occu[nlayer][nruns], x_current[nlayer][nruns];

  double x_max_neg = 0;
  double x_max_pos = 0;

  TString rootdir = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/GEM_luminosity/";
  
  //Loop over all runs
  for( int irun = 0; irun < nruns; irun++){
    
    TString Rootfile = Form(rootdir + "gen_replayed_%i_all.root",runs[irun]);

    //This bit of code here adds all root files together to analyze the total histograms
    gSystem->Exec("rm -f " + Rootfile);
    gSystem->Exec(Form("hadd -k -f %s %s/*%i*",Rootfile.Data(),rootdir.Data(),runs[irun]));
    
    TFile *hfile = new TFile(Rootfile,"read");  

    //Loop over all layers
    for( int ilayer = 0; ilayer < nlayer; ilayer++){
    
      x_current[ilayer][irun] = currents[irun]; //Get current for this run
        
      //Get all the histograms from the root files
      hADC_track_neg[ilayer][irun] = (TH1F*)hfile->Get(Form("hbb_gem_m%i_ADCmaxU_neg_good",ilayer));
      hADC_all_neg[ilayer][irun] = (TH1F*)hfile->Get(Form("hbb_gem_m%i_ADCmaxU_neg_all",ilayer));
      
      hADC_track_pos[ilayer][irun] = (TH1F*)hfile->Get(Form("hbb_gem_m%i_ADCmaxU_good",ilayer));
      hADC_all_pos[ilayer][irun] = (TH1F*)hfile->Get(Form("hbb_gem_m%i_ADCmaxU_all",ilayer));
      
      hstrip_track_neg[ilayer][irun] = (TH1F*)hfile->Get(Form("hbb_gem_m%i_stripU_neg_good",ilayer));
      hstrip_all_neg[ilayer][irun] = (TH1F*)hfile->Get(Form("hbb_gem_m%i_stripU_neg_all",ilayer));
      
      hstrip_track_pos[ilayer][irun] = (TH1F*)hfile->Get(Form("hbb_gem_m%i_stripU_good",ilayer));
      hstrip_all_pos[ilayer][irun] = (TH1F*)hfile->Get(Form("hbb_gem_m%i_stripU_all",ilayer));
      
      hstrip_SampMax_neg_good[ilayer][irun] = (TH1F*)hfile->Get("hbb_gem_m0_iSampMaxU_neg_good");
      hstrip_SampMax_neg_all[ilayer][irun] = (TH1F*)hfile->Get("hbb_gem_m0_iSampMaxU_neg_all");
      hstrip_SampMax_pos_good[ilayer][irun] = (TH1F*)hfile->Get("hbb_gem_m0_iSampMaxU_good");
      hstrip_SampMax_pos_all[ilayer][irun] = (TH1F*)hfile->Get("hbb_gem_m0_iSampMaxU_all");

      TH1F *hshould_neg = (TH1F*)hfile->Get(Form("hdidnothit_x_bb_gem_layer%i",ilayer));
      TH1F *hdid_neg = (TH1F*)hfile->Get(Form("hneghit_x_bb_gem_layer%i",ilayer));
      TH1F *hdid_neg_good = (TH1F*)hfile->Get(Form("hneghit_good_x_bb_gem_layer%i",ilayer));
      TH2F *hoccu_neg = (TH2F*)hfile->Get("hbb_gem_NstripsU_layer_neg");
      
      TH1F *hshould_pos = (TH1F*)hfile->Get(Form("hshouldhit_x_bb_gem_layer%i",ilayer));
      TH1F *hdid_pos = (TH1F*)hfile->Get(Form("hdidhit_x_bb_gem_layer%i",ilayer));
      TH1F *hdid_fullreadout = (TH1F*)hfile->Get(Form("hdidhit_fullreadout_x_bb_gem_layer%i",ilayer));
      TH2F *hoccu_pos = (TH2F*)hfile->Get("hbb_gem_NstripsU_layer");
      
      TH2F *htrack_nstrip_miss = (TH2F*)hfile->Get("hbb_gem_NstripsU_layer_neg_miss");
      TH2F *htrack_nstrip_hit = (TH2F*)hfile->Get("hbb_gem_NstripsU_layer_neg_hit");

      TH1D * htemp = hoccu_neg->ProjectionY("",0 + ilayer,1 + ilayer);
      htemp->GetXaxis()->SetRangeUser(3,2000);
      
      int nAPVs = 30;
      if(ilayer == 4) nAPVs = 10*4;

      //Get the occupancy numbers
      x_neg_occu[ilayer][irun] = htemp->GetMean() / (128*nAPVs);
      x_pos_occu[ilayer][irun] = hoccu_pos->ProjectionY("",0 + ilayer,1 + ilayer)->GetMean() / (128*30);
      
      if(x_neg_occu[ilayer][irun] > x_max_neg) x_max_neg = x_neg_occu[ilayer][irun];
      if(x_pos_occu[ilayer][irun] > x_max_pos) x_max_pos = x_pos_occu[ilayer][irun];

      //Get the efficiency numbers
      y_pos_eff[ilayer][irun] = hdid_pos->GetEntries()/hshould_pos->GetEntries();
      y_neg_frac[ilayer][irun] = hdid_neg->GetEntries()/hshould_neg->GetEntries();
      y_pos_frac[ilayer][irun] = hdid_neg_good->GetEntries()/hdid_fullreadout->GetEntries();
      y_neg_diff[ilayer][irun] = y_neg_frac[ilayer][irun] - y_pos_frac[ilayer][irun];
      
      
      htrack_nstrip_miss->GetYaxis()->SetRangeUser(0.8,15);
      htrack_nstrip_hit->GetYaxis()->SetRangeUser(0.8,15);
      
      y_track_miss[ilayer][irun] = htrack_nstrip_miss->ProjectionY("",0 + ilayer,1 + ilayer)->GetMean();
      y_track_hit[ilayer][irun] = htrack_nstrip_hit->ProjectionY("",0 + ilayer,1 + ilayer)->GetMean();
      

      
    }


 }

  TCanvas *c0 = new TCanvas("c0","",1600,600);
  c0->Divide(2,1);

  TCanvas *c1 = new TCanvas("c1","",1600,1000);
  c1->Divide(2,2);

  TCanvas *c2 = new TCanvas("c2","",1600,1000);
  c2->Divide(2,2);

  TCanvas *c3 = new TCanvas("c3","",1600,1000);
  c3->Divide(2,2);

  TCanvas *c4 = new TCanvas("c4","",1600,1000);
  c4->Divide(2,2);

  TLegend *legend = new TLegend(0.10,0.62,0.49,0.9);
  
  int icolor = 0;

  for(int irun=0; irun < nruns; irun++){
    icolor++;
    if(icolor == 5 || icolor == 10) icolor ++;

    c1->cd(1);
    hADC_all_neg[0][irun]->SetTitle("UV Layer 0 All Negative Strips");
    hADC_all_neg[0][irun]->SetLineColor(icolor);
    hADC_all_neg[0][irun]->Scale(1/hADC_all_neg[0][irun]->GetMaximum());
    hADC_all_neg[0][irun]->Draw("same hist");
    gPad->SetLogy();

    legend->AddEntry(hADC_all_neg[0][irun],titles[irun]);
    //legend->Draw("same");

    c1->cd(2);
    hADC_track_neg[0][irun]->SetTitle("UV Layer 0 Negative Strips On Tracks");
    hADC_track_neg[0][irun]->SetLineColor(icolor);
    hADC_track_neg[0][irun]->Scale(1/hADC_track_neg[0][irun]->GetMaximum());
    hADC_track_neg[0][irun]->Draw("same hist");
    gPad->SetLogy();
    

    legend->Draw("same");

    c1->cd(3);
    hADC_all_pos[0][irun]->SetTitle("UV Layer 0 All Positive Strips");
    hADC_all_pos[0][irun]->SetLineColor(icolor);
    hADC_all_pos[0][irun]->Scale(1/hADC_all_pos[0][irun]->GetMaximum());
    hADC_all_pos[0][irun]->Draw("same hist");
    gPad->SetLogy();

    c1->cd(4);
    hADC_track_pos[0][irun]->SetTitle("UV Layer 0 Positive Strips On Tracks");
    hADC_track_pos[0][irun]->SetLineColor(icolor);
    hADC_track_pos[0][irun]->Scale(1/hADC_track_pos[0][irun]->GetMaximum());
    hADC_track_pos[0][irun]->Draw("same hist");
    gPad->SetLogy();
   

  }

  c2->cd(1);
  hstrip_all_neg[0][2]->SetTitle("UV Layer 0 All Negative Strips");
  hstrip_all_neg[0][2]->Draw("same hist");

  c2->cd(2);
  hstrip_track_neg[0][2]->SetTitle("UV Layer 0 Negative Strips On Tracks");
  hstrip_track_neg[0][2]->Draw("same hist");
  
  c2->cd(3);
  hstrip_all_pos[0][2]->SetTitle("UV Layer 0 All Positive Strips");
  hstrip_all_pos[0][2]->Draw("same hist");

  c2->cd(4);
  hstrip_track_pos[0][2]->SetTitle("UV Layer 0 Positive Strips On Tracks");
  hstrip_track_pos[0][2]->Draw("same hist");


  c4->cd(1);
  hstrip_SampMax_neg_all[0][2]->SetTitle("Layer 0 All Negative Hits");
  hstrip_SampMax_neg_all[0][2]->Draw();

  c4->cd(2);
  hstrip_SampMax_neg_good[0][2]->SetTitle("Layer 0 Negative Hits on Tracks");
  hstrip_SampMax_neg_good[0][2]->Draw();

  c4->cd(3);
  hstrip_SampMax_pos_all[0][2]->SetTitle("Layer 0 All Positive Hits");
  hstrip_SampMax_pos_all[0][2]->Draw();

  c4->cd(4);
  hstrip_SampMax_pos_good[0][2]->SetTitle("Layer 0 Positive Hits on Tracks");
  hstrip_SampMax_pos_good[0][2]->Draw();
  

  TGraph *g_neg[nlayer];
  TGraph *g_neg_good[nlayer];
  TGraph *g_neg_diff[nlayer];
  TGraph *g_pos[nlayer];
  TGraph *g_occu_pos[nlayer];
  TGraph *g_occu_neg[nlayer];

  TGraph *g_track_miss[nlayer];
  TGraph *g_track_hit[nlayer];


  TLegend *legend2 = new TLegend(0.10,0.62,0.40,0.9);
  icolor = 0;

  for(int ilayer = 0; ilayer < nlayer; ilayer++){
    
    icolor++;
    if(icolor == 5) icolor++;

    
    g_neg[ilayer] = new TGraph(nruns,x_neg_occu[ilayer],y_neg_frac[ilayer]);
    g_neg_good[ilayer] = new TGraph(nruns,x_neg_occu[ilayer],y_pos_frac[ilayer]);
    g_neg_diff[ilayer] = new TGraph(nruns,x_neg_occu[ilayer],y_neg_diff[ilayer]);
    g_pos[ilayer] = new TGraph(nruns,x_pos_occu[ilayer],y_pos_eff[ilayer]);
    g_occu_pos[ilayer] = new TGraph(nruns,x_current[ilayer],x_pos_occu[ilayer]);
    g_occu_neg[ilayer] = new TGraph(nruns,x_current[ilayer],x_neg_occu[ilayer]);
    
    g_track_miss[ilayer] = new TGraph(nruns,x_neg_occu[ilayer],y_track_miss[ilayer]);
    g_track_hit[ilayer] = new TGraph(nruns,x_neg_occu[ilayer],y_track_hit[ilayer]);
    
    g_neg[ilayer]->SetTitle("Negative Strip Fraction on Missing Hits;Negative Occupancy;Fraction");    
    g_neg_good[ilayer]->SetTitle("Negative Strip Fraction on Found Hits;Negative Occupancy;Fraction");
    g_neg_diff[ilayer]->SetTitle("Negative Strip Fraction Difference;Negative Occupancy;% of missing hits - % of found hits");
    g_pos[ilayer]->SetTitle("Positive Strip Efficiency;Positive Occupancy;Efficiency");
    g_occu_pos[ilayer]->SetTitle("Positive Strip Occupancy;Beam Current (#mu A);Occupancy");
    g_occu_neg[ilayer]->SetTitle("Negative Strip Occupancy;Beam Current (#mu A);Occupancy");
    
    g_neg[ilayer]->SetMarkerStyle(20 + ilayer);
    g_neg_good[ilayer]->SetMarkerStyle(20 + ilayer);
    g_pos[ilayer]->SetMarkerStyle(20 + ilayer);
    g_neg_diff[ilayer]->SetMarkerStyle(20 + ilayer);
    g_occu_pos[ilayer]->SetMarkerStyle(20 + ilayer);
    g_occu_neg[ilayer]->SetMarkerStyle(20 + ilayer);
    
    g_neg[ilayer]->SetMarkerColor(icolor);
    g_neg_good[ilayer]->SetMarkerColor(icolor);
    g_neg_diff[ilayer]->SetMarkerColor(icolor);
    g_pos[ilayer]->SetMarkerColor(icolor);
    g_occu_pos[ilayer]->SetMarkerColor(icolor);
    g_occu_neg[ilayer]->SetMarkerColor(icolor);
   
    g_neg[ilayer]->GetYaxis()->SetRangeUser(0,1);
    g_neg_good[ilayer]->GetYaxis()->SetRangeUser(0,1);
    g_neg_diff[ilayer]->GetYaxis()->SetRangeUser(-0.2,0.2);
    g_pos[ilayer]->GetYaxis()->SetRangeUser(0,1);
    g_occu_pos[ilayer]->GetYaxis()->SetRangeUser(0,1);
    g_occu_neg[ilayer]->GetYaxis()->SetRangeUser(0,1);

    g_neg[ilayer]->GetXaxis()->SetLimits(0,x_max_neg*1.1);
    g_neg_good[ilayer]->GetXaxis()->SetLimits(0,x_max_neg*1.1);
    g_neg_diff[ilayer]->GetXaxis()->SetLimits(0,x_max_neg*1.1);
    g_pos[ilayer]->GetXaxis()->SetLimits(0,x_max_pos*1.1);
    
    

    legend2->AddEntry(g_neg[ilayer],Form("Layer %i",ilayer),"p");


    c3->cd(1);
    g_neg[ilayer]->GetYaxis()->SetTitleOffset(1.4);
    if(ilayer == 0) g_neg[ilayer]->Draw("AP");
    else g_neg[ilayer]->Draw("P same");

    c3->cd(2);
    if(ilayer == 0) g_neg_good[ilayer]->Draw("AP");
    else g_neg_good[ilayer]->Draw("P same");  
    
    c3->cd(3);
    if(ilayer == 0) g_neg_diff[ilayer]->Draw("AP");
    else g_neg_diff[ilayer]->Draw("P same");    

    c3->cd(4);
    if(ilayer == 0) g_pos[ilayer]->Draw("AP");
    else g_pos[ilayer]->Draw("P same");

    c0->cd(1);
    if(ilayer == 0) g_occu_pos[ilayer]->Draw("AP");
    else g_occu_pos[ilayer]->Draw("P same");
    
    c0->cd(2);
    if(ilayer == 0) g_occu_neg[ilayer]->Draw("AP");
    else g_occu_neg[ilayer]->Draw("P same");

 }

  c3->cd(1);
  legend2->Draw("same");

  c0->Print(output + "(");
  c1->Print(output);
  c2->Print(output);
  c3->Print(output);
  c4->Print(output + ")");

}
