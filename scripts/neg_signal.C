#include "pedestal_analysis.C"
#include "GEM_map_decode.C"
#include "Get_event_data.C"

//const int nruns = 4;
//int runs[nruns] = {2813,2815,2817,2820};
//int runs[nruns] = {2813};
//int currents[nruns] = {5,15,30,45};
TString run_title[nruns] = {"5","15","30","45"};

TString name_map[nmodules] = {"UV layer 0", "UV layer 1", "UV layer 2", "UV layer 3", "XY module 0", "XY module 1", "XY module 2", "XY module 3"}; 


void make_adc_plots(int irun, TCanvas *c, int color, TLegend *legend[8], TFile *f){
  
  int run = runs[irun];
   
  for(int imod=0; imod < nmodules; imod++){
    
    TH1F *h = (TH1F*)f->Get(Form("hADCpedsubU_allstrips_bb_gem_m%i",imod));
    h->SetTitle(name_map[imod]);

    c->cd(imod + 1);

    
    double xmin=350,xmax=3500;
    double ent=0;
    int minbin=h->GetXaxis()->FindBin(xmin);
    int maxbin=h->GetXaxis()->FindBin(xmax);
    for(int i=minbin;i<=maxbin;i++){
      ent+=h->GetBinContent(i);
    }

    double n_pos = ent;

    xmin=-500;
    xmax=-350;
    ent=0;
    minbin=h->GetXaxis()->FindBin(xmin);
    maxbin=h->GetXaxis()->FindBin(xmax);
    for(int i=minbin;i<=maxbin;i++){  
      ent+=h->GetBinContent(i);
    }

    double n_neg = ent;

    legend[imod]->AddEntry(h,Form(run_title[irun] + " %.2f",n_neg/n_pos*100));
    

    //h->GetXaxis()->SetRangeUser(-500,-100);
    h->Scale(1/h->GetMaximum());
    h->SetLineColor(color);
    h->Draw("same hist");

    gPad->SetLogy();
    legend[imod]->Draw("same");

    
    //cout<<currents[irun]<<" neg/pos = "<<n_neg/n_pos<<" "<<n_neg<<" "<<n_pos<<endl;
  
  }  


}


void make_CM_plots(int irun, TCanvas *c, int color, TLegend *legend[8], TFile *f){
  
  int run = runs[irun];

  for(int imod=0; imod < nmodules; imod++){
    
    TH2F *h2 = (TH2F*)f->Get(Form("hcommonmodeU_sorting_bb_gem_m%i",imod));

    double x[nAPV[imod][0]];
    double y[nAPV[imod][0]];
    double ex[nAPV[imod][0]];
    double ey[nAPV[imod][0]];

    for(int iapv=0; iapv < nAPV[imod][0]; iapv++){
      TH1D *h1 = h2->ProjectionY(Form("h%i",imod),iapv,iapv + 1);
      
      x[iapv] = iapv;
      y[iapv] = h1->GetMean();
      ex[iapv] = 0.5;
      ey[iapv] = h1->GetStdDev();
    
    }  
    
    TGraphErrors *g = new TGraphErrors(nAPV[imod][0],x,y,ex,ey);
    
    c->cd(imod+1);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->GetYaxis()->SetRangeUser(-60,60);
    g->SetTitle(name_map[imod] + ";APV;Common Mode Sorting - Common Mode Mean");
    
    if(color == 1) g->Draw("ap");
    else g->Draw("p");

    legend[imod]->Draw("same");
  }  


}



////// This is the main function ///////
void neg_signal(int run = 2813){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);



  TCanvas *c = new TCanvas("c","",1600,1200);
  c->Divide(4,2);

  TCanvas *c2 = new TCanvas("c2","",1600,1200);
  c2->Divide(4,2);

  TLegend *legend[nmodules];

  int color = 0;

  

  for(int imod=0; imod<nmodules; imod++){
    legend[imod] = new TLegend(0.48,0.7,0.9,0.9);
    legend[imod]->SetTextSize(0.04);
  }

  for(int irun=0; irun<nruns; irun++){
    color++;
    if(color == 5) color++;

    TString filename = Form("/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/GEM_luminosity/gen_replayed_%i_all.root",runs[irun]);

    TFile *f = new TFile(filename,"read");

    make_adc_plots(irun, c, color, legend, f); 
    make_CM_plots(irun, c2, color, legend, f); 
  }


  c->SaveAs("../plots/neg_signal_ADC.png");
  c2->SaveAs("../plots/neg_signal_CM.png");

}
