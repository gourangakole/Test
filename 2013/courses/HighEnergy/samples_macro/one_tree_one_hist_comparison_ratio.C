//=====Last used ======
// root -l one_tree_one_hist_comparison_ratio.C
//====================

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;
void one_tree_one_hist_comparison_ratio()
{
  //  gROOT->ProcessLine(".L tdrstyle.C"); // use when available
  //  setTDRStyle(); // just to set tdrstyle
  gStyle->SetOptStat(0);
  //TString copyname = plotname;
  //  TCanvas *c3 = new TCanvas("c3","A Simple Graph Example",200,10,700,500);
  TString Initial("base/");
  //plotname = Initial + plotname;

  TLegend *leg = new TLegend(0.60,0.50,0.75,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  // h1 means inclusive-loose & h2: loose but not tight(lNt)

  TString loosePrefix = "./";
  TFile *mine = new TFile(loosePrefix+"tree_13.root");

  TString lNtPrefix = "./";
  TFile *official = new TFile(lNtPrefix+"mcPileupUL2016.root");

  // loose
  //TCanvas *c1 = (TCanvas*)loose->Get("TauPt_AfterSelections");
  //TH1F *h1 = (TH1F*)c1->GetPrimitive("FakeTau");
  //  TH1F *h1 = (TH1F*)(wh_m_100->Get(plotname))->Clone("h1");
  TH1F *h1 = new TH1F("h1", "h1", 100, 0, 100.);
  TTree *tree = ((TTree*)mine->Get("Events"));
  tree->Draw("Pileup_nTrueInt>>h1","","goff");
  
  h1->SetLineColor(1);
  h1->SetLineWidth(2);


  //h1->Draw("hist");
  h1->Scale(1.0/h1->Integral());
  //h1->SetTitle(copyname);
  //h1->GetXaxis()->SetTitle("dijet mass");
  //h1->GetYaxis()->SetTitle("Normalized to unity/5 GeV");
  
  // lNt
  //TCanvas *c2 = (TCanvas*)lNt->Get("TauPt_AfterSelections");
  //TH1F *h2 = (TH1F*)c2->GetPrimitive("FakeTau");
  
  //  TH1F *h2 = (TH1F*)(ttbar->Get(plotname))->Clone("h2");
  TH1F *h2 = (TH1F*)(official->Get("pu_mc"))->Clone("h2");
  h2->SetLineColor(2);
  h2->SetLineWidth(2);
  //h2->Draw("hist");
  h2->Scale(1.0/h2->Integral());
  //  h2->SetTitle(copyname);
  //  h2->GetXaxis()->SetTitle("dijet mass");
  //  h2->GetYaxis()->SetTitle("Normalized to unity/5 GeV");
  
  //cout << "lNt: Mean->" << h2->GetMean() <<" Integral-> " << h2->Integral() << endl;
  //cout << "loose: Mean->" << h1->GetMean() <<" Integral-> " << h1->Integral() << endl;

  TCanvas *c3 = new TCanvas("c3","A Simple Graph Example",200,10,700,700);
  const float xpad[2] = {0.,1.0};
  const float ypad[4] = {0.,0.30,0.30,1.0};
  c3->Divide(1, 2); c3->cd(1);
  gPad->SetRightMargin(0.03);
  gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);

  gPad->SetLogy(1);
  h1->SetTitle("Pile up comparison");
  h1->GetYaxis()->SetTitle("Normalized to Unity");
  h1->GetXaxis()->SetTitle("NPU");

  h1->Draw("hist");
  h2->Draw("hist same");
  
  leg->AddEntry(h1,"CRAB DY","l");
  leg->AddEntry(h2,"mcPUUL2016","l");
  
  leg->Draw();
  
  // ratio
  c3->cd(2);
  c3->Update();
  c3->RedrawAxis();
  c3->GetFrame()->Draw();
  gPad->SetTopMargin(0); 
  gPad->SetBottomMargin(0.30); //gPad->SetGridy();
  gPad->SetRightMargin(0.03);
  gPad->SetTickx(0);
  //if(histDir=="") gPad->SetBottomMargin(0.55);
  gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
  TH1F *hRatio = (TH1F*)h2->Clone("hRatio");
  hRatio->Reset();
  hRatio->Add(h1);
  if (h2->Integral() > 0.0)
    hRatio->Divide(h2);
  hRatio->GetYaxis()->SetRangeUser(0.5, 2.5);
  hRatio->GetYaxis()->SetTitle("Ratio");
  hRatio->Draw("E");


  
  // TCanvas::Update() draws the frame, after which one can change it
  c3->Update();
  //   c1->GetFrame()->SetFillColor(21);
  c3->GetFrame()->SetBorderSize(12);
  c3->Modified();
  TString outFile("");
  //outFile += plotname;
  outFile += "PU_comparison.png";
  c3->SaveAs(outFile);
}

//void example_stack_all(){
//  example_stack("mjj_kfit");
//}


