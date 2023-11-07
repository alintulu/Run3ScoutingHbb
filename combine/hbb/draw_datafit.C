/************************************************
 * Jennet Dickinson 
 * Nov 19, 2020
 * Draw Roofit plots
 ************************************************/
#include <iostream>

using namespace RooFit;
using namespace RooStats;

bool blind = true;

void draw(int index, bool pass, bool is_ggf, bool log=false){

  double rZbb = 1;
  string year = "2022";
  string year_string = "1.1/fb (13.6 TeV)";
  string thisbin = "pt"+to_string(index+1);
  string thisbin_fit = "ptbin"+to_string(index) + "ggf";

  /* DATA */
  TFile* dataf = new TFile("signalregion.root");
  TH1D* data_obs;
  if( !pass && is_ggf )
    data_obs = (TH1D*)dataf->Get(("ggf_fail_pt"+to_string(index+1)+"_data_nominal").c_str());
  else if( pass && is_ggf)
    data_obs = (TH1D*)dataf->Get(("ggf_pass_pt"+to_string(index+1)+"_data_nominal").c_str());

  cout << data_obs->Integral() << endl;

  // blind data!
  if( blind && pass ){                                                                                        
    for(int i=10; i<15; i++){
      data_obs->SetBinContent(i,0);
      data_obs->SetBinError(i,0);
    }                            
  }                                                                                                          

  data_obs->SetLineColor(kBlack);
  data_obs->SetMarkerColor(kBlack);
  data_obs->SetMarkerStyle(20);

  string filename = "fitDiagnosticsTest.root";
  string name = thisbin_fit+"fail"+year;
  if( pass )
    name = thisbin_fit + "pass"+year;
  
  string histdirname = "shapes_fit_s/" + name+ "/";

  cout << histdirname << endl;

  TFile *f = new TFile(filename.c_str());

  /* ggF */
  TH1D* ggF = (TH1D*)f->Get((histdirname+"ggF").c_str());
  ggF->Scale(7.0);
  ggF->SetLineColor(kGreen+1);
  ggF->SetMarkerColor(kGreen+1);
  //ggF->SetLineColor(kRed+1);
  //ggF->SetMarkerColor(kRed+1);
  //ggF->SetLineStyle(2);
  ggF->SetLineWidth(3);

  /* bkg Higgs */
  TH1D* bkgHiggs = (TH1D*)f->Get((histdirname+"VBF").c_str());
  bkgHiggs->Reset();
  bkgHiggs->Add((TH1D*)f->Get((histdirname+"WH").c_str()));
  bkgHiggs->Add((TH1D*)f->Get((histdirname+"ZH").c_str()));
  bkgHiggs->Add((TH1D*)f->Get((histdirname+"ttH").c_str()));
  bkgHiggs->Add((TH1D*)f->Get((histdirname+"VBF").c_str()));
  bkgHiggs->Scale(7.0);
  bkgHiggs->SetLineWidth(1);
  bkgHiggs->SetLineColor(kBlack);
  bkgHiggs->SetFillColor(kOrange);

  /* VV */
  TH1D* VV = (TH1D*)f->Get((histdirname+"VV").c_str());;
  VV->Scale(7.0);
  VV->SetLineWidth(1);
  VV->SetLineColor(kBlack);
  VV->SetFillColor(kOrange-3);

  /* ttbar */
  TH1D* ttbar = (TH1D*)f->Get((histdirname+"TTbar").c_str());
  ttbar->Scale(7.0);
  ttbar->SetLineColor(kBlack);
  ttbar->SetFillColor(kViolet-5);

  /* Z + jets */
  TH1D* Zjets = (TH1D*)f->Get((histdirname+"ZJetsqq").c_str());
  Zjets->Scale(7.0);
  Zjets->SetLineColor(kBlack);
  Zjets->SetFillColor(kAzure+8);

  /* Z(bb) + jets */
  TH1D* Zjetsbb = (TH1D*)f->Get((histdirname+"ZJetsbb").c_str());
  Zjetsbb->Scale(7.0);
  Zjetsbb->Scale(rZbb);
  Zjetsbb->SetLineColor(kBlack);
  Zjetsbb->SetFillColor(kAzure-1);

  /* W + jets */
  TH1D* Wjets = (TH1D*)f->Get((histdirname+"W").c_str());
  Wjets->Scale(7.0);
  Wjets->SetLineColor(kBlack);
  Wjets->SetFillColor(kGray);
  
  /* QCD */
  TH1D* qcd = (TH1D*)f->Get((histdirname+"qcd").c_str());
  qcd->Scale(7.0);
  qcd->SetLineColor(kBlack);
  qcd->SetFillColor(kWhite);

  /* total background */
  TH1D* TotalBkg = (TH1D*)f->Get((histdirname+"/total_background").c_str());
  TotalBkg->Scale(7.0);
  TotalBkg->SetMarkerColor(kRed);
  TotalBkg->SetLineColor(kRed);
  TotalBkg->SetFillColor(kRed);
  TotalBkg->SetFillStyle(3001);

  double max = TotalBkg->GetMaximum();
  double xmax = 201;
  double xmin = 40;

  if(thisbin == "pt1") xmax = 166;
  if(thisbin == "pt2") xmax = 180;
  if(thisbin == "pt6") xmin = 47;

  TotalBkg->GetYaxis()->SetRangeUser(0.1,1000*max);
  TotalBkg->GetXaxis()->SetRangeUser(xmin,xmax);
  if( !log ) TotalBkg->GetYaxis()->SetRangeUser(0,1.3*max);
  TotalBkg->GetYaxis()->SetTitle("Events / 7 GeV");
  TotalBkg->GetXaxis()->SetTitle("m_{reg} [GeV]");

  THStack *bkg = new THStack("bkg","");
  if( log ){
    bkg->Add(bkgHiggs);
    bkg->Add(VV);
    bkg->Add(ttbar);
    bkg->Add(Zjets);
    bkg->Add(Zjetsbb);
    bkg->Add(Wjets);
    bkg->Add(qcd);
  }
  else{
    bkg->Add(qcd);
    bkg->Add(Wjets);
    bkg->Add(Zjetsbb);
    bkg->Add(Zjets);
    bkg->Add(ttbar);
    bkg->Add(VV);
    bkg->Add(bkgHiggs);
  }

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas(name.c_str(),name.c_str(),600,600);
  TPad *pad1 = new TPad("pad1","pad1",0,.33,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,.33);

  pad1->SetBottomMargin(0.00001);
  pad1->SetTopMargin(0.1);
  pad1->SetBorderMode(0);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.3);
  pad2->SetBorderMode(0);

  pad1->SetLeftMargin(0.15);
  pad2->SetLeftMargin(0.15);
  pad1->Draw();
  pad2->Draw();

  float textsize1 = 16/(pad1->GetWh()*pad1->GetAbsHNDC());
  float textsize2 = 16/(pad2->GetWh()*pad2->GetAbsHNDC());

  TotalBkg->GetYaxis()->SetTitleSize(textsize1);
  TotalBkg->GetYaxis()->SetLabelSize(textsize1);
  TotalBkg->GetYaxis()->SetTitleOffset(2*pad1->GetAbsHNDC());

  pad1->cd();
  if( log ) pad1->SetLogy();

  cout << "QCD: "     << qcd->Integral()     << endl;
  cout << "Wjets: "   << Wjets->Integral()   << endl;
  cout << "Zjets: "   << Zjets->Integral()   << endl;
  cout << "ttbar: "   << ttbar->Integral()   << endl;
  cout << "VV: "      << VV->Integral()      << endl;
  cout << "bkgHiggs: " << bkgHiggs->Integral() << endl;

  TotalBkg->Draw("e2");
  bkg->Draw("histsame");
  TotalBkg->Draw("e2same");
  ggF->Draw("histsame");
  data_obs->Draw("pesame");
  data_obs->Draw("axissame");
  
  double x1=.6, y1=.86;
  TLegend* leg = new TLegend(x1,y1,x1+.3,y1-.32);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(2);
  leg->SetTextSize(textsize1);

  leg->AddEntry(data_obs,"Data","p");
  leg->AddEntry(TotalBkg,"Bkg. Unc.","f");
  leg->AddEntry(qcd,"QCD","f");
  leg->AddEntry(Wjets,"W","f");
  leg->AddEntry(Zjets,"Z(qq)","f");
  leg->AddEntry(Zjetsbb,"Z(bb)","f");
  leg->AddEntry(ttbar,"t#bar{t}","f");
  leg->AddEntry(VV,"VV","f");
  leg->AddEntry(bkgHiggs,"Bkg. H","f");
  leg->AddEntry(ggF,"ggF","l");

  leg->Draw();

  TLatex l1;
  l1.SetNDC();
  l1.SetTextFont(42);
  l1.SetTextSize(textsize1);
  l1.DrawLatex(0.2,.82,"#bf{CMS} Preliminary");

  TLatex l2;
  l2.SetNDC();
  l2.SetTextFont(42);
  l2.SetTextSize(textsize1);
  l2.DrawLatex(0.7,.92,year_string.c_str());

  TLatex l3;
  l3.SetNDC();
  l3.SetTextFont(42);
  l3.SetTextSize(textsize1);
  string text = "DeepDoubleB Fail Region";
  if( pass )
    text = "DeepDoubleB Pass Region";
  l3.DrawLatex(0.2,.77,text.c_str());

  TLatex l4;
  l4.SetNDC();
  l4.SetTextFont(42);
  l4.SetTextSize(textsize1);
  string text2 = "";
  if( is_ggf )
    text2 += "ggF category p_{T} bin "+to_string(index+1);

  pad2->cd();
  
  TH1D* TotalBkg_sub = (TH1D*)TotalBkg->Clone("TotalBkg_sub");
  TotalBkg_sub->Reset();
  TH1D* data_obs_sub = (TH1D*)data_obs->Clone("data_obs_ratio");
  data_obs_sub->Reset();

  TH1D* ggF_sub = (TH1D*)ggF->Clone("ggF_sub");
  ggF_sub->Reset();

  for(int i=1; i<TotalBkg_sub->GetNbinsX()+1; i++){
    TotalBkg_sub->SetBinError(i,TotalBkg->GetBinError(i)/data_obs->GetBinError(i));


    data_obs_sub->SetBinContent(i,(data_obs->GetBinContent(i)-TotalBkg->GetBinContent(i))/data_obs->GetBinError(i));
    data_obs_sub->SetBinError(i,data_obs->GetBinError(i)/data_obs->GetBinError(i));

    ggF_sub->SetBinContent(i,ggF->GetBinContent(i)/data_obs->GetBinError(i));
  }

  TotalBkg_sub->GetYaxis()->SetTitleSize(textsize2);
  TotalBkg_sub->GetYaxis()->SetLabelSize(textsize2);
  TotalBkg_sub->GetXaxis()->SetTitleSize(textsize2);
  TotalBkg_sub->GetXaxis()->SetLabelSize(textsize2);
  TotalBkg_sub->GetYaxis()->SetTitleOffset(2*pad2->GetAbsHNDC());
  TotalBkg_sub->GetYaxis()->SetTitle("(Data - Bkg)/#sigma_{Data}");
  TotalBkg_sub->SetMarkerSize(0);

  // blind data!                                                                                                                                                              
  if( blind  ){
    for(int i=10; i<15; i++){
      TotalBkg_sub->SetBinError(i,0);
      data_obs_sub->SetBinContent(i,0);
      data_obs_sub->SetBinError(i,0);
    }
  }

  double min2 = data_obs_sub->GetMinimum();
  double max2 = data_obs_sub->GetMaximum();

  if ( ! pass) {
     max2 += 1;
     min2 -= 1;
  }

  TotalBkg_sub->GetYaxis()->SetRangeUser(min2 * 1.3, max2 * 1.3);

  TotalBkg_sub->Draw("e2");
  ggF_sub->Draw("histsame");
  data_obs_sub->Draw("pesame");

  if( !log ) name += "_lin";

  c->SaveAs(("plots/"+name+".png").c_str());
  c->SaveAs(("plots/"+name+".pdf").c_str());

  return;

}

void draw_datafit(){

  for(int i=0; i<5; i++){
    draw(i,0,1);
    draw(i,1,1);
  }

  return 0;

}
