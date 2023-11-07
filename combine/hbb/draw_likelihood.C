// Draw output of combine 2d contour
TGraph* get_contour(string poi, string filestring){

  TFile *f = new TFile(("higgsCombine"+filestring+".MultiDimFit.mH125.root").c_str());

  // Variables
  float x = 0;
  float y = 0;

  TTree* t = (TTree*)f->Get("limit");
  t->SetBranchAddress(poi.c_str(),&x);
  t->SetBranchAddress("deltaNLL",&y);
  int npoints = t->GetEntries();

  double r[npoints-1];
  double dnll[npoints-1];

  for(int i=0; i<npoints; i++){                                                                            

    t->GetEntry(i);

    if(i>0){
      r[i-1] = x;
      dnll[i-1] = 2*y;
    }
  }

  TGraph* g = new TGraph(npoints-1,r,dnll);

  g->SetLineWidth(3);
  g->SetLineColor(kBlack);

  return g;
}

void draw_likelihood(){

  string year = "all";
  string year_string = "137/fb, Run 2";
  double rZbb = 1;

  gStyle->SetOptTitle(0);

  TCanvas* c0 = new TCanvas();

  double x[2] = {-10,10};
  double y1sigma[2] = {1,1};
  double y2sigma[2] = {4,4};
  double y3sigma[2] = {9,9};
  double y5sigma[2] = {25,25};

  TGraph* g1sigma = new TGraph(2,x,y1sigma);
  TGraph* g2sigma = new TGraph(2,x,y2sigma);
  TGraph* g3sigma = new TGraph(2,x,y3sigma);
  TGraph* g5sigma = new TGraph(2,x,y5sigma);

  g1sigma->SetLineColor(kGray);
  g2sigma->SetLineColor(kGray);
  g3sigma->SetLineColor(kGray);
  g5sigma->SetLineColor(kGray);

  g1sigma->SetLineWidth(3);
  g2sigma->SetLineWidth(3);
  g3sigma->SetLineWidth(3);
  g5sigma->SetLineWidth(3);

  g1sigma->SetLineStyle(3);
  g2sigma->SetLineStyle(3);
  g3sigma->SetLineStyle(3);
  g5sigma->SetLineStyle(3);

  TGraph* gggf = get_contour("rggF","rggF");

  gggf->Draw("AC");
  gggf->GetYaxis()->SetRangeUser(0,2);
  gggf->GetXaxis()->SetRangeUser(-3,4);
  gggf->GetXaxis()->SetTitle("#mu_{ggF}");
  gggf->GetYaxis()->SetTitle("q_{#mu}");

  g1sigma->Draw("same");
  g2sigma->Draw("same");
  gggf->Draw("Csame");

  return;

}
