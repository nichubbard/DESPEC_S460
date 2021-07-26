void alpha_219Rn()
{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Fri Apr 16 19:41:53 2021) by ROOT version 6.18/04
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",258,74,538,323);
   Canvas_1->Range(2.53061,85.39121,2.542061,86.65186);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   
   TCutG *cutg = new TCutG("cID_Z1_AoQ0",8);
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,2.53878,86.41227);
   cutg->SetPoint(1,2.540152,85.95519);
   cutg->SetPoint(2,2.540152,85.71927);
   cutg->SetPoint(3,2.536636,85.60132);
   cutg->SetPoint(4,2.532862,85.73402);
   cutg->SetPoint(5,2.532518,86.22059);
   cutg->SetPoint(6,2.534834,86.44175);
   cutg->SetPoint(7,2.53878,86.41227);
   cutg->Draw("alp");
   
   TPaveText *pt = new TPaveText(0.4461567,0.9351869,0.5538433,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("Graph");
   pt->Draw();
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
