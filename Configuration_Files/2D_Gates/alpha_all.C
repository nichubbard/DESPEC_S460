void alpha_all()
{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Fri Apr 16 21:04:09 2021) by ROOT version 6.18/04
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",258,952,538,323);
   Canvas_1->Range(2.325579,66.28988,2.723435,99.25883);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   
   TCutG *cutg = new TCutG("cID_Z1_AoQ1",8);
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,2.596187,93.76401);
   cutg->SetPoint(1,2.650111,92.75498);
   cutg->SetPoint(2,2.657126,71.78471);
   cutg->SetPoint(3,2.59794,72.17954);
   cutg->SetPoint(4,2.593995,72.31116);
   cutg->SetPoint(5,2.391888,80.33952);
   cutg->SetPoint(6,2.494038,91.78982);
   cutg->SetPoint(7,2.596187,93.76401);
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
