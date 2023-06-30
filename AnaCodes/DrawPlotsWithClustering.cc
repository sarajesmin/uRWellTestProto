/* 
 * File:   DrawPlotsWithClustering.cc
 * Author: rafopar
 *
 * Created on June 7, 2023, 5:14 PM
 */

#include <cstdlib>

using namespace std;

void DrawActiveArea();

/*
 * 
 */
int DrawPlotsWithClustering(int run) {
    
    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetTopMargin(0.04);
    c1->SetRightMargin(0.04);
    TLatex *lat1 = new TLatex();
    lat1->SetNDC();
    lat1->SetTextFont(42);

    TLine *line1 = new TLine();
    line1->SetLineColor(2);
    line1->SetLineWidth(2);
    
    TFile *file_in = new TFile(Form("AnaClustering_%d.root", run), "Read");
    
    TH1D *h_U_HitStrip1 = (TH1D*)file_in->Get("h_U_HitStrip1");
    h_U_HitStrip1->SetTitle("; U strip number");
    TH1D *h_U_HitStrip2 = (TH1D*)file_in->Get("h_U_HitStrip2");
    h_U_HitStrip2->SetTitle("; U strip number");
    h_U_HitStrip2->SetLineColor(2);
    h_U_HitStrip1->Draw();
    h_U_HitStrip2->Draw("Same");
    lat1->SetTextColor(4);
    lat1->DrawLatex(0.6, 0.8, "3#sigma cut");
    lat1->SetTextColor(2    );
    lat1->DrawLatex(0.6, 0.7, "5#sigma cut");
    c1->Print(Form("Figs/U_Strip_Distr1_%d.pdf", run));
    c1->Print(Form("Figs/U_Strip_Distr1_%d.png", run));
    c1->Print(Form("Figs/U_Strip_Distr1_%d.root", run));

    TH1D *h_V_HitStrip1 = (TH1D*)file_in->Get("h_V_HitStrip1");
    h_V_HitStrip1->SetTitle("; V strip number");
    TH1D *h_V_HitStrip2 = (TH1D*)file_in->Get("h_V_HitStrip2");
    h_V_HitStrip2->SetTitle("; V strip number");
    h_V_HitStrip2->SetLineColor(2);
    h_V_HitStrip1->Draw();
    h_V_HitStrip2->Draw("Same");
    lat1->SetTextColor(4);
    lat1->DrawLatex(0.6, 0.8, "3#sigma cut");
    lat1->SetTextColor(2    );
    lat1->DrawLatex(0.6, 0.7, "5#sigma cut");
    c1->Print(Form("Figs/V_Strip_Distr1_%d.pdf", run));
    c1->Print(Form("Figs/V_Strip_Distr1_%d.png", run));
    c1->Print(Form("Figs/V_Strip_Distr1_%d.root", run));

    
    TH1D *h_n_UclHits1 = (TH1D*)file_in->Get("h_n_UclHits1");
    h_n_UclHits1->SetLineWidth(2);
    h_n_UclHits1->SetTitle("; The # of Hits in U clusters");
    int bin_2 = h_n_UclHits1->FindBin(2);
    h_n_UclHits1->SetMaximum(1.05* h_n_UclHits1->GetBinContent(bin_2) );
    h_n_UclHits1->SetAxisRange(2, h_n_UclHits1->GetBinCenter( h_n_UclHits1->GetNbinsX() ) );
    h_n_UclHits1->GetXaxis()->SetNdivisions(716);
    h_n_UclHits1->Draw();
    c1->Print(Form("Figs/N_Ucl_Hits_%d.pdf", run));
    c1->Print(Form("Figs/N_Ucl_Hits_%d.png", run));
    c1->Print(Form("Figs/N_Ucl_Hits_%d.root", run));
    
    TH1D *h_n_VclHits1 = (TH1D*)file_in->Get("h_n_VclHits1");
    h_n_VclHits1->SetLineWidth(2);
    h_n_VclHits1->SetTitle("; The # of Hits in V clusters");
    bin_2 = h_n_VclHits1->FindBin(2);
    h_n_VclHits1->SetMaximum(1.05* h_n_VclHits1->GetBinContent(bin_2) );
    h_n_VclHits1->SetAxisRange(2, h_n_VclHits1->GetBinCenter( h_n_VclHits1->GetNbinsX() ) );
    h_n_VclHits1->GetXaxis()->SetNdivisions(716);
    h_n_VclHits1->Draw();
    c1->Print(Form("Figs/N_Vcl_Hits_%d.pdf", run));
    c1->Print(Form("Figs/N_Vcl_Hits_%d.png", run));
    c1->Print(Form("Figs/N_Vcl_Hits_%d.root", run));
    
    TF1 *f_Landau = new TF1("f_Landau", "[0]*TMath::Landau(x, [1], [2])", 0., 1000.);
    f_Landau->SetNpx(4500);
            
    double adcMaxAxis = 500.;
    TH1D *h_U_PeakADC_MultiCl1 = (TH1D*)file_in->Get("h_U_PeakADC_MultiCl1");
    h_U_PeakADC_MultiCl1->SetLineWidth(2);
    h_U_PeakADC_MultiCl1->SetTitle("; ADC of U cluster peak");
    h_U_PeakADC_MultiCl1->SetAxisRange(0., adcMaxAxis, "X");
    f_Landau->SetParameters(5.* h_U_PeakADC_MultiCl1->GetMaximum(), h_U_PeakADC_MultiCl1->GetBinCenter( h_U_PeakADC_MultiCl1->GetMaximumBin() ) , 10);
    h_U_PeakADC_MultiCl1->Fit(f_Landau, "MeV", "", 0., adcMaxAxis);
    lat1->DrawLatex(0.5, 0.7, Form("MPV = %1.1f", f_Landau->GetParameter(1)));
    c1->Print(Form("Figs/Peak_ADC_UStrips_%d.pdf", run));
    c1->Print(Form("Figs/Peak_ADC_UStrips_%d.png", run));
    c1->Print(Form("Figs/Peak_ADC_UStrips_%d.root", run));
    
    TH1D *h_V_PeakADC_MultiCl1 = (TH1D*)file_in->Get("h_V_PeakADC_MultiCl1");
    h_V_PeakADC_MultiCl1->SetLineWidth(2);
    h_V_PeakADC_MultiCl1->SetTitle("; ADC of V cluster peak");
    h_V_PeakADC_MultiCl1->SetAxisRange(0., adcMaxAxis, "X");
    f_Landau->SetParameters(5.* h_V_PeakADC_MultiCl1->GetMaximum(), h_V_PeakADC_MultiCl1->GetBinCenter( h_V_PeakADC_MultiCl1->GetMaximumBin() ) , 2);
    h_V_PeakADC_MultiCl1->Fit(f_Landau, "MeV", "", 0., adcMaxAxis);
    lat1->DrawLatex(0.5, 0.7, Form("MPV = %1.1f", f_Landau->GetParameter(1)));
    c1->Print(Form("Figs/Peak_ADC_VStrips_%d.pdf", run));
    c1->Print(Form("Figs/Peak_ADC_VStrips_%d.png", run));
    c1->Print(Form("Figs/Peak_ADC_VStrips_%d.root", run));
    
    TH1D *h_U_Coord_MultrCl1 = (TH1D*)file_in->Get("h_U_Coord_MultrCl1");
    h_U_Coord_MultrCl1->SetTitle("; Cluster U coordinate [strip]");
    h_U_Coord_MultrCl1->SetLineWidth(2);
    h_U_Coord_MultrCl1->Draw();
    c1->Print(Form("Figs/cl_U_coordinate1_%d.pdf", run));
    c1->Print(Form("Figs/cl_U_coordinate1_%d.png", run));
    c1->Print(Form("Figs/cl_U_coordinate1_%d.root", run));
    
    TH1D *h_V_Coord_MultrCl1 = (TH1D*)file_in->Get("h_V_Coord_MultrCl1");
    h_V_Coord_MultrCl1->SetTitle("; Cluster V coordinate [strip]");
    h_V_Coord_MultrCl1->SetLineWidth(2);
    h_V_Coord_MultrCl1->Draw();
    c1->Print(Form("Figs/cl_V_coordinate1_%d.pdf", run));
    c1->Print(Form("Figs/cl_V_coordinate1_%d.png", run));
    c1->Print(Form("Figs/cl_V_coordinate1_%d.root", run));
    
    c1->SetLogz();
    TH2D *h_n_GEM_Y_vs_X_Hits1 = (TH2D*)file_in->Get("h_n_GEM_Y_vs_X_Hits1");
    h_n_GEM_Y_vs_X_Hits1->SetTitle("; The number of X hits in GEM; The number of Y hits in GEM");
    h_n_GEM_Y_vs_X_Hits1->Draw("colz");
    line1->DrawLine(1.5, 1.5, 1.5, 20.);
    line1->DrawLine(1.5, 1.5, 20., 1.5);
    c1->Print(Form("Figs/GEM_N_Y_vs_X_Hits1_%d.pdf", run));
    c1->Print(Form("Figs/GEM_N_Y_vs_X_Hits1_%d.png", run));
    c1->Print(Form("Figs/GEM_N_Y_vs_X_Hits1_%d.root", run));
    
    TH2D *h_GEM_XY1 = (TH2D*)file_in->Get("h_GEM_XY1");
    h_GEM_XY1->SetTitle("; GEM X Strip; GEM Y Strip");
    h_GEM_XY1->Draw("colz");
    c1->Print(Form("Figs/GEM_Y_vs_X_Strips_%d.pdf", run));
    c1->Print(Form("Figs/GEM_Y_vs_X_Strips_%d.png", run));
    c1->Print(Form("Figs/GEM_Y_vs_X_Strips_%d.root", run));
    
    TH2D *h_n_uRwell_V_vs_U_MultiHitCl = (TH2D*)file_in->Get("h_n_uRwell_V_vs_U_MultiHitCl");
    h_n_uRwell_V_vs_U_MultiHitCl->SetTitle("; Number of U clusters; umber of V clusters");
    int binx1 = h_n_uRwell_V_vs_U_MultiHitCl->GetXaxis()->FindBin(1);
    int binx2 = h_n_uRwell_V_vs_U_MultiHitCl->GetNbinsX() + 1;
    int biny1 = h_n_uRwell_V_vs_U_MultiHitCl->GetYaxis()->FindBin(1);
    int biny2 = h_n_uRwell_V_vs_U_MultiHitCl->GetNbinsY() + 1;
    double counts_integral = h_n_uRwell_V_vs_U_MultiHitCl->Integral();
    int counts_has_U_Cluster = h_n_uRwell_V_vs_U_MultiHitCl->Integral(binx1, binx2, 1, biny2);
    int counts_has_V_Cluster = h_n_uRwell_V_vs_U_MultiHitCl->Integral(1, binx2, biny1, biny2);
    int counts_has_AnyCluster = h_n_uRwell_V_vs_U_MultiHitCl->Integral() - h_n_uRwell_V_vs_U_MultiHitCl->GetBinContent(1, 1);
    int counts_has_U_AND_V_Cluster = h_n_uRwell_V_vs_U_MultiHitCl->Integral(binx1, binx2, biny1, biny2);
    h_n_uRwell_V_vs_U_MultiHitCl->SetAxisRange(0, 10., "Y");
    h_n_uRwell_V_vs_U_MultiHitCl->SetAxisRange(0, 10., "X");
    h_n_uRwell_V_vs_U_MultiHitCl->Draw("colz text");
    lat1->SetTextSize(0.03);
    lat1->DrawLatex(0.15, 0.9, Form("Has U cluster = %d #rightarrow %1.2f %%", counts_has_U_Cluster, 100.*double(counts_has_U_Cluster)/counts_integral ));
    lat1->DrawLatex(0.15, 0.85, Form("Has V cluster = %d #rightarrow %1.2f %% ", counts_has_V_Cluster, 100.*double(counts_has_V_Cluster)/counts_integral ));
    lat1->DrawLatex(0.15, 0.8, Form("Has any cluster = %d #rightarrow %1.2f %% ", counts_has_AnyCluster, 100.*double(counts_has_AnyCluster)/counts_integral ));
    lat1->DrawLatex(0.15, 0.75, Form("Has U and V cluster = %d #rightarrow %1.2f %% ", counts_has_U_AND_V_Cluster, 100.*double(counts_has_U_AND_V_Cluster)/counts_integral ));    
    c1->Print(Form("Figs/Number_OF_V_vs_U_clusters_%d.pdf", run));
    c1->Print(Form("Figs/Number_OF_V_vs_U_clusters_%d.png", run));
    c1->Print(Form("Figs/Number_OF_V_vs_U_clusters_%d.root", run));
    
    TCanvas *c2 = new TCanvas("c2", "", 1800., 1000.);
    c2->SetTopMargin(0.02);
    c2->SetRightMargin(0.02);
    TH2D *h_Cross_YXc1 = (TH2D*)file_in->Get("h_Cross_YXc1");
    h_Cross_YXc1->SetStats(0);
    h_Cross_YXc1->SetTitle("; Cross X coordinate [mm]; Cross Y coordinate [mm]");
    h_Cross_YXc1->Draw("col");
    DrawActiveArea();
    c2->Print(Form("Figs/Cross_YXc1_%d.pdf", run));
    c2->Print(Form("Figs/Cross_YXc1_%d.png", run));
    c2->Print(Form("Figs/Cross_YXc1_%d.root", run));
    
    TH2D *h_Cross_YXc2 = (TH2D*)file_in->Get("h_Cross_YXc2");
    h_Cross_YXc2->SetStats(0);
    h_Cross_YXc2->SetTitle("; Cross X coordinate [mm]; Cross Y coordinate [mm]");
    h_Cross_YXc2->Draw("col");
    DrawActiveArea();
    c2->Print(Form("Figs/Cross_YXc2_%d.pdf", run));
    c2->Print(Form("Figs/Cross_YXc2_%d.png", run));
    c2->Print(Form("Figs/Cross_YXc2_%d.root", run));
    
    return 0;
}

void DrawActiveArea() {

    TLine *line1 = new TLine();
    line1->SetLineColor(2);
    line1->SetLineWidth(3);
    line1->DrawLine(-723, 250., 723., 250.);
    line1->DrawLine(-506.14, -250., 506.14, -250.);
    line1->DrawLine(-506.14, -250., -723, 250.);
    line1->DrawLine(506.14, -250., 723, 250.);
}
