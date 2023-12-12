/* 
 * File:   DrawPlotsWithClustering.cc
 * Author: rafopar
 *
 * Created on June 7, 2023, 5:14 PM
 */

#include <cstdlib>
#include <cxxopts.hpp>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>

#include <uRwellTools.h>
using namespace std;


void DrawActiveArea();

/*
 * 
 */
int main(int argc, char **argv) {

    cxxopts::Options options("AnaClustering", "Performs clustering and also does analysis on cosmic data");

    options.add_options()
            ("r,Run", "Run number", cxxopts::value<int>())
            ("t,Threshold", "Hit Threshold in terms of sigma", cxxopts::value<double>()->default_value("5"))
            ("m,MinHits", "Number of minimum hits in the cluster", cxxopts::value<int>()->default_value("1"))
            //("f,File", "File number of the given run", cxxopts::value<int>())
            ;

    auto parsed_options = options.parse(argc, argv);

    int run = 0;
    int fnum = -1;


    if (parsed_options.count("Run")) {
        run = parsed_options["Run"].as<int>();
    } else {
        cout << "The run number is nor provided. Exiting..." << endl;
        exit(1);
    }


    const double threshold = parsed_options["Threshold"].as<double>();
    if (!parsed_options.count("Threshold")) {
        cout << "* You didn't provide the hit threshold. Exitint" << endl;
        exit(1);
    }

    const int MinClSize = parsed_options["MinHits"].as<int>();
    if (!parsed_options.count("MinHits")) {
        cout << "* You didn't provide the Minimum hits int the cluster. Exiting " << endl;
        exit(1);
    }

    cout << "The hit threshold is " << threshold << "\\sigma" << endl;
    cout << "The Minimum cluster size is " << MinClSize << "hits" << endl;


    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetTopMargin(0.04);
    c1->SetRightMargin(0.04);
    TLatex *lat1 = new TLatex();
    lat1->SetNDC();
    lat1->SetTextFont(42);

    TLine *line1 = new TLine();
    line1->SetLineColor(2);
    line1->SetLineWidth(2);
    TFile *file_in = new TFile(Form("AnaClustering_%d_Thr_%1.1f_MinHits_%d.root", run, threshold, MinClSize), "Read");

    TH1D *h_U_HitStrip1 = (TH1D*) file_in->Get("h_U_HitStrip1");
    h_U_HitStrip1->SetTitle("; U strip number");
    TH1D *h_U_HitStrip2 = (TH1D*) file_in->Get("h_U_HitStrip2");
    h_U_HitStrip2->SetTitle("; U strip number");
    h_U_HitStrip2->SetLineColor(2);
    h_U_HitStrip1->Draw();
    h_U_HitStrip2->Draw("Same");
    lat1->SetTextColor(4);
    lat1->DrawLatex(0.6, 0.8, "3#sigma cut");
    lat1->SetTextColor(2);
    lat1->DrawLatex(0.6, 0.7, "5#sigma cut");
    c1->Print(Form("Figs/U_Strip_Distr1_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/U_Strip_Distr1_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/U_Strip_Distr1_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    TH1D *h_V_HitStrip1 = (TH1D*) file_in->Get("h_V_HitStrip1");
    h_V_HitStrip1->SetTitle("; V strip number");
    TH1D *h_V_HitStrip2 = (TH1D*) file_in->Get("h_V_HitStrip2");
    h_V_HitStrip2->SetTitle("; V strip number");
    h_V_HitStrip2->SetLineColor(2);
    h_V_HitStrip1->Draw();
    h_V_HitStrip2->Draw("Same");
    lat1->SetTextColor(4);
    lat1->DrawLatex(0.6, 0.8, "3#sigma cut");
    lat1->SetTextColor(2);
    lat1->DrawLatex(0.6, 0.7, "5#sigma cut");
    c1->Print(Form("Figs/V_Strip_Distr1_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/V_Strip_Distr1_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/V_Strip_Distr1_%d_t%1.1f_m%d.root", run, threshold, MinClSize));


    TH1D *h_n_UclHits1 = (TH1D*) file_in->Get("h_n_UclHits1");
    h_n_UclHits1->SetLineWidth(2);
    h_n_UclHits1->SetTitle("; The # of Hits in U clusters");
    int bin_2 = h_n_UclHits1->FindBin(2);
    h_n_UclHits1->SetMaximum(1.05 * h_n_UclHits1->GetBinContent(bin_2));
    h_n_UclHits1->SetAxisRange(2, h_n_UclHits1->GetBinCenter(h_n_UclHits1->GetNbinsX()));
    h_n_UclHits1->GetXaxis()->SetNdivisions(716);
    h_n_UclHits1->Draw();
    c1->Print(Form("Figs/N_Ucl_Hits_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/N_Ucl_Hits_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/N_Ucl_Hits_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    TH1D *h_n_VclHits1 = (TH1D*) file_in->Get("h_n_VclHits1");
    h_n_VclHits1->SetLineWidth(2);
    h_n_VclHits1->SetTitle("; The # of Hits in V clusters");
    bin_2 = h_n_VclHits1->FindBin(2);
    h_n_VclHits1->SetMaximum(1.05 * h_n_VclHits1->GetBinContent(bin_2));
    h_n_VclHits1->SetAxisRange(2, h_n_VclHits1->GetBinCenter(h_n_VclHits1->GetNbinsX()));
    h_n_VclHits1->GetXaxis()->SetNdivisions(716);
    h_n_VclHits1->Draw();
    c1->Print(Form("Figs/N_Vcl_Hits_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/N_Vcl_Hits_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/N_Vcl_Hits_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    TF1 *f_Landau = new TF1("f_Landau", "[0]*TMath::Landau(x, [1], [2])", 0., 1000.);
    f_Landau->SetNpx(4500);

    double adcMaxAxis = 500.;
    TH1D *h_U_PeakADC_MultiCl1 = (TH1D*) file_in->Get("h_U_PeakADC_MultiCl1");
    h_U_PeakADC_MultiCl1->SetLineWidth(2);
    h_U_PeakADC_MultiCl1->SetTitle("; ADC of U cluster peak");
    h_U_PeakADC_MultiCl1->SetAxisRange(0., adcMaxAxis, "X");
    f_Landau->SetParameters(5. * h_U_PeakADC_MultiCl1->GetMaximum(), h_U_PeakADC_MultiCl1->GetBinCenter(h_U_PeakADC_MultiCl1->GetMaximumBin()), 10);
    h_U_PeakADC_MultiCl1->Fit(f_Landau, "MeV", "", 0., adcMaxAxis);
    lat1->DrawLatex(0.5, 0.7, Form("MPV = %1.1f", f_Landau->GetParameter(1)));
    c1->Print(Form("Figs/Peak_ADC_UStrips_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/Peak_ADC_UStrips_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/Peak_ADC_UStrips_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    TH1D *h_V_PeakADC_MultiCl1 = (TH1D*) file_in->Get("h_V_PeakADC_MultiCl1");
    h_V_PeakADC_MultiCl1->SetLineWidth(2);
    h_V_PeakADC_MultiCl1->SetTitle("; ADC of V cluster peak");
    h_V_PeakADC_MultiCl1->SetAxisRange(0., adcMaxAxis, "X");
    f_Landau->SetParameters(5. * h_V_PeakADC_MultiCl1->GetMaximum(), h_V_PeakADC_MultiCl1->GetBinCenter(h_V_PeakADC_MultiCl1->GetMaximumBin()), 2);
    h_V_PeakADC_MultiCl1->Fit(f_Landau, "MeV", "", 0., adcMaxAxis);
    lat1->DrawLatex(0.5, 0.7, Form("MPV = %1.1f", f_Landau->GetParameter(1)));
    c1->Print(Form("Figs/Peak_ADC_VStrips_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/Peak_ADC_VStrips_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/Peak_ADC_VStrips_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    TH1D *h_U_Coord_MultrCl1 = (TH1D*) file_in->Get("h_U_Coord_MultrCl1");
    h_U_Coord_MultrCl1->SetTitle("; Cluster U coordinate [strip]");
    h_U_Coord_MultrCl1->SetLineWidth(2);
    h_U_Coord_MultrCl1->Draw();
    c1->Print(Form("Figs/cl_U_coordinate1_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/cl_U_coordinate1_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/cl_U_coordinate1_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    TH1D *h_V_Coord_MultrCl1 = (TH1D*) file_in->Get("h_V_Coord_MultrCl1");
    h_V_Coord_MultrCl1->SetTitle("; Cluster V coordinate [strip]");
    h_V_Coord_MultrCl1->SetLineWidth(2);
    h_V_Coord_MultrCl1->Draw();
    c1->Print(Form("Figs/cl_V_coordinate1_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/cl_V_coordinate1_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/cl_V_coordinate1_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    TH1D *h_U_Coord_MultrCl2 = (TH1D*) file_in->Get("h_U_Coord_MultrCl2");
    h_U_Coord_MultrCl2->SetTitle("; Cluster U coordinate [strip]");
    h_U_Coord_MultrCl2->SetLineWidth(2);
    h_U_Coord_MultrCl2->Draw();
    c1->Print(Form("Figs/cl_U_coordinate2_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/cl_U_coordinate2_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/cl_U_coordinate2_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    TH1D *h_V_Coord_MultrCl2 = (TH1D*) file_in->Get("h_V_Coord_MultrCl2");
    h_V_Coord_MultrCl2->SetTitle("; Cluster V coordinate [strip]");
    h_V_Coord_MultrCl2->SetLineWidth(2);
    h_V_Coord_MultrCl2->Draw();
    c1->Print(Form("Figs/cl_V_coordinate2_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/cl_V_coordinate2_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/cl_V_coordinate2_%d_t%1.1f_m%d.root", run, threshold, MinClSize));


    c1->SetLogz();
    TH2D *h_n_GEM_Y_vs_X_Hits1 = (TH2D*) file_in->Get("h_n_GEM_Y_vs_X_Hits1");
    h_n_GEM_Y_vs_X_Hits1->SetTitle("; The number of X hits in GEM; The number of Y hits in GEM");
    h_n_GEM_Y_vs_X_Hits1->Draw("colz");
    line1->DrawLine(1.5, 1.5, 1.5, 20.);
    line1->DrawLine(1.5, 1.5, 20., 1.5);
    c1->Print(Form("Figs/GEM_N_Y_vs_X_Hits1_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/GEM_N_Y_vs_X_Hits1_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/GEM_N_Y_vs_X_Hits1_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    TH2D *h_GEM_XY1 = (TH2D*) file_in->Get("h_GEM_XY1");
    h_GEM_XY1->SetTitle("; GEM X Strip; GEM Y Strip");
    h_GEM_XY1->Draw("colz");
    c1->Print(Form("Figs/GEM_Y_vs_X_Strips_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/GEM_Y_vs_X_Strips_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/GEM_Y_vs_X_Strips_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    TH2D *h_n_uRwell_V_vs_U_MultiHitCl = (TH2D*) file_in->Get("h_n_uRwell_V_vs_U_MultiHitCl");
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
    lat1->DrawLatex(0.15, 0.9, Form("Has U cluster = %d #rightarrow %1.2f %%", counts_has_U_Cluster, 100. * double(counts_has_U_Cluster) / counts_integral));
    lat1->DrawLatex(0.15, 0.85, Form("Has V cluster = %d #rightarrow %1.2f %% ", counts_has_V_Cluster, 100. * double(counts_has_V_Cluster) / counts_integral));
    lat1->DrawLatex(0.15, 0.8, Form("Has any cluster = %d #rightarrow %1.2f %% ", counts_has_AnyCluster, 100. * double(counts_has_AnyCluster) / counts_integral));
    lat1->DrawLatex(0.15, 0.75, Form("Has U and V cluster = %d #rightarrow %1.2f %% ", counts_has_U_AND_V_Cluster, 100. * double(counts_has_U_AND_V_Cluster) / counts_integral));

    TCanvas *c2 = new TCanvas("c2", "", 1800., 1000.);
    c2->SetTopMargin(0.02);
    c2->SetRightMargin(0.02);

    TH2D *h_Cross_YXc2 = (TH2D*) file_in->Get("h_Cross_YXc2");
    h_Cross_YXc2->SetStats(0);
    h_Cross_YXc2->SetTitle("; Cross X coordinate [mm]; Cross Y coordinate [mm]");
    h_Cross_YXc2->Draw("col");
    DrawActiveArea();
    c2->Print(Form("Figs/Cross_YXc2_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_YXc2_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_YXc2_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    TF1 *f_Gaus = new TF1("f_Gaus", "[0]*TMath::Gaus(x, [1], [2])", -900., 900.);
    f_Gaus->SetNpx(4500);
    TF1 *f_GPol4 = new TF1("f_GPol4", "[0]*TMath::Gaus(x, [1], [2]) + [3] + x*( [4] + x*( [5] + x*( [6] + x*[7] ) )  )", -900., 900.);
    f_GPol4->SetNpx(4500);
    TF1 *f_Pol4 = new TF1("f_Pol4", "[0] + x*( [1] + x*( [2] + x*( [3] + x*[4] ) )  )", -900., 900.);
    f_Pol4->SetNpx(4500);
    f_Pol4->SetLineColor(4);

    TH1D *h_Cross_X2 = (TH1D*) h_Cross_YXc2->ProjectionX("h_Cross_X2", 1, h_Cross_YXc2->GetNbinsY());
    h_Cross_X2->Draw();
    f_GPol4->SetParameters(h_Cross_X2->GetMaximum(), h_Cross_X2->GetBinCenter(h_Cross_X2->GetMaximumBin()), 75.);
    h_Cross_X2->Fit(f_GPol4, "MeV", "", -800, 800.);
    double pars[8];
    f_GPol4->GetParameters(pars);
    f_Gaus->SetParameters(pars);
    f_Pol4->SetParameters(&pars[3]);
    f_Pol4->Draw("Same");
    double N_Cross_BgrSubtr = f_Gaus->Integral(-800., 800.) / h_Cross_X2->GetBinWidth(15);
    c2->Print(Form("Figs/Cross_X2_Fit_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_X2_Fit_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_X2_Fit_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    c1->cd();
    lat1->DrawLatex(0.15, 0.7, Form("Has Cross BgrSubtr = %d #rightarrow %1.2f %% ", int(N_Cross_BgrSubtr), 100. * double(N_Cross_BgrSubtr) / counts_integral));
    c1->Print(Form("Figs/Number_OF_V_vs_U_clusters_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c1->Print(Form("Figs/Number_OF_V_vs_U_clusters_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c1->Print(Form("Figs/Number_OF_V_vs_U_clusters_%d_t%1.1f_m%d.root", run, threshold, MinClSize));



    c2->cd();
    TH2D *h_Cross_YXc1 = (TH2D*) file_in->Get("h_Cross_YXc1");
    h_Cross_YXc1->SetStats(0);
    h_Cross_YXc1->SetTitle("; Cross X coordinate [mm]; Cross Y coordinate [mm]");
    h_Cross_YXc1->SetTitleSize(0.05, "Y");
    h_Cross_YXc1->SetLabelSize(0.05, "Y");
    h_Cross_YXc1->SetTitleOffset(0.9, "Y");
    h_Cross_YXc1->SetTitleSize(0.05, "X");
    h_Cross_YXc1->SetLabelSize(0.05, "X");
    h_Cross_YXc1->SetMaximum(h_Cross_YXc1->GetMaximum() / 5.);
    h_Cross_YXc1->SetMinimum(3);
    h_Cross_YXc1->Draw("col");
    DrawActiveArea();
    c2->Print(Form("Figs/Cross_YXc1_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_YXc1_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_YXc1_%d_t%1.1f_m%d.root", run, threshold, MinClSize));


    TH2D *h_Cross_YXc3 = (TH2D*) file_in->Get("h_Cross_YXc3");
    h_Cross_YXc3->SetStats(0);
    h_Cross_YXc3->SetTitle("; Cross X coordinate [mm]; Cross Y coordinate [mm]");
    h_Cross_YXc3->SetTitleSize(0.05, "Y");
    h_Cross_YXc3->SetLabelSize(0.05, "Y");
    h_Cross_YXc3->SetTitleOffset(0.9, "Y");
    h_Cross_YXc3->SetTitleSize(0.05, "X");
    h_Cross_YXc3->SetLabelSize(0.05, "X");
    h_Cross_YXc3->Draw("col");
    h_Cross_YXc3->SetMaximum(h_Cross_YXc3->GetMaximum() / 5.);
    uRwellTools::DrawGroupStripBiundaries();
    c2->Print(Form("Figs/Cross_YXc3_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_YXc3_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_YXc3_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    int Y_bin2 = h_Cross_YXc3->GetYaxis()->FindBin(245);
    int Y_bin1 = h_Cross_YXc3->GetYaxis()->FindBin(-245);
    TH1D *h_Cross_X3 = (TH1D*) h_Cross_YXc3->ProjectionX("h_Cross_X3", Y_bin1, Y_bin2);
    h_Cross_X3->SetLineColor(95);
    h_Cross_X3->SetFillColor(95);
    h_Cross_X3->Draw();

    c2->Print(Form("Figs/Cross_X3_%d_%1.1f_%d.pdf", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_X3_%d_%1.1f_%d.png", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_X3_%d_%1.1f_%d.root", run, threshold, MinClSize));

    h_Cross_X3->SetAxisRange(-90., 90.);
    h_Cross_X3->SetMinimum(0.);
    h_Cross_X3->SetStats(0);
    h_Cross_X3->SetTitleSize(0.05, "X");
    h_Cross_X3->SetLabelSize(0.05, "X");
    h_Cross_X3->SetLabelSize(0.05, "Y");
    c2->SetGridy();
    c2->SetGridx();
    c2->Modified();
    c2->Update();
    c2->Print(Form("Figs/Cross_X3_%d_%1.1f_%d_ZoomedOnX.pdf", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_X3_%d_%1.1f_%d_ZoomedOnX.png", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_X3_%d_%1.1f_%d_ZoomedOnX.root", run, threshold, MinClSize));


    TH2D *h_Cross_YXc_Weighted3 = (TH2D*) file_in->Get("h_Cross_YXc_Weighted3");
    h_Cross_YXc_Weighted3->SetStats(0);
    h_Cross_YXc_Weighted3->SetTitle("; Cross X coordinate [mm]; Cross Y coordinate [mm]");
    h_Cross_YXc_Weighted3->SetTitleSize(0.05, "Y");
    h_Cross_YXc_Weighted3->SetLabelSize(0.05, "Y");
    h_Cross_YXc_Weighted3->SetTitleOffset(0.9, "Y");
    h_Cross_YXc_Weighted3->SetTitleSize(0.05, "X");
    h_Cross_YXc_Weighted3->SetLabelSize(0.05, "X");
    h_Cross_YXc_Weighted3->Draw("col");
    h_Cross_YXc_Weighted3->SetMaximum(h_Cross_YXc_Weighted3->GetMaximum() / 5.);
    uRwellTools::DrawGroupStripBiundaries();
    c2->Print(Form("Figs/Cross_YXc3_Weighted_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_YXc3_Weighted_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_YXc3_Weighted_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

    Y_bin2 = h_Cross_YXc_Weighted3->GetYaxis()->FindBin(245);
    Y_bin1 = h_Cross_YXc_Weighted3->GetYaxis()->FindBin(-245);

    
    TH1D *h_Cross_X_Weighted3 = (TH1D*) h_Cross_YXc_Weighted3->ProjectionX("h_Cross_X_Weighted3", Y_bin1, Y_bin2);
    h_Cross_X_Weighted3->SetLineColor(4);
    h_Cross_X_Weighted3->SetFillColor(4);
    int bin2Nomalize = h_Cross_X_Weighted3->FindBin(10.);
    h_Cross_X3->Scale(h_Cross_X_Weighted3->GetBinContent(bin2Nomalize) / h_Cross_X3->GetBinContent(bin2Nomalize));
    h_Cross_X3->GetXaxis()->UnZoom();

    TLegend *leg2 = new TLegend(0.15, 0.7, 0.3, 0.97);
    leg2->SetBorderSize(0);
    leg2->AddEntry(h_Cross_X_Weighted3, "Weighted");
    leg2->AddEntry(h_Cross_X3, "Not Weighted");
    h_Cross_X_Weighted3->SetStats(0);
    h_Cross_X_Weighted3->Draw("hist");
    h_Cross_X3->Draw("hist same");
    leg2->Draw();
    c2->Print(Form("Figs/Cross_X3_Weighted_Unweighted_%d_t%1.1f_m%d.pdf", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_X3_Weighted_Unweighted_%d_t%1.1f_m%d.png", run, threshold, MinClSize));
    c2->Print(Form("Figs/Cross_X3_Weighted_Unweighted_%d_t%1.1f_m%d.root", run, threshold, MinClSize));

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
