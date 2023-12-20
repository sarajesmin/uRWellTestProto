/* 
 * File:   DrawHVDependecePlots.cc
 * Author: rafopar
 *
 * Created on December 1, 2023, 8:59 PM
 */

#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <cxxopts.hpp>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>

#include <uRwellTools.h>

using namespace std;

/*
 * 
 */
int main(int argc, char **argv) {

    cxxopts::Options options("DrawHVDependencePlots", "Draws plots of efficiencies (or others staff) as a function of HV");

    options.add_options()
            ("s,Series", "HV Scan series", cxxopts::value<int>())
            ("t,Threshold", "Hit Threshold in terms of sigma", cxxopts::value<double>())
            ("m,MinHits", "Number of minimum hits in the cluster", cxxopts::value<int>())
            ;

    auto parsed_options = options.parse(argc, argv);
    
    if (!parsed_options.count("Series")) {        
        cout << "The Series of runs is nor provided. Exiting..." << endl;
        exit(1);
    }
    const int series = parsed_options["Series"].as<int>();


    if (!parsed_options.count("Threshold")) {
        cout << "* You didn't provide the hit threshold. Exiting" << endl;
        exit(1);
    }
    const double threshold = parsed_options["Threshold"].as<double>();

    if (!parsed_options.count("MinHits")) {
        cout << "* You didn't provide the Minimum hits int the cluster. Exiting " << endl;
        exit(1);
    }
    const int MinHits = parsed_options["MinHits"].as<int>();

    cout << "The hit threshold is " << threshold << "\\sigma" << endl;
    cout << "The Minimum cluster size is " << MinHits << "hits" << endl;


    std::map<int, std::vector<int> > mv_runs;
    mv_runs[1] = {1731, 1745, 1750, 1753, 1790, 1761, 1779};
    mv_runs[2] = {1824, 1826, 1828, 1833, 1835, 1837};
    std::map<int, double> m_MESH_HV; // The key is the run number, the value is the MESH_HV
    std::map<int, double> m_Cathode_HV; // The key is the run number, the value is the Cathode_HV

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.02);
    c1->SetGridy();
    c1->SetGridx();

    TGraphAsymmErrors *gr_Eff_U = new TGraphAsymmErrors();
    gr_Eff_U->SetMarkerStyle(20);
    gr_Eff_U->SetMarkerColor(2);
    gr_Eff_U->SetMarkerSize(2);
    TGraphAsymmErrors *gr_Eff_V = new TGraphAsymmErrors();
    gr_Eff_V->SetMarkerStyle(21);
    gr_Eff_V->SetMarkerColor(4);
    gr_Eff_V->SetMarkerSize(2);
    TGraphAsymmErrors *gr_Eff_OR = new TGraphAsymmErrors();
    gr_Eff_OR->SetMarkerStyle(22);
    gr_Eff_OR->SetMarkerColor(6);
    gr_Eff_OR->SetMarkerSize(2);
    TGraphAsymmErrors *gr_Eff_AND = new TGraphAsymmErrors();
    gr_Eff_AND->SetMarkerStyle(23);
    gr_Eff_AND->SetMarkerColor(7);
    gr_Eff_AND->SetMarkerSize(2);
    TGraphAsymmErrors *gr_Eff_CrsBgrSubtr = new TGraphAsymmErrors();
    gr_Eff_CrsBgrSubtr->SetMarkerStyle(24);
    gr_Eff_CrsBgrSubtr->SetMarkerColor(8);
    gr_Eff_CrsBgrSubtr->SetMarkerSize(2);


    std::string hvTablefileName = Form("HV_Table_%d.dat", series);
    ifstream inp_HVTable(hvTablefileName.c_str());

    if (inp_HVTable.is_open()) {


        while (!inp_HVTable.eof()) {
            int run, HV_MESH, HV_Cathode, HV_GEM, HV_Drift;

            inp_HVTable >> run >> HV_MESH >> HV_Cathode>>HV_GEM;

            cout << run << "  " << HV_MESH << "   " << HV_Cathode << "   " << HV_GEM << endl;
            m_MESH_HV[run] = HV_MESH;
            m_Cathode_HV[run] = HV_Cathode;
        }

    } else {
        cout << "Can not opern the file" << hvTablefileName.c_str() << endl;
    }

    for (int i = 0; i < mv_runs[series].size(); i++) {

        int run = mv_runs[series].at(i);


        double eff_U, eff_V, eff_OR, eff_AND, eff_BgrSubtr;
        uRwellTools::uRwellEff eff;

        TFile *file_in = new TFile(Form("AnaClustering_%d_Thr_%1.1f_MinHits_%d.root", run, threshold, MinHits));

        TH2D *h_n_uRwell_V_vs_U_MultiHitCl = (TH2D*) file_in->Get("h_n_uRwell_V_vs_U_MultiHitCl");
        double All_GoodEventts = h_n_uRwell_V_vs_U_MultiHitCl->Integral();

        uRwellTools::CalcEfficiencies(h_n_uRwell_V_vs_U_MultiHitCl, eff);

        TH2D *h_Cross_YXc2 = (TH2D*) file_in->Get("h_Cross_YXc2");

        double n_CrossBgrSubtracted = uRwellTools::getNofBgrSbtrCrosses(h_Cross_YXc2);

        eff.eff_crs_BgrSubtr = 100. * n_CrossBgrSubtracted / All_GoodEventts;

        gr_Eff_U->SetPoint(i, m_MESH_HV[run], eff.eff_U);
        gr_Eff_U->SetPointError(i, 0., 0., eff.errLow_eff_U, eff.errUp_eff_U);
        gr_Eff_V->SetPoint(i, m_MESH_HV[run], eff.eff_V);
        gr_Eff_U->SetPointError(i, 0., 0., eff.errLow_eff_V, eff.errUp_eff_V);
        gr_Eff_OR->SetPoint(i, m_MESH_HV[run], eff.eff_OR);
        gr_Eff_U->SetPointError(i, 0., 0., eff.errLow_eff_OR, eff.errUp_eff_OR);
        gr_Eff_AND->SetPoint(i, m_MESH_HV[run], eff.eff_AND);
        gr_Eff_U->SetPointError(i, 0., 0., eff.errLow_eff_AND, eff.errUp_eff_AND);

        eff.errUp_eff_crs_BgrSubtr = 100 * TEfficiency::ClopperPearson(All_GoodEventts, n_CrossBgrSubtracted, uRwellTools::OneSigma, true) - eff.eff_crs_BgrSubtr;
        eff.errLow_eff_crs_BgrSubtr = eff.eff_crs_BgrSubtr - 100 * TEfficiency::ClopperPearson(All_GoodEventts, n_CrossBgrSubtracted, uRwellTools::OneSigma, false);
        gr_Eff_CrsBgrSubtr->SetPoint(i, m_MESH_HV[run], eff.eff_crs_BgrSubtr);
        gr_Eff_CrsBgrSubtr->SetPointError(i, 0., 0., eff.errLow_eff_crs_BgrSubtr, eff.errUp_eff_crs_BgrSubtr);

        //cout<<eff.errLow_eff_U<<"  "<<eff.errUp_eff_U<<endl;

        delete h_Cross_YXc2;
        delete h_n_uRwell_V_vs_U_MultiHitCl;
        delete file_in;
    }

    TLegend *leg1 = new TLegend(0.12, 0.75, 0.4, 0.97);
    leg1->SetBorderSize(0);
    leg1->AddEntry(gr_Eff_U, "U Cluster");
    leg1->AddEntry(gr_Eff_V, "V Cluster");
    leg1->AddEntry(gr_Eff_OR, "Any Cluster");
    leg1->AddEntry(gr_Eff_AND, "U and V Cluster");
    leg1->AddEntry(gr_Eff_CrsBgrSubtr, "Bgr Subtracted");

    TMultiGraph *mtgr1 = new TMultiGraph();
    mtgr1->Add(gr_Eff_U);
    mtgr1->Add(gr_Eff_V);
    mtgr1->Add(gr_Eff_OR);
    mtgr1->Add(gr_Eff_AND);
    mtgr1->Add(gr_Eff_CrsBgrSubtr);

    mtgr1->Draw("APL");
    mtgr1->SetTitle("; MESH HV [V]; Efficiency");
    mtgr1->SetMaximum(100.);
    mtgr1->SetMinimum(0.);
    leg1->Draw();
    c1->Print(Form("Figs/HV_Eff_Dependence_Thr_%1.1f_MinHits_%d_Series_%d.pdf", threshold, MinHits, series));
    c1->Print(Form("Figs/HV_Eff_Dependence_Thr_%1.1f_MinHits_%d_Series_%d.png", threshold, MinHits, series));
    c1->Print(Form("Figs/HV_Eff_Dependence_Thr_%1.1f_MinHits_%d_Series_%d.root", threshold, MinHits, series));

    return 0;
}