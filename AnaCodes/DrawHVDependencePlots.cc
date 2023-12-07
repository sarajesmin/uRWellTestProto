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

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMultiGraph.h>

#include <uRwellTools.h>

using namespace std;

/*
 * 
 */
int main() {

    std::vector<int> v_runs = {1731, 1745, 1750, 1753, 1789, 1761, 1779};
    std::map<int, double> m_MESH_HV; // The key is the run number, the value is the MESH_HV
    std::map<int, double> m_Cathode_HV; // The key is the run number, the value is the Cathode_HV

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.02);
    c1->SetGridy();
    c1->SetGridx();

    TGraph *gr_Eff_U = new TGraph();
    gr_Eff_U->SetMarkerStyle(20);
    gr_Eff_U->SetMarkerColor(2);
    gr_Eff_U->SetMarkerSize(2);
    TGraph *gr_Eff_V = new TGraph();
    gr_Eff_V->SetMarkerStyle(21);
    gr_Eff_V->SetMarkerColor(4);
    gr_Eff_V->SetMarkerSize(2);
    TGraph *gr_Eff_OR = new TGraph();
    gr_Eff_OR->SetMarkerStyle(22);
    gr_Eff_OR->SetMarkerColor(6);
    gr_Eff_OR->SetMarkerSize(2);
    TGraph *gr_Eff_AND = new TGraph();
    gr_Eff_AND->SetMarkerStyle(23);
    gr_Eff_AND->SetMarkerColor(7);
    gr_Eff_AND->SetMarkerSize(2);
    TGraph *gr_Eff_CrsBgrSubtr = new TGraph();
    gr_Eff_CrsBgrSubtr->SetMarkerStyle(24);
    gr_Eff_CrsBgrSubtr->SetMarkerColor(8);
    gr_Eff_CrsBgrSubtr->SetMarkerSize(2);

    int series = 1; //just to identify series of runs for HV studies.

    const double threshold = 4;
    const int MinHits = 1;

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

    for (int i = 0; i < v_runs.size(); i++) {

        int run = v_runs.at(i);


        double eff_U, eff_V, eff_OR, eff_AND, eff_BgrSubtr;
        TFile *file_in = new TFile(Form("AnaClustering_%d_Thr_%1.1f_MinHits_%d.root", run, threshold, MinHits));

        TH2D *h_n_uRwell_V_vs_U_MultiHitCl = (TH2D*) file_in->Get("h_n_uRwell_V_vs_U_MultiHitCl");
        double All_GoodEventts = h_n_uRwell_V_vs_U_MultiHitCl->Integral();

        uRwellTools::CalcEfficiencies(h_n_uRwell_V_vs_U_MultiHitCl, eff_U, eff_V, eff_OR, eff_AND);

        TH2D *h_Cross_YXc2 = (TH2D*) file_in->Get("h_Cross_YXc2");

        double n_CrossBgrSubtracted = uRwellTools::getNofBgrSbtrCrosses(h_Cross_YXc2);

        eff_BgrSubtr = 100. * n_CrossBgrSubtracted / All_GoodEventts;

        gr_Eff_U->SetPoint(i, m_MESH_HV[run], eff_U);
        gr_Eff_V->SetPoint(i, m_MESH_HV[run], eff_V);
        gr_Eff_OR->SetPoint(i, m_MESH_HV[run], eff_OR);
        gr_Eff_AND->SetPoint(i, m_MESH_HV[run], eff_AND);
        gr_Eff_CrsBgrSubtr->SetPoint(i, m_MESH_HV[run], eff_BgrSubtr);

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