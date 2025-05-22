#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "hipo4/reader.h"
#include "clas12reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "func.h"

int Rad_analysis(string nameFile = "Rad_lepton", int version = -19, int model = 6, int lepton_ID = 11, bool shift = true)
{
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();

    gStyle->SetPaintTextFormat("4.1f");
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetLabelSize(.03, "xyz");
    gStyle->SetTitleSize(.03, "xyz");
    gStyle->SetTitleSize(.07, "t");
    gStyle->SetFrameLineWidth(1);
    gStyle->SetLineWidth(1);
    gStyle->SetHistLineWidth(1);
    gStyle->SetMarkerStyle(13);
    gStyle->SetTitleW(0.8); // per cent of the pad width
    gStyle->SetTitleH(0.1);

    //********************************
    // DECLARATION OF VARIABLES
    //********************************
    double_t lepton_theta, lepton_phi, lepton_p, lepton_sfpcal, lepton_sfecout, lepton_sfecin;
    double_t lepton_m2pcal, lepton_m2ecin, lepton_m2ecout, lepton_nphe, lepton_gamma_angle, lepton_gamma_angle_phi, lepton_score, lepton_score_upd;
    int lepton_photon;

    TString root_file;
    
    root_file = "/work/clas12/mtenorio/LeptonID/6_var/" + nameFile + ".root";

    TFile *file = new TFile(root_file, "READ");
    TTree *tree = (TTree *)file->Get("analysis");

    tree->SetMakeClass(1);

    tree->SetBranchAddress("lepton_theta", &lepton_theta);
    tree->SetBranchAddress("lepton_phi", &lepton_phi);
    tree->SetBranchAddress("lepton_p", &lepton_p);
    tree->SetBranchAddress("lepton_sfpcal", &lepton_sfpcal);
    tree->SetBranchAddress("lepton_sfecin", &lepton_sfecin);
    tree->SetBranchAddress("lepton_sfecout", &lepton_sfecout);
    tree->SetBranchAddress("lepton_m2pcal", &lepton_m2pcal);
    tree->SetBranchAddress("lepton_m2ecin", &lepton_m2ecin);
    tree->SetBranchAddress("lepton_m2ecout", &lepton_m2ecout);
    tree->SetBranchAddress("lepton_nphe", &lepton_nphe);
    tree->SetBranchAddress("lepton_gamma_angle", &lepton_gamma_angle);
    // tree->SetBranchAddress("lepton_score",&lepton_score);
    tree->SetBranchAddress("lepton_photon", &lepton_photon);
    tree->SetBranchAddress("lepton_gamma_angle_phi", &lepton_gamma_angle_phi);
    // if(model==6&&Lepton_ID==-11)
    // tree->SetBranchAddress("lepton_score_upd",&lepton_score_upd);

    // Create a set of variables and declare them to the reader
    Float_t P, Theta, Phi, PCAL, ECIN, ECOUT;
    Float_t m2PCAL = -1;
    Float_t m2ECIN = -1;
    Float_t m2ECOUT = -1;
    Float_t Nphe;

    //-********************************************
    //:::::::::::::::::::TMVA::::::::::::::::::::::
    //-********************************************
    TMVA::Reader *readerTMVA = new TMVA::Reader("!Color:Silent");
    Double_t score;
    Double_t score6;
    Double_t score9;

    if (model == 9)
    {
        readerTMVA->AddVariable("P", &P);
        readerTMVA->AddVariable("Theta", &Theta);
        readerTMVA->AddVariable("Phi", &Phi);
    }
    // readerTMVA->AddVariable( "Nphe",&Nphe);
    readerTMVA->AddVariable("SFPCAL", &PCAL);
    readerTMVA->AddVariable("SFECIN", &ECIN);
    readerTMVA->AddVariable("SFECOUT", &ECOUT);
    readerTMVA->AddVariable("m2PCAL", &m2PCAL);
    readerTMVA->AddVariable("m2ECIN", &m2ECIN);
    readerTMVA->AddVariable("m2ECOUT", &m2ECOUT);

    // Book Methods
    TString weightfile;
    TString dir_weights = "/work/clas12/mtenorio/ML_weights_pass2/";
    TString name_model;
    if (model == 6)
        name_model = "/TMVAClassification_BDT_6.weights.xml";
    if (model == 9)
        name_model = "/TMVAClassification_BDT.weights.xml";

    TString name_configuration;

    int from, to;
    if (version == -19)
    {
        if (lepton_ID == 11)
            name_configuration = "S19neg";
        else if (lepton_ID == -11)
            name_configuration = "S19pos";
    }
    if (version == -18)
    {
        if (lepton_ID == 11)
            name_configuration = "F18inneg";
        else if (lepton_ID == -11)
            name_configuration = "F18inpos";
    }
    if (version == +18)
    {
        if (lepton_ID == 11)
            name_configuration = "F18outneg";
        else if (lepton_ID == -11)
            name_configuration = "F18outpos";
    }
    weightfile = dir_weights + name_configuration + name_model;
    readerTMVA->BookMVA("BDT method", weightfile);

    float_t cuts[15] = {-0.6, -0.5, -0.4, -0.3, -0.2, -0.15, -0.1, -0.06, -0.04, -0.02, 0.0, 0.05, 0.1, 0.2, 0.3};
    int k_events[10][15] = {0};

    //********************************
    // HISTOGRAMAS
    //********************************
    TH1F *h_e_gamma[10];
    TH1F *h_e_gamma_cut[10][15];

    for (int i = 0; i < 10; i++)
    {
        ostringstream name_ss;
        name_ss << "h_e_gamma_" << i;
        ostringstream sstr_ss;
        sstr_ss << "e -> e#gamma " << i << " ;#Delta #Theta, degrees; ";
        h_e_gamma[i] = new TH1F(name_ss.str().c_str(), sstr_ss.str().c_str(), 400, -2, 2);

        for (int j = 0; j < 15; j++)
        {
            ostringstream name;
            name << "Cut:" << i << "." << j;
            ostringstream sstr;
            sstr << "Cut " << cuts[j] << " - Section " << i << " ;#Delta #Theta, degrees; ";
            h_e_gamma_cut[i][j] = new TH1F(name.str().c_str(), sstr.str().c_str(), 400, -2, 2);
            ;
            name.str("");
            sstr.str("");
        }
        name_ss.str("");
        sstr_ss.str("");
    }

    // Start
    int counter = 0;
    for (int fc = 0; fc < tree->GetEntries(); fc++)
    { // Run 5032 to 5038 // 6616 6618. //5423 5424
        tree->GetEntry(fc);
        if (lepton_photon != 1)
            continue;

        if (lepton_ID == -11 && (lepton_gamma_angle_phi > 8 || lepton_gamma_angle_phi < 2))
            continue;

        if (lepton_nphe <= 3)
            continue;

        Double_t SF_corr, m2_corr;
        if (lepton_ID == -11)
        {
            if (version == -19)
            {
                SF_corr = 0.01;
                m2_corr = 0.8;
            }
            else if (version == 18)
            {
                SF_corr = 0.03;
                m2_corr = 1.0;
            }
            else
            {
                SF_corr = 0.01;
                m2_corr = 0.8;
            }
        }
        else
        {
            if (version == -19)
            {
                SF_corr = 0.03;
                m2_corr = 0.8;
            }
            else if (version == 18)
            {
                SF_corr = 0.05;
                m2_corr = 1.1;
            }
            else
            {
                SF_corr = 0.02;
                m2_corr = 0.8;
            }
        }

        if (!shift)
        {
            SF_corr = 0.0;
            m2_corr = 1.0;
        }

        P = lepton_p;
        Theta = lepton_theta;
        Phi = lepton_phi;
        Nphe = lepton_nphe;
        PCAL = lepton_sfpcal;
        ECIN = lepton_sfecin + SF_corr;
        ECOUT = lepton_sfecout;
        m2PCAL = lepton_m2pcal;
        m2ECIN = lepton_m2ecin * m2_corr;
        m2ECOUT = lepton_m2ecout;

        lepton_score = readerTMVA->EvaluateMVA("BDT method");

        counter++;
        int events = 7950;

        if (version == -19)
            events = 18000;
        else if (version == -18)
            events = 65000;

        int sect;

        if (counter < events)
            sect = 0;
        else
            sect = 1;
        /*else if(counter<events*3)
            sect=2;
        else if(counter<events*4)
            sect=3;
        else if(counter<events*5)
            sect=4;
        else if(counter<events*6)
            sect=5;
        else if(counter<events*7)
            sect=6;
        else if(counter<events*8)
            sect=7;
        else if(counter<events*9)
            sect=8;
        else
            sect=9;*/

        h_e_gamma[sect]->Fill(lepton_gamma_angle);
        for (int j = 0; j < 15; j++)
        {
            if (lepton_score > cuts[j])
            {
                h_e_gamma_cut[sect][j]->Fill(lepton_gamma_angle);
                if (abs(lepton_gamma_angle) < 0.5)
                    k_events[sect][j]++;
            }
        }
    } // End for "Runs"

    stringstream ss;
    ss << model;
    string mod = ss.str();

    string pdf_original = nameFile + "_" + mod + "TEMP.pdf";
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0001);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.15);

    double_t n[10];
    double_t k;
    double_t acc, err;
    float_t acceptance[10][15];
    float_t erry_Bayessian[10][15];
    float_t erry_Classic[10][15];
    float_t erry_Poisson[10][15];
    float_t errx[15] = {0};

    double_t error_n;

    float min, max;
    if (version == 18)
    {
        min = -1.6;
        max = 1.6;
    }
    else
    {
        min = -1.5;
        max = 2;
    }

    float_t neve, mean, sigma, a, b, c;
    TCanvas *can2 = new TCanvas("can2", "canvas", 200, 10, 700, 500);

    int max_sec;
    if (version == 18)
        max_sec = 2;
    else if (version == -19)
        max_sec = 2;
    else
        max_sec = 2;

    for (int s = 0; s < max_sec; s++)
    {
        cout << "---------SECTION---------: " << s + 1;
        if (version == -18)
        {
            neve = 2E4;
            mean = 0.03;
            sigma = 0.2;
        }
        else
        {
            neve = 5E3;
            mean = 0.04;
            sigma = 0.1;
        }
        a = 7;
        b = 3;
        c = 7;

        can2->Divide(4, 4);

        can2->cd(1);
        h_e_gamma[s]->Draw();
        /*TF1 *g = g43_free(h_e_gamma[s], neve, mean, sigma, a, b, c, min, max);

        n[s] = g->GetParameter(0);
        error_n = g->GetParError(0);
        mean = g->GetParameter(1);
        sigma = g->GetParameter(2);
        a = g->GetParameter(3);
        b = g->GetParameter(4);
        c = g->GetParameter(5);

        k = n[s];
        */

        for (int i = 0; i < 15; i++)
        {
            can2->cd(i + 2);
            cout << "---------FITTING CUT: " << cuts[i];
            h_e_gamma_cut[s][i]->Draw();
/*
            TF1 *g = g43(h_e_gamma_cut[s][i], k, mean, sigma, a, b, c, min, max); // h_e_gamma->GetMaximum()
            k = g->GetParameter(0);
            // if(i==0)
            //   n=k;

            if (h_e_gamma_cut[s][i]->GetEntries() == 0 || k_events[s][i] == 0)
            {
                acceptance[s][i] = 0.0;
            }
            else
            {
                acceptance[s][i] = k / n[s];
            }

            // Calculate error
            err = sqrt((((k + 1) * (k + 2)) / ((n[s] + 2) * (n[s] + 3))) - (pow(k + 1, 2) / pow(n[s] + 2, 2)));
            erry_Classic[s][i] = acceptance[s][i] * sqrt(pow(g->GetParError(0) / k, 2) + pow(error_n / n[s], 2));
            erry_Poisson[s][i] = acceptance[s][i] * sqrt((1 / k) + (1 / n[s]));
            erry_Bayessian[s][i] = err;

            // We use the background parameters for the next iteration
            a = g->GetParameter(3);
            b = g->GetParameter(4);
            c = g->GetParameter(5);
            */
        }
        // if(s<9)
        can2->Print((pdf_original + "(").c_str());
        // else
        //   can2->Print( (pdf_original+")").c_str());

        can2->Clear();
    }
    can2->Print((pdf_original + ")").c_str());

    ////-----------Print and average results---------------
    double_t average_k[15] = {0};
    double_t average_error[15] = {0};

    int ave_sec;
    // if(version!=18)

    // else
    //   ave_sec=max_sec;

    for (int j = 0; j < 15; j++)
    {
        ave_sec = max_sec;
        for (int i = 0; i < max_sec; i++)
        {
            average_k[j] = average_k[j] + (acceptance[i][j]);
            if (isnan(erry_Bayessian[i][j]))
                ave_sec = ave_sec - 1;
            else
                average_error[j] = average_error[j] + erry_Bayessian[i][j];
        }
        average_k[j] = average_k[j] / max_sec;
        average_error[j] = average_error[j] / ave_sec;
    }

   TString Sver;
    if (version == -19)
        Sver = "S19";
    else if (version == -18)
        Sver = "F18in";
    else
        Sver = "F18out";

    TString Slep;
    if (lepton_ID == -11)
        Slep = "pos";
    else
        Slep = "neg";

    FILE *f_readme;

    f_readme = fopen("/lustre24/expphy/volatile/clas12/mtenorio/Electron_AI_PID/Validation/DATA.txt", "a");
    fprintf(f_readme, "\n//Model %d Version: %s%s\n", model, Sver.Data(), Slep.Data());
    for (int i = 0; i < max_sec; i++)
    {
        fprintf(f_readme, "n[%d]=%f\n", i, n[i]);
         fprintf(f_readme, "Double_t x_%dbdt_%s%s[15]={", model, Sver.Data(), Slep.Data());
        for (int j = 0; j < 15; j++)
            fprintf(f_readme, "%f, ", acceptance[i][j]);
        fprintf(f_readme, "} \n");
       fprintf(f_readme, "Double_t xerrB_%dbdt_%s%s[15]={", model, Sver.Data(), Slep.Data());
        for (int j = 0; j < 15; j++)
        {
            if (isnan(erry_Bayessian[i][j]))
                fprintf(f_readme, "---nan--, ");
            else
                fprintf(f_readme, "%f, ", erry_Bayessian[i][j]);
        }
        fprintf(f_readme, "} \n");
    }
     fprintf(f_readme, "AVERAGE\n");
    fprintf(f_readme, "Double_t x_%dbdt_%s%s[15]={", model, Sver.Data(), Slep.Data());
    for (int i = 0; i < 15; i++)
        fprintf(f_readme, "%f, ", average_k[i]);
    fprintf(f_readme, "} \n");
    fprintf(f_readme, "Double_t xerrB_%dbdt_%s%s[15]={", model, Sver.Data(), Slep.Data());
    for (int i = 0; i < 15; i++)
        fprintf(f_readme, "%f, ", average_error[i]);
    fprintf(f_readme, "} \n");

    float_t Eff[15] = {0};
    float_t Eff_err[15] = {0};
    for (int i = 0; i < 15; i++)
    {
        Eff[i] = average_k[i];
        Eff_err[i] = average_error[i];
    }

    TCanvas *can = new TCanvas("can", "canvas", 200, 10, 700, 700);
    TGraphErrors *gr_Acceptance = new TGraphErrors(15, cuts, Eff, errx, Eff_err);
    gr_Acceptance->SetMarkerStyle(8);
    gr_Acceptance->GetYaxis()->SetRangeUser(0, 1.3);
    gr_Acceptance->SetTitle(" ; Cut value ;Efficiency ");
    gr_Acceptance->Draw("ALEP");

    


    /*auto legend = new TLegend(0.1, 0.50, 0.87, 0.67);
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    legend->SetTextSize(0.02);
    legend->AddEntry(gr_Acceptance, Form("Efficiency at -0.01 cut =%.5f",Eff[5]), "p");
    legend->Draw("same");*/
    //can->Print((pdf_original + ")").c_str());

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
    cout << "Counter: " << counter << endl;

    return 0;
}
