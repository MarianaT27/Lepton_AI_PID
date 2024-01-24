#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "THStack.h"
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;



void setupCanvas(TCanvas* canvas){
    
    canvas->SetFillColor(0);
    canvas->SetBorderMode(0);
    canvas->SetBorderSize(0);
    canvas->SetFrameFillColor(0);
    canvas->SetFrameBorderMode(0);
    //canvas->SetGrid();
    
}

void Variable_Plots_TMVA(TString infile,TString name_model,Float_t score_cut, TString outfile)
{
    gROOT->SetBatch(kTRUE);
    
    gStyle->SetPaintTextFormat("4.1f");
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(.03, "xyz");
    gStyle->SetTitleSize(.03, "xyz");
    gStyle->SetTitleSize(.07, "t");
    gStyle->SetFrameLineWidth(1);
    gStyle->SetLineWidth(1);
    gStyle->SetHistLineWidth(1);
    gStyle->SetMarkerStyle(13);
    gStyle->SetTitleW(0.8); // per cent of the pad width
    gStyle->SetTitleH(0.1); // per cent of the pad height
    
    
    std::vector<vector<TString>> labels1D{
        {"P", "P", "custom_range", "0.", "14.",
            "custom_binning", "100", "", "legend", "P"},
        
        {"Theta", "#Theta", "custom_range", "0.", "1",
            "custom_binning", "100", "", "legend", "electron_SF"},
        
        {"Phi", "#Phi", "custom_range", "-3.2", "3.2",
            "custom_binning", "100", "", "legend", "electron_SF"},
        
        {"SFPCAL", "SF PCAL", "custom_range", "0.", "0.3",
            "custom_binning", "100", "", "legend", "positron_SF_4"},
        
        {"SFECIN", "SF ECIN", "custom_range", "0.", "0.25",
            "custom_binning", "100", "", "legend", "positron_SF_4"},
        
        {"SFECOUT", "SF ECOUT", "custom_range", "0.", "0.25",
            "custom_binning", "100", "", "legend", "positron_SF_4"},
        
        {"m2PCAL", "m2 PCAL", "custom_range", "0.", "200",
            "custom_binning", "100", "", "legend", "positron_SF_4"},
        
        {"m2ECIN", "m2 ECIN", "custom_range", "0.", "500",
            "custom_binning", "100", "", "legend", "positron_SF_4"},
        
        {"m2ECOUT", "m2 ECOUT", "custom_range", "0.", "500",
            "custom_binning", "100", "", "legend", "positron_SF_4"}
    };
    
    //TFile *Data_file = new TFile("Documents/Lepton_ID_root/"+name+"_Pion.root");
    
    TFile *Data_file = new TFile(infile+".root");
    TTree *Data_tree = (TTree *)Data_file->Get("resP");
    TTree *Data_tree_1 = (TTree *)Data_file->Get("resN");
    int entries_tree = Data_tree->GetEntries();
    int entries_tree_1 = Data_tree_1->GetEntries();
    TCut cut_0 = "";
    TCut cut_1 = "";
    
    
    Float_t P_P, Theta_P, Phi_P,SFPCAL_P,SFECIN_P,SFECOUT_P, m2PCAL_P, m2ECIN_P, m2ECOUT_P;
    Float_t P_N, Theta_N, Phi_N,SFPCAL_N,SFECIN_N,SFECOUT_N, m2PCAL_N, m2ECIN_N, m2ECOUT_N;
    
    Float_t scoreBDT_P,scoreBDT_N,score,scoreMLP_P,scoreMLP_N;
    
    //---------------Tree 1-------------------------------
    //TMVAPos.root, has true positives and false negatives
    //--------------------------------------------------
    //Data: #positives, #True positives, #False Negatives #score
    Data_tree->SetBranchAddress("scoreBDT_P",&scoreBDT_P);
    Data_tree->SetBranchAddress("scoreMLP_P",&scoreMLP_P);
    //positives sample variables
    //Data_tree->SetBranchAddress("P_P", &P_P);
    Data_tree->SetBranchAddress("P_P", &P_P);
    Data_tree->SetBranchAddress("Theta_P", &Theta_P);
    Data_tree->SetBranchAddress("Phi_P", &Phi_P);
    Data_tree->SetBranchAddress("SFPCAL_P", &SFPCAL_P);
    Data_tree->SetBranchAddress("SFECIN_P", &SFECIN_P);
    Data_tree->SetBranchAddress("SFECOUT_P", &SFECOUT_P);
    Data_tree->SetBranchAddress("m2PCAL_P", &m2PCAL_P);
    Data_tree->SetBranchAddress("m2ECIN_P", &m2ECIN_P);
    Data_tree->SetBranchAddress("m2ECOUT_P", &m2ECOUT_P);
    
    
    //---------------Tree 1-------------------------------
    //TMVAPos.root, has true positives and false negatives
    //--------------------------------------------------
    //Data: #negatives, #True negatives, #False positives #score
    Data_tree_1->SetBranchAddress("scoreBDT_N",&scoreBDT_N);
    Data_tree_1->SetBranchAddress("scoreMLP_N",&scoreMLP_N);
    //negatives sample variables
    //Data_tree->SetBranchAddress("P_N", &P_N);
    Data_tree_1->SetBranchAddress("P_N", &P_N);
    Data_tree_1->SetBranchAddress("Theta_N", &Theta_N);
    Data_tree_1->SetBranchAddress("Phi_N", &Phi_N);
    Data_tree_1->SetBranchAddress("SFPCAL_N", &SFPCAL_N);
    Data_tree_1->SetBranchAddress("SFECIN_N", &SFECIN_N);
    Data_tree_1->SetBranchAddress("SFECOUT_N", &SFECOUT_N);
    Data_tree_1->SetBranchAddress("m2PCAL_N", &m2PCAL_N);
    Data_tree_1->SetBranchAddress("m2ECIN_N", &m2ECIN_N);
    Data_tree_1->SetBranchAddress("m2ECOUT_N", &m2ECOUT_N);
    
    
    TString name_pdf=outfile+name_model;
    TString n1="_P";
    TString n2="_N";
    TString label = "True Positives";
    TString label_1 = "Negative Sample";
    TString label_2 = "True Negatives";
    TString label_3 = "Predicted Negative";
    
    int PS=entries_tree;
    int NS=entries_tree_1;
    int PP=0;
    int PN=0;
    
    TCanvas *cancG0 = new TCanvas("cancG0", "cancG0", 1500, 1200);
    cancG0->Divide(3,3,0.01,0.01,0);
    
    TCanvas *cancG1 = new TCanvas("cancG1", "cancG1", 1500, 1200);
    cancG1->Divide(3,3,0.01,0.01,0);
    
    TCanvas *cancG2 = new TCanvas("cancG2", "cancG2", 1500, 1200);
    cancG2->Divide(3,3,0.01,0.01,0);
    
    TCanvas *cancG3 = new TCanvas("cancG3", "cancG3", 1500, 1200);
    cancG3->Divide(3,3,0.01,0.01,0);
    
    
    cout << "-------------------------------------------------" << endl;
    cout << "Plotting TP and TN" << endl;
    cout << "-------------------------------------------------" << endl;
    
    for (int i = 0; i < labels1D.size(); i++){
        int PP=0;
        int PN=0;
        cancG0->cd(i+1);
        
        //////////////////////////////////////////////////
        // Set options for each label
        TString label1 = labels1D[i][0];
        TString xAxis_label = labels1D[i][1];
        TString range_x_option = labels1D[i][2];
        TString min_x_option = labels1D[i][3];
        TString max_x_option = labels1D[i][4];
        TString binning_x_option = labels1D[i][5];
        TString nb_x_bins = labels1D[i][6];
        TString string_cut = labels1D[i][7];
        TString legend_option = labels1D[i][8];
        TString output_string = labels1D[i][9];
        //////////////////////////////////////////////////
        cout << "Doing 1D " << label1 << " plot" << endl;
        
        
        TCut cut = string_cut.Data();
        
        TString cut_string = cut.GetTitle();
        
        float min_X_histo_ini;
        float max_X_histo_ini;
        
        int nBins_X = 20;
        //Extract custom range options
        if (range_x_option == "custom_range"){
            //cout << "Reading range " << endl;
            min_X_histo_ini = stof((string)min_x_option.Data());
            max_X_histo_ini = stof((string)max_x_option.Data());
        }
        
        //Extract custom bins options
        if (binning_x_option == "custom_binning"){
            //cout << "Reading # of data bins: " << nb_x_bins.Data() << endl;
            nBins_X = stoi((string)nb_x_bins.Data());
        }
        
        //True Positives:Positives from the positive sample
        TH1F *Data_hist = new TH1F("Data_hist", "Data_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
        if(name_model=="MLP")
            Data_tree->Draw(label1 +n1+ ">>Data_hist", Form(cut_1 * cut&&"scoreMLP_P>=%g",score_cut));
        if(name_model=="BDT")
            Data_tree->Draw(label1 +n1+ ">>Data_hist",Form( cut_1 * cut&&"scoreBDT_P>=%g",score_cut));
        
        
        
        
        //True Negatives: Negatives from the negative sample
        TH1F *Data_hist_2 = new TH1F("Data_hist_2", "Data_hist_2", nBins_X, min_X_histo_ini, max_X_histo_ini);
        if(name_model=="MLP")
            Data_tree_1->Draw(label1 +n2+ ">>Data_hist_2", Form(cut_1 * cut&&"scoreMLP_N<%g",score_cut));
        if(name_model=="BDT")
            Data_tree_1->Draw(label1 +n2+ ">>Data_hist_2",Form( cut_1 * cut&&"scoreBDT_N<%g",score_cut));
        
        
        Data_hist->Scale(1./Data_hist->Integral());
        Data_hist_2->Scale(1./Data_hist_2->Integral());
        
        
        Data_hist->SetLineColorAlpha(kCyan+2,0.35);
        Data_hist->SetMarkerColor(kCyan-3);
        Data_hist->SetFillColor(0);
        Data_hist->SetMarkerStyle(20);
        Data_hist->SetFillStyle(0);
        Data_hist->SetTitle(";" + xAxis_label+";counts");
        
        
        //Data_hist_1->SetLineWidth(2);
        //Data_hist_2->SetLineColor(kAzure-9);
        Data_hist_2->SetMarkerColor(kRed-7);
        //Data_hist_2->SetMarkerStyle(5);
        Data_hist_2->SetFillStyle(0);
        Data_hist_2->SetFillColor(0);
        Data_hist_2->SetTitle(";" + xAxis_label+";counts");
        if(i==2)
            Data_hist->SetMaximum(std::max(Data_hist_2->GetMaximum(),Data_hist->GetMaximum()) * 1.5);
        else
            Data_hist->SetMaximum(std::max(Data_hist_2->GetMaximum(),Data_hist->GetMaximum()) * 1.2);
        
        
        //Data_hist->Draw();
        //Data_hist_2->Draw("same");
        
        TRatioPlot *rp = new TRatioPlot(Data_hist, Data_hist_2);
        //cancG1->SetTicks(0, 1);
        rp->Draw();
        rp->GetLowerRefGraph()->SetMinimum(0);
        rp->GetLowerRefGraph()->SetMaximum(2);
        rp->GetLowYaxis()->SetNdivisions(505);
        
        rp->GetUpperPad()->cd();
        Data_hist->Draw("P same");
        Data_hist_2->SetMarkerStyle(20);
        Data_hist_2->Draw("PE same");
        
        int positive=Data_hist->GetEntries();
        int negative=Data_hist_2->GetEntries();
        
        
        
        auto legend = new TLegend(0.65, 0.87, 0.87, 0.67);
        legend->AddEntry(Data_hist, Form("True Positives(%d)",positive), "p");
        legend->AddEntry(Data_hist_2, Form("True Negatives(%d)",negative), "p");
        legend->SetFillStyle(0);
        legend->SetLineWidth(0);
        
        
        legend->Draw("same ");
        
    }
    
    
    
    cout<<"Saving in "<<name_pdf<<".pdf"<<endl;
    cancG0->SaveAs(name_pdf + ".pdf(");
    
    cout << "-------------------------------------------------" << endl;
    cout << "Plotting TP and FP" << endl;
    cout << "-------------------------------------------------" << endl;
    
    for (int i = 0; i < labels1D.size(); i++){
        int PP=0;
        int PN=0;
        cancG1->cd(i+1);
        
        //////////////////////////////////////////////////
        // Set options for each label
        TString label1 = labels1D[i][0];
        TString xAxis_label = labels1D[i][1];
        TString range_x_option = labels1D[i][2];
        TString min_x_option = labels1D[i][3];
        TString max_x_option = labels1D[i][4];
        TString binning_x_option = labels1D[i][5];
        TString nb_x_bins = labels1D[i][6];
        TString string_cut = labels1D[i][7];
        TString legend_option = labels1D[i][8];
        TString output_string = labels1D[i][9];
        //////////////////////////////////////////////////
        cout << "Doing 1D " << label1 << " plot" << endl;
        
        
        TCut cut = string_cut.Data();
        
        TString cut_string = cut.GetTitle();
        
        float min_X_histo_ini;
        float max_X_histo_ini;
        
        int nBins_X = 20;
        //Extract custom range options
        if (range_x_option == "custom_range"){
            //cout << "Reading range " << endl;
            min_X_histo_ini = stof((string)min_x_option.Data());
            max_X_histo_ini = stof((string)max_x_option.Data());
        }
        
        //Extract custom bins options
        if (binning_x_option == "custom_binning"){
            //cout << "Reading # of data bins: " << nb_x_bins.Data() << endl;
            nBins_X = stoi((string)nb_x_bins.Data());
        }
        
        
        //True Positives
        TH1F *Data_hist_1 = new TH1F("Data_hist_1", "Data_hist_1", nBins_X, min_X_histo_ini, max_X_histo_ini);
        if(name_model=="MLP")
            Data_tree->Draw(label1 +n1+ ">>Data_hist_1", Form(cut_1 * cut&&"scoreMLP_P>=%g",score_cut));
        if(name_model=="BDT")
            Data_tree->Draw(label1 +n1+ ">>Data_hist_1",Form( cut_1 * cut&&"scoreBDT_P>=%g",score_cut));
        
        //False Positives
        TH1F *Data_hist_3 = new TH1F("Data_hist_3", "Data_hist_3", nBins_X, min_X_histo_ini, max_X_histo_ini);
        if(name_model=="MLP")
            Data_tree_1->Draw(label1 +n2+ ">>Data_hist_3", Form(cut_1 * cut&&"scoreMLP_N>=%g",score_cut));
        if(name_model=="BDT")
            Data_tree_1->Draw(label1 +n2+ ">>Data_hist_3",Form( cut_1 * cut&&"scoreBDT_N>=%g",score_cut));
        
        
        
        Data_hist_1->Scale(1./Data_hist_1->Integral());
        Data_hist_3->Scale(1./Data_hist_3->Integral());
        
        
        
        Data_hist_1->SetMarkerStyle(20);
        Data_hist_1->SetLineColorAlpha(kCyan+2,0.35);
        Data_hist_1->SetMarkerColor(kCyan-3);
        //Data_hist_1->SetLineColor(kAzure-6);
        Data_hist_1->SetTitle(";" + xAxis_label+";counts");
        
        
        //Data_hist_1->SetLineWidth(2);
        Data_hist_3->SetLineColor(0);
        Data_hist_3->SetFillStyle(0);
        Data_hist_3->SetFillColor(0);
        Data_hist_3->SetMarkerColor(kOrange-3);
        Data_hist_3->SetTitle(";" + xAxis_label+";counts");
        
        Data_hist_1->SetMaximum(std::max(Data_hist_3->GetMaximum(),Data_hist_1->GetMaximum()) * 1.2);
        
        //Data_hist_1->Draw();
        
        int TP=Data_hist_1->GetEntries();
        int FP=Data_hist_3->GetEntries();
        
        TRatioPlot *rp = new TRatioPlot(Data_hist_1, Data_hist_3);
        cancG1->SetTicks(0, 1);
        rp->Draw();
        rp->GetLowerRefGraph()->SetMinimum(0);
        rp->GetLowerRefGraph()->SetMaximum(2);
        rp->GetLowYaxis()->SetNdivisions(505);
        //rp->SetSeparationMargin(0.0);
        rp->GetUpperPad()->cd();
        Data_hist_1->Draw("P same");
        
        Data_hist_3->SetMarkerStyle(20);
        Data_hist_3->Draw("PE same");
        
        auto legend = new TLegend(0.65, 0.87, 0.87, 0.67);
        legend->AddEntry(Data_hist_1, Form("True Positives (%d)", TP), "p");
        legend->AddEntry(Data_hist_3, Form("False Positives (%d)", FP), "p");
        legend->SetFillStyle(0);
        legend->SetLineWidth(0);
        
        
        legend->Draw("same ");
        
    }
    
    cout<<"Saving in "<<name_pdf<<".pdf"<<endl;
    cancG1->SaveAs(name_pdf + ".pdf");
    
    cout << "-------------------------------------------------" << endl;
    cout << "Plotting NS and FP" << endl;
    cout << "-------------------------------------------------" << endl;
    
    for (int i = 0; i < labels1D.size(); i++){
        cancG2->cd(i+1);
        
        //////////////////////////////////////////////////
        // Set options for each label
        TString label1 = labels1D[i][0];
        TString xAxis_label = labels1D[i][1];
        TString range_x_option = labels1D[i][2];
        TString min_x_option = labels1D[i][3];
        TString max_x_option = labels1D[i][4];
        TString binning_x_option = labels1D[i][5];
        TString nb_x_bins = labels1D[i][6];
        TString string_cut = labels1D[i][7];
        TString legend_option = labels1D[i][8];
        TString output_string = labels1D[i][9];
        //////////////////////////////////////////////////
        cout << "Doing 1D " << label1 << " plot" << endl;
        
        
        TCut cut = string_cut.Data();
        
        TString cut_string = cut.GetTitle();
        
        float min_X_histo_ini;
        float max_X_histo_ini;
        
        int nBins_X = 20;
        //Extract custom range options
        if (range_x_option == "custom_range"){
            //cout << "Reading range " << endl;
            min_X_histo_ini = stof((string)min_x_option.Data());
            max_X_histo_ini = stof((string)max_x_option.Data());
        }
        
        //Extract custom bins options
        if (binning_x_option == "custom_binning"){
            //cout << "Reading # of data bins: " << nb_x_bins.Data() << endl;
            nBins_X = stoi((string)nb_x_bins.Data());
        }
        
        
        //Negative Sample
        TH1F *Data_hist_1 = new TH1F("Data_hist_1", "Data_hist_1", nBins_X, min_X_histo_ini, max_X_histo_ini);
        Data_tree_1->Draw(label1 +n2+ ">>Data_hist_1", cut_1 * cut);
        
        //False Positives
        TH1F *Data_hist_3 = new TH1F("Data_hist_3", "Data_hist_3", nBins_X, min_X_histo_ini, max_X_histo_ini);
        if(name_model=="MLP")
            Data_tree_1->Draw(label1 +n2+ ">>Data_hist_3", Form(cut_1 * cut&&"scoreMLP_N>=%g",score_cut));
        if(name_model=="BDT")
            Data_tree_1->Draw(label1 +n2+ ">>Data_hist_3",Form( cut_1 * cut&&"scoreBDT_N>=%g",score_cut));
        
        
        
        //Data_hist_1->Scale(1./Data_hist_1->Integral());
        //Data_hist_3->Scale(1./Data_hist_3->Integral());
        
        
        
        //Data_hist->SetLineWidth(2);
        Data_hist_1->SetLineColor(kRed);
        Data_hist_1->SetFillColor(0);
        Data_hist_1->SetFillStyle(0);
        Data_hist_1->SetTitle(";" + xAxis_label+";counts");
        
        
        //Data_hist_1->SetLineWidth(2);
        Data_hist_3->SetLineColor(0);
        Data_hist_3->SetFillStyle(0);
        Data_hist_3->SetFillColor(0);
        Data_hist_3->SetMarkerColor(kOrange-3);
        Data_hist_3->SetTitle(";" + xAxis_label+";counts");
        
        Data_hist_3->SetMaximum(std::max(Data_hist_3->GetMaximum(),Data_hist_1->GetMaximum()) * 1.2);
        
        
        
        int FP=Data_hist_3->GetEntries();
        
        TRatioPlot *rp = new TRatioPlot(Data_hist_3, Data_hist_1);
        cancG2->SetTicks(0, 1);
        rp->Draw();
        rp->GetLowerRefGraph()->SetMinimum(0.0);
        rp->GetLowerRefGraph()->SetMaximum(0.2);
        rp->GetLowYaxis()->SetNdivisions(505);
        //rp->SetSeparationMargin(0.0);
        rp->GetUpperPad()->cd();
        
        //Data_hist_1->SetLineColor(20);
        Data_hist_1->Draw("same");
        Data_hist_3->SetMarkerStyle(20);
        Data_hist_3->Draw("PE same");
        
        auto legend = new TLegend(0.65, 0.87, 0.87, 0.67);
        legend->AddEntry(Data_hist_1, Form("Negative Sample (%d)", NS), "l");
        legend->AddEntry(Data_hist_3, Form("False Positives (%d)", FP), "p");
        legend->SetFillStyle(0);
        legend->SetLineWidth(0);
        
        
        legend->Draw("same ");
        
    }
    
    cout<<"Saving in "<<name_pdf<<".pdf"<<endl;
    cancG2->SaveAs(name_pdf + ".pdf");
    
    
    
    cout << "-------------------------------------------------" << endl;
    cout << "Plotting Positive Sample and TP" << endl;
    cout << "-------------------------------------------------" << endl;
    
    for (int i = 0; i < labels1D.size(); i++){
        int PP=0;
        int PN=0;
        cancG3->cd(i+1);
        
        //////////////////////////////////////////////////
        // Set options for each label
        TString label1 = labels1D[i][0];
        TString xAxis_label = labels1D[i][1];
        TString range_x_option = labels1D[i][2];
        TString min_x_option = labels1D[i][3];
        TString max_x_option = labels1D[i][4];
        TString binning_x_option = labels1D[i][5];
        TString nb_x_bins = labels1D[i][6];
        TString string_cut = labels1D[i][7];
        TString legend_option = labels1D[i][8];
        TString output_string = labels1D[i][9];
        //////////////////////////////////////////////////
        cout << "Doing 1D " << label1 << " plot" << endl;
        
        
        TCut cut = string_cut.Data();
        
        TString cut_string = cut.GetTitle();
        
        float min_X_histo_ini;
        float max_X_histo_ini;
        
        int nBins_X = 20;
        //Extract custom range options
        if (range_x_option == "custom_range"){
            //cout << "Reading range " << endl;
            min_X_histo_ini = stof((string)min_x_option.Data());
            max_X_histo_ini = stof((string)max_x_option.Data());
        }
        
        //Extract custom bins options
        if (binning_x_option == "custom_binning"){
            //cout << "Reading # of data bins: " << nb_x_bins.Data() << endl;
            nBins_X = stoi((string)nb_x_bins.Data());
        }
        
        //Positive Sample
        TH1F *Data_hist = new TH1F("Data_hist", "Data_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
        Data_tree->Draw(label1 +n1+">>Data_hist", cut_0 * cut);
        
        int positive=Data_tree->GetEntries();
        
        
        
        //True Positive
        TH1F *Data_hist_2 = new TH1F("Data_hist_2", "Data_hist_2", nBins_X, min_X_histo_ini, max_X_histo_ini);
        if(name_model=="MLP")
            Data_tree->Draw(label1 +n1+ ">>Data_hist_2", Form(cut_1 * cut&&"scoreMLP_P>=%g",score_cut));
        if(name_model=="BDT")
            Data_tree->Draw(label1 +n1+ ">>Data_hist_2",Form( cut_1 * cut&&"scoreBDT_P>=%g",score_cut));
        
        
        //Data_hist->Scale(1./Data_hist->Integral());
        //Data_hist_2->Scale(1./Data_hist_2->Integral());
        
        
        //Data_hist->SetLineWidth(2);
        Data_hist->SetLineColor(kAzure-6);
        Data_hist->SetFillColor(0);
        Data_hist->SetFillStyle(0);
        Data_hist->SetTitle(";" + xAxis_label+";counts");
        
        
        //Data_hist_1->SetLineWidth(2);
        //Data_hist_2->SetLineColor(kAzure-9);
        Data_hist_2->SetMarkerColor(kCyan-3);
        //Data_hist_2->SetMarkerStyle(5);
        Data_hist_2->SetFillStyle(0);
        Data_hist_2->SetFillColor(0);
        Data_hist_2->SetTitle(";" + xAxis_label+";counts");
        if(i==2)
            Data_hist->SetMaximum(std::max(Data_hist_2->GetMaximum(),Data_hist->GetMaximum()) * 1.5);
        else
            Data_hist->SetMaximum(std::max(Data_hist_2->GetMaximum(),Data_hist->GetMaximum()) * 1.2);
        
        
        Data_hist->Draw("same");
        //Data_hist_2->Draw("same");
        
        TRatioPlot *rp = new TRatioPlot(Data_hist_2, Data_hist);
        //cancG1->SetTicks(0, 1);
        rp->Draw();
        rp->GetLowerRefGraph()->SetMinimum(0.8);
        rp->GetLowerRefGraph()->SetMaximum(1.2);
        rp->GetLowYaxis()->SetNdivisions(505);
        
        rp->GetUpperPad()->cd();
        
        Data_hist_2->SetMarkerStyle(20);
        Data_hist_2->Draw("PE same");
        
        int TP=Data_hist_2->GetEntries();
        
        auto legend = new TLegend(0.65, 0.87, 0.87, 0.67);
        legend->AddEntry(Data_hist, Form("Positives Sample(%d)",PS), "l");
        legend->AddEntry(Data_hist_2, Form("True Positives(%d)",TP), "p");
        legend->SetFillStyle(0);
        legend->SetLineWidth(0);
        
        
        legend->Draw("same ");
        
    }
    
    
    
    cout<<"Saving in "<<name_pdf<<".pdf"<<endl;
    cancG3->SaveAs(name_pdf + ".pdf)");
    
}

void True_False(TString infile, TString name_model,Float_t score_cut, TString outfile)
{
    gROOT->SetBatch(kTRUE);
    
    gStyle->SetPaintTextFormat("4.1f");
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(.03, "xyz");
    gStyle->SetTitleSize(.03, "xyz");
    gStyle->SetTitleSize(.07, "t");
    gStyle->SetFrameLineWidth(1);
    gStyle->SetLineWidth(1);
    gStyle->SetHistLineWidth(1);
    gStyle->SetMarkerStyle(13);
    gStyle->SetTitleW(0.8); // per cent of the pad width
    gStyle->SetTitleH(0.1); // per cent of the pad height
    
    
    std::vector<vector<TString>> labels1D{
        {"P", "P", "custom_range", "0.", "14.",
            "custom_binning", "100", "", "legend", "P"},
        
        {"Theta", "#Theta", "custom_range", "0.", "1",
            "custom_binning", "100", "", "legend", "electron_SF"},
        
        {"Phi", "#Phi", "custom_range", "-3.2", "3.2",
            "custom_binning", "100", "", "legend", "electron_SF"},
        
        {"SFPCAL", "SF PCAL", "custom_range", "0.", "0.3",
            "custom_binning", "100", "", "legend", "positron_SF_4"},
        
        {"SFECIN", "SF ECIN", "custom_range", "0.", "0.25",
            "custom_binning", "100", "", "legend", "positron_SF_4"},
        
        {"SFECOUT", "SF ECOUT", "custom_range", "0.", "0.25",
            "custom_binning", "100", "", "legend", "positron_SF_4"},
        
        {"m2PCAL", "m2 PCAL", "custom_range", "0.", "200",
            "custom_binning", "100", "", "legend", "positron_SF_4"},
        
        {"m2ECIN", "m2 ECIN", "custom_range", "0.", "500",
            "custom_binning", "100", "", "legend", "positron_SF_4"},
        
        {"m2ECOUT", "m2 ECOUT", "custom_range", "0.", "500",
            "custom_binning", "100", "", "legend", "positron_SF_4"}
    };
    
    //TFile *Data_file = new TFile("Documents/Lepton_ID_root/"+name+"_Pion.root");
    
    TFile *Data_file = new TFile(infile+".root");
    TTree *Data_tree = (TTree *)Data_file->Get("resP");
    TTree *Data_tree_1 = (TTree *)Data_file->Get("resN");
    int entries_tree = Data_tree->GetEntries();
    int entries_tree_1 = Data_tree_1->GetEntries();
    TCut cut_0 = "";
    TCut cut_1 = "";
    
    int TP,TN,FP,FN;
    
    Float_t P_P, Theta_P, Phi_P,SFPCAL_P,SFECIN_P,SFECOUT_P, m2PCAL_P, m2ECIN_P, m2ECOUT_P;
    Float_t P_N, Theta_N, Phi_N,SFPCAL_N,SFECIN_N,SFECOUT_N, m2PCAL_N, m2ECIN_N, m2ECOUT_N;
    
    Float_t scoreBDT_P,scoreBDT_N,score,scoreMLP_P,scoreMLP_N;
    
    //---------------Tree 1-------------------------------
    //TMVAPos.root, has true positives and false negatives
    //--------------------------------------------------
    //Data: #positives, #True positives, #False Negatives #score
    Data_tree->SetBranchAddress("scoreBDT_P",&scoreBDT_P);
    Data_tree->SetBranchAddress("scoreMLP_P",&scoreMLP_P);
    //positives sample variables
    //Data_tree->SetBranchAddress("P_P", &P_P);
    Data_tree->SetBranchAddress("P_P", &P_P);
    Data_tree->SetBranchAddress("Theta_P", &Theta_P);
    Data_tree->SetBranchAddress("Phi_P", &Phi_P);
    Data_tree->SetBranchAddress("SFPCAL_P", &SFPCAL_P);
    Data_tree->SetBranchAddress("SFECIN_P", &SFECIN_P);
    Data_tree->SetBranchAddress("SFECOUT_P", &SFECOUT_P);
    Data_tree->SetBranchAddress("m2PCAL_P", &m2PCAL_P);
    Data_tree->SetBranchAddress("m2ECIN_P", &m2ECIN_P);
    Data_tree->SetBranchAddress("m2ECOUT_P", &m2ECOUT_P);
    
    
    //---------------Tree 1-------------------------------
    //TMVAPos.root, has true positives and false negatives
    //--------------------------------------------------
    //Data: #negatives, #True negatives, #False positives #score
    Data_tree_1->SetBranchAddress("scoreBDT_N",&scoreBDT_N);
    Data_tree_1->SetBranchAddress("scoreMLP_N",&scoreMLP_N);
    //negatives sample variables
    //Data_tree->SetBranchAddress("P_N", &P_N);
    Data_tree_1->SetBranchAddress("P_N", &P_N);
    Data_tree_1->SetBranchAddress("Theta_N", &Theta_N);
    Data_tree_1->SetBranchAddress("Phi_N", &Phi_N);
    Data_tree_1->SetBranchAddress("SFPCAL_N", &SFPCAL_N);
    Data_tree_1->SetBranchAddress("SFECIN_N", &SFECIN_N);
    Data_tree_1->SetBranchAddress("SFECOUT_N", &SFECOUT_N);
    Data_tree_1->SetBranchAddress("m2PCAL_N", &m2PCAL_N);
    Data_tree_1->SetBranchAddress("m2ECIN_N", &m2ECIN_N);
    Data_tree_1->SetBranchAddress("m2ECOUT_N", &m2ECOUT_N);
    
    
    TString name1="TPvsFP";
    TString name2="FNvsTN";
    TString name_pdf=outfile+name_model;
    //Float_t score_cut=0.86;
    
    //From file
    TString n1="_P";
    TString n2="_N";
    //TString n3="_PP";
    //TString n4="_PN";
    
    
    TString label = "True Positive";
    TString label_1 = "False Positive";
    TString label_2 = "True Negative";
    TString label_3 = "False Negative";
    
    int PS=entries_tree;
    int NS=entries_tree_1;
    int T=0;
    int F=0;
    
    TCanvas *cancG0 = new TCanvas("cancG0", "cancG0", 1500, 1200);
    cancG0->Divide(3,3,0.01,0.01,0);
    
    cout << "-------------------------------------------------" << endl;
    cout << "Plotting True Positive vs False Positive" << endl;
    cout << "-------------------------------------------------" << endl;
    for (int i = 0; i < labels1D.size(); i++){
        int PP=0;
        int PN=0;
        cancG0->cd(i+1);
        
        //////////////////////////////////////////////////
        // Set options for each label
        TString label1 = labels1D[i][0];
        TString xAxis_label = labels1D[i][1];
        TString range_x_option = labels1D[i][2];
        TString min_x_option = labels1D[i][3];
        TString max_x_option = labels1D[i][4];
        TString binning_x_option = labels1D[i][5];
        TString nb_x_bins = labels1D[i][6];
        TString string_cut = labels1D[i][7];
        TString legend_option = labels1D[i][8];
        TString output_string = labels1D[i][9];
        //////////////////////////////////////////////////
        cout << "Doing 1D " << label1 << " plot" << endl;
        
        TCut cut = string_cut.Data();
        
        TString cut_string = cut.GetTitle();
        
        float min_X_histo_ini;
        float max_X_histo_ini;
        
        int nBins_X = 20;
        //Extract custom range options
        if (range_x_option == "custom_range"){
            //cout << "Reading range " << endl;
            min_X_histo_ini = stof((string)min_x_option.Data());
            max_X_histo_ini = stof((string)max_x_option.Data());
        }
        
        //Extract custom bins options
        if (binning_x_option == "custom_binning"){
            //cout << "Reading # of data bins: " << nb_x_bins.Data() << endl;
            nBins_X = stoi((string)nb_x_bins.Data());
        }
        
        //TCanvas *c1 = new TCanvas("c1","c1");
        //setupCanvas(c1);
        
        
        //True Positives: Positives (score>=0.86) from positive sample (n1=_P)
        TH1F *Data_hist = new TH1F("Data_hist", "Data_hist", nBins_X, min_X_histo_ini, max_X_histo_ini);
        //Float_t test_cut=0.01;
        if(name_model=="MLP")
            Data_tree->Draw(label1 +n1+">>Data_hist", Form(cut_0 * cut&&"scoreMLP_P>=%g",score_cut));
        if(name_model=="BDT")
            Data_tree->Draw(label1 +n1+">>Data_hist", Form(cut_0 * cut&&"scoreBDT_P>=%g",score_cut));
        
        //TCut weight = Form("%i/%i", entries_tree, entries_tree_1);
        TCut weight = "";
        
        //False Positives: Positives (score>=0.86) from negative sample (n2=_N)
        TH1F *Data_hist_1 = new TH1F("Data_hist_1", "Data_hist_1", nBins_X, min_X_histo_ini, max_X_histo_ini);
        if(name_model=="MLP")
            Data_tree_1->Draw(label1 +n2+ ">>Data_hist_1", Form(weight*cut_1 * cut&&"scoreMLP_N>=%g",score_cut));
        if(name_model=="BDT")
            Data_tree_1->Draw(label1 +n2+ ">>Data_hist_1",Form(weight* cut_1 * cut&&"scoreBDT_N>=%g",score_cut));
        
        Data_hist->Scale(1./Data_hist->Integral(),"nosw2");
        Data_hist_1->Scale(1./Data_hist_1->Integral(),"nosw2");
        
        
        Data_hist->SetMaximum(std::max(Data_hist->GetMaximum(),Data_hist_1->GetMaximum())*1.2);
        
        //Data_hist->SetLineWidth(2);
        Data_hist->SetLineColor(kBlue);
        Data_hist->SetFillStyle(0);
        Data_hist->SetFillColor(0);
        Data_hist->SetTitle(";" + xAxis_label+";counts");
        
        //Data_hist_1->SetLineWidth(2);
        Data_hist_1->SetLineColor(0);
        Data_hist_1->SetFillStyle(0);
        Data_hist_1->SetFillColor(0);
        Data_hist_1->SetMarkerColor(kMagenta);
        Data_hist_1->SetTitle(";" + xAxis_label+";counts");
        
        Data_hist->Draw();
        Data_hist_1->Draw("same");
        
        TRatioPlot *rp = new TRatioPlot(Data_hist, Data_hist_1);
        cancG0->SetTicks(0, 1);
        rp->Draw();
        rp->GetLowerRefGraph()->SetMinimum(0);
        rp->GetLowerRefGraph()->SetMaximum(2);
        rp->GetLowYaxis()->SetNdivisions(505);
        
        
        
        auto legend = new TLegend(0.65, 0.87, 0.87, 0.67);
        legend->SetFillStyle(0);
        legend->SetLineWidth(0);
        TP=Data_hist->GetEntries();
        FP=Data_hist_1->GetEntries();
        legend->AddEntry(Data_hist, Form("%s(%d)", label.Data(),TP), "l");
        legend->AddEntry(Data_hist_1, Form("%s(%d)", label_1.Data(),FP), "lp");
        legend->Draw("same ");
        
        
        
    }
    //cancG0->cd(0);
    cout<<"Saving in "<<name_pdf<<".pdf"<<endl;
    cancG0->SaveAs(name_pdf + ".pdf");
    
    
    TCanvas *cancG1 = new TCanvas("cancG1", "cancG1", 1500, 1200);
    cancG1->Divide(3,3,0.01,0.01,0);
    
    cout << "-------------------------------------------------" << endl;
    cout << "Plotting True Negative vs False Negative" << endl;
    cout << "-------------------------------------------------" << endl;
    for (int i = 0; i < labels1D.size(); i++){
        int PP=0;
        int PN=0;
        cancG1->cd(i+1);
        
        //////////////////////////////////////////////////
        // Set options for each label
        TString label1 = labels1D[i][0];
        TString xAxis_label = labels1D[i][1];
        TString range_x_option = labels1D[i][2];
        TString min_x_option = labels1D[i][3];
        TString max_x_option = labels1D[i][4];
        TString binning_x_option = labels1D[i][5];
        TString nb_x_bins = labels1D[i][6];
        TString string_cut = labels1D[i][7];
        TString legend_option = labels1D[i][8];
        TString output_string = labels1D[i][9];
        //////////////////////////////////////////////////
        cout << "Doing 1D " << label1 << " plot" << endl;
        
        TCut cut = string_cut.Data();
        
        TString cut_string = cut.GetTitle();
        
        float min_X_histo_ini;
        float max_X_histo_ini;
        
        int nBins_X = 20;
        //Extract custom range options
        if (range_x_option == "custom_range"){
            //cout << "Reading range " << endl;
            min_X_histo_ini = stof((string)min_x_option.Data());
            max_X_histo_ini = stof((string)max_x_option.Data());
        }
        
        //Extract custom bins options
        if (binning_x_option == "custom_binning"){
            //cout << "Reading # of data bins: " << nb_x_bins.Data() << endl;
            nBins_X = stoi((string)nb_x_bins.Data());
        }
        
        TCut weight = "";
        
        //True Negatives: Negatives (score<0.86) from a negative sample (n2=_N)
        TH1F *Data_hist_2 = new TH1F("Data_hist_2", "Data_hist_2", nBins_X, min_X_histo_ini, max_X_histo_ini);
        if(name_model=="MLP")
            Data_tree_1->Draw(label1 +n2+">>Data_hist_2", Form(weight *cut_1 * cut&&"scoreMLP_N<%g",score_cut));
        if(name_model=="BDT")
            Data_tree_1->Draw(label1 +n2+">>Data_hist_2",Form(weight *cut_1 * cut&&"scoreBDT_N<%g",score_cut));
        
        //False Negatives: Negatives (score<0.86) from a positive sample (n1=_P)
        TH1F *Data_hist_3 = new TH1F("Data_hist_3", "Data_hist_3", nBins_X, min_X_histo_ini, max_X_histo_ini);
        if(name_model=="MLP")
            Data_tree->Draw(label1 +n1+">>Data_hist_3", Form(cut_0 * cut&&"scoreMLP_P<%g",score_cut));
        if(name_model=="BDT")
            Data_tree->Draw(label1 +n1+">>Data_hist_3", Form(cut_0 * cut&&"scoreBDT_P<%g",score_cut));
        
        
        
        
        //Double_t scale = Data_hist_3->GetXaxis()->GetBinWidth(1)/(Data_hist_3->Integral());
        
        Data_hist_2->Scale(1./Data_hist_2->Integral(),"nosw2");
        Data_hist_3->Scale(1./Data_hist_3->Integral(),"nosw2");
        
        //Data_hist_1->SetLineWidth(2);
        Data_hist_2->SetLineColor(kRed);
        Data_hist_2->SetFillStyle(3004);
        Data_hist_2->SetFillColor(kRed);
        Data_hist_2->SetTitle(";" + xAxis_label+";counts");
        
        
        //Data_hist_1->SetLineWidth(2);
        Data_hist_3->SetLineColor(kYellow+1);
        Data_hist_3->SetFillStyle(3004);
        Data_hist_3->SetFillColor(kYellow+1);
        Data_hist_3->SetTitle(";" + xAxis_label+";counts");
        
        
        Data_hist_3->SetMaximum(std::max(Data_hist_3->GetMaximum(),Data_hist_2->GetMaximum()) * 1.2);
        
        
        
        auto legend = new TLegend(0.65, 0.87, 0.87, 0.67);
        legend->SetFillStyle(0);
        legend->SetLineWidth(0);
        
        FN=Data_hist_2->GetEntries();
        TN=Data_hist_3->GetEntries();
        
        
        legend->AddEntry(Data_hist_2, Form("%s(%d)", label_2.Data(),FN), "l");
        legend->AddEntry(Data_hist_3, Form("%s(%d)", label_3.Data(),TN), "l");
        Data_hist_2->Draw("same");
        Data_hist_3->Draw("same");
        legend->Draw("same ");
        
    }
    //cancG0->cd(0);
    cout<<"Saving in "<<name_pdf<<".pdf"<<endl;
    
    cancG1->SaveAs(name_pdf+ ".pdf)");
    
    //gApplication->Terminate();
    
    
}



void ROC(TString infile, TString name_model="All"){
    
    gROOT->SetBatch(kTRUE);
    
    gStyle->SetPaintTextFormat("4.1f");
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(.03, "xyz");
    gStyle->SetTitleSize(.03, "xyz");
    gStyle->SetTitleSize(.07, "t");
    gStyle->SetFrameLineWidth(1);
    gStyle->SetLineWidth(1);
    gStyle->SetHistLineWidth(1);
    gStyle->SetMarkerStyle(13);
    gStyle->SetTitleW(0.8); // per cent of the pad width
    gStyle->SetTitleH(0.1); // per cent of the pad height
    
    TFile *Data_file = new TFile(infile+".root");
    TTree *Data_tree = (TTree *)Data_file->Get("resP");
    TTree *Data_tree_1 = (TTree *)Data_file->Get("resN");
    int entries_tree = Data_tree->GetEntries();
    int entries_tree_1 = Data_tree_1->GetEntries();
    
    
    int TP,TN,FP,FN;
    
    Float_t scoreBDT_P,scoreBDT_N,score,scoreMLP_P,scoreMLP_N;
    
    //---------------Tree 1-------------------------------
    //TMVAPos.root, has true positives and false negatives
    //--------------------------------------------------
    //Data: #positives, #True positives, #False Negatives #score
    Data_tree->SetBranchAddress("scoreBDT_P",&scoreBDT_P);
    Data_tree->SetBranchAddress("scoreMLP_P",&scoreMLP_P);
    
    
    //---------------Tree 1-------------------------------
    //TMVAPos.root, has true positives and false negatives
    //--------------------------------------------------
    //Data: #negatives, #True negatives, #False positives #score
    Data_tree_1->SetBranchAddress("scoreBDT_N",&scoreBDT_N);
    Data_tree_1->SetBranchAddress("scoreMLP_N",&scoreMLP_N);
    
    Float_t segments, start, end;
    Float_t  start_MLP, end_MLP;
    Float_t  start_BDT, end_BDT;
    Float_t TPR,FPR;
    Float_t score_cut;
    Int_t max=500;
    int i=1;
    Double_t x[max], y[max];
    Double_t x2[max], y2[max];
    
    
    start_MLP=-1.25;
    end_MLP=1.5;
    start_BDT=-0.8;
    end_BDT=0.8;
    
    
    
    
    
    if(name_model=="MLP"||name_model=="All"){
        cout<<"ROC--MLP---Start"<<endl;
        start=start_MLP;
        end=end_MLP;
        segments=(end-start)/max;
        
        score_cut=start;
        
        
        segments=(end-start)/max;
        while (score_cut<=end){
            TP=Data_tree->GetEntries(Form("scoreMLP_P>=%g",score_cut));
            TN=Data_tree_1->GetEntries(Form("scoreMLP_N<%g",score_cut));
            FP=Data_tree_1->GetEntries(Form("scoreMLP_N>=%g",score_cut));
            FN=Data_tree->GetEntries(Form("scoreMLP_P<%g",score_cut));
            
            TPR=(1.0*TP)/(TP+FN);
            FPR=(1.0*FP)/(FP+TN);
            
            
            x[i-1]=TPR;
            y[i-1]=1-FPR;
            
            score_cut=start+i*segments;
            i++;
        }
        cout<<"ROC--MLP---End"<<endl;
    }
    i=1;
    
    if(name_model=="BDT"||name_model=="All"){
        cout<<"ROC--BDT---Start"<<endl;
        
        start=start_BDT;
        end=end_BDT;
        segments=(end-start)/max;
        
        score_cut=start;
        
        while (score_cut<=end){
            TP=Data_tree->GetEntries(Form("scoreBDT_P>=%g",score_cut));
            TN=Data_tree_1->GetEntries(Form("scoreBDT_N<%g",score_cut));
            FP=Data_tree_1->GetEntries(Form("scoreBDT_N>=%g",score_cut));
            FN=Data_tree->GetEntries(Form("scoreBDT_P<%g",score_cut));
            
            TPR=(1.0*TP)/(TP+FN);
            FPR=(1.0*FP)/(FP+TN);
            
            
            x2[i-1]=TPR;
            y2[i-1]=1-FPR;
            
            score_cut=start+i*segments;
            i++;
        }
        cout<<"ROC--BDT---End"<<endl;
        
    }
    
    TCanvas *c1 = new TCanvas("c1","c1",1000, 800);
    
    TGraph *gr_MLP  = new TGraph(max,x,y);
    TGraph *gr_BDT  = new TGraph(max,x2,y2);
    
    
    
    if(name_model=="MLP"){
        gr_MLP->SetLineColor(1);
        gr_MLP->SetLineWidth(3);
        gr_MLP->SetMarkerStyle(20);
        gr_MLP->Draw("ALP");
        c1->SaveAs(name_model+"_ROC.png");
    }
    else if(name_model=="BDT"){
        gr_BDT->SetLineWidth(3);
        gr_BDT->SetMarkerStyle(21);
        gr_BDT->SetLineColor(2);
        gr_BDT->Draw("ALP");
        c1->SaveAs(name_model+"_ROC.png");
    }
    else{
        cout<<"PLOT ROC--MLP AND BDT"<<endl;
        Float_t areaMLP, areaBDT;
        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle("ROC curve for both methods; True Positives Rate; 1- False Positives Rate");
        
        auto legend = new TLegend();
        legend->SetFillStyle(0);
        legend->SetLineWidth(0);
        mg->GetXaxis()->SetRangeUser(0., 1.);
        mg->GetYaxis()->SetRangeUser(0., 1.01);
        
        
        gr_MLP->SetName("gr1");
        gr_MLP->SetLineColorAlpha(kBlue,0.95);
        gr_MLP->SetLineWidth(4);
        
        
        //gr_BDT->SetLineStyle(10);
        gr_BDT->SetName("gr2");
        gr_BDT->SetLineWidth(4);
        gr_BDT->SetLineColorAlpha(kRed, 0.45);
        
        
        areaMLP=gr_MLP->Integral();
        areaBDT=gr_BDT->Integral();
        
        
        mg->Add(gr_MLP);
        mg->Add(gr_BDT);
        mg->Draw("AL");
        
        legend->AddEntry("gr1",Form("MLP"),"l");
        legend->AddEntry("gr2",Form("BDT"),"l");
        legend->Draw();
        
        
        c1->SaveAs("ROC.png");
    }
    
    
    
    
    
    
    
    //TCanvas *cancG1 = new TCanvas("cancG1", "cancG1", 1500, 1200);
    //cancG1->Divide(3,3,0.01,0.01,0);
    
    //TH1F *Data_hist_3 = new TH1F("Data_hist_3", "Data_hist_3", nBins_X, min_X_histo_ini, max_X_histo_ini);
    
    
    
}

void Print_Table(int TP, int FP, int TN, int FN){
    cout << "              "    << " | "
         << "Actual e+  "    << " | "
         << "Actual pi+"<< "\n";
    
    cout << std::string(10*3 + 2*3, '-') << "\n";
    
    
    cout << "Predicted e^+ "   << " | "<< TP << " | "<< FP << "\n";
    
    cout << "Predicted pi^+"   << " | "<< FN << " | "<< TN << "\n";
}

void ROC_All(TString infile, TString infile2){
    
    gROOT->SetBatch(kTRUE);
    
    gStyle->SetPaintTextFormat("4.1f");
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(.03, "xyz");
    gStyle->SetTitleSize(.03, "xyz");
    gStyle->SetTitleSize(.07, "t");
    gStyle->SetFrameLineWidth(1);
    gStyle->SetLineWidth(1);
    gStyle->SetHistLineWidth(1);
    gStyle->SetMarkerStyle(13);
    gStyle->SetTitleW(0.8); // per cent of the pad width
    gStyle->SetTitleH(0.1); // per cent of the pad height
    
    TFile *Data_file = new TFile(infile+".root");
    TTree *Data_tree_P = (TTree *)Data_file->Get("resP");
    TTree *Data_tree_N = (TTree *)Data_file->Get("resN");
    
    TFile *Data_file2 = new TFile(infile2+".root");
    TTree *Data_tree2_P = (TTree *)Data_file2->Get("resP");
    TTree *Data_tree2_N = (TTree *)Data_file2->Get("resN");
    
    
    Float_t scoreBDT_P,scoreBDT_N,score,scoreMLP_P,scoreMLP_N;
    
    //---------------File 1-------------------------------
    
    Data_tree_P->SetBranchAddress("scoreBDT_P",&scoreBDT_P);
    Data_tree_P->SetBranchAddress("scoreMLP_P",&scoreMLP_P);
    Data_tree_N->SetBranchAddress("scoreBDT_N",&scoreBDT_N);
    Data_tree_N->SetBranchAddress("scoreMLP_N",&scoreMLP_N);
    
    
    //---------------File 2-------------------------------
    
    Data_tree2_P->SetBranchAddress("scoreBDT_P",&scoreBDT_P);
    Data_tree2_P->SetBranchAddress("scoreMLP_P",&scoreMLP_P);
    Data_tree2_N->SetBranchAddress("scoreBDT_N",&scoreBDT_N);
    Data_tree2_N->SetBranchAddress("scoreMLP_N",&scoreMLP_N);
    
    Float_t segments, start, end;
    Float_t  start_MLP, end_MLP;
    Float_t  start_BDT, end_BDT;
    Float_t TPR,FPR;
    Float_t score_cut;
    Int_t max=500;
    int TP=0;
    int FP=0;
    int TN=0;
    int FN=0;
    int i=1;
    Double_t x_mlp[max], y_mlp[max];
    Double_t x_bdt[max], y_bdt[max];
    
    Double_t x2_mlp[max], y2_mlp[max];
    Double_t x2_bdt[max], y2_bdt[max];
    
    
    start_MLP=-1.25;
    end_MLP=1.5;
    start_BDT=-0.8;
    end_BDT=0.8;
    
    cout<<"ROC "<<infile<<"--MLP"<<endl;
    
    start=start_MLP;
    end=end_MLP;
    segments=(end-start)/max;
    score_cut=start;
    
    while (score_cut<=end){
        TP=Data_tree_P->GetEntries(Form("scoreMLP_P>=%g",score_cut));
        TN=Data_tree_N->GetEntries(Form("scoreMLP_N<%g",score_cut));
        FP=Data_tree_N->GetEntries(Form("scoreMLP_N>=%g",score_cut));
        FN=Data_tree_P->GetEntries(Form("scoreMLP_P<%g",score_cut));
        
        TPR=(1.0*TP)/(TP+FN);
        FPR=(1.0*FP)/(FP+TN);
        
        
        x_mlp[i-1]=TPR;
        y_mlp[i-1]=1-FPR;
        
        score_cut=start+i*segments;
        i++;
    }
    
    
    i=1;
    
    cout<<"ROC "<<infile<<"--BDT"<<endl;
    
    start=start_BDT;
    end=end_BDT;
    segments=(end-start)/max;
    score_cut=start;
    
    while (score_cut<=end){
        TP=Data_tree_P->GetEntries(Form("scoreBDT_P>=%g",score_cut));
        TN=Data_tree_N->GetEntries(Form("scoreBDT_N<%g",score_cut));
        FP=Data_tree_N->GetEntries(Form("scoreBDT_N>=%g",score_cut));
        FN=Data_tree_P->GetEntries(Form("scoreBDT_P<%g",score_cut));
        
        TPR=(1.0*TP)/(TP+FN);
        FPR=(1.0*FP)/(FP+TN);
        
        
        x_bdt[i-1]=TPR;
        y_bdt[i-1]=1-FPR;
        
        score_cut=start+i*segments;
        i++;
    }
    
    //Print_Table(TP,FP,TN,FN);
    i=1;
    
    cout<<"ROC "<<infile2<<"--MLP"<<endl;
    
    start=start_MLP;
    end=end_MLP;
    segments=(end-start)/max;
    score_cut=start;
    
    while (score_cut<=end){
        TP=Data_tree2_P->GetEntries(Form("scoreMLP_P>=%g",score_cut));
        TN=Data_tree2_N->GetEntries(Form("scoreMLP_N<%g",score_cut));
        FP=Data_tree2_N->GetEntries(Form("scoreMLP_N>=%g",score_cut));
        FN=Data_tree2_P->GetEntries(Form("scoreMLP_P<%g",score_cut));
        
        TPR=(1.0*TP)/(TP+FN);
        FPR=(1.0*FP)/(FP+TN);
        
        
        x2_mlp[i-1]=TPR;
        y2_mlp[i-1]=1-FPR;
        
        score_cut=start+i*segments;
        i++;
    }
    
    //Print_Table(TP,FP,TN,FN);
    
    i=1;
    
    cout<<"ROC "<<infile2<<"--BDT"<<endl;
    
    start=start_BDT;
    end=end_BDT;
    segments=(end-start)/max;
    score_cut=start;
    
    while (score_cut<=end){
        TP=Data_tree2_P->GetEntries(Form("scoreBDT_P>=%g",score_cut));
        TN=Data_tree2_N->GetEntries(Form("scoreBDT_N<%g",score_cut));
        FP=Data_tree2_N->GetEntries(Form("scoreBDT_N>=%g",score_cut));
        FN=Data_tree2_P->GetEntries(Form("scoreBDT_P<%g",score_cut));
        
        TPR=(1.0*TP)/(TP+FN);
        FPR=(1.0*FP)/(FP+TN);
        
        
        x2_bdt[i-1]=TPR;
        y2_bdt[i-1]=1-FPR;
        
        score_cut=start+i*segments;
        i++;
    }
    //Print_Table(TP,FP,TN,FN);
    
    
   TCanvas *c1 = new TCanvas("c1","c1",1000, 800);
    
    TGraph *gr_MLP  = new TGraph(max,x_mlp,y_mlp);
    TGraph *gr_BDT  = new TGraph(max,x_bdt,y_bdt);
    TGraph *gr2_MLP  = new TGraph(max,x2_mlp,y2_mlp);
    TGraph *gr2_BDT  = new TGraph(max,x2_bdt,y2_bdt);
    
    
    cout<<"PLOT ROC--MLP AND BDT"<<endl;
    Float_t areaMLP, areaBDT;
    Float_t areaMLP2, areaBDT2;
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("ROC; True Positives Rate; 1- False Positives Rate");
    
    auto legend = new TLegend();
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    mg->GetXaxis()->SetRangeUser(0., 1.);
    mg->GetYaxis()->SetRangeUser(0., 1.01);
    
    
    gr_MLP->SetName("gr1");
    gr_MLP->SetLineColorAlpha(kBlue,0.95);
    gr_MLP->SetLineWidth(4);
    
    gr_BDT->SetName("gr2");
    gr_BDT->SetLineWidth(4);
    gr_BDT->SetLineColorAlpha(kRed, 0.45);
    
    gr2_MLP->SetName("gr3");
    gr2_MLP->SetLineColorAlpha(kViolet,0.95);
    gr2_MLP->SetLineWidth(4);
    
    gr2_BDT->SetName("gr4");
    gr2_BDT->SetLineWidth(4);
    gr2_BDT->SetLineColorAlpha(kOrange, 0.45);
    
    double x[] = { 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
    double x1[] = { 0.0,0.2,0.4,0.6,0.8,1};
    TF1 f1("f1",[&](double *x, double *){ return gr_MLP->Eval(x[0]); },0,1,0);
    TF1 f2("f2",[&](double *x1, double *){ return gr_MLP->Eval(x1[0]); },0,1,0);
    TF1 f3("f3",[&](double *x, double *){ return gr_MLP->Eval(x[0]); },0,1,0);
    TF1 f4("f4",[&](double *x, double *){ return gr_MLP->Eval(x[0]); },0,1,0);
    
    //areaMLP = f1.Integral(0.0,1.0);
    
    //TF1 f1("f",[&](double *x, double *){ return g.Eval(x[0]); },0,1);
    
    
    areaMLP=gr_MLP->Integral()+0.5;
    areaBDT=gr_BDT->Integral()+0.5;
    areaMLP2=gr2_MLP->Integral()+0.5;
    areaBDT2=gr2_BDT->Integral()+0.5;
    
    
    mg->Add(gr_MLP);
    mg->Add(gr_BDT);
    mg->Add(gr2_MLP);
    mg->Add(gr2_BDT);
    mg->Draw("AL");
    
    legend->AddEntry("gr1",Form("MLP 9 variables (%f)",areaMLP),"l");
    legend->AddEntry("gr2",Form("BDT 9 variables (%f)",areaBDT),"l");
    legend->AddEntry("gr3",Form("MLP 6 variables (%f)",areaMLP2),"l");
    legend->AddEntry("gr4",Form("BDT 6 variables (%f)",areaBDT2),"l");
    legend->Draw();
    
    
    c1->SaveAs("ROC.png");
    
    
    
}

int plot(){
    TString model="BDT";
    Float_t cut=(-0.06);//5336
    TString file="Test";
    TString file2="Test_mod";
    TString outfile="Ratios_6Var";
    Variable_Plots_TMVA(file2,model,cut, outfile);
    //True_False(file,model,cut, outfile);
    
    //ROC(file);
    //ROC_All(file,file2);
    
    gApplication->Terminate();
    return 0;
}
