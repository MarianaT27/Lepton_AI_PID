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

struct cartesian
{

    double x;
    double y;
    double z;
};

struct response
{

    cartesian pos;
    cartesian tpos;
    double time;
    double energy;
    double path;
    int sector;
    int layer;
    int index;
    int component;
    double du;
    double dv;
    double dw;
    double m2u;
    double m2v;
    double m2w;
    double m3u;
    double m3v;
    double m3w;
    double u;
    double v;
    double w;
    double widthu;
    double widthv;
    double widthw;
    double x;
    double y;
    double z;
    double quality;
    int degree;
    int detector;
};

struct particl
{

    TLorentzVector lorentz;
    cartesian vertexinfo;
    map<int, response> responses;
    int index = -1;
    double beta;
    double chi2pid;
    int status;
    int pid;
    double vtime;
    double E;
};

Double_t fitf(Double_t *x, Double_t *par)
{

    Double_t arg = 0;
    if (par[2] != 0)
        arg = (x[0] - par[1]) / par[2];
    Double_t fitval = par[0] * TMath::Exp(-0.5 * arg * arg);
    return fitval;
}

int Background_updated(string nameFile = "Rad_lepton", int version = -19, int model = 6, string train = "jpsitcs", int lepton_ID = 11)
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

    // CONFIG VARIABLES
    double Beam_E;
    if (version == -19)
        Beam_E = 10.2;
    else
        Beam_E = 10.6;
    TLorentzVector beam(0, 0, Beam_E, Beam_E);
    TLorentzVector target(0, 0, 0, 0.938);
    TLorentzVector miss;
    Int_t run_number;
    Int_t event_number;

    // lepton variables
    TLorentzVector lepton;
    Double_t lepton_vx, lepton_vy, lepton_vz, lepton_vt;
    int number_of_leptons;

    Float_t n;
    Float_t k;
    Float_t acceptance[15];

    // Float_t cuts[11]={-0.6,-0.4,-0.2, -0.1,-0.05, -0.01,0.05,0.1,0.2,0.3,0.4};
    Float_t cuts[15] = {-0.6, -0.5, -0.4, -0.3, -0.2, -0.15, -0.1, -0.06, -0.04, -0.02, 0.0, 0.05, 0.1, 0.2, 0.3};
    // Double_t cuts[7]={0,0.2,0.4, 0.5,0.6, 0.8,1};

    //********************************
    // HISTOGRAMAS
    //********************************
    TH1F *h_X = new TH1F("h_X", "Missing Mass P>4.5; MM [GeV] ; ", 100, 0.5, 1.5);
    TH1F *h_X_1 = new TH1F("h_X_1", "Missing Mass P<4.5; MM [GeV]; ", 100, 0.5, 1.5);
    TH1F *h_X_0 = new TH1F("h_X_0", "No cuts; MM [GeV]; ", 100, 0.5, 1.5);

    TH1F *h_electronR_P = new TH1F("h_electronR_P", "Real e^{-}; P; ", 100, 0, 11);
    TH1F *h_electronF_P = new TH1F("h_electronF_P", "False e^{-}; P;", 100, 0, 11);

    TH1F *h_X_cut[15] = {0};

    TH1F *h_e_P = new TH1F("h_e_P", "P", 100, 0, 14);
    TH1F *h_e_theta = new TH1F("h_e_theta", "#Theta", 100, 0, 1);
    TH1F *h_e_phi = new TH1F("h_e_phi", "#Phi", 100, -3.2, 3.2);
    TH1F *h_e_SFPCAL = new TH1F("h_e_SFPCAL", "SFPCAL", 100, 0, 0.3);
    TH1F *h_e_SFECIN = new TH1F("h_e_SFECIN", "SFECIN", 100, 0, 0.25);
    TH1F *h_e_SFECOUT = new TH1F("h_e_SFECOUT", "SFECOUT", 100, 0, 0.25);
    TH1F *h_e_m2PCAL = new TH1F("h_e_m2PCAL", "m2PCAL", 100, 0, 200);
    TH1F *h_e_m2ECIN = new TH1F("h_e_m2ECIN", "m2ECIN", 100, 0, 500);
    TH1F *h_e_m2ECOUT = new TH1F("h_e_m2ECOUT", "m2ECOUT", 100, 0, 500);

    TH1F *h_e_P_cut = new TH1F("h_e_P_cut", "P", 100, 0, 14);
    TH1F *h_e_theta_cut = new TH1F("h_e_theta_cut", "#Theta", 100, 0, 1);
    TH1F *h_e_phi_cut = new TH1F("h_e_phi_cut", "#Phi", 100, -3.2, 3.2);
    TH1F *h_e_SFPCAL_cut = new TH1F("h_e_SFPCAL_cut", "SFPCAL", 100, 0, 0.3);
    TH1F *h_e_SFECIN_cut = new TH1F("h_e_SFECIN_cut", "SFECIN", 100, 0, 0.25);
    TH1F *h_e_SFECOUT_cut = new TH1F("h_e_SFECOUT_cut", "SFECOUT", 100, 0, 0.25);
    TH1F *h_e_m2PCAL_cut = new TH1F("h_e_m2PCAL_cut", "m2PCAL", 100, 0, 200);
    TH1F *h_e_m2ECIN_cut = new TH1F("h_e_m2ECIN_cut", "m2ECIN", 100, 0, 500);
    TH1F *h_e_m2ECOUT_cut = new TH1F("h_e_m2ECOUT_cut", "m2ECOUT", 100, 0, 500);

    for (int j = 0; j < 15; j++)
    {
        ostringstream name;
        name << "Cut " << j;
        ostringstream sstr;
        sstr << "Cut " << cuts[j] << " ;Missing Mass; ";
        h_X_cut[j] = new TH1F(name.str().c_str(), sstr.str().c_str(), 100, 0.5, 1.5);
        ;
        name.str("");
        sstr.str("");
    }

    //********************************************
    //:::::::::::::::::::TMVA::::::::::::::::::::::
    //********************************************
    TMVA::Reader *readerTMVA = new TMVA::Reader("Silent");
    Double_t score;
    // Create a set of variables and declare them to the reader
    Float_t P, Theta, Phi, PCAL, ECIN, ECOUT;
    Float_t m2PCAL = -1;
    Float_t m2ECIN = -1;
    Float_t m2ECOUT = -1;
    Float_t Nphe;

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
        from = 0;
        to = 10;
    }
    if (version == -18)
    {
        if (lepton_ID == 11)
            name_configuration = "F18inneg";
        else if (lepton_ID == -11)
            name_configuration = "F18inpos";
        from = 0;
        to = 40;
    }
    if (version == +18)
    {
        if (lepton_ID == 11)
            name_configuration = "F18outneg";
        else if (lepton_ID == -11)
            name_configuration = "F18outpos";

        from = 0;
        to = 40;
    }
    weightfile = dir_weights + name_configuration + name_model;
    readerTMVA->BookMVA("BDT method", weightfile);

    // TString weightfile= "/work/clas12/mtenorio/ML_weights_pass2/S19neg/TMVAClassification_MLP.weights.xml";
    // readerTMVA->BookMVA( "MLP method", weightfile );

    // Start
    for (int fc = from; fc < to; fc++)
    { // Run 5032 to 5038 // 6616 6618. //5423 5424

        char filename1[500];

        if (version == -19)
            sprintf(filename1, "/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/%s/%s_00%d.hipo", train.c_str(), train.c_str(), IP_S19[fc]);
        else if (version == -18)
            sprintf(filename1, "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/%s/%s_00%d.hipo", train.c_str(), train.c_str(), IP_F18in[fc]);
        else
            sprintf(filename1, "/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/%s/%s_00%d.hipo", train.c_str(), train.c_str(), IP_F18out[fc]);

        // sprintf(filename1,"/volatile/clas12/osg/marianat/job_6992/output/S19_neg_mc-Lund_BG_%d-6992-%d.hipo",fc,fc);

        hipo::reader reader;
        reader.open(filename1);

        hipo::dictionary factory;

        reader.readDictionary(factory);

        // factory.show();
        hipo::structure particles;
        hipo::structure detectors;
        hipo::event event;

        hipo::bank dataPART;
        hipo::bank dataEVENT;
        hipo::bank dataHEADER;
        hipo::bank dataCALO;
        hipo::bank dataCHE;

        hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
        hipo::bank EVENT(factory.getSchema("REC::Event"));
        hipo::bank PART(factory.getSchema("REC::Particle"));
        hipo::bank CHE(factory.getSchema("REC::Cherenkov"));
        hipo::bank HEADER(factory.getSchema("RUN::config"));

        int counter = 0;
        while (reader.next() == true)
        { // Loops all events

            counter++;

            particl lepton;
            particl photon;
            particl electronR;
            particl electronF;
            particl proton;
            particl piplus;
            reader.read(event);

            event.getStructure(EVENT);
            event.getStructure(PART);
            event.getStructure(CALO);
            event.getStructure(CHE);

            int rn = 0;
            int en = 0;

            if (PART.getSize() < 1)
                continue;

            if (HEADER.getRows() == 1)
            {
                for (int i = 0; i < HEADER.getRows(); i++)
                {
                    rn = HEADER.getInt("run", i);
                    en = HEADER.getInt("event", i);
                }
            }

            int event_start_time;

            for (int i = 0; i < EVENT.getRows(); i++)
            {
                event_start_time = EVENT.getFloat("startTime", i);
            }

            number_of_leptons = 0;
            int number_of_photons = 0;
            int number_eR = 0;
            int number_eF = 0;
            int number_pip = 0;
            int number_proton = 0;
            int number_other = 0;

            int lepton_photon = 0;

            int nrows = PART.getRows();

            for (int i = 0; i < PART.getRows(); i++)
            {
                float px = PART.getFloat("px", i);
                float py = PART.getFloat("py", i);
                float pz = PART.getFloat("pz", i);
                // PID 11 Means lepton, Status Between 2000 and 4000 Means Forward Detector
                if (PART.getInt("pid", i) == 11 && abs(PART.getInt("status", i)) >= 2000 && abs(PART.getInt("status", i)) < 4000)
                {
                    number_of_leptons = number_of_leptons + 1;
                    TVector3 lepton3(px, py, pz);
                    electronR.lorentz.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.0005 * 0.0005));
                    electronR.index = i;
                    number_eR++;
                }
                else if (PART.getInt("pid", i) == -11 && abs(PART.getInt("status", i)) >= 2000 && abs(PART.getInt("status", i)) < 4000)
                {
                    number_pip = number_pip + 1;
                    piplus.lorentz.SetPxPyPzE(px, py, pz, sqrt(px * px + py * py + pz * pz + 0.13957 * 0.13957));
                    piplus.index = i;
                }
                else
                {
                    number_other++;
                }
            }

            // REC::Calorimeter
            for (int i1 = 0; i1 < CALO.getRows(); i1++)
            {
                response ecalresponse;
                ecalresponse.energy = CALO.getFloat("energy", i1);
                ecalresponse.index = CALO.getFloat("index", i1);
                ecalresponse.m2u = CALO.getFloat("m2u", i1);
                ecalresponse.m2v = CALO.getFloat("m2v", i1);
                ecalresponse.m2w = CALO.getFloat("m2w", i1);
                if (CALO.getShort("pindex", i1) == piplus.index && CALO.getByte("layer", i1) == 1) // PCAL
                    piplus.responses.insert(pair<int, response>(109, ecalresponse));
                if (CALO.getShort("pindex", i1) == piplus.index && CALO.getByte("layer", i1) == 4) // ECIN
                    piplus.responses.insert(pair<int, response>(110, ecalresponse));
                if (CALO.getShort("pindex", i1) == piplus.index && CALO.getByte("layer", i1) == 7) // ECOUT
                    piplus.responses.insert(pair<int, response>(111, ecalresponse));

                if (CALO.getShort("pindex", i1) == electronR.index && CALO.getByte("layer", i1) == 1) // PCAL
                    electronR.responses.insert(pair<int, response>(109, ecalresponse));
                if (CALO.getShort("pindex", i1) == electronR.index && CALO.getByte("layer", i1) == 4) // ECIN
                    electronR.responses.insert(pair<int, response>(110, ecalresponse));
                if (CALO.getShort("pindex", i1) == electronR.index && CALO.getByte("layer", i1) == 7) // ECOUT
                    electronR.responses.insert(pair<int, response>(111, ecalresponse));
            }
            for (int i = 0; i < CHE.getRows(); i++)
            {
                response htccresponse;
                htccresponse.detector = CHE.getByte("detector", i);
                htccresponse.energy = CHE.getFloat("nphe", i);
                htccresponse.time = CHE.getFloat("time", i);
                if (CHE.getShort("pindex", i) == electronR.index && CHE.getByte("detector", i) == 15)
                    electronR.responses.insert(pair<int, response>(6, htccresponse));
                if (CHE.getShort("pindex", i) == piplus.index && CHE.getByte("detector", i) == 15)
                    piplus.responses.insert(pair<int, response>(6, htccresponse));
            }

            if (number_of_leptons == 1 && number_pip == 1)
            {
                miss = (beam + target) - (electronR.lorentz + piplus.lorentz);
                if (piplus.lorentz.P() < 4.5)
                {
                    h_X_1->Fill(abs(miss.M()));
                    continue;
                }
                h_X->Fill(abs(miss.M()));

                // if(number_other==0){
                // miss=(beam+target)-(electronR.lorentz+piplus.lorentz);
                h_X_0->Fill(miss.M());

                P = piplus.lorentz.P();
                Theta = piplus.lorentz.Theta();
                Phi = piplus.lorentz.Phi();
                Nphe = piplus.responses[6].energy;
                PCAL = piplus.responses[109].energy / P;
                ECIN = piplus.responses[110].energy / P;
                ECOUT = piplus.responses[111].energy / P;
                m2PCAL = (piplus.responses[109].m2u + piplus.responses[109].m2v + piplus.responses[109].m2w) / 3;
                m2ECIN = (piplus.responses[110].m2u + piplus.responses[110].m2v + piplus.responses[110].m2w) / 3;
                m2ECOUT = (piplus.responses[111].m2u + piplus.responses[111].m2v + piplus.responses[111].m2w) / 3;

                score = readerTMVA->EvaluateMVA("BDT method");

                if (miss.M() > 0.86 && miss.M() < 1.05)
                {
                    h_e_P->Fill(P);
                    h_e_theta->Fill(Theta);
                    h_e_phi->Fill(Phi);
                    h_e_SFPCAL->Fill(PCAL);
                    h_e_SFECIN->Fill(ECIN);
                    h_e_SFECOUT->Fill(ECOUT);
                    h_e_m2PCAL->Fill(m2PCAL);
                    h_e_m2ECIN->Fill(m2ECIN);
                    h_e_m2ECOUT->Fill(m2ECOUT);
                }

                for (int j = 0; j < 15; j++)
                {
                    if (score > cuts[j])
                    {
                        h_X_cut[j]->Fill(miss.M());
                        if (miss.M() > 0.86 && miss.M() < 1.05)
                        {
                            if (j == 10)
                            {
                                h_e_P_cut->Fill(P);
                                h_e_theta_cut->Fill(Theta);
                                h_e_phi_cut->Fill(Phi);
                                h_e_SFPCAL_cut->Fill(PCAL);
                                h_e_SFECIN_cut->Fill(ECIN);
                                h_e_SFECOUT_cut->Fill(ECOUT);
                                h_e_m2PCAL_cut->Fill(m2PCAL);
                                h_e_m2ECIN_cut->Fill(m2ECIN);
                                h_e_m2ECOUT_cut->Fill(m2ECOUT);
                            }
                        }
                    }

                } // for j
                //}//
            }

        } // end while
        printf("processed events = %d\n", counter);
        printf("run = %d\n", fc);
    } // End for "Runs"

    string pdf_original = nameFile + "_" + train + ".pdf";
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0001);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.15);

    TCanvas *cansq = new TCanvas("cansq", "canvas", 1000, 900);

    // Define custom colors if desired (ROOT has limited built-ins, but we can define new ones)
    TColor *teal = new TColor(2001, 0.4, 0.7, 0.7);   // soft teal
    TColor *coral = new TColor(2002, 0.98, 0.6, 0.5); // soft coral

    TColor *deepTeal = new TColor(2003, 0.0, 0.5, 0.5);   // deeper teal
    TColor *warmCoral = new TColor(2004, 0.95, 0.4, 0.3); // deeper coral

    // Apply fill colors and styles
    h_X_0->SetFillColor(2003); // deeper teal
    h_X_0->SetFillStyle(1001); // solid fill

    h_X_cut[10]->SetFillColor(2004); // deeper coral
    h_X_cut[10]->SetFillStyle(1001); // solid fill

    // Set line colors slightly darker for contrast
    h_X_0->SetLineColor(2003);
    h_X_0->SetLineWidth(2);

    h_X_cut[10]->SetLineColor(2004);
    h_X_cut[10]->SetLineWidth(2);

    // Titles and axis labels
    h_X_0->SetTitle("");
    h_X_0->GetXaxis()->SetTitle("Missing Mass, GeV");
    h_X_0->GetXaxis()->SetTitleSize(0.045);
    h_X_0->GetXaxis()->SetLabelSize(0.04);
    // h_e_gamma_cut[0]->GetXaxis()->SetTitleOffset(1.2);

    h_X_0->GetYaxis()->SetTitle("");
    h_X_0->GetYaxis()->SetTitleSize(0.044);
    h_X_0->GetYaxis()->SetLabelSize(0.04);
    // h_e_gamma_cut[0]->GetYaxis()->SetTitleOffset(1.4);

    // Draw both histograms
    h_X_0->Draw("HIST");
    h_X_cut[10]->Draw("HIST SAME");

    // Legend
    TLegend *legend = new TLegend(0.65, 0.7, 0.88, 0.88);
    legend->AddEntry(h_X_0, "Before Cut", "f");
    legend->AddEntry(h_X_cut[10], "After Cut", "f");

    legend->Draw();

    // can2->Print((pdf_original + ")").c_str());

    cansq->Print((nameFile + "_" + train + ".png").c_str());

    Float_t erry[15];
    Float_t errx[15] = {0};
    Float_t erry_Classic[15];
    Float_t erry_Poisson[15];
    Float_t erry_Bayessian[15];

    Float_t error_n, err;

    float min, max;
    min = 0.6;
    max = 1.5;

    TCanvas *can2 = new TCanvas("can2", "canvas", 200, 10, 700, 500);

    can2->Divide(4, 4);
    int can_num = 1;
    can2->cd(can_num);
    h_X_0->Draw();
    TF1 *g = g43_free_bg(h_X_0, 3000, 0.93, 0.10, 3, 30, 10, min, max);
    n = g->GetParameter(0);
    error_n = g->GetParError(0);
    auto mean = g->GetParameter(1);
    auto sigma = g->GetParameter(2);
    auto a = g->GetParameter(3);
    auto b = g->GetParameter(4);
    auto c = g->GetParameter(5);
    can_num++;
    for (int i = 0; i < 15; i++)
    {
        can2->cd(can_num);
        h_X_cut[i]->Draw();
        TF1 *g = g43_bg(h_X_cut[i], n, mean, sigma, a, b, c, min, max);
        k = g->GetParameter(0);

        if (h_X_cut[i]->GetEntries() == 0 || k < 0)
        {
            acceptance[i] = 0.0;
        }
        else
        {
            acceptance[i] = (k / n);
        }
        erry_Classic[i] = acceptance[i] * sqrt(pow(g->GetParError(0) / k, 2) + pow(error_n / n, 2));
        erry_Poisson[i] = acceptance[i] * sqrt((1 / k) + (1 / n));
        err = sqrt((((k + 1) * (k + 2)) / ((n + 2) * (n + 3))) - (pow(k + 1, 2) / pow(n + 2, 2)));
        erry_Bayessian[i] = err;
        can_num++;
    }
    can2->Print((pdf_original + "(").c_str());
    can2->Clear();

    FILE *f_readme;

    f_readme = fopen("For_Ratios_BS.txt", "a");
    fprintf(f_readme, "------------------------- \n");
    fprintf(f_readme, "Model: %d Version: %d Lepton: %d \n", model, version, lepton_ID);
    fprintf(f_readme, "Initial events: %f \n", n);
    fprintf(f_readme, "Cut   Suppression Error\n");
    fprintf(f_readme, "Cuts={");
    for (int i = 0; i < 15; i++)
    {
        fprintf(f_readme, "%.2f, ", cuts[i]);
    }
    fprintf(f_readme, "} \nDouble_t y_%dbdt_F18outpos[%d]={", model, 15);
    for (int i = 0; i < 15; i++)
    {
        fprintf(f_readme, "%f, ", acceptance[i]);
    }
    fprintf(f_readme, "}; \nDouble_t yerrC_%dbdt_F18outpos[%d]={", model, 15);
    for (int i = 0; i < 15; i++)
    {
        fprintf(f_readme, "%f, ", erry_Classic[i]);
    }
    fprintf(f_readme, "}; \nDouble_t yerrP_%dbdt_F18outpos[%d]={", model, 15);
    for (int i = 0; i < 15; i++)
    {
        fprintf(f_readme, "%f, ", erry_Poisson[i]);
    }
    fprintf(f_readme, "}; \nDouble_t yerrB_%dbdt_F18outpos[%d]={", model, 15);
    for (int i = 0; i < 15; i++)
    {
        fprintf(f_readme, "%f, ", erry_Bayessian[i]);
    }
    fprintf(f_readme, "};\n------------------------- \n");
    fclose(f_readme);
    TString nameroot = nameFile + ".root";
    TFile *outputfile = new TFile(nameroot, "RECREATE");
    outputfile->cd();

    h_e_P->Write();
    h_e_theta->Write();
    h_e_phi->Write();
    h_e_SFPCAL->Write();
    h_e_SFECIN->Write();
    h_e_SFECOUT->Write();
    h_e_m2PCAL->Write();
    h_e_m2ECIN->Write();
    h_e_m2ECOUT->Write();
    h_e_P_cut->Write();
    h_e_theta_cut->Write();
    h_e_phi_cut->Write();
    h_e_SFPCAL_cut->Write();
    h_e_SFECIN_cut->Write();
    h_e_SFECOUT_cut->Write();
    h_e_m2PCAL_cut->Write();
    h_e_m2ECIN_cut->Write();
    h_e_m2ECOUT_cut->Write();
    outputfile->Close();

    cout << "---------RESULTS----------" << endl;
    cout << "Initial events: " << n << endl;
    for (int i = 0; i < 15; i++)
    {
        cout << "Cut :" << cuts[i] << " -> " << acceptance[i] << " " << erry_Bayessian[i] << endl;
    }

    TCanvas *can = new TCanvas("can", "canvas", 200, 10, 700, 700);
    TGraphErrors *gr_Acceptance = new TGraphErrors(15, cuts, acceptance, errx, erry_Bayessian);
    gr_Acceptance->SetMarkerStyle(8);
    gr_Acceptance->GetYaxis()->SetRangeUser(0, 1.05);
    // gr_Acceptance->GetAttLine(0)->SetLineColor(kRed);
    gr_Acceptance->SetTitle(" ; Cut value ;Background suppression");
    gr_Acceptance->Draw("ALEP");
    // auto legend = new TLegend(0.1, 0.50, 0.87, 0.67);
    // legend->SetFillStyle(0);
    // legend->SetLineWidth(0);
    // legend->SetTextSize(0.02);
    // legend->AddEntry(gr_Acceptance, Form("at -0.01 cut =%.5f",acceptance[8]), "p");
    // legend->Draw("same");
    can->Print((pdf_original + "(").c_str());

    can->Clear();

    gStyle->SetOptStat(0);

    if (h_X->GetEntries() != 0)
    {
        can->Divide(1, 2);
        can->cd(1);
        h_X_1->Draw();
        can->cd(2);
        h_X->Draw();
        can->Print((pdf_original + "(").c_str());
    }

    can->Clear();
    gStyle->SetOptStat(1);
    can->Divide(3, 3);
    can->cd(1);

    h_e_P->SetLineColor(kRed);
    h_e_P->SetTitle("P");
    h_e_P->Draw();
    h_e_P_cut->Draw("same");

    can->cd(2);
    h_e_theta->SetLineColor(kRed);
    h_e_theta->SetTitle("Theta");
    h_e_theta->Draw();
    h_e_theta_cut->Draw("same");

    can->cd(3);
    h_e_phi->SetLineColor(kRed);
    h_e_phi->SetTitle("Phi");
    h_e_phi->Draw();
    h_e_phi_cut->Draw("same");

    can->cd(4);
    h_e_SFPCAL->SetLineColor(kRed);
    h_e_SFPCAL->SetTitle("SFPCAL");
    h_e_SFPCAL->Draw();
    h_e_SFPCAL_cut->Draw("same");

    can->cd(5);
    h_e_SFECIN->SetLineColor(kRed);
    h_e_SFECIN->SetTitle("SFECIN");
    h_e_SFECIN->Draw();
    h_e_SFECIN_cut->Draw("same");

    can->cd(6);
    h_e_SFECOUT->SetLineColor(kRed);
    h_e_SFECOUT->SetTitle("SFECOUT");
    h_e_SFECOUT->Draw();
    h_e_SFECOUT_cut->Draw("same");

    can->cd(7);
    h_e_m2PCAL->SetLineColor(kRed);
    h_e_m2PCAL->SetTitle("m2PCAL");
    h_e_m2PCAL->Draw();
    h_e_m2PCAL_cut->Draw("same");

    can->cd(8);
    h_e_m2ECIN->SetLineColor(kRed);
    h_e_m2ECIN->SetTitle("m2ECIN");
    h_e_m2ECIN->Draw();
    h_e_m2ECIN_cut->Draw("same");

    can->cd(9);
    h_e_m2ECOUT->SetLineColor(kRed);
    h_e_m2ECOUT->SetTitle("m2ECOUT");
    h_e_m2ECOUT->Draw();
    h_e_m2ECOUT_cut->Draw("same");

    can->Print((pdf_original + ")").c_str());
    can->Clear();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    return 0;
}
