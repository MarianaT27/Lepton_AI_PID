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
#include "reader.h"
#include "clas12reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

struct cartesian {

  double x;
  double y;
  double z;

};


struct response {

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


struct particl {

  TLorentzVector lorentz;
  cartesian vertexinfo;
  map<int,response> responses;
  int index = -1;
  double beta;
  double chi2pid;
  int status;
  int pid;
  double vtime;
  double E;

};


Double_t fitf(Double_t *x,Double_t *par) {

    Double_t arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}


TF1* g43(TH1* h, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double low, double high){

   TF1* f = new TF1(((TString)"background_fit") + h->GetName(),"[0]*0.398942*(0.01)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2]) + [3]*(x-[1])*(x-[1]) - [4]*(x-[1]) + [5] + [6]*(x-[1])*(x-[1])*(x-[1])",
    low, high);
   TF1* back = new TF1(((TString)"background_fit")," [1]*(x-[0])*(x-[0]) - [2]*(x-[0])+ [3]+[4]*(x-[0])*(x-[0])*(x-[0])",
    low, high);
   TF1* gauss = new TF1(((TString)"peak"),"[0]*0.398942*(0.01)*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])",
    low, high);

   f->SetParameters(x1, x2, x3, x4, x5, x6, x7);
   f->SetParNames("X1", "X2", "X3", "X4","X5", "X6", "X7");
   h->Fit(f, "MEL", "same", low, high);

   double par[6];
   f->GetParameters(par);
   back->SetParameters(par[1], par[3], par[4], par[5], par[6]);
   gauss->SetParameters(par[0], par[1], par[2]);

   back->SetLineColor(kBlue);
   gauss->SetLineColor(kGreen);

   back->Draw("same");
   gauss->Draw("same");

   return f;
}




int Rad_electron(string nameFile="Rad_lepton", int version=-19, int model=6, string train="jpsitcs", int lepton_ID=11) {
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
    //DECLARATION OF VARIABLES
    //********************************

    //CONFIG VARIABLES
    //TLorentzVector beam(0,0,Beam_E,Beam_E);
    TLorentzVector target(0,0,0,0.938);
    Int_t run_number;
    Int_t event_number;

    //lepton variables
    TLorentzVector lepton;
    Double_t lepton_vx, lepton_vy, lepton_vz, lepton_vt;
    int number_of_leptons;

    Float_t cuts[11]={-0.6,-0.4,-0.2, -0.1,-0.05, -0.01,0.05,0.1,0.2,0.3,0.4};
    //Double_t cuts[7]={0,0.2,0.4, 0.5,0.6, 0.8,1};


    //********************************
    //HISTOGRAMAS
    //********************************
    TH1F* h_e_gamma = new TH1F("h_e_gamma","e -> e#gamma; #Delta #phi, degrees; ",85,-4,4);

    TH1F* h_e_gamma_cut[11];

    TH1F* h_e_P=       new TH1F("h_e_P","No Cut; P;",100,0,11);
    TH1F* h_e_theta =  new TH1F("h_e_theta","No Cut; #Theta;",100,0,1);
    TH1F* h_e_phi=     new TH1F("h_e_phi","No Cut; #Phi;",100,-3,3);
    TH1F* h_e_SFPCAL = new TH1F("h_e_SFPCAL","No Cut; SFPCAL;",150,0,0.3);
    TH1F* h_e_SFECIN = new TH1F("h_e_SFECIN","No Cut; SFECIN;",150,0,0.2);
    TH1F* h_e_SFECOUT =new TH1F("h_e_SFECOUT","No Cut; SFECOUT;",150,0,0.05);
    TH1F* h_e_m2PCAL = new TH1F("h_e_m2PCAL","No Cut; m2PCAL;",150,0,80);
    TH1F* h_e_m2ECIN = new TH1F("h_e_m2ECIN","No Cut; m2ECIN;",150,0,300);
    TH1F* h_e_m2ECOUT= new TH1F("h_e_m2ECOUT","No Cut; m2ECOUT;",150,0,130);

    TH1F* h_e_P_cut=       new TH1F("h_e_P_cut","-0.01 Cut; P;",100,0,11);
    TH1F* h_e_theta_cut =  new TH1F("h_e_theta_cut","-0.01 Cut; #Theta;",100,0,1);
    TH1F* h_e_phi_cut=     new TH1F("h_e_phi_cut","-0.01 Cut; #Phi;",100,-3,3);
    TH1F* h_e_SFPCAL_cut = new TH1F("h_e_SFPCAL_cut","-0.01 Cut; SFPCAL;",150,0,0.3);
    TH1F* h_e_SFECIN_cut = new TH1F("h_e_SFECIN_cut","-0.01 Cut; SFECIN;",150,0,0.2);
    TH1F* h_e_SFECOUT_cut =new TH1F("h_e_SFECOUT_cut","-0.01 Cut; SFECOUT;",150,0,0.05);
    TH1F* h_e_m2PCAL_cut = new TH1F("h_e_m2PCAL_cut","-0.01 Cut; m2PCAL;",150,0,80);
    TH1F* h_e_m2ECIN_cut = new TH1F("h_e_m2ECIN_cut","-0.01 Cut; m2ECIN;",150,0,300);
    TH1F* h_e_m2ECOUT_cut= new TH1F("h_e_m2ECOUT_cut","-0.01 Cut; m2ECOUT;",150,0,130);


    for(int j = 0; j < 11; j++){
      ostringstream name;
      name << "Cut " << j;
      ostringstream sstr;
      sstr<<"Cut "<< cuts[j]<<" ;#Delta #phi, degrees; ";
      h_e_gamma_cut[j] = new TH1F(name.str().c_str(), sstr.str().c_str(), 85, -4,4);;
      name.str("");
      sstr.str("");
      
    }

    //********************************************
    //:::::::::::::::::::TMVA::::::::::::::::::::::
    //********************************************
     TMVA::Reader *readerTMVA = new TMVA::Reader( "!Color:!Silent" );
     Double_t score;
    // Create a set of variables and declare them to the reader
     Float_t P, Theta, Phi, PCAL,ECIN,ECOUT, m2PCAL, m2ECIN, m2ECOUT,Nphe;

     if(model==9){
        readerTMVA->AddVariable( "P",&P );
        readerTMVA->AddVariable( "Theta",&Theta);
        readerTMVA->AddVariable( "Phi",&Phi);
     }
     //readerTMVA->AddVariable( "Nphe",&Nphe);
     readerTMVA->AddVariable( "SFPCAL",&PCAL);
     readerTMVA->AddVariable( "SFECIN",&ECIN);
    readerTMVA->AddVariable( "SFECOUT",&ECOUT );
    readerTMVA->AddVariable( "m2PCAL",&m2PCAL);
    readerTMVA->AddVariable( "m2ECIN",&m2ECIN);
    readerTMVA->AddVariable( "m2ECOUT",&m2ECOUT);

    //Book Methods
    TString weightfile;
    int from, to;
    if(version==-19){
        if(lepton_ID==11){
            if(model==6)
                weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/S19neg/TMVAClassification_BDT_6.weights.xml";
            if(model==9)
                weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/S19neg/TMVAClassification_BDT.weights.xml";
        }
        else if(lepton_ID==-11){
            if(model==6)
                weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/S19pos/TMVAClassification_BDT_6.weights.xml";
            if(model==9)
                 weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/S19pos/TMVAClassification_BDT.weights.xml";
        }
        from=6616;
        to=6620;
    }
    if(version==-18){
        if(lepton_ID==11){
            if(model==6)
                weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/F18inneg/TMVAClassification_BDT_6.weights.xml";
            if(model==9)
                weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/F18inneg/TMVAClassification_BDT.weights.xml";
        }
        else if(lepton_ID==-11){
            if(model==6)
                weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/F18inpos/TMVAClassification_BDT_6.weights.xml";
            if(model==9)
                weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/F18inpos/TMVAClassification_BDT.weights.xml";
        }
        from=5032;
        to=5038;
    }
    if(version==+18){
        if(lepton_ID==11){
            if(model==6)
                weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/F18outneg/TMVAClassification_BDT_6.weights.xml";
            if(model==9)
                weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/F18outneg/TMVAClassification_BDT.weights.xml";
        }
        else if(lepton_ID==-11){
            if(model==6)
                weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/F18outpos/TMVAClassification_BDT_6.weights.xml";
            if(model==9)
                weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/F18outpos/TMVAClassification_BDT.weights.xml";
        }

        from=5423;
        to=5426;
    }

    readerTMVA->BookMVA( "BDT method", weightfile );


    //TString weightfile= "/lustre19/expphy/volatile/clas12/mtenorio/weights/S19neg/TMVAClassification_MLP.weights.xml";
    //readerTMVA->BookMVA( "MLP method", weightfile );


    //Start
    for(int fc =from; fc<=to; fc++) {//Run 5032 to 5038 // 6616 6618. //5423 5424

        char filename1[500];

        if(version==-19)
            sprintf(filename1,"/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/%s/%s_00%d.hipo",train.c_str(),train.c_str(),fc);
        else if(version ==-18)
            sprintf(filename1,"/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/%s/%s_00%d.hipo",train.c_str(),train.c_str(),fc);
        else
           sprintf(filename1,"/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass2/train/%s/%s_00%d.hipo",train.c_str(),train.c_str(),fc);

        //sprintf(filename1,"/volatile/clas12/osg/marianat/job_6992/output/S19_neg_mc-Lund_BG_%d-6992-%d.hipo",fc,fc);




       hipo::reader  reader;
       reader.open(filename1);

       hipo::dictionary  factory;

       reader.readDictionary(factory);

       factory.show();
       hipo::structure  particles;
       hipo::structure  detectors;
       hipo::event      event;

       hipo::bank  dataPART;
       hipo::bank  dataEVENT;
       hipo::bank  dataHEADER;
       hipo::bank  dataCALO;
       hipo::bank  dataCHE;


       hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
       hipo::bank EVENT(factory.getSchema("REC::Event"));
       hipo::bank PART(factory.getSchema("REC::Particle"));
       hipo::bank CHE(factory.getSchema("REC::Cherenkov"));
       hipo::bank HEADER(factory.getSchema("RUN::config"));

       int counter = 0;
        while(reader.next()==true ){//Loops all events

            counter++;

            particl lepton;
            particl photon;
            reader.read(event);

            event.getStructure(EVENT);
            event.getStructure(PART);
            event.getStructure(CALO);
            event.getStructure(CHE);

            int rn = 0;
            int en = 0;


            if(PART.getSize()<1) 
              continue;


          if(HEADER.getRows()==1) {
            for(int i = 0; i < HEADER.getRows(); i++) {
                rn = HEADER.getInt("run",i);
                en = HEADER.getInt("event",i);
            }
        }



        int event_start_time;

        for(int i = 0; i < EVENT.getRows(); i++) {
            event_start_time = EVENT.getFloat("startTime",i);

        }

        number_of_leptons = 0;
        int number_of_photons = 0;

        int lepton_photon=0;

        int nrows = PART.getRows();


        for(int i = 0; i < PART.getRows(); i++){
            float  px = PART.getFloat("px",i);
            float  py = PART.getFloat("py",i);
            float  pz = PART.getFloat("pz",i);
                //PID 11 Means lepton, Status Between 2000 and 4000 Means Forward Detector
            if(PART.getInt("pid",i)==lepton_ID && abs(PART.getInt("status",i))>=2000 && abs(PART.getInt("status",i))<4000) {
                number_of_leptons = number_of_leptons + 1;
                lepton.lorentz.SetPxPyPzE(px, py, pz, sqrt(px*px +py*py + pz*pz + 0.0005*0.0005));
                lepton.index = i;
            }
        }

            //REC::Calorimeter
        for (int i1 = 0; i1 < CALO.getRows(); i1++) {
            response ecalresponse;
            ecalresponse.energy = CALO.getFloat("energy",i1);
            ecalresponse.index = CALO.getFloat("index",i1);
            ecalresponse.m2u = CALO.getFloat("m2u",i1);
            ecalresponse.m2v = CALO.getFloat("m2v",i1);
            ecalresponse.m2w = CALO.getFloat("m2w",i1);
                if(CALO.getShort("pindex",i1)==lepton.index && CALO.getByte("layer",i1)==1) //PCAL
                  lepton.responses.insert(pair<int,response>(109,ecalresponse));
                if(CALO.getShort("pindex",i1)==lepton.index && CALO.getByte("layer",i1)==4)//ECIN
                  lepton.responses.insert(pair<int, response>(110,ecalresponse));
                if(CALO.getShort("pindex",i1)==lepton.index && CALO.getByte("layer",i1)==7)//ECOUT
                  lepton.responses.insert(pair<int, response>(111,ecalresponse));
          }
          for(int i = 0; i < CHE.getRows(); i++){
            response htccresponse;
            htccresponse.detector= CHE.getByte("detector",i);
            htccresponse.energy = CHE.getFloat("nphe",i);
            htccresponse.time = CHE.getFloat("time",i);
            if(CHE.getShort("pindex",i)==lepton.index && CHE.getByte("detector",i)==15) 
                lepton.responses.insert(pair<int, response>(6,htccresponse));
        }



        if(number_of_leptons==1 &&lepton.lorentz.P()>4.5 ){
            //if(lepton_ID==-11&&lepton.lorentz.P()<4.0)
              //  continue;


            P=lepton.lorentz.P();
            Theta=lepton.lorentz.Theta();
            Phi=lepton.lorentz.Phi();
            Nphe=lepton.responses[6].energy;
            PCAL=lepton.responses[109].energy/P;
            ECIN=lepton.responses[110].energy/P;
            ECOUT=lepton.responses[111].energy/P;
            m2PCAL=(lepton.responses[109].m2u+lepton.responses[109].m2v+lepton.responses[109].m2w)/3;
            m2ECIN=(lepton.responses[110].m2u+lepton.responses[110].m2v+lepton.responses[110].m2w)/3;
            m2ECOUT=(lepton.responses[111].m2u+lepton.responses[111].m2v+lepton.responses[111].m2w)/3;

            score=readerTMVA->EvaluateMVA("BDT method");




            for(int i = 0; i < PART.getRows(); i++){
                int   pid = PART.getInt("pid",i);
                int charge = PART.getByte("charge",i);
                float  px = PART.getFloat("px",i);
                float  py = PART.getFloat("py",i);
                float  pz = PART.getFloat("pz",i);
                int    status = PART.getInt("status",i);
                float beta = PART.getFloat("beta",i);
                if(pid==22 && abs(status)>=2000 && abs(status)<4000 ) { 
                    TLorentzVector photon_vector;
                    number_of_photons=number_of_photons+1;
                    photon_vector.SetPxPyPzE(px, py, pz, sqrt(px*px +py*py + pz*pz + 0.0*0.0));
                    photon.lorentz=photon_vector;
                //Plot e peak without cut
                    h_e_gamma->Fill(photon.lorentz.Theta()*57.2958-lepton.lorentz.Theta()*57.2958);
                    if(abs(photon.lorentz.Theta()*57.2958-lepton.lorentz.Theta()*57.2958)<1 ){
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

                    for(int j=0;j<11;j++){
                        if(score>cuts[j]){
                            h_e_gamma_cut[j]->Fill(photon.lorentz.Theta()*57.2958-lepton.lorentz.Theta()*57.2958);
                            if(abs(photon.lorentz.Theta()*57.2958-lepton.lorentz.Theta()*57.2958)<1 && j==5){
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


                }//if
            }//for part
        }


        }//end while
        printf("processed events = %d\n",counter);
        printf("run = %d\n",fc);
    }//End for "Runs"


    string pdf_original=nameFile+"_"+train+".pdf";
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0001);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.15);


    int e_Begin,e_Cut;
    Float_t acceptance[11];

    float min, max;
    if(train=="jpsitcs"){
        min=-2.5;
        max=2.5;
    }
    else{
        min=-3;//befoe -2 to 2
        max=3;
    }

    TCanvas *can2 = new TCanvas("can2","canvas",200,10,700,500);
    can2->Divide(4,3);

    can2->cd(1);
    h_e_gamma->Draw();
    TF1 *g= g43(h_e_gamma,1E5,0.0,0.5,0.5,30,500,0.5,min,max);
    e_Begin=g->GetParameter(0);
    //can->Print( (pdf_original+"(").c_str());

    for(int i=0;i<11;i++){
        can2->cd(i+2);
        h_e_gamma_cut[i]->Draw();

        TF1 *g= g43(h_e_gamma_cut[i],1E5,0.0,0.5,0.5,30,500,0.5,min,max);//e-e+ untagged
        e_Cut=g->GetParameter(0);

        if(h_e_gamma_cut[i]->GetEntries()==0||e_Cut<0)
            acceptance[i]=0.0;
        else
            acceptance[i]=(e_Cut*1.0/e_Begin*1.0);
        //if(i<10)
            //can->Print( (pdf_original+"(").c_str());
        //can->Clear();

    }
    can2->Print( (pdf_original+"(").c_str());
    can2->Clear();


    TCanvas *can = new TCanvas("can","canvas",200,10,700,700);
    TGraph *gr_Acceptance = new TGraph(11,cuts,acceptance);
    gr_Acceptance->SetMarkerStyle(21);
    gr_Acceptance->GetYaxis()->SetRangeUser(0, 1.1);
    gr_Acceptance->SetTitle("Signal Efficiency; Efficiency ; Cut value");
    gr_Acceptance->Draw("ALP");
    auto legend = new TLegend(0.1, 0.50, 0.87, 0.67);
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    legend->SetTextSize(0.02);
    legend->AddEntry(gr_Acceptance, Form("Efficiency at -0.01 cut =%.5f",acceptance[5]), "l");
    legend->Draw("same");
    can->Print( (pdf_original+"(").c_str());

    can->Clear();

    gStyle->SetOptStat(0);

    can->Divide(3,3);
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

    can->Print( (pdf_original+")").c_str());
    can->Clear();


    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<<" s\n";

    return 0;

}



