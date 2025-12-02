#include <cstdlib>
#include <iostream>
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
//These are the libraries to have for using the TMVA methods
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "h_LeptonID.h"

using namespace TMVA;

int ExampleGetScore(){
    //Adding root file of Fall 2018 inbending, it's important to keep track of the "version" which is -18 (+ for outbending - for inbending, the number is the year)
    TString root_file = "/work/clas12/mtenorio/Analysis/Latest/F18in_All.root";
    TChain *tree = new TChain("analysis"); // "analysis" is the tree name in the files
    tree->Add(root_file);
    int version=-18;

    tree->SetMakeClass(1);
    //For lepton ID
    Double_t electron_m2pcal, electron_m2ecin, electron_m2ecout;
    Double_t electron_sfpcal, electron_sfecin, electron_sfecout;
    Double_t positron_m2pcal, positron_m2ecin, positron_m2ecout;
    Double_t positron_sfpcal, positron_sfecin, positron_sfecout;
    int number_of_electrons, number_of_positrons;

    tree->SetBranchAddress("electron_m2ecin", &electron_m2ecin);
    tree->SetBranchAddress("electron_m2ecout", &electron_m2ecout);
    tree->SetBranchAddress("electron_m2pcal", &electron_m2pcal);
    tree->SetBranchAddress("electron_sfecin", &electron_sfecin);
    tree->SetBranchAddress("electron_sfpcal", &electron_sfpcal);
    tree->SetBranchAddress("electron_sfecout", &electron_sfecout);
    tree->SetBranchAddress("positron_m2ecin", &positron_m2ecin);
    tree->SetBranchAddress("positron_m2ecout", &positron_m2ecout);
    tree->SetBranchAddress("positron_m2pcal", &positron_m2pcal);
    tree->SetBranchAddress("positron_sfecin", &positron_sfecin);
    tree->SetBranchAddress("positron_sfpcal", &positron_sfpcal);
    tree->SetBranchAddress("positron_sfecout", &positron_sfecout);

    tree->SetBranchAddress("number_of_positrons", &number_of_positrons);
    tree->SetBranchAddress("number_of_electrons", &number_of_electrons);


    //Add TMVA reader
    TMVA::Reader *readerTMVA = new TMVA::Reader("!Color:Silent");

    // Create a set of variables and declare them to the reader
    // TMVA variables for binding (different names, the first is the name of the variables in the weights file, the second is the name of the variables in the .h file)
    readerTMVA->AddVariable("SFPCAL", &PCAL_TMVA);
    readerTMVA->AddVariable("SFECIN", &ECIN_TMVA);
    readerTMVA->AddVariable("SFECOUT", &ECOUT_TMVA);
    readerTMVA->AddVariable("m2PCAL", &m2PCAL_TMVA);
    readerTMVA->AddVariable("m2ECIN", &m2ECIN_TMVA);
    readerTMVA->AddVariable("m2ECOUT", &m2ECOUT_TMVA);

    // Book Methods
    //You can change the location, all weight files are on /work/clas12/mtenorio/   
    readerTMVA->BookMVA("BDT_pos_F18in", "/work/clas12/mtenorio/ML_weights_pass2/F18inpos/TMVAClassification_BDT_6.weights.xml");
    readerTMVA->BookMVA("BDT_ele_F18in", "/work/clas12/mtenorio/ML_weights_pass2/F18inneg/TMVAClassification_BDT_6.weights.xml");
    readerTMVA->BookMVA("BDT_pos_F18out", "/work/clas12/mtenorio/ML_weights_pass2/F18outpos/TMVAClassification_BDT_6.weights.xml");
    readerTMVA->BookMVA("BDT_ele_F18out", "/work/clas12/mtenorio/ML_weights_pass2/F18outneg/TMVAClassification_BDT_6.weights.xml");
    readerTMVA->BookMVA("BDT_pos_S19", "/work/clas12/mtenorio/ML_weights_pass2/S19pos/TMVAClassification_BDT_6.weights.xml");
    readerTMVA->BookMVA("BDT_ele_S19", "/work/clas12/mtenorio/ML_weights_pass2/S19neg/TMVAClassification_BDT_6.weights.xml");
    int passedevents=0;
    int percentageStep = 5;
    int step = tree->GetEntries() * percentageStep / 100;

    for (Long64_t fc = 0; fc < tree->GetEntries(); fc++)
    { 
        tree->GetEntry(fc);
           
       if (fc % step == 0)
        {
            double percentage = (fc * 100.0) / tree->GetEntries();
            std::cout << "Progress: " << percentage << "%" << std::endl;
            //if(percentage>5)
              //break;
        }
        
        //Selectec 1 e- and 1 e+ 
        if(number_of_positrons != 1 || number_of_electrons != 1)
            continue;

      //Initialize the variables
      double score_pos=1;
      double score_ele=1;

      //Select the appropiate method depending on the version
      auto method_ele = (version==-19 ? "BDT_ele_S19" :(version==-18 ||version==-28) ? "BDT_ele_F18in" : "BDT_ele_F18out");
      auto method_pos = (version==-19 ? "BDT_pos_S19" :(version==-18 ||version==-28) ? "BDT_pos_F18in" : "BDT_pos_F18out");

      //Set the variables with the correction for positron ID
      VarTMVA varpos = getTMVAInput(positron_sfpcal, positron_sfecin, positron_sfecout,
                                                    positron_m2pcal, positron_m2ecin, positron_m2ecout,
                                                    -11, version);
      //Calculate positron score
      score_pos = calculateTMVAScore(varpos, readerTMVA, method_pos);

       //Set the variables with the correction for electron ID
        VarTMVA varele = getTMVAInput(electron_sfpcal, electron_sfecin, electron_sfecout,
                                                    electron_m2pcal, electron_m2ecin, electron_m2ecout,
                                                    11, version);
       //Calculate electron score
        score_ele = calculateTMVAScore(varele, readerTMVA, method_ele);

        if(score_ele<0.0||score_pos<0.05)
            continue;
        else
          passedevents++;
  
    }
    cout<<"Number of events that passed the score cut: "<<passedevents<<endl;

return 0;

}
