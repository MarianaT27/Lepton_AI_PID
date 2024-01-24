

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

void TMVAClassificationApplication( TString myMethodList = "" )
{
    
    gStyle->SetPalette(55);
    // This loads the library
    TMVA::Tools::Instance();
    
    
    std::cout << std::endl;
    std::cout << "==> Start TMVAClassificationApplication" << std::endl;
    
    
    // Create the Reader object
    
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
    
    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    Float_t P, Theta, Phi, PCAL,ECIN,ECOUT, m2PCAL, m2ECIN, m2ECOUT;
    
    //reader->AddVariable( "P",&P );
    //reader->AddVariable( "Theta",&Theta);
    //reader->AddVariable( "Phi",&Phi);
    reader->AddVariable( "SFPCAL",&PCAL);
    reader->AddVariable( "SFECIN",&ECIN);
    reader->AddVariable( "SFECOUT",&ECOUT );
    reader->AddVariable( "m2PCAL",&m2PCAL);
    reader->AddVariable( "m2ECIN",&m2ECIN);
    reader->AddVariable( "m2ECOUT",&m2ECOUT);
    
    //Create output root file
    TFile *target  = new TFile( "Test_mod.root","RECREATE" );
    TTree *resP = new TTree("resP","Test_mod.root");
    TTree *resN = new TTree("resN","Test_mod.root");
    
    Float_t P_P, Theta_P, Phi_P,SFPCAL_P,SFECIN_P,SFECOUT_P, m2PCAL_P, m2ECIN_P, m2ECOUT_P;
    Float_t P_N, Theta_N, Phi_N,SFPCAL_N,SFECIN_N,SFECOUT_N, m2PCAL_N, m2ECIN_N, m2ECOUT_N;

    Float_t scoreBDT_P,scoreMLP_P,scoreBDT_N,scoreMLP_N;

    resP->Branch("scoreBDT_P", &scoreBDT_P, "scoreBDT_P/F");
    resP->Branch("scoreMLP_P", &scoreMLP_P, "scoreMLP_P/F");
    
    resN->Branch("scoreBDT_N", &scoreBDT_N, "scoreBDT_N/F");
    resN->Branch("scoreMLP_N", &scoreMLP_N, "scoreMLP_N/F");
    
    resP->Branch("P_P", &P_P, "P_P/F");
    resP->Branch("Theta_P", &Theta_P, "Theta_P/F");
    resP->Branch("Phi_P", &Phi_P, "Phi_P/F");
    resP->Branch("SFPCAL_P", &SFPCAL_P, "SFPCAL_P/F");
    resP->Branch("SFECIN_P", &SFECIN_P, "SFECIN_P/F");
    resP->Branch("SFECOUT_P", &SFECOUT_P, "SFECOUT_P/F");
    resP->Branch("m2PCAL_P", &m2PCAL_P, "m2PCAL_P/F");
    resP->Branch("m2ECIN_P", &m2ECIN_P, "m2ECIN_P/F");
    resP->Branch("m2ECOUT_P", &m2ECOUT_P, "m2ECOUT_P/F");
    
    resN->Branch("P_N", &P_N, "P_N/F");
    resN->Branch("Theta_N", &Theta_N, "Theta_N/F");
    resN->Branch("Phi_N", &Phi_N, "Phi_N/F");
    resN->Branch("SFPCAL_N", &SFPCAL_N, "SFPCAL_N/F");
    resN->Branch("SFECIN_N", &SFECIN_N, "SFECIN_N/F");
    resN->Branch("SFECOUT_N", &SFECOUT_N, "SFECOUT_N/F");
    resN->Branch("m2PCAL_N", &m2PCAL_N, "m2PCAL_N/F");
    resN->Branch("m2ECIN_N", &m2ECIN_N, "m2ECIN_N/F");
    resN->Branch("m2ECOUT_N", &m2ECOUT_N, "m2ECOUT_N/F");
    
    
   
    
    // Book the MVA methods
    
    TString name="S19_positives_mod";
    
    TString dir    = "dataset_"+name+"/weights/";
    TString prefix = "TMVAClassification";
    
    //Book Methods
    TString weightfile = dir + prefix +"_MLP" +".weights.xml";
    reader->BookMVA( "MLP method", weightfile );
    
    weightfile = dir + prefix +"_BDT" +".weights.xml";
    reader->BookMVA( "BDT method", weightfile );
    

    
    //Change name files
    TFile *input1 = new TFile("Test_pos_Lepton.root");//e+ or e-
    //TFile *input2 = new TFile("outFile_Pion_positron.root");//pi+ or pi-
    //TFile *input2 = new TFile("Documents/Lepton_ID_root/"+name+"_Pion.root");  //Pions
    
    std::cout << "--- TMVAClassificationApp    : Using input file: " << input1->GetName() << std::endl;
    
    std::cout << "--- Select signal sample" << std::endl;
    TTree* theTree = (TTree*)input1->Get("tree");
    //Float_t P, Theta, Phi, PCAL,ECIN,ECOUT, m2PCAL, m2ECIN, m2ECOUT;
    theTree->SetBranchAddress( "P",&P );
    theTree->SetBranchAddress( "Theta",&Theta);
    theTree->SetBranchAddress( "Phi",&Phi);
    theTree->SetBranchAddress( "SFPCAL",&PCAL);
    theTree->SetBranchAddress( "SFECIN",&ECIN);
    theTree->SetBranchAddress( "SFECOUT",&ECOUT );
    theTree->SetBranchAddress( "m2PCAL",&m2PCAL);
    theTree->SetBranchAddress( "m2ECIN",&m2ECIN);
    theTree->SetBranchAddress( "m2ECOUT",&m2ECOUT);
    
    std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests
    
    std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();
    for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
        
        if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
        
        theTree->GetEntry(ievt);
        P_P=P;
        Theta_P=Theta;
        Phi_P=Phi;
        SFPCAL_P=PCAL;
        SFECIN_P=ECIN;
        SFECOUT_P=ECOUT;
        m2PCAL_P=m2PCAL;
        m2ECIN_P=m2ECIN;
        m2ECOUT_P=m2ECOUT;
        
        
        scoreBDT_P=reader->EvaluateMVA("BDT method");
        scoreMLP_P=reader->EvaluateMVA("MLP method");
       
        resP->Fill();
    }
    
    input1->Close();
    
    //Change name files
    //TFile *input1 = new TFile("outFile_Lepton_positron.root");//e+ or e-
    TFile *input2 = new TFile("Test_pos_Pion.root");//pi+ or pi-

    std::cout << "--- TMVAClassificationApp    : Using input file: " << input2->GetName() << std::endl;
    
    std::cout << "--- Select signal sample" << std::endl;
    TTree* theTree2 = (TTree*)input2->Get("tree");
    //Float_t P, Theta, Phi, PCAL,ECIN,ECOUT, m2PCAL, m2ECIN, m2ECOUT;
    theTree2->SetBranchAddress( "P",&P );
    theTree2->SetBranchAddress( "Theta",&Theta);
    theTree2->SetBranchAddress( "Phi",&Phi);
    theTree2->SetBranchAddress( "SFPCAL",&PCAL);
    theTree2->SetBranchAddress( "SFECIN",&ECIN);
    theTree2->SetBranchAddress( "SFECOUT",&ECOUT );
    theTree2->SetBranchAddress( "m2PCAL",&m2PCAL);
    theTree2->SetBranchAddress( "m2ECIN",&m2ECIN);
    theTree2->SetBranchAddress( "m2ECOUT",&m2ECOUT);
    
    //std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests
    
    std::cout << "--- Processing: " << theTree2->GetEntries() << " events" << std::endl;
    for (Long64_t ievt=0; ievt<theTree2->GetEntries();ievt++) {
        
        if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
        
        theTree2->GetEntry(ievt);
        P_N=P;
        Theta_N=Theta;
        Phi_N=Phi;
        SFPCAL_N=PCAL;
        SFECIN_N=ECIN;
        SFECOUT_N=ECOUT;
        m2PCAL_N=m2PCAL;
        m2ECIN_N=m2ECIN;
        m2ECOUT_N=m2ECOUT;
        
        
        scoreBDT_N=reader->EvaluateMVA("BDT method");
        scoreMLP_N=reader->EvaluateMVA("MLP method");
       
        resN->Fill();
    }
    
    
    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();
    
    target->Write();
    
    target->Close();
    
    std::cout << "--- Created root file: \"test.root\" containing the MVA output histograms" << std::endl;
    
    delete reader;
    
    std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}

int main( int argc, char** argv )
{
    TString methodList;
    for (int i=1; i<argc; i++) {
        TString regMethod(argv[i]);
        if(regMethod=="-b" || regMethod=="--batch") continue;
        if (!methodList.IsNull()) methodList += TString(",");
        methodList += regMethod;
    }
    TMVAClassificationApplication(methodList);
    return 0;
}
