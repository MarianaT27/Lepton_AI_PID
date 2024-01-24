

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/MethodCategory.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int TMVAClassification( TString myMethodList = "" )
{

   gStyle->SetPalette(55);	
   TMVA::Tools::Instance();


   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;
    
    TString name="S19_negatives";

   //Change name files
   TFile *input1 = new TFile("Documents/Lepton_ID_root/"+name+"_Lepton.root");//e+ or e-
   TFile *input2 = new TFile("Documents/Lepton_ID_root/"+name+"_Pion.root");  //Pions
  
   TTree *signalTree     = (TTree*)input1->Get("tree");
   TTree *background     = (TTree*)input2->Get("tree");

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "Documents/"+name+"_TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset_"+name+"_mod");
   

   dataloader->AddVariable( "P","P","", 'F' );
   dataloader->AddVariable( "Theta","Theta","", 'F' );
   dataloader->AddVariable( "Phi","Phi","", 'F' );
   dataloader->AddVariable( "SFPCAL","PCAL","", 'F' );
   dataloader->AddVariable( "SFECIN","ECIN", "", 'F' );
   dataloader->AddVariable( "SFECOUT","ECOUT", "", 'F' );
   dataloader->AddVariable( "m2PCAL","m2PCAL","", 'F' );
   dataloader->AddVariable( "m2ECIN","m2ECIN", "", 'F' );
   dataloader->AddVariable( "m2ECOUT","m2ECOUT", "", 'F' );


   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree    ( signalTree,     signalWeight );
   dataloader->AddBackgroundTree( background, backgroundWeight );

   
   
   dataloader->PrepareTrainingAndTestTree( "", "",
   					"SplitMode=Random:NormMode=NumEvents:!V" );
   


   //"MLP"
   factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   
   //"BDT"  // Adaptive Boost
   factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

   // "SVM" Support Vector Machine
   //factory->BookMethod( dataloader, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001" );
    
    //factory->BookMethod(dataloader, TMVA::Types::kFisher, "Fisher","H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=60:NsmoothMVAPdf=10" );
    
    //factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "Likelihood", "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

  
   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();

   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;
   delete dataloader;
   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

   return 0;
}

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   return TMVAClassification(methodList);
}
