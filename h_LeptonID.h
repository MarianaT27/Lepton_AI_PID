
#ifndef h_LeptonID
#define h_LeptonID

static Float_t PCAL_TMVA, ECIN_TMVA, ECOUT_TMVA;
    static Float_t m2PCAL_TMVA, m2ECIN_TMVA, m2ECOUT_TMVA;



struct VarTMVA {
    Double_t PCAL;     // lepton_sfpcal
    Double_t ECIN;     // lepton_sfecin + correction
    Double_t ECOUT;    // lepton_sfecout
    Double_t m2PCAL;   // lepton_m2pcal
    Double_t m2ECIN;   // lepton_m2ecin multiplied by correction
    Double_t m2ECOUT;  // lepton_m2ecout
};


VarTMVA getTMVAInput(Double_t lepton_sfpcal,Double_t lepton_sfecin, Double_t lepton_sfecout,Double_t lepton_m2pcal,Double_t lepton_m2ecin,Double_t lepton_m2ecout,int lepton, int version){
    VarTMVA vars;
    Double_t SF_corr, m2_corr;
    if(lepton==-11){
        if (version == -19) {
            SF_corr = 0.01;
            m2_corr = 0.8;
        } else if (version == 18) {
            SF_corr = 0.03;
            m2_corr = 1.0;
        } else {
            SF_corr = 0.01;
            m2_corr = 0.8;
        }
    }
    else{
        if (version == -19){
            SF_corr = 0.03;
            m2_corr = 0.8;
        } else if (version == 18){
            SF_corr = 0.05;
            m2_corr = 1.1;
        } else{
            SF_corr = 0.02;
            m2_corr = 0.8;
        }

    }

    vars.PCAL     = lepton_sfpcal;
    vars.ECIN     = lepton_sfecin + SF_corr;
    vars.ECOUT    = lepton_sfecout;
    vars.m2PCAL   = lepton_m2pcal;
    vars.m2ECIN   = lepton_m2ecin * m2_corr;
    vars.m2ECOUT  = lepton_m2ecout;

    return vars;
}

Double_t calculateTMVAScore(const VarTMVA &vars, TMVA::Reader* reader, const char* methodName)
{
    PCAL_TMVA   = vars.PCAL;
    ECIN_TMVA   = vars.ECIN;
    ECOUT_TMVA  = vars.ECOUT;
    m2PCAL_TMVA = vars.m2PCAL;
    m2ECIN_TMVA = vars.m2ECIN;
    m2ECOUT_TMVA= vars.m2ECOUT;
    // Evaluate the MVA using the configured method name (adjust "MyMethod" as needed).
    return reader->EvaluateMVA(methodName);
}

#endif
