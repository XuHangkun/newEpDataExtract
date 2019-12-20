#ifndef MYFUNCTION
#define MYFUNCTION
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
class myFunction
{
    public:
        myFunction();
        ~myFunction();
        //double Gaus(double *x,double *par);

        double MCShape(double *x,double *par);


        /*
    private:
        double m_mcAmplitude;
        double m_mcEnergyScale;
        double m_mcResolution;
        double m_coAmplitude;
        double m_coEnergyScale;
        TGraph *m_gausShape;
        TGraph *m_lossShape;
        TGraph *m_highShape;
        TFile* f1;
        TH1F* lossH;
        TH1F* highH;
        TH1F* gausH;*/
};
#endif
