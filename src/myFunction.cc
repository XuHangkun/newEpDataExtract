#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "myFunction.h"



myFunction::myFunction()
{
        /*
        m_gausShape = new TGraph();
        m_lossShape = new TGraph();
        m_highShape = new TGraph();
        f1=new TFile("../MCShape/co60_specs_guwq.root");
        lossH = (TH1F*)f1->Get  ("specLoss");
        gausH = (TH1F*)f1->Get  ("specGaus");
        highH = (TH1F*)f1->Get  ("specHigh");
        m_mcAmplitude  = 7239;
        m_mcEnergyScale = 2.482;
        m_mcResolution  = 0.1296/m_mcEnergyScale*100.;
        m_coAmplitude  = 691.3;
        m_coEnergyScale = 2.512;
        gausH->Smooth(500);
        highH->Smooth(500);
        lossH->Smooth(500);
        int n = highH->GetNbinsX();
        for(int i = 0; i < n-1; i++)
        {
                double energy  = highH->GetBinCenter (i+1);
                double conGaus = gausH->GetBinContent(i+1);
                double conHigh = highH->GetBinContent(i+1);
                double conLoss = lossH->GetBinContent(i+1);
                m_gausShape->SetPoint(i,energy,conGaus);
                m_lossShape->SetPoint(i,energy,conLoss);
                m_highShape->SetPoint(i,energy,conHigh);
        }*/
}
myFunction::~myFunction()
{
        /*
        delete m_gausShape;
        delete m_lossShape;
        delete m_highShape;
        delete f1;
        delete lossH;
        delete gausH;
        delete highH;*/
}
/*
double myFunction::Gaus(double *x, double *par)
{
        double energy = par[1];
        double res    = par[2]* par[1]*0.01;
        return par[0]*exp(-0.5*pow((x[0]-energy)/res,2));

}*/
double myFunction::MCShape(double *x,double *par)
{
        /*
        double energy    = x  [0];
        double eScale    = par[1]/m_mcEnergyScale;
        double lossAmp   = par[3]*par[0]/m_mcAmplitude;
        double highAmp   = par[0]/m_mcAmplitude;
        double gausShape;
        gausShape = par[0]*exp(-0.5*pow((x[0]-par[1])/(par[2]* par[1]*0.01),2));
        double highShape = m_highShape->Eval(energy/eScale);
        double lossShape = m_lossShape->Eval(energy/eScale);
        return gausShape + highAmp*highShape + lossAmp*lossShape;
        */
       double energy = par[1];
       double res    = par[2]* par[1]*0.01;
       return par[0]*exp(-0.5*pow((x[0]-energy)/res,2));

}
