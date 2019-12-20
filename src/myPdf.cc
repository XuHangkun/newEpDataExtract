#include "myPdf.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"


//ClassImp(myPdf)

myPdf::myPdf()
{
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
        }

}

myPdf::myPdf(const char* name,const char* title,
                RooAbsReal& x_p,
                RooAbsReal& mean_p,
                RooAbsReal& sigma_p,
                RooAbsReal& d_p):
    RooAbsPdf(name,title),
    x("x","x",this,x_p),
    mean("mean","mean",this,mean_p),
    sigma("sigma","sigam",this,sigma_p),
    d("d","d",this,d_p)
{
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
        }
}

myPdf::myPdf(const myPdf& other,const char* name):
    RooAbsPdf(other,name),
    x("x",this,other.x),
    mean("mean",this,other.mean),
    sigma("sigma",this,other.sigma),
    d("d",this,other.d)
{
    this->m_mcAmplitude=other.m_mcAmplitude;
    this->m_mcEnergyScale=other.m_mcEnergyScale;
    this->m_mcResolution=other.m_mcResolution;
    this->m_coAmplitude=other.m_coAmplitude;
    this->m_coEnergyScale=other.m_coEnergyScale;
    this->f1=other.f1;
    this->m_gausShape=other.m_gausShape;
    this->m_lossShape=other.m_lossShape;
    this->m_highShape=other.m_highShape;
    this->lossH=other.lossH;
    this->highH=other.highH;
    this->gausH=other.gausH;
}

myPdf::~myPdf()
{

}


Double_t myPdf::evaluate() const
{
        double energy    = x;
        double eScale    = mean/m_mcEnergyScale;
        double lossAmp   = d/m_mcAmplitude;
        double highAmp   = 1./m_mcAmplitude;
        double gausShape;
        gausShape = exp(-0.5*pow((x-mean)/(sigma*mean*0.01),2));
        double highShape = m_highShape->Eval(energy/eScale);
        double lossShape = m_lossShape->Eval(energy/eScale);
        return gausShape + highAmp*highShape + lossAmp*lossShape;
 
}
