#ifndef MYPDF
#define MYPDF



#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
class myPdf:public RooAbsPdf
{
    public:
        myPdf();
        myPdf(const char* name,const char* title,
                RooAbsReal& x_p,
                RooAbsReal& mean_p,
                RooAbsReal& sigma_p,
                RooAbsReal& d_p);
        virtual ~myPdf();
        myPdf(const myPdf& other,const char* name=0);
        virtual TObject* clone(const char* newname) const
        {
            return new myPdf(*this,newname);
        }

    protected:
        RooRealProxy x;
        RooRealProxy mean;
        RooRealProxy sigma;
        RooRealProxy d;
        Double_t evaluate() const;

    private:
        //ClassDef(myPdf,1)
        double m_mcAmplitude;
        double m_mcEnergyScale;
        double m_mcResolution;
        double m_coAmplitude;
        double m_coEnergyScale;
        TGraph* m_gausShape;
        TGraph* m_lossShape;
        TGraph* m_highShape;
        TFile* f1;
        TH1F* lossH;
        TH1F* highH;
        TH1F* gausH;
};


#endif
