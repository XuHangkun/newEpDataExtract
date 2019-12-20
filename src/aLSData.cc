#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TH2F.h"
#include "aLSData.h"
#include "myFunction.h"
#include "MCShape.h"
#include "myPdf.h"
#include <fstream>
#include <iostream>
#include "TStyle.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooDataHist.h"
#include "RooClassFactory.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include <fstream>
using namespace std;

double MCShape(double* x,double* par);



aLSData::aLSData(float ppo, float bisMSB)
{
        concentrationOfPPO=ppo;
        concentrationOfBisMSB=bisMSB;


        ifstream f(Form("/dybfs/users/xuhangkun/Dayabay/newEpDataExtract/runList/runList_%.1fPPO_%.1fbisMSB.dat",concentrationOfPPO,concentrationOfBisMSB));
        for(int i=0; i<12; i++){
                f>>adCalibRunNumber[i]>>height[i];
        }
        f>>physicsRunNumber;
        f.close();


        for(int i=0;i<12;i++){
                totalData[i]=new aRunData("ADCalib",adCalibRunNumber[i],concentrationOfPPO,concentrationOfBisMSB,height[i]);
                
                signalData[i]=new aRunData("Signal",adCalibRunNumber[i],concentrationOfPPO,concentrationOfBisMSB,height[i]);
                lightYield[i]=0;
                fitErrorOfLightYield[i]=0;
                for(int j=0;j<8;j++){
                        ringChargeRatio[i][j]=0;
                        fitErrorOfRingChargeRatio[i][j]=0;
                }
                for(int j=0;j<192;j++){
                        pmtChargeRatio[i][j]=0;
                        fitErrorOfPmtChargeRatio[i][j]=0;
                }

        }
        bkgData=new aRunData("Physics",physicsRunNumber,concentrationOfPPO,concentrationOfBisMSB);
}

aLSData::~aLSData()
{
        for(int i=0;i<12;i++){
                delete totalData[i];
                delete signalData[i];
        }
        delete bkgData;
}

void aLSData::getData()
{
        for(int i=0; i<12;i++){
                cout<<Form("Start getting AdCalib data at height %.1f mZ (LS: %.1f g/L PPO && %.1f mg/L bisMSB)",height[i],concentrationOfPPO,concentrationOfBisMSB)<<endl;
                totalData[i]->getData();
                cout<<"Finish --------------------------<<<<<<<<"<<endl;
        }
        
        cout<<Form("Start getting Physics data (LS: %.1f g/L PPO && %.1f mg/L bisMSB)",concentrationOfPPO,concentrationOfBisMSB)<<endl;
        bkgData->getData();
        cout<<"Finish --------------------------<<<<<<<<"<<endl;
}

void aLSData::calculateSignal()
{
        for(int i=0;i<12;i++){
                long double ratio=-(totalData[i]->getRunTime())/(bkgData->getRunTime());
                //calculate light yield 
                signalData[i]->getLightYield()->Sumw2();
                signalData[i]->getLightYield()->Add(totalData[i]->getLightYield(),bkgData->getLightYield(),1,ratio);
                for(int j=0;j<8;j++){
                        signalData[i]->getRingChargeRatio(j)->Sumw2();
                        signalData[i]->getRingChargeRatio(j)->Add(totalData[i]->getRingChargeRatio(j),bkgData->getRingChargeRatio(j),1,ratio);
                }
                for(int j=0;j<192;j++){
                        signalData[i]->getPmtChargeRatio(j)->Sumw2();
                        signalData[i]->getPmtChargeRatio(j)->Add(totalData[i]->getPmtChargeRatio(j),bkgData->getPmtChargeRatio(j),1,ratio);
                }
                
                signalData[i]->getFirstHitTime()->Sumw2();
                signalData[i]->getFirstHitTime()->Add(totalData[i]->getFirstHitTime(),bkgData->getFirstHitTime(),1,ratio);
                
        }
}
void aLSData::fitLightYield(int index_height)
{






         /*
         TH1F* h=(TH1F*)signalData[index_height]->getLightYield()->Clone();
         //use rooFit
         int maxBin=h->GetMaximumBin();
         double low=h->GetBinLowEdge(maxBin);
         double initpar_1=low+h->GetBinWidth(maxBin)/2;
         RooRealVar x("x","x",200,initpar_1+80);
         RooRealVar mean("mean","mean",initpar_1,initpar_1-50,initpar_1+80);
         RooRealVar sigma("sigma","sigma",8,1,30);
         RooRealVar d("d","d",6,0,30);
         RooDataHist data("data","dataset with x",x,h);

         RooClassFactory::makePdf("myPdf","x,A,mean,sigma,d");
         myPdf MCShape("MCShape","MCShape",x,mean,sigma,d);
         MCShape.fitTo(data);
         */
    



        
        TF1* function_mc=new TF1("function_mc",MCShape,200,560,4);
        function_mc->SetParameter(0,300);
        function_mc->SetParLimits(0,10,2000);
        int maxBin=signalData[index_height]->getLightYield()->GetMaximumBin();
        double low_1=signalData[index_height]->getLightYield()->GetBinLowEdge(maxBin);
        double initpar1=low_1+signalData[index_height]->getLightYield()->GetBinWidth(maxBin)/2;
        function_mc->SetParameter(1,initpar1);
        function_mc->SetParLimits(1,250,500);
        function_mc->SetParameter(2,5);
        function_mc->SetParLimits(2,1,20);
        function_mc->SetParameter(3,1);
        function_mc->SetParLimits(3,0,40);
        gStyle->SetOptFit(1111); 
        TFitResultPtr fitResultLightYield=signalData[index_height]->getLightYield()->Fit("function_mc","QRS");
        lightYield[index_height]=function_mc->GetParameter(1);
        fitErrorOfLightYield[index_height]=function_mc->GetParError(1);
        double Chi2LightYield=fitResultLightYield->Chi2()/(fitResultLightYield->Ndf());
        if(Chi2LightYield<5 && Chi2LightYield>2){
                cout<<"WARNING: BIG FIT CHISQUARE FOR LIGHT YIELD! "<<endl;
                cout<<"PPO concentration: "<<concentrationOfPPO<<"  BisMSB concentration: "<<concentrationOfBisMSB<<"  Height: "<<height[index_height]<<endl;
        }else if(Chi2LightYield>=5){
                cout<<"ERROR: BAD FIT CHISQUARE FOR LIGHT YIELD! "<<endl;
                cout<<"PPO concentration: "<<concentrationOfPPO<<"  BisMSB concentration: "<<concentrationOfBisMSB<<"  Height: "<<height[index_height]<<endl;
        }
        delete function_mc;

}

void aLSData::fitAllLightYield()
{
    for(int i=0;i<12;i++){
        fitLightYield(i);
    }

}

void aLSData::fitRingChargeRatio(int index_height, int index_ring)
{
                TFitResultPtr fitResultRingChargeRatio;
                double initpar1=signalData[index_height]->getRingChargeRatio(index_ring)->GetEntries();
                int maxBin=signalData[index_height]->getRingChargeRatio(index_ring)->GetMaximumBin();
                double low_1=signalData[index_height]->getRingChargeRatio(index_ring)->GetBinLowEdge(maxBin);
                double initpar2=low_1+signalData[index_height]->getRingChargeRatio(index_ring)->GetBinWidth(maxBin)/2;
                double initpar3=0.02;
                gStyle->SetOptFit(1111);
                if(height[index_height]<-0.8 && index_ring<1){
                        double low=initpar2-5*initpar3;
                        if(low<0.01){ low=0.01; }
                        double high=initpar2+1.0*initpar3;
                        if(high>0.3) {high=0.3; }
                        signalData[index_height]->getRingChargeRatio(index_ring)->Fit("gaus","QP","",low,high);
                        TF1* fit_1=signalData[index_height]->getRingChargeRatio(index_ring)->GetFunction("gaus");
                        initpar1=fit_1->GetParameter(0);
                        initpar2=fit_1->GetParameter(1);
                        initpar3=fit_1->GetParameter(2);
                        low=initpar2-5*initpar3;
                        if(low<0.01){ low=0.01; }
                        high=initpar2+1.0*initpar3;
                        if(high>0.3) {high=0.3; }
                        
                        fitResultRingChargeRatio=signalData[index_height]->getRingChargeRatio(index_ring)->Fit("gaus", "QS","",low,high);

                }else if(height[index_height]<-0.8 && index_ring<3){
                        double low=initpar2-1*initpar3;
                        if(low<0.01){ low=0.01; }
                        double high=initpar2+5.0*initpar3;
                        if(high>0.3) {high=0.3; }
                        
                        signalData[index_height]->getRingChargeRatio(index_ring)->Fit("gaus","QP",  "",low,high);
                        TF1* fit_2=signalData[index_height]->getRingChargeRatio(index_ring)->GetFunction("gaus");
                        initpar1=fit_2->GetParameter(0);
                        initpar2=fit_2->GetParameter(1);
                        initpar3=fit_2->GetParameter(2);
                        low=initpar2-1*initpar3;
                        if(low<0.01){ low=0.01; }
                        high=initpar2+5*initpar3;
                        if(high>0.3) {high=0.3; }
                       
                        fitResultRingChargeRatio=signalData[index_height]->getRingChargeRatio(index_ring)->Fit("gaus","QS","",low,high);
                }else{
                        fitResultRingChargeRatio=signalData[index_height]->getRingChargeRatio(index_ring)->Fit("gaus","QS");
                }
                TF1* func=signalData[index_height]->getRingChargeRatio(index_ring)->GetFunction("gaus");
                ringChargeRatio[index_height][index_ring]=func->GetParameter(1);
                fitErrorOfRingChargeRatio[index_height][index_ring]=func->GetParError(1);
                double Chi2RingChargeRatio=(fitResultRingChargeRatio->Chi2())/(fitResultRingChargeRatio->Ndf());
               if(Chi2RingChargeRatio<5 && Chi2RingChargeRatio>2){
                cout<<"WARNING: BIG FIT CHISQUARE FOR RING CHARGE RATIO!      "<<Chi2RingChargeRatio<<endl;
                cout<<"PPO concentration: "<<concentrationOfPPO<<"  BisMSB concentration: "<<concentrationOfBisMSB<<"  Height: "<<height[index_height]<<"   Ring: "<<(index_ring+1)<<endl;
        }else if(Chi2RingChargeRatio>=5){
                cout<<"ERROR: BAD FIT CHISQUARE FOR RING CHARGE RATIO!      "<<Chi2RingChargeRatio<<endl;
                cout<<"PPO concentration: "<<concentrationOfPPO<<"  BisMSB concentration: "<<concentrationOfBisMSB<<"  Height: "<<height[index_height]<<"   Ring: "<<(index_ring+1)<<endl;
        }
 

}

void aLSData::fitAllRingChargeRatio()
{
        for(int i=0;i<12;i++){
            for(int j=0;j<8;j++){
                fitRingChargeRatio(i,j);
            }
        }
}

void aLSData::fitPmtChargeRatio(int index_height, int index_pmt)
{


                pmtChargeRatio[index_height][index_pmt]=signalData[index_height]->getPmtChargeRatio(index_pmt)->GetMean();
                fitErrorOfPmtChargeRatio[index_height][index_pmt]=signalData[index_height]->getPmtChargeRatio(index_pmt)->GetMeanError();
                /*
                TFitResultPtr fitResultPmtChargeRatio;
                double initpar1=signalData[index_height]->getPmtChargeRatio(index_pmt)->GetEntries();
                int maxBin=signalData[index_height]->getPmtChargeRatio(index_pmt)->GetMaximumBin();
                double low_1=signalData[index_height]->getPmtChargeRatio(index_pmt)->GetBinLowEdge(maxBin);
                double initpar2=low_1+signalData[index_height]->getPmtChargeRatio(index_pmt)->GetBinWidth(maxBin)/2;
                double initpar3=0.004;
                gStyle->SetOptFit(1111);
                
                int NUMBERS=signalData[index_height]->getPmtChargeRatio(index_pmt)->GetNbinsX();
                if(NUMBERS<100){
                    return;
                }
                if(height[index_height]<-0.8 && index_pmt<1){
                        double low=initpar2-5*initpar3;
                        if(low<0.0005){ low=0.0005; }
                        double high=initpar2+1.0*initpar3;
                        if(high>0.025) {high=0.025; } 
                        signalData[index_height]->getPmtChargeRatio(index_pmt)->Fit("gaus","QP","",low,high);
                        TF1* fit_1=signalData[index_height]->getPmtChargeRatio(index_pmt)->GetFunction("gaus");
                        initpar1=fit_1->GetParameter(0);
                        initpar2=fit_1->GetParameter(1);
                        initpar3=fit_1->GetParameter(2);
                        low=initpar2-5*initpar3;
                        if(low<0.0005){ low=0.0005; }
                        high=initpar2+1.0*initpar3;
                        if(high>0.025) {high=0.025; }
                        fitResultPmtChargeRatio=signalData[index_height]->getPmtChargeRatio(index_pmt)->Fit("gaus", "QS","",low,high);

                }else if(height[index_height]<-0.8 && index_pmt<3){
                        double low=initpar2-5*initpar3;
                        if(low<0.0005){ low=0.0005; }
                        double high=initpar2+1.0*initpar3;
                        if(high>0.025) {high=0.025; }
                        signalData[index_height]->getPmtChargeRatio(index_pmt)->Fit("gaus","QP","",low,high);
                        TF1* fit_2=signalData[index_height]->getPmtChargeRatio(index_pmt)->GetFunction("gaus");
                        initpar1=fit_2->GetParameter(0);
                        initpar2=fit_2->GetParameter(1);
                        initpar3=fit_2->GetParameter(2);
                        low=initpar2-5*initpar3;
                        if(low<0.0005){ low=0.0005; }
                        high=initpar2+1.0*initpar3;
                        if(high>0.025) {high=0.025; }
                        fitResultPmtChargeRatio=signalData[index_height]->getPmtChargeRatio(index_pmt)->Fit("gaus","QS","",low,high);
                }else{
                        fitResultPmtChargeRatio=signalData[index_height]->getPmtChargeRatio(index_pmt)->Fit("gaus","QS");
                }
                TF1* func=signalData[index_height]->getPmtChargeRatio(index_pmt)->GetFunction("gaus");
                pmtChargeRatio[index_height][index_pmt]=func->GetParameter(1);
                fitErrorOfPmtChargeRatio[index_height][index_pmt]=func->GetParError(1);
                double Chi2PmtChargeRatio=(fitResultPmtChargeRatio->Chi2())/(fitResultPmtChargeRatio->Ndf());
                if(Chi2PmtChargeRatio<5 && Chi2PmtChargeRatio>2){
                    cout<<"WARNING: BIG FIT CHISQUARE FOR PMT CHARGE RATIO!      "<<Chi2PmtChargeRatio<<endl;
                    cout<<"PPO concentration: "<<concentrationOfPPO<<"  BisMSB concentration: "<<concentrationOfBisMSB<<"  Height: "<<height[index_height]<<"   Pmt: "<<(index_pmt+1)<<endl;
                }else if(Chi2PmtChargeRatio>=5){
                    cout<<"ERROR: BAD FIT CHISQUARE FOR PMT CHARGE RATIO!      "<<Chi2PmtChargeRatio<<endl;
                    cout<<"PPO concentration: "<<concentrationOfPPO<<"  BisMSB concentration: "<<concentrationOfBisMSB<<"  Height: "<<height[index_height]<<"   Pmt: "<<(index_pmt+1)<<endl;
                }
               */      

}

void aLSData::fitAllPmtChargeRatio()
{
        for(int i=0;i<12;i++){
            for(int j=0;j<192;j++){
                fitPmtChargeRatio(i,j);
            }
        }
}

void aLSData::save()
{


    //save all histogram in a root file
    TString filename=Form("/dybfs/users/xuhangkun/Dayabay/newEpDataExtract/result/%.1fgPPO_%.1fmgBisMSB.root",concentrationOfPPO,concentrationOfBisMSB);
    TFile f(filename,"recreate");
    for(int i=0;i<12;i++){
        f.mkdir(Form("%.1fmZ",height[i]));
        f.mkdir(Form("%.1fmZ/lightYield",height[i]));
        f.mkdir(Form("%.1fmZ/ringChargeRatio",height[i]));
        f.mkdir(Form("%.1fmZ/pmtChargeRatio",height[i]));
    }
    for(int i=0;i<12;i++){
        f.cd(Form("%.1fmZ/lightYield",height[i]));
        signalData[i]->getLightYield()->Write();
        for(int j=0;j<8;j++){
            f.cd(Form("%.1fmZ/ringChargeRatio",height[i]));
            signalData[i]->getRingChargeRatio(j)->Write();
        }
        for(int j=0;j<192;j++){
            f.cd(Form("%.1fmZ/pmtChargeRatio",height[i]));
            signalData[i]->getPmtChargeRatio(j)->Write();
        }
    }

    //save all fit result in a txt file
    for(int i=0;i<12;i++){
        ofstream outdat(Form("/dybfs/users/xuhangkun/Dayabay/newEpDataExtract/result/%.1fgPPO_%.1fmgBisMSB_%.1fmZ.txt",concentrationOfPPO,concentrationOfBisMSB,height[i]));
        outdat<<lightYield[i]<<"    "<<fitErrorOfLightYield[i]<<endl;
        for(int j=0;j<8;j++){
            outdat<<ringChargeRatio[i][j]<<"    "<<fitErrorOfRingChargeRatio[i][j]<<endl;
            for(int k=0;k<24;k++){
                if(k>=1){ outdat<<"      "; }
                outdat<<pmtChargeRatio[i][j*24+k]<<"    "<<fitErrorOfPmtChargeRatio[i][j*24+k];
            }
            outdat<<endl;
        }
        outdat.close();
    }

}

void aLSData::checkData()
{
    //we check pmt charge ratio here
    TCanvas c1("c1","check data");
    c1.Print(Form("/dybfs/users/xuhangkun/Dayabay/newEpDataExtract/dataCheck/checkDataOf_%.1fgPPO_%.1fmgBisMSB.pdf[",concentrationOfPPO,concentrationOfBisMSB));
    for(int k=0;k<12;k++){
        TString histoName=Form("checkData_%.1fmZ",height[k]);
        TH2F *h_checkData=new TH2F(histoName,histoName,8,0,8,24,0,24);
        for(int i=0;i<8;i++){
            for(int j=0;j<24;j++){
                h_checkData->SetBinContent(i+1,j+1,pmtChargeRatio[k][i*24+j]);
            }
        }
        h_checkData->Draw("colz");
        c1.Print(Form("/dybfs/users/xuhangkun/Dayabay/newEpDataExtract/dataCheck/checkDataOf_%.1fgPPO_%.1fmgBisMSB.pdf",concentrationOfPPO,concentrationOfBisMSB));
        delete h_checkData;
    }
    c1.Print(Form("/dybfs/users/xuhangkun/Dayabay/newEpDataExtract/dataCheck/checkDataOf_%.1fgPPO_%.1fmgBisMSB.pdf]",concentrationOfPPO,concentrationOfBisMSB));


    //add pmt charge ratios of same ring
    /*
    for(int i=0;i<12;i++){
        for(int j=0;j<8;j++){
            double r=0;
            for(int k=0;k<24;k++){
                r+=pmtChargeRatio[i][j*24+k];
            }
            double ratio=100*(signalData[i]->getRingChargeRatio(j)->GetMean()-r)/r;
            if(fabs(ratio)>2){
                cout<<"ERROR: PMT CHARGE RATIO MUST BE WRONG!!!      <<<----------------"<<endl;
            }
        }
    }*/
}
