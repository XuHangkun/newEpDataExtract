#ifndef ARUNDATA
#define ARUNDATA


#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#define ringNumber 8
#define pmtNumber 192
using namespace std;

class aRunData
{
public:
        aRunData(TString in_type,int in_runNumber,float in_concentrationOfPPO,float in_concentrationOfBisMSB);
        aRunData(TString in_type,int in_runNumber,float in_concentrationOfPPO,float in_concentrationOfBisMSB,float height);
        ~aRunData();

        //get and set function
        TH1F* getLightYield() { return lightYield; }
        void setLightYield(TH1F* value) { lightYield=value; }
        TH1F* getRingChargeRatio(int i) { return ringChargeRatio[i]; }
        void setRingChargeRatio(int i,TH1F* value) { ringChargeRatio[i]=value; }
        TH1F* getPmtChargeRatio(int i) { return pmtChargeRatio[i]; }
        void setPmtChargeRatio(int i,TH1F* value) { pmtChargeRatio[i]=value; }
        
        TH1F* getFirstHitTime() { return firstHitTime; }
        void setFirstHitTime(TH1F* value) { firstHitTime=value; } 

        float getConcentrationOfPPO() {return concentrationOfPPO; }
        void setConcentrationOfPPO(float value) { concentrationOfPPO=value; }
        float getConcentrationOfBisMSB() {return concentrationOfBisMSB; }
        void setConcentrationOfBisMSB(float value) { concentrationOfBisMSB=value; }
        float getHeightOfCalib() {return heightOfCalib; }
        void setHeightOfCalib(float value) { heightOfCalib=value; }
        TString getType() { return type; }
        void setType(TString value) { type=value; }
        int getRunNumber() { return runNumber; }
        void setRunNumber(int value) { runNumber=value; }
        long double getRunTime() { return runTime; }
        void setRunTime(long double value) { runTime=value; }
        
        void getData();





private:
         
        int   nBinsOfLightYield;
        float   lowOfLightYield;
        float   highOfLightYield;
        TString nameOfLightYield;
        TH1F* lightYield;

        int   nBinsOfRingChargeRatio;
        float   lowOfRingChargeRatio;
        float   highOfRingChargeRatio;
        TString nameOfRingChargeRatio[ringNumber];
        TH1F* ringChargeRatio[ringNumber];

        int   nBinsOfPmtChargeRatio;
        float   lowOfPmtChargeRatio;
        float   highOfPmtChargeRatio;
        TString nameOfPmtChargeRatio[pmtNumber];
        TH1F* pmtChargeRatio[pmtNumber];



  
        int   nBinsOfFirstHitTime;
        float   lowOfFirstHitTime;
        float   highOfFirstHitTime;   
        TString nameOfFirstHitTime;
        TH1F* firstHitTime;

        float concentrationOfPPO;
        float concentrationOfBisMSB;
        float heightOfCalib;
        TString type;
        int runNumber;
        long double runTime;

};

#endif
