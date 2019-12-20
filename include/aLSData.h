#ifndef ALSDATA
#define ALSDATA


#include "aRunData.h"
#include "myFunction.h"

using namespace std;

class aLSData
{
public:
        aLSData(float ppo, float bisMSB);
        ~aLSData();

        //set and get function
        float getConcentrationOfPPO() { return concentrationOfPPO; }
        void  setConcentrationOfPPO(float value ) { concentrationOfPPO=value; }
        float getConcentrationOfBisMSB() { return concentrationOfBisMSB; }
        void  setConcentrationOfBisMSB(float value) { concentrationOfBisMSB=value; }
        float getHeight(int i) { return height[i]; }
        void  setHeight(int i,float value) { height[i]=value; }
        float getLightYield(int i) { return lightYield[i]; }
        void  setLightYield(int i, float value) { lightYield[i]=value; }
        float getFitErrorOfLightYieldError(int i) { return fitErrorOfLightYield[i]; }
        void  setFitErrorOfLightYieldError(int i, float value) { fitErrorOfLightYield[i]=value; }
        float getRingChargeRatio(int i,int j) { return ringChargeRatio[i][j]; }
        void  setRingChargeRatio(int i,int j,float value) { ringChargeRatio[i][j]=value; }
        float getFitErrorOfRingChargeRatio(int i,int j)   { return fitErrorOfRingChargeRatio[i][j]; }
        void  setFitErrorOfRingChargeRatio(int i,int j,float value) { fitErrorOfRingChargeRatio[i][j]=value; }
        float getPmtChargeRatio(int i,int j) { return pmtChargeRatio[i][j]; }
        void  setPmtChargeRatio(int i,int j,float value) { pmtChargeRatio[i][j]=value; }
        float getFitErrorOfPmtChargeRatio(int i,int j)   { return fitErrorOfPmtChargeRatio[i][j]; }
        void  setFitErrorOfPmtChargeRatio(int i,int j,float value) { fitErrorOfPmtChargeRatio[i][j]=value; }


        //get total data and bkgdata from run list
        void getData();
        //calculate background data from total data and signal data
        void calculateSignal();


        //fit light yield by MC shape
        void fitLightYield(int index_height);
        void fitAllLightYield();
        void fitRingChargeRatio(int index_height,int index_ring);
        void fitAllRingChargeRatio();
        void fitPmtChargeRatio(int index_height,int index_ring);
        void fitAllPmtChargeRatio();


        //save all result
        void save();

        //check data
        void checkData();

private:
        float concentrationOfPPO;
        float concentrationOfBisMSB;
        float height[12];

        int adCalibRunNumber[12];
        int physicsRunNumber;

        aRunData* totalData[12];
        aRunData* bkgData;
        aRunData* signalData[12];

        float lightYield[12];
        float fitErrorOfLightYield[12];
        float ringChargeRatio[12][8];
        float pmtChargeRatio[12][192];
        float fitErrorOfRingChargeRatio[12][8];
        float fitErrorOfPmtChargeRatio[12][192];
};


#endif
