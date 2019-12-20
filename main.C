#include "aRunData.h"
#include "aLSData.h"
#include <stdio.h>
#include <iostream>
using namespace std;

int main(int argc,char** argv)
{
        cout<<" @ ~_~ @   -----> let's start! "<<endl;


        char* pEnd_1;
        float PPO=0;
        PPO=strtof(argv[1],&pEnd_1);
        char* pEnd_2;
        float BisMSB=0;
        BisMSB=strtof(argv[2],&pEnd_2);

        
        aLSData testData(PPO,BisMSB);
        testData.getData();
        testData.calculateSignal();
        testData.fitAllLightYield();
        testData.fitAllRingChargeRatio();
        testData.fitAllPmtChargeRatio();
        testData.save();
        testData.checkData();





        cout<<" @ ~u~ @   -----> We finish here! "<<endl;
}
