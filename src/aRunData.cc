#include "TROOT.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TH1F.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "aRunData.h"
#define VERBOSE 0
using namespace std;

aRunData::aRunData(TString in_type, int in_runNumber,float in_concentrationOfPPO, float in_concentrationOfBisMSB)
{
        type=in_type;
        runNumber=in_runNumber;
        concentrationOfPPO=in_concentrationOfPPO;
        concentrationOfBisMSB=in_concentrationOfBisMSB;
        heightOfCalib=0;
        runTime=0;
        
        nBinsOfLightYield=500;
        lowOfLightYield=200;
        highOfLightYield=700;
        nameOfLightYield=Form(" light yield at %.1f mZ (LS : %.1f g/L ppo %.1f mg/L bisMSB)",heightOfCalib,concentrationOfPPO,concentrationOfBisMSB);
        lightYield=new TH1F(type+nameOfLightYield,nameOfLightYield,nBinsOfLightYield,lowOfLightYield,highOfLightYield);

        nBinsOfRingChargeRatio=600;
        lowOfRingChargeRatio=0.01;
        highOfRingChargeRatio=0.3;
        for(int i=0;i<ringNumber;i++){
                nameOfRingChargeRatio[i]=Form(" charge ratio of ring %d at %.1f mZ (LS : %.1f g/L ppo %.1f mg/L bisMSB)",i+1,heightOfCalib,concentrationOfPPO,concentrationOfBisMSB);
                ringChargeRatio[i]=new TH1F(type+nameOfRingChargeRatio[i],nameOfRingChargeRatio[i],nBinsOfRingChargeRatio,lowOfRingChargeRatio,highOfRingChargeRatio);
        }

        nBinsOfPmtChargeRatio=500;
        lowOfPmtChargeRatio=0.0005;
        highOfPmtChargeRatio=0.025;
        for(int i=0;i<pmtNumber;i++){
                nameOfPmtChargeRatio[i]=Form(" charge ratio of pmt %d at %.1f mZ (LS : %.1f g/L ppo %.1f mg/L bisMSB)",i+1,heightOfCalib,concentrationOfPPO,concentrationOfBisMSB);
                pmtChargeRatio[i]=new TH1F(type+nameOfPmtChargeRatio[i],nameOfPmtChargeRatio[i],nBinsOfPmtChargeRatio,lowOfPmtChargeRatio,highOfPmtChargeRatio);
        }


        nBinsOfFirstHitTime=800;
        lowOfFirstHitTime=-1650;
        highOfFirstHitTime=-1250;
        nameOfFirstHitTime=Form(" first hit time at %.1f mZ (LS : %.1f g/L ppo %.1f mg/L bisMSB)",heightOfCalib,concentrationOfPPO,concentrationOfBisMSB);
        firstHitTime=new TH1F(type+nameOfFirstHitTime,nameOfFirstHitTime,nBinsOfFirstHitTime,lowOfFirstHitTime,highOfFirstHitTime);


}


aRunData::aRunData(TString in_type, int in_runNumber,float in_concentrationOfPPO, float in_concentrationOfBisMSB,float height)
{
        type=in_type;
        runNumber=in_runNumber;
        concentrationOfPPO=in_concentrationOfPPO;
        concentrationOfBisMSB=in_concentrationOfBisMSB;
        heightOfCalib=height;
        runTime=0;
        
        nBinsOfLightYield=500;
        lowOfLightYield=200;
        highOfLightYield=700;
        nameOfLightYield=Form(" light yield at %.1f mZ (LS : %.1f g/L ppo %.1f mg/L bisMSB)",heightOfCalib,concentrationOfPPO,concentrationOfBisMSB);
        lightYield=new TH1F(type+nameOfLightYield,nameOfLightYield,nBinsOfLightYield,lowOfLightYield,highOfLightYield);

        nBinsOfRingChargeRatio=600;
        lowOfRingChargeRatio=0.01;
        highOfRingChargeRatio=0.3;
        for(int i=0;i<ringNumber;i++){
                nameOfRingChargeRatio[i]=Form(" charge ratio of ring %d at %.1f mZ (LS : %.1f g/L ppo %.1f mg/L bisMSB)",i,heightOfCalib,concentrationOfPPO,concentrationOfBisMSB);
                ringChargeRatio[i]=new TH1F(type+nameOfRingChargeRatio[i],nameOfRingChargeRatio[i],nBinsOfRingChargeRatio,lowOfRingChargeRatio,highOfRingChargeRatio);
        }
        nBinsOfPmtChargeRatio=500;
        lowOfPmtChargeRatio=0.0005;
        highOfPmtChargeRatio=0.025;
        for(int i=0;i<pmtNumber;i++){
                nameOfPmtChargeRatio[i]=Form(" charge ratio of pmt %d at %.1f mZ (LS : %.1f g/L ppo %.1f mg/L bisMSB)",i,heightOfCalib,concentrationOfPPO,concentrationOfBisMSB);
                pmtChargeRatio[i]=new TH1F(type+nameOfPmtChargeRatio[i],nameOfPmtChargeRatio[i],nBinsOfPmtChargeRatio,lowOfPmtChargeRatio,highOfPmtChargeRatio);
        }


        nBinsOfFirstHitTime=800;
        lowOfFirstHitTime=-1650;
        highOfFirstHitTime=-1250;
        nameOfFirstHitTime=Form(" first hit time at %.1f mZ (LS : %.1f g/L ppo %.1f mg/L bisMSB)",heightOfCalib,concentrationOfPPO,concentrationOfBisMSB);
        firstHitTime=new TH1F(type+nameOfFirstHitTime,nameOfFirstHitTime,nBinsOfFirstHitTime,lowOfFirstHitTime,highOfFirstHitTime);


}

aRunData::~aRunData()
{
        delete lightYield;
        for(int i=0; i<ringNumber; i++){
                delete ringChargeRatio[i];
        }
        delete firstHitTime;
}

void aRunData::getData()
{
        int a=(runNumber/1000)*1000;
        int b=(runNumber/100)*100;
        TString fileName=TString(Form("/dybfs/rec/P18A/rec/runs_00%d/runs_00%d/run_00%d/recon.NoTag.00%d.",a,b,runNumber,runNumber))+type+TString(Form(".EH1-AD1.P18A-I._0001.root"));
        cout<<"File name: "<<fileName<<endl;
        TChain data("/Event/CalibReadout/CalibReadoutHeader");
        TChain adScaled("/Event/Rec/AdScaled");

        int status=data.Add(fileName);
        if(status==0)
        {
                cout<<"Bad filename: "<<fileName<<endl;
                return;
        }
        data.SetMakeClass(1);

        status=adScaled.Add(fileName);
        if(status==0)
        {
                cout<<"Bad filename: "<<fileName<<endl;
                return;
        }
        adScaled.SetMakeClass(1);
        data.AddFriend(&adScaled,"adScaled");

        data.SetBranchStatus("*",0);
        adScaled.SetBranchStatus("*",0);
       

        
        short detector=0;
        data.SetBranchStatus("adScaled.detector",kTRUE);
        data.SetBranchAddress("adScaled.detector",&detector);
        unsigned int triggerTimeSec=0;
        data.SetBranchStatus("adScaled.triggerTimeSec",kTRUE);
        data.SetBranchAddress("adScaled.triggerTimeSec",&triggerTimeSec);
        unsigned int triggerTimeNanoSec=0;
        data.SetBranchStatus("adScaled.triggerTimeNanoSec",kTRUE);
        data.SetBranchAddress("adScaled.triggerTimeNanoSec",&triggerTimeNanoSec);
        vector<float> chargeAD;
        data.SetBranchStatus("chargeAD",kTRUE);
        data.SetBranchAddress("chargeAD",&chargeAD);
        vector<unsigned int> ring;
        data.SetBranchStatus("ring",kTRUE);
        data.SetBranchAddress("ring",&ring);
        vector<unsigned int> column;
        data.SetBranchStatus("column",kTRUE);
        data.SetBranchAddress("column",&column);
        vector<float> timeAD;
        data.SetBranchStatus("timeAD",kTRUE);
        data.SetBranchAddress("timeAD",&timeAD);
        float energy=0;
        data.SetBranchStatus("adScaled.rawEvis",kTRUE);
        data.SetBranchAddress("adScaled.rawEvis",&energy);
 

        //select events here and save the information
        int maxEntries=data.GetEntries();
        long double beginTime=0;
        long double lastMuonTime=0;
        long double endTime=0;
        bool muonVeto=false;
        for(int entry=0;entry<maxEntries;entry++)
        {
                data.GetEntry(entry);
                
                if(VERBOSE>3){
                        cout<<"read tree "<<endl;
                        cout<<"Energy: "<<energy<<endl;
                }
                
                long double time=triggerTimeSec+1.e-9*triggerTimeNanoSec;
                if(entry==0)
                {
                        beginTime=time;
                }
                endTime=time;
                if((detector !=1 )) { continue; }
                if( energy>20 )
                {
                        lastMuonTime=time;
                        muonVeto=true;
                        continue;
                }
                if(muonVeto && (time-lastMuonTime)<1.e-3)
                {
                        continue;
                }else{
                        muonVeto=false;
                }
                if( energy<0.5 ) {  continue; }
                if( timeAD.size()==0 ) { continue; }
                float firstHit[192]={0};
                float pmtCharge[192]={0};
                float ly=0;
                float rcr[ringNumber]={0};
                for(int i=0;i<timeAD.size();i++){
                        int pmtIndex=(ring[i]-1)*24+column[i]-1;
                        if(timeAD[i]<(-1650) || timeAD[i]>(-1250)) { continue; }
                        if(firstHit[pmtIndex]>-1)
                        {
                                firstHit[pmtIndex]=timeAD[i];
                                if(chargeAD[i])
                                pmtCharge[pmtIndex]=chargeAD[i];
                                ly+=chargeAD[i];
                                rcr[ring[i]-1]+=chargeAD[i];
                        }
                        
                }
                //modify 177 pmt of LS: 
                //              2.0 g/L ppo and 7 mg/L bisMSB
                //              2.5 g/L ppo and 7 mg/L bisMSB
                if(concentrationOfPPO<2.6 && concentrationOfPPO>1.9){
                    if(concentrationOfBisMSB>6.9 && concentrationOfBisMSB<7.1){
                        ly-=pmtCharge[177];
                        rcr[7]-=pmtCharge[177];
                        pmtCharge[177]=rcr[7]/23.;
                        rcr[7]+=pmtCharge[177];
                        ly+=pmtCharge[177];
                    }
                }

                for(int j=0;j<ringNumber;j++){
                        rcr[j]/=ly;
                        ringChargeRatio[j]->Fill(rcr[j]);
                }
                lightYield->Fill(ly);
                for(int j=0;j<pmtNumber;j++){
                        float a=pmtCharge[j]/ly;
                        pmtChargeRatio[j]->Fill(a);
                        firstHitTime->Fill(firstHit[j]);
                }
                
                //manage 177 pmt
                //ofstream outdat;
                //outdat.open(Form("177PmtCharge_%.1fgPPO_%.1fmgBisMSB_%.1fmZ.txt",concentrationOfPPO,concentrationOfBisMSB,heightOfCalib),ios::app);
                //cout<<pmtCharge[176]<<endl;
                ///////////////



        }
        runTime=endTime-beginTime;                
}
