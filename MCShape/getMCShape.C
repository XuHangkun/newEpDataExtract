{
    ofstream outdat_2("lossShape.txt");
    ofstream outdat_1("highShape.txt");
    TGraph* m_lossShape = new TGraph();
    TGraph* m_highShape = new TGraph();
    TFile* f1=new TFile("../MCShape/co60_specs_guwq.root");
    TH1F* lossH = (TH1F*)f1->Get  ("specLoss");
    TH1F* highH = (TH1F*)f1->Get  ("specHigh");
    highH->Smooth(500);
    lossH->Smooth(500);
    int n = highH->GetNbinsX();
    for(int i = 0; i < n-1; i++)
    {
        double energy  = highH->GetBinCenter (i+1);
        double conHigh = highH->GetBinContent(i+1);
        double conLoss = lossH->GetBinContent(i+1);
        outdat_1<<conHigh<<" ,"<<endl;
        outdat_2<<conLoss<<" ,"<<endl;        
    }
    outdat_1.close();
    outdat_2.close();
}
