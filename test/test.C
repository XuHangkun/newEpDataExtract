class myFunction
{
public:
    myFunction(){}
    ~myFunction(){}
    double Gaus(double *x, double *par)
    {
            double energy = par[1];
            double res    = par[2]* par[1]*0.01;
            return par[0]*exp(-0.5*pow((x[0]-energy)/res,2));

    }

};

void test()
{
    myFunction* f=new myFunction();
    TF1* func=new TF1("func",f,&myFunction::Gaus,-5,5,3,"myFunction","Gaus");
    func->SetParameter(0,10);
    func->SetParameter(1,1);
    func->SetParameter(2,1);
    TH1F* h=new TH1F("h","h",100,-5,5);
    h->FillRandom("gaus",10000);
    h->Fit(func,"LQR");




}
