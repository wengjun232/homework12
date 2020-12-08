const int npoints=5000;
Double_t x[npoints];

int get_input_data();
void fcn(Int_t &npar, Double_t *gin, Double_t &chi2, Double_t *par, Int_t iflag);
void minuit();

void MinuitFit()
{
  // Get data points
   get_input_data();
   minuit();
}

void minuit() {
   const int npar=2;
   TMinuit *gMinuit = new TMinuit(npar); 
   gMinuit->SetFCN(fcn);

   Double_t arglist[10];
   Int_t ierflg = 0;
   
   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

// Set starting values and step sizes for parameters
   Double_t vstart[npar] = { 0.5 };
   Double_t step[npar] = { 0.1 };
   gMinuit->mnparm(0, "tau", vstart[0], step[0], 0,0,ierflg);

// Now ready for minimization step
   arglist[0] = 500;
   arglist[1] = 1.;
   gMinuit->mnexcm("MIGRAD", arglist ,0,ierflg);
   gMinuit->mnexcm("HESSE", arglist ,0,ierflg);
   gMinuit->mnexcm("MINOS", arglist ,0,ierflg);   
// Print results
   Double_t fmin,fedm,errdef;
   Double_t tau,tau_err;
   Int_t nvpar,nparx,icstat;
   gMinuit->mnstat(fmin,fedm,errdef,nvpar,nparx,icstat);
   gMinuit->GetParameter(0, tau, tau_err );
   cout <<"tau           = "<<tau<<" +/- "<<tau_err<<endl;
   cout <<"log_l           = "<<-0.5*fmin<<endl; 
}

void fcn(Int_t &npar, Double_t *gin, Double_t &chi2, 
        Double_t *par, Int_t iflag)
{
   // Calculate log-likelihood
   Double_t log_l = 0;
   Double_t beta = par[0];
   Double_t alpha=par[1];
   for (Int_t i=0;i<npoints; i++) {
     log_l += log(3)-log(2*(3+bata))+log(1+alpha*x[i]+beta*x[i]*x[i];
   }
   // factor 2 to get errors right
   chi2 = -2.0*log_l;  
}

int get_input_data()
{
   // Generate data points according to the user defined pdf
   Double_t tauTrue = 1.0;
   TH1F *h1 = new TH1F("h1","x",100,-1,1);
   double sum=0.0;
    float p1=1;
        float p2=1;
        fnc->SetParameters(p1,p2);
        float maxy=fnc->GetMaximum();
        float miny=fnc->GetMinimum();
        int i=1;
        while(i<=npoints){
                float a=gRandom->Uniform(-1,1);
                float b=gRandom->Uniform(0,maxy);
                h1->Fill(b);
                float c=3/(2*(3+p1))*(1+p2*a+p1*a*a);
                if(b<=c){
                        hist->Fill(a);
                        x[i-1]=a;
			sum+=a;
                        i++;
                }
        }

/*   for (int i=0;i<npoints;i++) {
       Double_t r = gRandom->Exp( tauTrue );
       x[i]=r;
       sum += r;
       h1->Fill(r);
   }*/
   double average = sum/(double)npoints;
   cout << average << endl;
   h1->Draw();
   return 0;
}

