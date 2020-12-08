void homeex(){
	TF1 *fnc=new TF1("fnc","[n]*3/(2*(3+[beta]))*(1+[alpha]*x+[beta]*x*x)",-0.95,0.95);
	float p1=1;
	float p2=1;
	float x[5000];
	fnc->SetParameters(1,p1,p2);
	float maxy=fnc->GetMaximum();
	float miny=fnc->GetMinimum();
	TH1F *hist=new TH1F("hist","hist",60,-0.95,0.95);
	int i=1;
	while(i<=5000){
		float a=gRandom->Uniform(-0.95,0.95);
		float b=gRandom->Uniform(0,maxy);
		float c=3/(2*(3+p1))*(1+p2*a+p1*a*a);
		if(b<=c){
			hist->Fill(a);
			x[i-1]=a;
			i++;
		}
	}
	TCanvas *c1=new TCanvas("c1","");
//	fnc->Draw();
//	hist->Draw();
	hist->Fit(fnc,"l");
//	TCanvas *c2=new TCanvas("c2","");
//	hist->Draw();
}


