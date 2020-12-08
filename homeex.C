void homeex(){
	TF1 *fnc=new TF1("fnc","3/(2*(3+[beta]))*(1+[alpha]*x+[beta]*x*x)",-1,1);
	float p1=1;
	float p2=1;
	float x[5000];
	fnc->SetParameters(p1,p2);
	float maxy=fnc->GetMaximum();
	float miny=fnc->GetMinimum();
	TH1F *hist=new TH1F("hist","hist",60,-1,1);
	int i=1;
	while(i<=5000){
		float a=gRandom->Uniform(-1,1);
		float b=gRandom->Uniform(0,maxy);
		hist0->Fill(b);
		float c=3/(2*(3+p1))*(1+p2*a+p1*a*a);
		if(b<=c){
			hist->Fill(a);
			x[i-1]=a;
			i++;
		}
	}
	TCanvas *c1=new TCanvas("c1","");
	hist->Draw();
}


