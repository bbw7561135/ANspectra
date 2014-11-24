#include "TStyle.h"

void SetStyle(){
	gStyle->SetFillColor(0);
	gStyle->SetFrameLineWidth(1);
	gStyle->SetPadLeftMargin(0.1);
	gStyle->SetPadRightMargin(0.01);
	gStyle->SetPadTopMargin(0.005);
	gStyle->SetPadBottomMargin(0.08);
	gStyle->SetTitleOffset(1.35,"Y");
	gStyle->SetTitleOffset(1.15,"X");
}

void SIspectra(int nf, int nt){
	const int stre=10;
	const int cht=2, chf=2, che=92, cha=12, ched=(che-1)/stre+1;
	const int col[cha-1]={40,2,4,6,8,46,1,9,11,28,30}, sh[9]={2,7,4,3,1,6,9,8,5};
	const char t[cht]={'n','a'}, f[chf]={'e','m'}, fn[chf][5]={"e","#mu"}, tn[cht][8]={"#nu","#bar#nu"};
	const char *prefix="./data/", *postfix="E2.dat", *siname="HGm_KM", *hname="Honda";
	char sifile[256],sifile2[256],pict[256];
	FILE *fsi=101,*fsi2=101;
	double sp[che][cha],sp2[che][cha],y;
	int ed, e, a, ndr;
	TCanvas csp("csp","Spectra",200,200), cza("cza","Zenith Angle Distribution",200,200);
	TGraph grsp[cha-1],grsp2[cha-1],grza[ched];

	SetStyle();
//	sprintf(sifile,"%s%s_%c%c_%s",prefix,siname,t[nt],f[nf],postfix);
	sprintf(sifile,"./output/%c%c.dat",f[nf],t[nt]);
	fsi=fopen(sifile,"r");
//	sp[0][0]=0;
	for(a=0; a<cha; a++){
//		sp[0][a]=1-0.1*(a-1);
		fscanf(fsi,"%le",&(sp[0][a]));
	}
	for(e=1; e<che; e++){
		for(a=0; a<cha; a++){
			fscanf(fsi,"%le",&(sp[e][a]));
		}
	}
	fclose(fsi);
	for(a=1; a<cha; a++){
		for(e=1; e<che; e++){
			/*if((sp[e][a]==0)&&(sp[e][1]==0)){
				y=0;
			}else{
				y=sp[e][a]/sp[e][1];
			}*/
			y=sp[e][a];
			grsp[a-1].SetPoint(e-1,sp[e][0],y);
		}
		grsp[a-1].SetLineColor(col[a-1]);
	}

	TLegend* legname=new TLegend(0.15,0.15,0.2,0.2);
	legname->AddEntry((TObject*)0,Form("#font[132]{%s_{%s}}",tn[nt],fn[nf]),"");
	legname->SetBorderSize(0);

	for(ndr=0; ndr<2; ndr++){
		csp.cd();
		grsp[0].GetXaxis()->SetTitle("#font[132]{E_{#nu}, GeV}");
		//grsp[0].GetYaxis()->SetTitle("#font[132]{#Phi_{#nu}(E,#theta)/#Phi_{#nu}(E,0)}");
		//grsp[0].GetYaxis()->SetRangeUser(0.05,1.05);
		grsp[0].GetYaxis()->SetTitle("#font[132]{#Phi_{#nu}(E,#theta) E^{3}, GeV^{2}/(cm^{2} s ster)}");
		//grsp[0].GetYaxis()->SetRangeUser(1e-12,2e-1);
		//grsp[0].GetYaxis()->SetRangeUser(1e-8,2e-1);
		//grsp[0].GetYaxis()->SetRangeUser(2e-3,6e-2);//m
		grsp[0].GetYaxis()->SetRangeUser(7e-5,3e-2);//e
		grsp[0].GetXaxis()->SetLimits(7e0,2e4);
		grsp[0].Draw("AL");
		for(a=1; a<cha-1; a++){
			grsp[a].Draw("LSAME");
		}
		legname->Draw();
		//TLegend* legsp=new TLegend(0.8,0.3,0.95,0.75);
		TLegend* legsp=new TLegend(0.2,0.1,0.4,0.55);
		for(a=cha-1; a>0; a--){
			legsp->AddEntry(&grsp[a-1],Form("#font[132]{%5.2f}",sp[0][a]),"L");
		}
		legsp->SetBorderSize(0);
		//legsp->Draw();
		csp.SetLogx();
		csp.SetLogy();
		sprintf(pict,"./pictures/%s_%s_%c%c_spectra.eps",hname,siname,f[nf],t[nt]);
		csp.Print(pict);
		delete legsp;
	}
	for(e=1; e<che; e+=stre){
		ed=(e-1)/stre;
		for(a=1; a<cha; a++){
			if((sp[e][a]==0)&&(sp[e][1]==0)){
				y=0;
			}else{
				y=sp[e][a]/sp[e][1];
			}
			grza[ed].SetPoint(a-1,sp[0][a],y);
		}
		grza[ed].SetLineColor(col[ed]);
		grza[ed].SetLineStyle(sh[ed]);
	}
	cza.cd();
	grza[0].GetXaxis()->SetTitle("#font[132]{cos(#theta)}");
	grza[0].GetXaxis()->SetRangeUser(0,1);
	grza[0].GetYaxis()->SetTitle("#font[132]{#Phi_{#nu}(E,#theta)/#Phi_{#nu}(E,0)}");
	grza[0].GetYaxis()->SetRangeUser(0,1.05);
	grza[0].Draw("AC");
	for(ed=1; ed<ched; ed++){
		grza[ed].Draw("CSAME");
	}
	legname->Draw();
//	TLegend* legza=new TLegend(0.5,0.3,0.85,0.75);
	TLegend* legza=new TLegend(0.65,0.7,0.9,0.95);
	for(e=che-1; e>0; e-=stre){
		ed=(e-1)/stre;
		legza->AddEntry(&grza[ed],Form("#font[132]{10^{%1.0f} GeV}",log10(sp[e][0])),"L");
	}
	legza->SetBorderSize(0);
	legza->Draw();
	sprintf(pict,"./pictures/%s_%s_%c%c_ZAD.eps",hname,siname,f[nf],t[nt]);
	cza.Print(pict);
	delete legza;
	delete legname;
}
