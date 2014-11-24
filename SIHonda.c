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

void SIHonda(int nf, int nt){
	const int stre=10;
	const int cht=2, chf=2, che=92, cha=12, ched=(che-1)/stre+1;
	const int col[cha-1]={40,2,4,6,8,46,1,9,11,28,30}, sh[9]={2,7,4,3,1,6,9,8,5};
	const char t[cht]={'n','a'}, f[chf]={'e','m'}, fn[chf][5]={"e","#mu"}, tn[cht][8]={"#nu","#bar#nu"};
	char sifile[256], sifile2[256], pict[256];
	FILE *fsi=101, *fsi2=101;
	double sp[che][cha], sp2[che][cha];
	int ed, e, a, ndr;
	TCanvas csp("csp","Spectra",200,200);
	TGraph grsp[cha-1], grsp2[cha-1];

	SetStyle();
	sprintf(sifile,"./output/%c%c1.dat",f[nf],t[nt]);
	sprintf(sifile2,"./output/%c%c4.dat",f[nf],t[nt]);
	fsi=fopen(sifile,"r");
	fsi2=fopen(sifile2,"r");
	for(a=0; a<cha; a++){
		fscanf(fsi,"%le",&(sp[0][a]));
		fscanf(fsi2,"%le",&(sp2[0][a]));
	}
	for(e=1; e<che; e++){
		for(a=0; a<cha; a++){
			fscanf(fsi,"%le",&(sp[e][a]));
			fscanf(fsi2,"%le",&(sp2[e][a]));
		}
	}
	fclose(fsi);
	fclose(fsi2);

	for(a=1; a<cha; a++){
		for(e=1; e<che; e++){
			grsp[a-1].SetPoint(e-1,sp[e][0],sp[e][a]);
			grsp2[a-1].SetPoint(e-1,sp2[e][0],sp2[e][a]);
		}
		grsp[a-1].SetLineColor(col[a-1]);
		grsp2[a-1].SetLineColor(col[a-1]);
		grsp2[a-1].SetLineStyle(2);
	}

	for(ndr=0; ndr<2; ndr++){
		csp.cd();
		grsp[0].GetXaxis()->SetTitle("#font[132]{E_{#nu}, GeV}");
		grsp[0].GetYaxis()->SetTitle("#font[132]{#Phi_{#nu}(E,#theta) E^{3}, GeV^{2}/(cm^{2} s ster)}");
		grsp[0].GetYaxis()->SetRangeUser(2e-3,6e-2);//m
		//grsp[0].GetYaxis()->SetRangeUser(7e-5,3e-2);//e
		grsp[0].GetXaxis()->SetLimits(7e0,2e4);
		grsp[0].Draw("AL");
		grsp2[0].Draw("LSAME");
		for(a=1; a<cha-1; a++){
			grsp[a].Draw("LSAME");
			grsp2[a].Draw("LSAME");
		}
		TLegend* legsp=new TLegend(0.2,0.15,0.4,0.4);
		legsp->AddEntry((TObject*)0,Form("#font[132]{%s_{%s}}",tn[nt],fn[nf]),"");
		legsp->AddEntry(&grsp2[0],Form("#font[132]{%5.2f, Honda}",sp2[0][1]),"L");
		legsp->AddEntry(&grsp2[1],Form("#font[132]{%5.2f, Honda}",sp2[0][2]),"L");
		legsp->AddEntry(&grsp[cha-3],Form("#font[132]{%5.2f, ISU14}",sp[0][cha-2]),"L");
		legsp->AddEntry(&grsp[cha-2],Form("#font[132]{%5.2f, ISU14}",sp[0][cha-1]),"L");
		legsp->SetBorderSize(0);
		legsp->Draw();
		csp.SetLogx();
		csp.SetLogy();
		sprintf(pict,"./pictures/%c%c_comparison.eps",f[nf],t[nt]);
		csp.Print(pict);
		delete legsp;
	}
}
