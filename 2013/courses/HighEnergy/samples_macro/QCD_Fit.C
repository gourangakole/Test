#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TString.h"
#include "TRandom3.h"
#include "TLine.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TStyle.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooGlobalFunc.h"
#include "RooBifurGauss.h"
#include "RooLandau.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "RooCBShape.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "RooPlot.h"

using namespace RooFit;
//using namespace RooStats;
const double pi=acos(-1.0);

RooFitResult* Fit(TH1D* tChannel, TH1D* TTbar, TH1D* WJets, TH1D* QCD, TH1D* data, double *qcdEv, double *nqcdEv){
	
	TH1D* nonQCD; TH1D* hQCD_;
	hQCD_=(TH1D*) QCD->Clone();
	nonQCD=(TH1D*)tChannel->Clone(); //nonQCD->Sumw2();
	nonQCD->Add(TTbar); nonQCD->Add(WJets); //nonQCD->Add(WW);
	double NnonQCD=nonQCD->Integral();
	double Nqcd=QCD->Integral();
	
	hQCD_->Scale(1.0/hQCD_->Integral());
	nonQCD->Scale(1.0/nonQCD->Integral());
	
	double sigma_qcd=sqrt(Nqcd); double max_qcd= Nqcd + (450.0 * sigma_qcd); double min_qcd= Nqcd - (450.0 * sigma_qcd);
	double sigma_nqcd= sqrt(NnonQCD); double max_nqcd= NnonQCD + (100.0 * sigma_nqcd); double min_nqcd= NnonQCD - (100.0 * sigma_nqcd);
	
	cout<<"Total QCD: "<<Nqcd<<"\tsigma_qcd: "<<sigma_qcd<<"\tmax_qcd: "<<max_qcd<<"\tmin_qcd: "<<min_qcd<<endl;
	cout<<"Total nonQCD: "<<NnonQCD<<"\tsigma_nqcd: "<<sigma_nqcd<<"\tmax_nqcd: "<<max_nqcd<<"\tmin_nqcd: "<<min_nqcd<<endl;	 

	TCanvas* can=new TCanvas("Canvas","Canvas"); //c->Divide(2,1);
	TCanvas* can1=new TCanvas("Corelation","Corelation");
	
/// For mT fit

	RooRealVar mtw("mtw","M_{T} of W",0.0,200.0);
	RooPlot* frame = mtw.frame(Title("m_{T}"));
	RooDataHist ds("ds","ds",mtw,Import(*data));

	RooDataHist qcd_hist("qcd_hist","qcd_hist",RooArgList(mtw),hQCD_);
	RooDataHist Nqcd_hist("Nqcd_hist","Nqcd_hist",RooArgList(mtw),nonQCD);
	
	RooHistPdf QCD_PDF("QCD_PDF","QCD_PDF",mtw,qcd_hist);	
	RooHistPdf NonQCD_PDF("NonQCD_PDF","NonQCD_PDF",mtw,Nqcd_hist);

/// For MET fit

/*	RooRealVar met("met","#slash{E}_{T} ",0.0,200.0);
	RooPlot* frame = met.frame(Title("#slash{E}_{T}"));
	RooDataHist ds("ds","ds",met,Import(*data));

	RooDataHist qcd_hist("qcd_hist","qcd_hist",RooArgList(met),hQCD_);
	RooDataHist Nqcd_hist("Nqcd_hist","Nqcd_hist",RooArgList(met),nonQCD);
	
	RooHistPdf QCD_PDF("QCD_PDF","QCD_PDF",met,qcd_hist);	
	RooHistPdf NonQCD_PDF("NonQCD_PDF","NonQCD_PDF",met,Nqcd_hist);
*/
	ds.Print();
	ds.plotOn(frame,Name("Data"));
	
	if(min_qcd<0.0) min_qcd=0.0;
	RooRealVar nonQCDEvents("nonQCDEvents","nonQCDEvents",NnonQCD,min_nqcd,max_nqcd);
	RooRealVar QCDEvents("QCDEvents","QCDEvents",Nqcd,min_qcd,max_qcd);
	
	RooAddPdf model("model","model",RooArgList(NonQCD_PDF,QCD_PDF),RooArgList(nonQCDEvents,QCDEvents));

	RooFitResult* res = model.fitTo(ds,Extended(kTRUE),Save());
	model.plotOn(frame,Name("Fit")); model.paramOn(frame);

	TH2D* corr_hist=(TH2D*)res->correlationHist();
	
//	RooPlot* pull_frame =mtw.frame(Title("Pull Distribution"));
//	RooPlot* pull_frame =met.frame(Title("Pull Distribution"));
//	RooHist* hpull=frame->pullHist();
//	pull_frame->addPlotable(hpull,"P");
	
	cout<<"chi square/NDF : "<<frame->chiSquare()<<endl;
	
//	TPaveText *box= new TPaveText(0.3, 0.85, 0.6, 0.9,"BRNDC");
//	box->SetFillColor(10);
//	box->SetBorderSize(1);
//    box->SetTextAlign(12);
//    box->SetTextSize(0.04F);
//    box->SetFillStyle(1001);
//    box->SetFillColor(10);
//    TText *text = 0;
//    Char_t buf[30];
//    sprintf( buf,  "#chi^{2}/ndf = %f", frame->chiSquare() );
//    text = box->AddText( buf );
//    frame->addObject(box) ;  
	
	double xmax=200.0; double xmin=0.0;
	TLine *line=new TLine(xmin,0.0,xmax,0.0);
	line->SetLineColor(kRed); line->SetLineWidth(2);
	
	TLegend* leg= new TLegend(0.81,0.27,0.93,0.90);
	leg->SetTextSize(0.037); leg->SetBorderSize(0); leg->SetLineStyle(0); leg->SetTextSize(0.027); leg->SetFillStyle(0); leg->SetFillColor(0);
		
//	can->Divide(2,1);
//	can->cd(1);
	can->cd();
	model.plotOn(frame, Name("qcd"), Components(RooArgList(QCD_PDF)), LineColor(kRed),LineStyle(1), LineWidth(4));
	model.plotOn(frame, Name("non-qcd"), Components(RooArgList(NonQCD_PDF)), LineColor(kGreen),LineStyle(1), LineWidth(4));

	frame->GetXaxis()->CenterTitle(1); frame->Draw();
	leg->AddEntry("Data","Data","ple1"); leg->AddEntry("Fit","Fit","l"); leg->AddEntry("qcd","QCD","l"); leg->AddEntry("non-qcd","NonQCD","l"); leg->Draw();
	
//	can->cd(2);
//	pull_frame->GetYaxis()->SetTitle("Pull ( #sigma )"); pull_frame->GetYaxis()->CenterTitle(1); pull_frame->GetYaxis()->SetRangeUser(-5.0,5.0); pull_frame->GetXaxis()->CenterTitle(1); pull_frame->Draw(); line->Draw("SAME");
	
	can->Draw();
	
	can1->cd();corr_hist->Draw("COLZ TEXT");
	can1->Draw();

	*qcdEv = (double) QCDEvents.getVal(); *nqcdEv = (double) nonQCDEvents.getVal();
		
	res->correlationMatrix().Print();
		
/// For mT fit	
	mtw.setRange("myRange",50.0,200.0);
	RooAbsReal* qcd_no=QCD_PDF.createIntegral(mtw,NormSet(mtw),Range("myRange"));
	RooAbsReal* nonqcd_no=NonQCD_PDF.createIntegral(mtw,NormSet(mtw),Range("myRange"));
	cout<<"QCD fraction in mT > 50: "<<qcd_no->getVal()<<"\tnonQCD fraction in mT > 50: "<<nonqcd_no->getVal()<<endl;
	cout<<"=========================================\nExtrapolated Yields above 50 GeV\n=========================================="<<endl;
	cout<<"Data: "<<frame->getFitRangeNEvt(50.0,200.0)<<"+/-"<<sqrt(frame->getFitRangeNEvt(50.0,200.0))<<"\t QCD: "<<(qcd_no->getVal())*(QCDEvents.getVal())<<"+/-"<<(qcd_no->getVal())*(QCDEvents.getError())<<"\t nonQCD: "<<(nonqcd_no->getVal())*(nonQCDEvents.getVal())<<"+/-"<<(nonqcd_no->getVal())*(nonQCDEvents.getError())<<endl; 


/// For met fit
/*	met.setRange("myRange",45.0,200.0);
	RooAbsReal* qcd_no=QCD_PDF.createIntegral(met,NormSet(met),Range("myRange"));
	RooAbsReal* nonqcd_no=NonQCD_PDF.createIntegral(met,NormSet(met),Range("myRange"));

	cout<<"=========================================\nExtrapolated Yields above 45 GeV\n=========================================="<<endl;
	cout<<"Data: "<<frame->getFitRangeNEvt(45.0,200.0)<<"+/-"<<sqrt(frame->getFitRangeNEvt(45.0,200.0))<<"\t QCD: "<<(qcd_no->getVal())*QCDEvents.getVal()<<"+/-"<<(qcd_no->getVal())*(QCDEvents.getError())<<"\t nonQCD: "<<(nonqcd_no->getVal())*(nonQCDEvents.getVal())<<"+/-"<<(nonqcd_no->getVal())*(nonQCDEvents.getError())<<endl; 
*/	
	
	return res;
}	

int QCD_Fit(){
	
	gStyle->SetOptStat(0);

	std::vector<std::string> channel;
	channel.push_back("TChannel_Powheg");
	channel.push_back("TbarChannel_Powheg"); 
	channel.push_back("TWChannel_Powheg");
	channel.push_back("TbarWChannel_Powheg");
	channel.push_back("TTbar");
	channel.push_back("WJets_aMCatNLO");
	channel.push_back("DYJets_aMCatNLO");
	channel.push_back("QCD");
	channel.push_back("Data");

	unsigned int ChannelSize=channel.size();
	
//	channel.push_back("SChannel");
//	channel.push_back("WW");
//	channel.push_back("WZ");

	std::vector<std::string> basicObs;
	basicObs.push_back("muPt");
	basicObs.push_back("muEta");
	basicObs.push_back("muIso");
	basicObs.push_back("muCharge");
	basicObs.push_back("bJetPt");
	basicObs.push_back("bJetEta");
	basicObs.push_back("lJetPt");
	basicObs.push_back("lJetEta");
	basicObs.push_back("met");
	basicObs.push_back("metPhi");
	basicObs.push_back("nPV");
	basicObs.push_back("nGoodPV");
	basicObs.push_back("mtwMass");
	basicObs.push_back("topMass");
	basicObs.push_back("cosThetaStar");
	basicObs.push_back("diJetMass");
	basicObs.push_back("jetpTSum");
	basicObs.push_back("dEta_mu_bJet");
	basicObs.push_back("dR_bJet_lJet");
//	basicObs.push_back("px_nu");
//	basicObs.push_back("py_nu");
//	basicObs.push_back("pz_nu");
//	TLorentzVector muon4v, neutrino4v, Wboson4v;
//	basicObs.push_back("BDT_response");
	
	TString Channel, Var, treeName, HistNameIso, HistNameAntiIso;
	TString cutMC="Xsec_wgt*LHEWeightSign*bWeight*puWeight*(muPt>24.0)";
//	TString cutMC="";
	TString cutData="PassesTrig*(muPt>24.0)";
//	TString cutData="";


	TCanvas* cIso[19]; 
//	TCanvas* cAntiIso[15];
	TLegend* leg[19];
//	TLegend* legIso[15];
//	TLegend* legAntiIso[15];
	THStack* hsIso[19]; 
//	THStack* hsAntiIso[15];
	 

	for(unsigned int ii=0;ii<basicObs.size();ii++){
		Var=basicObs.at(ii);
		cIso[ii]= new TCanvas(Var,Var); cIso[ii]->Divide(1,2);
		hsIso[ii]=new THStack(Var,Var);

//		cAntiIso[ii]= new TCanvas(Var,Var); cAntiIso[ii]->Divide(1,2);
//		hsAntiIso[ii]=new THStack(Var,Var);

		leg[ii]= new TLegend(0.6,0.8,0.75,0.95);
		leg[ii]->SetFillColor(kWhite);
	}

	
	TTree* tr_Iso[8]; 
	TTree* tr_AntiIso[8];
	TH1D* hIso[9][19]; TH1D* hAntiIso[9][19];	
	TH1D* tChannel_Iso[19]; TH1D* tWChannel_Iso[19]; TH1D* WJets_Iso[19]; TH1D* DYJets_Iso[19]; TH1D* TTbar_Iso[19]; TH1D* QCD_Iso[19]; TH1D* AllMC_Iso[19]; // TH1D* WW_Iso[15]; TH1D* WZ_Iso[15];  TH1D* sChannel_Iso[15];
	TH1D* tChannel_AntiIso[19]; TH1D* tWChannel_AntiIso[19]; TH1D* WJets_AntiIso[19]; TH1D* DYJets_AntiIso[19]; TH1D* TTbar_AntiIso[19]; TH1D* QCD_AntiIso[19]; // TH1D* AllMC_AntiIso[15]; TH1D* WW_AntiIso[15]; TH1D* WZ_AntiIso[15]; TH1D* sChannel_AntiIso[15]; 

	TH1D* QCD_AntiIso_DD[19]; 

	TH1D* Data_Iso[19]; TH1D* Data_AntiIso[19];

	for(unsigned int ii=0;ii<ChannelSize;ii++){
		Channel=channel.at(ii);
		
		for(unsigned int kk=0;kk<basicObs.size();kk++){
			Var=basicObs.at(kk);
			HistNameIso=Channel+"_"+Var+"_Iso";
			HistNameAntiIso=Channel+"_"+Var+"_AntiIso";
			
			if(kk==0){ hIso[ii][kk]=new TH1D(HistNameIso,Var,30,0.0,150.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,30,0.0,150.0);} 
			if(kk==1){ hIso[ii][kk]=new TH1D(HistNameIso,Var,25,0.0,2.5); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,25,0.0,2.5);} 
			if(kk==2){ hIso[ii][kk]= new TH1D(HistNameIso,Var,50,0.0,0.5); hAntiIso[ii][kk]= new TH1D(HistNameAntiIso,Var,50,0.0,0.5);} 
			if(kk==3){ hIso[ii][kk]= new TH1D(HistNameIso,Var,3,-1.5,1.5); hAntiIso[ii][kk]= new TH1D(HistNameAntiIso,Var,3,-1.5,1.5);} 
			if(kk==4 || kk==6){ hIso[ii][kk]=new TH1D(HistNameIso,Var,40,0.0,200.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,40,0.0,200.0);}
			if(kk==5 || kk==7){ hIso[ii][kk]=new TH1D(HistNameIso,Var,10,0.0,5.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,10,0.0,5.0); }
			if(kk==8 || kk==12){ hIso[ii][kk]=new TH1D(HistNameIso,Var,20,0.0,200.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,20,0.0,200.0);}
			if(kk==9){ hIso[ii][kk]=new TH1D(HistNameIso,Var,20,-pi,pi); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,20,-pi,pi);} 
			if(kk==10 || kk==11){ hIso[ii][kk]=new TH1D(HistNameIso,Var,25,0.0,50.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,25,0.0,50.0);}
			if(kk==13){ hIso[ii][kk]= new TH1D(HistNameIso,Var,30,100.0,400.0); hAntiIso[ii][kk]= new TH1D(HistNameAntiIso,Var,30,100.0,400.0);}
			if(kk==14){ hIso[ii][kk]=new TH1D(HistNameIso,Var,10,-1.0,1.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,10,-1.0,1.0);}
			if(kk==15){ hIso[ii][kk]=new TH1D(HistNameIso,Var,45,50.0,500.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,45,50.0,500.0);}
			if(kk==16){ hIso[ii][kk]=new TH1D(HistNameIso,Var,45,50.0,500.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,45,50.0,500.0);}
			if(kk==17){ hIso[ii][kk]=new TH1D(HistNameIso,Var,100,0.0,10.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,100,0.0,10.0);}
			if(kk==18){ hIso[ii][kk]=new TH1D(HistNameIso,Var,100,0.0,10.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,100,0.0,10.0);}
//			if(kk==19){ hIso[ii][kk]=new TH1D(HistNameIso,Var,20,-1.0,1.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,20,-1.0,1.0);}
			

/*			if(kk==15){ hIso[ii][kk]=new TH1D(HistNameIso,Var,40,-200,0, 200.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,40,-200,0, 200.0);}
			if(kk==16){ hIso[ii][kk]=new TH1D(HistNameIso,Var,40,-200,0, 200.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,40,-200,0, 200.0);}
			if(kk==17){ hIso[ii][kk]=new TH1D(HistNameIso,Var,40,-200,0, 200.0); hAntiIso[ii][kk]=new TH1D(HistNameAntiIso,Var,40,-200,0, 200.0);}
*/				
//			if(kk==15) h[ii][kk]=new TH1D(HistName,Var,24,-0.12,0.12);
			hIso[ii][kk]->Reset(); hAntiIso[ii][kk]->Reset();
			hIso[ii][kk]->Sumw2(); hAntiIso[ii][kk]->Sumw2();
		}
	}

	TFile* inFileMC[2];
	TFile* inFileData[2];
	
	inFileMC[0]=new TFile("Trees_MC_2J1T_Iso.root","Read");
	inFileMC[1]=new TFile("Trees_MC_2J1T_AntiIso.root","Read");
	
	inFileData[0]=new TFile("Trees_Data_2J1T_Iso.root","Read");
	inFileData[1]=new TFile("Trees_Data_2J1T_AntiIso.root","Read");
	
	std::cout<<"Initialization Complete. Start Collecting Trees"<<std::endl;

	gROOT->cd();
	
	for(unsigned int ll=0;ll<ChannelSize;ll++){
		Channel=channel.at(ll);
		treeName=Channel;
		std::cout<<"Channel: "<<Channel<<"\ttreeName: "<<treeName<<std::endl;
		
		if(Channel!="Data"){ 
			tr_Iso[ll]=(TTree*)inFileMC[0]->Get(treeName);
			tr_AntiIso[ll]=(TTree*)inFileMC[1]->Get(treeName);
		}
		
		else{
			tr_Iso[ll]=(TTree*)inFileData[0]->Get(treeName);
			tr_AntiIso[ll]=(TTree*)inFileData[1]->Get(treeName);
		}				
			
//		if(Channel.Contains(TString("WJets"))) cutMC+="*0.324";
//		if(Channel.Contains(TString("DYJets"))) cutMC+="*0.3";

		for(unsigned int kk=0;kk<basicObs.size();kk++){
			Var=basicObs.at(kk);
			HistNameIso=Channel+"_"+Var+"_Iso";
			HistNameAntiIso=Channel+"_"+Var+"_AntiIso";
			std::cout<<"Channel: "<<Channel<<"\tVariable: "<<Var<<std::endl;

//			if(Var!="mtwMass"){ cutMC=cutMC+"*(mtwMass>50.0)"; cutData=cutData+"(mtwMass>50.0)";}

//			if(Var=="dPhi_mu_met") Var="abs("+Var+")";
			
			if(Channel!="Data"){ 
				tr_Iso[ll]->Project(HistNameIso,Var,cutMC);
				tr_AntiIso[ll]->Project(HistNameAntiIso,Var,cutMC);
			}		
			
			else{ 
				tr_Iso[ll]->Project(HistNameIso,Var,cutData);
				tr_AntiIso[ll]->Project(HistNameAntiIso,Var,cutData);
			}
								
			if(ll==0){ tChannel_Iso[kk]=(TH1D*)hIso[ll][kk]->Clone(); tChannel_AntiIso[kk]=(TH1D*)hAntiIso[ll][kk]->Clone();} 
			if(ll==1){ tChannel_Iso[kk]->Add(hIso[ll][kk]); tChannel_AntiIso[kk]->Add(hAntiIso[ll][kk]);}
			if(ll==2){ tWChannel_Iso[kk]=(TH1D*)hIso[ll][kk]->Clone(); tWChannel_AntiIso[kk]=(TH1D*)hAntiIso[ll][kk]->Clone();}
			if(ll==3){ tWChannel_Iso[kk]->Add(hIso[ll][kk]); tWChannel_AntiIso[kk]->Add(hAntiIso[ll][kk]);}
//			if(ll==4){ sChannel_Iso[kk]=(TH1D*)hIso[ll][kk]->Clone(); sChannel_AntiIso[kk]=(TH1D*)hAntiIso[ll][kk]->Clone();}
			if(ll==4){ TTbar_Iso[kk]=(TH1D*)hIso[ll][kk]->Clone(); TTbar_AntiIso[kk]=(TH1D*)hAntiIso[ll][kk]->Clone();} 
			if(ll==5){ WJets_Iso[kk]=(TH1D*)hIso[ll][kk]->Clone(); WJets_AntiIso[kk]=(TH1D*)hAntiIso[ll][kk]->Clone();} 
			if(ll==6){ DYJets_Iso[kk]=(TH1D*)hIso[ll][kk]->Clone(); DYJets_AntiIso[kk]=(TH1D*)hAntiIso[ll][kk]->Clone();}
//			if(ll==6) { WJets_[kk]=(TH1D*)h[kk][ll]->Clone(); WJets_[kk]->Scale(0.324);}
//			if(ll==7) { DYJets_[kk]=(TH1D*) h[kk][ll]->Clone();DYJets_[kk]->Scale(0.033);}
//			if(ll==8){ WW_Iso[kk]=(TH1D*)hIso[ll][kk]->Clone(); WW_AntiIso[kk]=(TH1D*)hAntiIso[ll][kk]->Clone();}
//			if(ll==9){ WZ_Iso[kk]=(TH1D*)hIso[ll][kk]->Clone(); WZ_AntiIso[kk]=(TH1D*)hAntiIso[ll][kk]->Clone();}
			if(ll==7){ QCD_Iso[kk]=(TH1D*)hIso[ll][kk]->Clone(); QCD_AntiIso[kk]=(TH1D*)hAntiIso[ll][kk]->Clone();}  
			if(ll==8){ Data_Iso[kk]=(TH1D*)hIso[ll][kk]->Clone(); Data_AntiIso[kk]=(TH1D*)hAntiIso[ll][kk]->Clone();}
		}

		std::cout<<"ll= "<<ll<<"\t"<<Channel<<" completed\tTotal: "<<ChannelSize<<std::endl;
//		std::cout<<Channel<<"\t Histograms collected"<<std::endl;
//		delete tr_Iso[ll]; delete tr_AntiIso[ll];
	}

	std::cout<<"All Histograms collected"<<std::endl;
	
	for(unsigned int zz=0;zz<basicObs.size();zz++){	
		Var=basicObs.at(zz);
		
//		if(Var=="metPt") std::cout<<"Data in Iso region: "<<Data_Iso[zz]->Integral()<<"\tQCD in Iso region: "<<QCD_Iso[zz]->Integral()<<"\nData in anti-Iso region: "<<Data_AntiIso[zz]->Integral()<<"\tQCD in anti-Iso region: "<<QCD_AntiIso[zz]->Integral()<<std::endl;
//		WJets_Iso[zz]->Scale(0.324); DYJets_Iso[zz]->Scale(0.11);
		WJets_Iso[zz]->Add(DYJets_Iso[zz]); WJets_AntiIso[zz]->Add(DYJets_AntiIso[zz]);

//		tWChannel_Iso[zz]->Add(sChannel_Iso[zz]); tWChannel_AntiIso[zz]->Add(sChannel_AntiIso[zz]); 

		TTbar_Iso[zz]->Add(tWChannel_Iso[zz]); TTbar_AntiIso[zz]->Add(tWChannel_AntiIso[zz]);

//		WW_Iso[zz]->Add(WZ_Iso[zz]); WW_AntiIso[zz]->Add(WZ_AntiIso[zz]);

		Data_Iso[zz]->SetMarkerStyle(20); Data_AntiIso[zz]->SetMarkerStyle(20);
		
		QCD_AntiIso_DD[zz]=(TH1D*)Data_AntiIso[zz]->Clone();

//		QCD_Iso[zz]->Scale(lumiScale); 	WJets_Iso[zz]->Scale(lumiScale);  TTbar_Iso[zz]->Scale(lumiScale); tChannel_Iso[zz]->Scale(lumiScale);
//		QCD_AntiIso[zz]->Scale(lumiScale); 	WJets_AntiIso[zz]->Scale(lumiScale);  TTbar_AntiIso[zz]->Scale(lumiScale); tChannel_AntiIso[zz]->Scale(lumiScale);

		QCD_AntiIso_DD[zz]->Add(tChannel_AntiIso[zz],-1); QCD_AntiIso_DD[zz]->Add(WJets_AntiIso[zz],-1); QCD_AntiIso_DD[zz]->Add(TTbar_AntiIso[zz],-1); //QCD_AntiIso_DD[zz]->Add(WW_AntiIso[zz],-1);  

	}	
	
	double QCDEvents, nonQCDEvents;
	
	double 	tot_nonQCD= WJets_Iso[12]->Integral() + tChannel_Iso[12]->Integral() + TTbar_Iso[12]->Integral();

	cout<<"Non QCD: "<<tot_nonQCD<<"\tQCD: "<<QCD_Iso[12]->Integral()<<endl;

	QCD_AntiIso_DD[12]->Scale(QCD_Iso[12]->Integral()/QCD_AntiIso_DD[12]->Integral());
	QCD_AntiIso[12]->Scale(QCD_Iso[12]->Integral()/QCD_AntiIso[12]->Integral());
	
//	RooFitResult* res=Fit(tChannel_Iso[12],TTbar_Iso[12],WJets_Iso[12],QCD_Iso[12], Data_Iso[12], &QCDEvents, &nonQCDEvents);
//	RooFitResult* res=Fit(tChannel_Iso[12],TTbar_Iso[12],WJets_Iso[12],QCD_AntiIso[12], Data_Iso[12], &QCDEvents, &nonQCDEvents);
	RooFitResult* res=Fit(tChannel_Iso[12],TTbar_Iso[12],WJets_Iso[12],QCD_AntiIso_DD[12], Data_Iso[12], &QCDEvents, &nonQCDEvents);


//	double 	tot_nonQCD= WJets_Iso[8]->Integral() + tChannel_Iso[8]->Integral() + TTbar_Iso[8]->Integral();
//	cout<<"Non QCD: "<<tot_nonQCD<<"\tQCD: "<<QCD_Iso[8]->Integral()<<endl;

//	QCD_AntiIso_DD[8]->Scale(QCD_Iso[8]->Integral()/QCD_AntiIso_DD[8]->Integral());
//	QCD_AntiIso[8]->Scale(QCD_Iso[8]->Integral()/QCD_AntiIso[8]->Integral());

//	RooFitResult* res=Fit(tChannel_Iso[8],TTbar_Iso[8],WJets_Iso[8],QCD_Iso[8], WW_Iso[8], Data_Iso[8], &QCDEvents, &nonQCDEvents);
//	RooFitResult* res=Fit(tChannel_Iso[8],TTbar_Iso[8],WJets_Iso[8],QCD_AntiIso[8], Data_Iso[8], &QCDEvents, &nonQCDEvents);
//	RooFitResult* res=Fit(tChannel_Iso[8],TTbar_Iso[8],WJets_Iso[8],QCD_AntiIso_DD[8], Data_Iso[8], &QCDEvents, &nonQCDEvents);

	res->floatParsFinal().Print("s");	


	TH1D* band[17]; int nBin;
	double lEdge[17]; double uEdge[17];
	double err;	

	TH1D* dummyData[17];
	const float xpad[2] = {0.f,1.f};
	const float ypad[4] = {0.,0.2351916,0.2351916,0.98};

//	double qcd_sf, nonqcd_sf; 	
//	double qcd_sf=0.92425; double nonqcd_sf=1.01253; 	

	
	for(unsigned int zz=0;zz<basicObs.size();zz++){

/// Scale according to the fit

		Var=basicObs.at(zz);	
		
/// Add histograms to the stack
		  
		
		QCD_Iso[zz]->Scale(QCDEvents/QCD_Iso[zz]->Integral()); WJets_Iso[zz]->Scale(nonQCDEvents/tot_nonQCD); TTbar_Iso[zz]->Scale(nonQCDEvents/tot_nonQCD); tChannel_Iso[zz]->Scale(nonQCDEvents/tot_nonQCD);
		QCD_AntiIso_DD[zz]->Scale(QCDEvents/QCD_AntiIso_DD[zz]->Integral());
//		QCD_AntiIso_DD[zz]->Scale(QCD_Iso[zz]->Integral()/QCD_AntiIso_DD[zz]->Integral());

				
		if(Var=="muIso"){ QCD_Iso[zz]->SetFillColor(kGray); QCD_Iso[zz]->SetLineColor(kBlack);  hsIso[zz]->Add(QCD_Iso[zz]);  AllMC_Iso[zz]=(TH1D*)QCD_Iso[zz]->Clone();}
		else{ QCD_AntiIso_DD[zz]->SetFillColor(kGray); QCD_AntiIso_DD[zz]->SetLineColor(kBlack);  hsIso[zz]->Add(QCD_AntiIso_DD[zz]);  AllMC_Iso[zz]=(TH1D*)QCD_AntiIso_DD[zz]->Clone();}
//		WW_Iso[zz]->SetFillColor(kBlue); WW_Iso[zz]->SetLineColor(kBlack); hsIso[zz]->Add(WW_Iso[zz]);  AllMC_Iso[zz]->Add(WW_Iso[zz]);
		WJets_Iso[zz]->SetFillColor(kGreen-2); WJets_Iso[zz]->SetLineColor(kBlack); hsIso[zz]->Add(WJets_Iso[zz]);  AllMC_Iso[zz]->Add(WJets_Iso[zz]);
		TTbar_Iso[zz]->SetFillColor(kOrange-3); TTbar_Iso[zz]->SetLineColor(kBlack); hsIso[zz]->Add(TTbar_Iso[zz]); AllMC_Iso[zz]->Add(TTbar_Iso[zz]);
		tChannel_Iso[zz]->SetFillColor(kRed); tChannel_Iso[zz]->SetLineColor(kBlack); hsIso[zz]->Add(tChannel_Iso[zz]);  AllMC_Iso[zz]->Add(tChannel_Iso[zz]);

/*		WW_AntiIso[zz]->SetFillColor(kBlue); WW_AntiIso[zz]->SetLineColor(kBlack); hsAntiIso[zz]->Add(WW_AntiIso[zz]);  AllMC_AntiIso[zz]=(TH1D*)WW_AntiIso[zz]->Clone();
		WJets_AntiIso[zz]->SetFillColor(kGreen-2); WJets_AntiIso[zz]->SetLineColor(kBlack); hsAntiIso[zz]->Add(WJets_AntiIso[zz]); AllMC_AntiIso[zz]=(TH1D*)WJets_AntiIso[zz]->Clone();// AllMC_AntiIso[zz]->Add(WJets_AntiIso[zz]);
		TTbar_AntiIso[zz]->SetFillColor(kOrange-3); TTbar_AntiIso[zz]->SetLineColor(kBlack); hsAntiIso[zz]->Add(TTbar_AntiIso[zz]); AllMC_AntiIso[zz]->Add(TTbar_AntiIso[zz]);
		tChannel_AntiIso[zz]->SetFillColor(kRed); tChannel_AntiIso[zz]->SetLineColor(kBlack); hsAntiIso[zz]->Add(tChannel_AntiIso[zz]);  AllMC_AntiIso[zz]->Add(tChannel_AntiIso[zz]);
		QCD_AntiIso[zz]->SetFillColor(kGray); QCD_AntiIso[zz]->SetLineColor(kBlack);  hsAntiIso[zz]->Add(QCD_AntiIso[zz]);  AllMC_AntiIso[zz]->Add(QCD_AntiIso[zz]);
*/

		if(AllMC_Iso[zz]->Integral()>0){
			dummyData[zz]=(TH1D*)Data_Iso[zz]->Clone();
			dummyData[zz]->Divide(AllMC_Iso[zz]);
			dummyData[zz]->GetYaxis()->SetTitle("Data/MC"); dummyData[zz]->GetYaxis()->CenterTitle(1);
			nBin=dummyData[zz]->GetXaxis()->GetNbins(); lEdge[zz]=dummyData[zz]->GetXaxis()->GetXmin(); uEdge[zz]=dummyData[zz]->GetXaxis()->GetXmax();
			TString bandTitle="Band_"+basicObs.at(zz);
			band[zz]=new TH1D(bandTitle,"",nBin,lEdge[zz],uEdge[zz]);
			for(int nn=0; nn<=nBin; nn++){
				band[zz]->SetBinContent(nn+1,1.0);
				if(AllMC_Iso[zz]->GetBinContent(nn+1)!=0) err=(AllMC_Iso[zz]->GetBinError(nn+1))*(dummyData[zz]->GetBinContent(nn+1))/AllMC_Iso[zz]->GetBinContent(nn+1);
				else err=0.0;	
				band[zz]->SetBinError(nn+1,err);
			}
			band[zz]->SetFillColor(kGray+2); band[zz]->SetFillStyle(3006);
			band[zz]->GetYaxis()->SetTitle("Data/MC"); band[zz]->GetYaxis()->CenterTitle(1); band[zz]->GetYaxis()->SetTitleOffset(0.19); band[zz]->GetYaxis()->SetTitleSize(0.05);	
		} 
		
		leg[zz]->AddEntry(Data_Iso[zz],"Data","ple1");
		leg[zz]->AddEntry(tChannel_Iso[zz],"t-Channel","f"); 
		leg[zz]->AddEntry(TTbar_Iso[zz],"t#bar{t}+tW","f"); 
		leg[zz]->AddEntry(WJets_Iso[zz],"W/Z+Jets","f");
//		leg[zz]->AddEntry(QCD_Iso[zz],"QCD","f");

		if(Var!="muIso")leg[zz]->AddEntry(QCD_AntiIso_DD[zz],"QCD (DD)","f");
		else leg[zz]->AddEntry(QCD_Iso[zz],"QCD","f");

//		if(Var!="mtwMass1" && Var!="muCharge")leg[zz]->AddEntry(QCD_Iso[zz],"QCD","f");
//		else leg[zz]->AddEntry(QCD_AntiIso_DD[zz],"QCD(DD)","f");
//		leg[zz]->AddEntry(WW_Iso[zz],"WW+WZ","f");

		leg[zz]->AddEntry(band[zz],"MC stat.","f"); 

		cIso[zz]->cd(1); //gPad->SetLogy(1);
		gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]); 
		hsIso[zz]->Draw("HIST");
		 //leg[kk]->Draw();
		Data_Iso[zz]->SetLineColor(kBlack);
		Data_Iso[zz]->Draw("E1 SAME"); leg[zz]->Draw(); //gPad->Update();
		
		cIso[zz]->cd(2);
		gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
		band[zz]->Draw("E2");	
		dummyData[zz]->SetMarkerStyle(20);
		dummyData[zz]->SetLineColor(kBlack); 
		dummyData[zz]->Draw("E1 SAME");
		 //dummyData[zz]->SetFillColor(kGray);
		gPad->SetGridy(1);  
		
		cIso[zz]->Draw();

/*	if(AllMC_AntiIso[zz]->Integral()>0){
			dummyData[zz]=(TH1D*)Data_AntiIso[zz]->Clone();
			dummyData[zz]->Divide(AllMC_AntiIso[zz]);
			dummyData[zz]->GetYaxis()->SetTitle("Data/MC"); dummyData[zz]->GetYaxis()->CenterTitle(1);
			nBin=dummyData[zz]->GetXaxis()->GetNbins(); lEdge[zz]=dummyData[zz]->GetXaxis()->GetXmin(); uEdge[zz]=dummyData[zz]->GetXaxis()->GetXmax();
			TString bandTitle="Band_"+basicObs.at(zz);
			band[zz]=new TH1D(bandTitle,"",nBin,lEdge[zz],uEdge[zz]);
			for(int nn=0; nn<=nBin; nn++){
				band[zz]->SetBinContent(nn+1,1.0);
				if(AllMC_AntiIso[zz]->GetBinContent(nn+1)!=0) err=(AllMC_AntiIso[zz]->GetBinError(nn+1))*(dummyData[zz]->GetBinContent(nn+1))/AllMC_AntiIso[zz]->GetBinContent(nn+1);
				else err=0.0;	
				band[zz]->SetBinError(nn+1,err);
			}
			band[zz]->SetFillColor(kGray+2); band[zz]->SetFillStyle(3006);
			band[zz]->GetYaxis()->SetTitle("Data/MC"); band[zz]->GetYaxis()->CenterTitle(1); band[zz]->GetYaxis()->SetTitleOffset(0.19); band[zz]->GetYaxis()->SetTitleSize(0.05);	
		} 
		
		leg[zz]->AddEntry(Data_AntiIso[zz],"Data","ple1");
		leg[zz]->AddEntry(tChannel_AntiIso[zz],"t-Channel","f"); 
		leg[zz]->AddEntry(TTbar_AntiIso[zz],"t#bar{t}+tW,s-channel single top","f"); 
		leg[zz]->AddEntry(WJets_AntiIso[zz],"W/Z+Jets","f");
		leg[zz]->AddEntry(QCD_AntiIso[zz],"QCD","f");
//		leg[zz]->AddEntry(WW_AntiIso[zz],"WW+WZ","f");
		leg[zz]->AddEntry(band[zz],"MC stat.","f"); 

		cAntiIso[zz]->cd(1); //gPad->SetLogy(1);
		gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]); 
		hsAntiIso[zz]->Draw("HIST");
		 //leg[kk]->Draw();
		Data_AntiIso[zz]->SetLineColor(kBlack);
		Data_AntiIso[zz]->Draw("E1 SAME"); leg[zz]->Draw(); //gPad->Update();
		
		cAntiIso[zz]->cd(2);
		gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
		band[zz]->Draw("E2");	
		dummyData[zz]->SetMarkerStyle(20);
		dummyData[zz]->SetLineColor(kBlack); 
		dummyData[zz]->Draw("E1 SAME");
		 //dummyData[zz]->SetFillColor(kGray);
		gPad->SetGridy(1);  
		
		cAntiIso[zz]->Draw();
*/ 
	}	

	
	hIso[7][12]->SetLineColor(kBlue); hAntiIso[7][12]->SetLineColor(kRed); QCD_AntiIso_DD[12]->SetLineColor(kBlack);	
	hIso[7][12]->Scale(1.0/hIso[7][12]->Integral()); hAntiIso[7][12]->Scale(1.0/hAntiIso[7][12]->Integral());  QCD_AntiIso_DD[12]->Scale(1.0/QCD_AntiIso_DD[12]->Integral());

//	hIso[6][8]->SetLineColor(kBlue); hAntiIso[6][8]->SetLineColor(kRed); QCD_AntiIso_DD[8]->SetLineColor(kBlack);	
//	hIso[6][8]->Scale(1.0/hIso[6][8]->Integral()); hAntiIso[6][8]->Scale(1.0/hAntiIso[6][8]->Integral());  QCD_AntiIso_DD[8]->Scale(1.0/QCD_AntiIso_DD[8]->Integral());

//	hIso[6][9]->SetLineColor(kBlue); hAntiIso[6][9]->SetLineColor(kRed); QCD_AntiIso_DD[9]->SetLineColor(kBlack);	
//	hIso[6][9]->Scale(1.0/hIso[6][9]->Integral()); hAntiIso[6][9]->Scale(1.0/hAntiIso[6][9]->Integral());  QCD_AntiIso_DD[9]->Scale(1.0/QCD_AntiIso_DD[9]->Integral());

//	hIso[6][0]->SetLineColor(kBlue); hAntiIso[6][0]->SetLineColor(kRed); QCD_AntiIso_DD[0]->SetLineColor(kBlack);	
//	hIso[6][0]->Scale(1.0/hIso[6][0]->Integral()); hAntiIso[6][0]->Scale(1.0/hAntiIso[6][0]->Integral());  QCD_AntiIso_DD[0]->Scale(1.0/QCD_AntiIso_DD[0]->Integral());
		
	TCanvas* shapeComp=new TCanvas("QCD Shapes","QCD shapes"); //shapeComp->Divide(2,1);
	TLegend* legShape=new TLegend(0.6,0.75,0.8,0.95);
	
	legShape->AddEntry(hIso[7][12],"MC QCD shape for I_{rel} < 0.06","ple1");
	legShape->AddEntry(hAntiIso[7][12],"MC QCD shape for I_{rel} > 0.12","ple1");
	legShape->AddEntry(QCD_AntiIso_DD[12],"Data-Driven QCD shape for I_{rel} > 0.12","ple1");
	
//	legShape->AddEntry(hIso[6][8],"MC QCD shape for I_{rel} < 0.06","ple1");
//	legShape->AddEntry(hAntiIso[6][8],"MC QCD shape for I_{rel} > 0.12","ple1");
//	legShape->AddEntry(QCD_AntiIso_DD[8],"Data-Driven QCD shape for I_{rel} > 0.12","ple1");

//	legShape->AddEntry(hIso[6][9],"MC QCD shape for I_{rel} < 0.06","ple1");
//	legShape->AddEntry(hAntiIso[6][9],"MC QCD shape for I_{rel} > 0.12","ple1");
//	legShape->AddEntry(QCD_AntiIso_DD[9],"Data-Driven QCD shape for I_{rel} > 0.12","ple1");

//	legShape->AddEntry(hIso[6][0],"MC QCD shape for I_{rel} < 0.06","ple1");
//	legShape->AddEntry(hAntiIso[6][0],"MC QCD shape for I_{rel} > 0.12","ple1");
//	legShape->AddEntry(QCD_AntiIso_DD[0],"Data-Driven QCD shape for I_{rel} > 0.12","ple1");

//	shapeComp->cd(1);

	shapeComp->cd();
	
	hIso[7][12]->Draw("PLE1"); hAntiIso[7][12]->Draw("PLE1 SAME"); QCD_AntiIso_DD[12]->Draw("PLE1 SAME");
	
//	shapeComp->cd(2);
//	hIso[6][8]->Draw("PLE1"); hAntiIso[6][8]->Draw("PLE1 SAME"); QCD_AntiIso_DD[8]->Draw("PLE1 SAME");

//	hIso[6][9]->Draw("PLE1"); hAntiIso[6][9]->Draw("PLE1 SAME"); QCD_AntiIso_DD[9]->Draw("PLE1 SAME");

//	hIso[6][0]->Draw("PLE1"); hAntiIso[6][0]->Draw("PLE1 SAME"); QCD_AntiIso_DD[0]->Draw("PLE1 SAME");

	legShape->Draw();
	
	shapeComp->Draw();
	
//	for(unsigned int mm=0;mm<ChannelSize;mm++){
//		delete tr_Iso[mm]; delete tr_AntiIso[mm];
//	}	
	return 1;
}				 	
