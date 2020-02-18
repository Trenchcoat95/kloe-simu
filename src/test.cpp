#include <TGeoManager.h>
#include <TString.h>
#include <TGeoNode.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TApplication.h> 
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TDirectoryFile.h>



#include <vector>
#include <map>
#include <iostream>

void test(const char* fIn, const char* fIn2)
{
   TFile* f= new TFile(fIn,"");
   TFile* f2= new TFile(fIn2,"");
   
   TH1D* hEreco =  new TH1D("hEreco","",100,0,16000);
   hEreco = (TH1D*)f->Get("hmuonErecoInner");
   TH1D* hEreco2 =  new TH1D("hEreco2","",100,0,16000);
   hEreco2 = (TH1D*)f2->Get("hmuonErecoInner");
   
   TH1D* hEreco10 =  new TH1D("hEreco10","",100,0,16000);
   TH1D* hEreco20 =  new TH1D("hEreco20","",100,0,16000);
   TH1D* hEreco30 =  new TH1D("hEreco30","",100,0,16000);
   TH1D* hEreco40 =  new TH1D("hEreco40","",100,0,16000);
   TH1D* hEreco50 =  new TH1D("hEreco50","",100,0,16000);
   TH1D* hEreco60 =  new TH1D("hEreco60","",100,0,16000);
   TH1D* hEreco70 =  new TH1D("hEreco70","",100,0,16000);
   TH1D* hEreco80 =  new TH1D("hEreco80","",100,0,16000);
   TH1D* hEreco90 =  new TH1D("hEreco90","",100,0,16000);
   
   TH1D* hEreco102 =  new TH1D("hEreco102","",100,0,16000);
   TH1D* hEreco202 =  new TH1D("hEreco202","",100,0,16000);
   TH1D* hEreco302 =  new TH1D("hEreco302","",100,0,16000);
   TH1D* hEreco402 =  new TH1D("hEreco402","",100,0,16000);
   TH1D* hEreco502 =  new TH1D("hEreco502","",100,0,16000);
   TH1D* hEreco602 =  new TH1D("hEreco602","",100,0,16000);
   TH1D* hEreco702 =  new TH1D("hEreco702","",100,0,16000);
   TH1D* hEreco802 =  new TH1D("hEreco802","",100,0,16000);
   TH1D* hEreco902 =  new TH1D("hEreco902","",100,0,16000);
   
    TCanvas *c2 = new TCanvas("c2","c2",200,10,900,800);
	gStyle->SetOptStat(0);
	//gPad->SetLogy();

	hEreco->SetTitle("");
	hEreco->GetXaxis()->SetTitle("Energy (MeV)");
	hEreco->GetYaxis()->SetTitle("events");
	hEreco->SetLineColor(kBlue);
	hEreco->SetFillColor(kBlue-10);
	//hEreco->GetXaxis()->SetRangeUser(0., 80.);
	hEreco->SetFillStyle(3001);

	hEreco2->SetLineColor(kOrange);
	hEreco2->SetFillColor(kOrange);
	//hEreco2->GetXaxis()->SetRangeUser(0., 80.);
	hEreco2->SetFillStyle(3001);

	hEreco->Draw(" hist");
	hEreco2->Draw(" hist  same");


	TLegend* legend2 = new TLegend(0.62,0.8,0.894,0.9);
	legend2->AddEntry( hEreco, "E_{reco} #nu_{#mu} standard", "f");
	legend2->AddEntry( hEreco2, "E_{reco} #nu_{#mu} modified", "f");
	legend2->SetFillColor(0);
	legend2->SetBorderSize(1);
    legend2->Draw("");
	
	c2->Modified();
	c2->Update();
	gSystem->ProcessEvents();
	
   /*
   hEreco->Chi2Test(hEreco2,"P");
   hEreco->KolmogorovTest(hEreco2,"D");
   hEreco->AndersonDarlingTest(hEreco2,"D");
   */
   TTree* tLay = (TTree*) f->Get("tLay");
   int MuonReco;
   double Ereco;
   tLay->SetBranchAddress("MuonReco",&MuonReco);
   tLay->SetBranchAddress("Ereco",&Ereco);
   
   
   TTree* tLay2 = (TTree*) f2->Get("tLay");
   int MuonReco2;
   double Ereco2;
   tLay2->SetBranchAddress("MuonReco",&MuonReco2);
   tLay2->SetBranchAddress("Ereco",&Ereco2);
   
   vector<double> vectorE;
   vector<double> vectorE2;
   
   for (int i = 0; i < tLay->GetEntries(); i++) 
	{
		tLay->GetEntry(i);
		if (MuonReco==1){vectorE.push_back(Ereco);}
	}
	
   for (int i = 0; i < tLay2->GetEntries(); i++) 
	{
		tLay2->GetEntry(i);
		if (MuonReco2==1){vectorE2.push_back(Ereco2);}
	}
   
   Double_t* sample1 = new Double_t[vectorE.size()];
   Double_t* sample2 = new Double_t[vectorE2.size()];
   
   
   int percent[9];
   for(int i=1; i < 10; i++) {percent[i-1]=vectorE.size() * 0.1 * i;}  
   
   int percent2[9];
   for(int i=1; i < 10; i++) {percent2[i-1]=vectorE2.size() * 0.1 * i;}  
   
   //std::cout<<"percent[0]"<<percent[1]<<std::endl;
   
   Double_t* sample101 = new Double_t[percent[0]];
   Double_t* sample102 = new Double_t[percent2[0]];
   Double_t* sample201 = new Double_t[percent[1]];
   Double_t* sample202 = new Double_t[percent2[1]];
   Double_t* sample301 = new Double_t[percent[2]];
   Double_t* sample302 = new Double_t[percent2[2]];
   Double_t* sample401 = new Double_t[percent[3]];
   Double_t* sample402 = new Double_t[percent2[3]];
   Double_t* sample501 = new Double_t[percent[4]];
   Double_t* sample502 = new Double_t[percent2[4]];
   Double_t* sample601 = new Double_t[percent[5]];
   Double_t* sample602 = new Double_t[percent2[5]];
   Double_t* sample701 = new Double_t[percent[6]];
   Double_t* sample702 = new Double_t[percent2[6]];
   Double_t* sample801 = new Double_t[percent[7]];
   Double_t* sample802 = new Double_t[percent2[7]];
   Double_t* sample901 = new Double_t[percent[8]];
   Double_t* sample902 = new Double_t[percent2[8]];
   
   
   //std::cout<<"check2";
   
   for (int i = 0; i < vectorE.size(); i++) 
    {
	   if ( i < percent[0] ){sample101[i]=vectorE[i];
	                         hEreco10->Fill(vectorE[i]);}
	   if ( i < percent[1] ){sample201[i]=vectorE[i];
	                         hEreco20->Fill(vectorE[i]);}
	   if ( i < percent[2] ){sample301[i]=vectorE[i];
	                         hEreco30->Fill(vectorE[i]);}
	   if ( i < percent[3] ){sample401[i]=vectorE[i];
	                         hEreco40->Fill(vectorE[i]);}
	   if ( i < percent[4] ){sample501[i]=vectorE[i];
	                         hEreco50->Fill(vectorE[i]);}
	   if ( i < percent[5] ){sample601[i]=vectorE[i];
	                         hEreco60->Fill(vectorE[i]);}
	   if ( i < percent[6] ){sample701[i]=vectorE[i];
	                         hEreco70->Fill(vectorE[i]);}
	   if ( i < percent[7] ){sample801[i]=vectorE[i];
	                         hEreco80->Fill(vectorE[i]);}
	   if ( i < percent[8] ){sample901[i]=vectorE[i];
	                         hEreco90->Fill(vectorE[i]);}
	   sample1[i]=vectorE[i];
	}
   //std::cout<<"check3";
   for (int i = 0; i < vectorE2.size(); i++) 
    {
	   if ( i < percent2[0] ){sample102[i]=vectorE2[i];
	                          hEreco102->Fill(vectorE2[i]);}
	   if ( i < percent2[1] ){sample202[i]=vectorE2[i];
	                          hEreco202->Fill(vectorE2[i]);}
	   if ( i < percent2[2] ){sample302[i]=vectorE2[i];
	                          hEreco302->Fill(vectorE2[i]);}
	   if ( i < percent2[3] ){sample402[i]=vectorE2[i];
	                          hEreco402->Fill(vectorE2[i]);}
	   if ( i < percent2[4] ){sample502[i]=vectorE2[i];
	                          hEreco502->Fill(vectorE2[i]);}
	   if ( i < percent2[5] ){sample602[i]=vectorE2[i];
	                          hEreco602->Fill(vectorE2[i]);}
	   if ( i < percent2[6] ){sample702[i]=vectorE2[i];
	                          hEreco702->Fill(vectorE2[i]);}
	   if ( i < percent2[7] ){sample802[i]=vectorE2[i];
	                          hEreco802->Fill(vectorE2[i]);}
	   if ( i < percent2[8] ){sample902[i]=vectorE2[i];
	                          hEreco902->Fill(vectorE2[i]);}
	   sample2[i]=vectorE2[i];
    }
	
   //std::cout<<"check4";
   double pvalueAD[10];
   double pvalueKS[10];
   double pvaluechi[10];
   
   ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(vectorE.size(),sample1, vectorE2.size(), sample2);
   pvalueAD[9] = goftest-> AndersonDarling2SamplesTest();
   pvalueKS[9] = goftest-> KolmogorovSmirnov2SamplesTest();
   
   ROOT::Math::GoFTest* goftest10 = new ROOT::Math::GoFTest(percent[0],sample101, percent2[0], sample102);
   pvalueAD[0] = goftest10-> AndersonDarling2SamplesTest();
   pvalueKS[0] = goftest10-> KolmogorovSmirnov2SamplesTest();
   
   ROOT::Math::GoFTest* goftest20 = new ROOT::Math::GoFTest(percent[1],sample201, percent2[1], sample202);
   pvalueAD[1] = goftest20-> AndersonDarling2SamplesTest();
   pvalueKS[1] = goftest20-> KolmogorovSmirnov2SamplesTest();
   
   ROOT::Math::GoFTest* goftest30 = new ROOT::Math::GoFTest(percent[2],sample301, percent2[2], sample302);
   pvalueAD[2] = goftest30-> AndersonDarling2SamplesTest();
   pvalueKS[2] = goftest30-> KolmogorovSmirnov2SamplesTest();
   
   ROOT::Math::GoFTest* goftest40 = new ROOT::Math::GoFTest(percent[3],sample401, percent2[3], sample402);
   pvalueAD[3] = goftest40-> AndersonDarling2SamplesTest();
   pvalueKS[3] = goftest40-> KolmogorovSmirnov2SamplesTest();
   
   ROOT::Math::GoFTest* goftest50 = new ROOT::Math::GoFTest(percent[4],sample501, percent2[4], sample502);
   pvalueAD[4] = goftest50-> AndersonDarling2SamplesTest();
   pvalueKS[4] = goftest50-> KolmogorovSmirnov2SamplesTest();
   
   ROOT::Math::GoFTest* goftest60 = new ROOT::Math::GoFTest(percent[5],sample601, percent2[5], sample602);
   pvalueAD[5] = goftest60-> AndersonDarling2SamplesTest();
   pvalueKS[5] = goftest60-> KolmogorovSmirnov2SamplesTest();
   
   ROOT::Math::GoFTest* goftest70 = new ROOT::Math::GoFTest(percent[6],sample701, percent2[6], sample702);
   pvalueAD[6] = goftest70-> AndersonDarling2SamplesTest();
   pvalueKS[6] = goftest70-> KolmogorovSmirnov2SamplesTest();
   
   ROOT::Math::GoFTest* goftest80 = new ROOT::Math::GoFTest(percent[7],sample801, percent2[7], sample802);
   pvalueAD[7] = goftest80-> AndersonDarling2SamplesTest();
   pvalueKS[7] = goftest80-> KolmogorovSmirnov2SamplesTest();
   
   ROOT::Math::GoFTest* goftest90 = new ROOT::Math::GoFTest(percent[8],sample901, percent2[8], sample902);
   pvalueAD[8] = goftest90-> AndersonDarling2SamplesTest();
   pvalueKS[8] = goftest90-> KolmogorovSmirnov2SamplesTest();
   
   pvaluechi[0] = hEreco10->Chi2Test(hEreco102);
   pvaluechi[1] = hEreco20->Chi2Test(hEreco202);
   pvaluechi[2] = hEreco30->Chi2Test(hEreco302);
   pvaluechi[3] = hEreco40->Chi2Test(hEreco402);
   pvaluechi[4] = hEreco50->Chi2Test(hEreco502);
   pvaluechi[5] = hEreco60->Chi2Test(hEreco602);
   pvaluechi[6] = hEreco70->Chi2Test(hEreco702);
   pvaluechi[7] = hEreco80->Chi2Test(hEreco802);
   pvaluechi[8] = hEreco90->Chi2Test(hEreco902);
   pvaluechi[9] = hEreco->Chi2Test(hEreco2);
   
   Double_t x[10] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
   
   TCanvas *c0 = new TCanvas("c0","c0",200,10,900,600);
   TGraph* grAD = new TGraph(10,x,pvalueAD);
   grAD->SetTitle("");
   grAD->GetYaxis()->SetTitle("p-value");
   grAD->GetXaxis()->SetTitle("Sample fraction");
   grAD->SetMarkerColor(kBlue);
   grAD->Draw("AC*");

   TCanvas *c1 = new TCanvas("c1","c1",200,10,900,600);
   TGraph* grKS = new TGraph(10,x,pvalueKS);
   grKS->SetTitle("");
   grKS->GetYaxis()->SetTitle("p-value");
   grKS->GetXaxis()->SetTitle("Sample fraction");
   grKS->SetMarkerColor(kBlue);
   grKS->Draw("AC*");
   
   TCanvas *c3 = new TCanvas("c3","c3",200,10,900,600);
   TGraph* grchi = new TGraph(10,x,pvaluechi);
   grchi->SetTitle("");
   grchi->GetYaxis()->SetTitle("p-value");
   grchi->GetXaxis()->SetTitle("Sample fraction");
   grchi->SetMarkerColor(kBlue);
   grchi->Draw("AC*");
   
   
   
   std::cout<< " pvalueAD = " <<  pvalueAD[9] << ";   pvalueKS = "<< pvalueKS[9]<< std::endl;
   
   hEreco->Chi2Test(hEreco2,"P");
   
   c2->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/ErecoBoth.png");
   c0->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/pvalueAD.png");
   c1->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/pvalueKD.png");
   c3->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/pvaluechi.png");
   
   
}