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

void test(const char* fIn, const char* fIn2, const char* fChain, const char* fChain2, const char* fOut)
{
   TFile* f= new TFile(fIn,"");
   TFile* f2= new TFile(fIn2,"");
   
   TH1D* hEreco = (TH1D*)f->Get("hpreco_muonreco_In");
   TH1D* hEreco2 = (TH1D*)f2->Get("hpreco_muonreco_In");
   
   TH1D* hEnu = (TH1D*)f->Get("hEtrue_nu");
   TH1D* hEnu2 = (TH1D*)f2->Get("hEtrue_nu");
   
   //TH1S *hErecopercent[10];
   TH1D* hEreco10 =  new TH1D("hEreco10","",100,0,16000);
   TH1D* hEreco20 =  new TH1D("hEreco20","",100,0,16000);
   TH1D* hEreco30 =  new TH1D("hEreco30","",100,0,16000);
   TH1D* hEreco40 =  new TH1D("hEreco40","",100,0,16000);
   TH1D* hEreco50 =  new TH1D("hEreco50","",100,0,16000);
   TH1D* hEreco50s =  new TH1D("hEreco50s","",100,0,16000);
   TH1D* hEreco60 =  new TH1D("hEreco60","",100,0,16000);
   TH1D* hEreco70 =  new TH1D("hEreco70","",100,0,16000);
   TH1D* hEreco80 =  new TH1D("hEreco80","",100,0,16000);
   TH1D* hEreco90 =  new TH1D("hEreco90","",100,0,16000);
   
   TH1D* hEreco102 =  new TH1D("hEreco102","",100,0,16000);
   TH1D* hEreco202 =  new TH1D("hEreco202","",100,0,16000);
   TH1D* hEreco302 =  new TH1D("hEreco302","",100,0,16000);
   TH1D* hEreco402 =  new TH1D("hEreco402","",100,0,16000);
   TH1D* hEreco502 =  new TH1D("hEreco502","",100,0,16000);
   TH1D* hEreco502s =  new TH1D("hEreco502s","",100,0,16000);
   TH1D* hEreco602 =  new TH1D("hEreco602","",100,0,16000);
   TH1D* hEreco702 =  new TH1D("hEreco702","",100,0,16000);
   TH1D* hEreco802 =  new TH1D("hEreco802","",100,0,16000);
   TH1D* hEreco902 =  new TH1D("hEreco902","",100,0,16000);
   
    TCanvas *c2 = new TCanvas("Eboth","Eboth",200,10,900,800);
	gStyle->SetOptStat(0);
	//gPad->SetLogy();

	hEreco->SetTitle("");
	hEreco->GetXaxis()->SetTitle("p_{#mu}[MeV/c]");
	hEreco->GetYaxis()->SetTitle("Events");
	hEreco->SetLineColor(kBlue);
	hEreco->SetFillColor(kBlue-10);
	//hEreco->GetXaxis()->SetRangeUser(0., 80.);
	hEreco->SetFillStyle(3001);
    hEreco2->SetLineColor(kOrange);
	hEreco2->SetFillColor(kOrange);
	//hEreco2->GetXaxis()->SetRangeUser(0., 80.);
	hEreco2->SetFillStyle(3001);
	hEreco2->Scale(hEreco->Integral()/hEreco2->Integral());

	hEreco->Draw(" hist");
	hEreco2->Draw(" hist  same");


	TLegend* legend2 = new TLegend(0.62,0.8,0.894,0.9);
	legend2->AddEntry( hEreco, "p^{reco}_{#mu} from nominal flux", "f");
	legend2->AddEntry( hEreco2, "p^{reco}_{#mu} from modified flux", "f");
	legend2->SetFillColor(0);
	legend2->SetBorderSize(1);
    legend2->Draw("");
	
	c2->Modified();
	c2->Update();
	gSystem->ProcessEvents();
	
	
    TCanvas *c2b = new TCanvas("Enuboth","Enuboth",200,10,900,800);
	gStyle->SetOptStat(0);
	//gPad->SetLogy();

	hEnu->SetTitle("");
	hEnu->GetXaxis()->SetTitle("p_{#mu}[MeV/c]");
	hEnu->GetYaxis()->SetTitle("Events");
	hEnu->SetLineColor(kBlue);
	hEnu->SetFillColor(kBlue-10);
	//hEreco->GetXaxis()->SetRangeUser(0., 80.);
	hEnu->SetFillStyle(3001);
    hEnu2->SetLineColor(kOrange);
	hEnu2->SetFillColor(kOrange);
	//hEreco2->GetXaxis()->SetRangeUser(0., 80.);
	hEnu2->SetFillStyle(3001);
	hEnu2->Scale(hEnu->Integral()/hEnu2->Integral());
	

	hEnu->Draw(" hist");
	hEnu2->Draw(" hist  same");


	TLegend* legend2b = new TLegend(0.62,0.8,0.894,0.9);
	legend2b->AddEntry( hEreco, "E_{#nu} from nominal flux", "f");
	legend2b->AddEntry( hEreco2, "E_{#nu} from modified flux", "f");
	legend2b->SetFillColor(0);
	legend2b->SetBorderSize(1);
    legend2b->Draw("");
	
	c2b->Modified();
	c2b->Update();
	gSystem->ProcessEvents();
	
   /*
   hEreco->Chi2Test(hEreco2,"P");
   hEreco->KolmogorovTest(hEreco2,"D");
   hEreco->AndersonDarlingTest(hEreco2,"D");
   */
   TChain* tout = new TChain("tout");
   tout->Add(fChain);
   int isPmuOK;
   double pxmu_reco,pymu_reco,pzmu_reco,maxde_frontoutlayer,xv_reco;
   tout->SetBranchAddress("isPmuOK",&isPmuOK);
   tout->SetBranchAddress("maxde_frontoutlayer",&maxde_frontoutlayer);
   tout->SetBranchAddress("xv_reco",&xv_reco);
   tout->SetBranchAddress("pxmu_reco",&pxmu_reco);
   tout->SetBranchAddress("pymu_reco",&pymu_reco);
   tout->SetBranchAddress("pzmu_reco",&pzmu_reco);
  
   
   
   TChain* tout2 = new TChain("tout");
   tout2->Add(fChain2);
   int isPmuOK2;
   double pxmu_reco2,pymu_reco2,pzmu_reco2,maxde_frontoutlayer2,xv_reco2;
   tout2->SetBranchAddress("isPmuOK",&isPmuOK2);
   tout2->SetBranchAddress("maxde_frontoutlayer",&maxde_frontoutlayer2);
   tout2->SetBranchAddress("xv_reco",&xv_reco2);
   tout2->SetBranchAddress("pxmu_reco",&pxmu_reco2);
   tout2->SetBranchAddress("pymu_reco",&pymu_reco2);
   tout2->SetBranchAddress("pzmu_reco",&pzmu_reco2);
   
   vector<double> vectorE;
   vector<double> vectorE2;
   
   for (int i = 0; i < tout->GetEntries(); i++) 
	{
		
		tout->GetEntry(i);
		if (isPmuOK>0 && maxde_frontoutlayer < 40 && xv_reco <= 1690 && xv_reco >= -1690)
		{
			double ptot = std::sqrt(pxmu_reco*pxmu_reco + pymu_reco*pymu_reco + pzmu_reco*pzmu_reco);
			vectorE.push_back(ptot);
		    //std::cout<<"vectorE="<<vectorE.at(j)<<std::endl;
			

		}
	}
	
   for (int i = 0; i < tout2->GetEntries(); i++) 
	{
		tout2->GetEntry(i);
		if (isPmuOK2 >0 && maxde_frontoutlayer2 < 40 && xv_reco2 <= 1690 && xv_reco2 >= -1690){vectorE2.push_back(std::sqrt(pxmu_reco2*pxmu_reco2 + pymu_reco2*pymu_reco2 + pzmu_reco2*pzmu_reco2));}
	}
   
   Double_t* sample1 = new Double_t[vectorE.size()];
   Double_t* sample2 = new Double_t[vectorE2.size()];
   
   
   int percent[9];
   for(int i=1; i < 10; i++) {percent[9-i]=vectorE.size() * std::pow(0.1,i);}  
   
   int percent2[9];
   for(int i=1; i < 10; i++) {percent2[9-i]=vectorE2.size() * std::pow(0.1,i);}  
   
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
   Double_t* sample501s = new Double_t[vectorE.size()-percent[4]];
   Double_t* sample502s = new Double_t[vectorE2.size()-percent2[4]];
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
	   if ( i < percent[0] ){sample101[i]=vectorE.at(i);
	                         hEreco10->Fill(vectorE.at(i));}
	   if ( i < percent[1] ){sample201[i]=vectorE.at(i);
	                         hEreco20->Fill(vectorE.at(i));}
	   if ( i < percent[2] ){sample301[i]=vectorE.at(i);
	                         hEreco30->Fill(vectorE.at(i));}
	   if ( i < percent[3] ){sample401[i]=vectorE.at(i);
	                         hEreco40->Fill(vectorE.at(i));}
	   if ( i < percent[4] ){sample501[i]=vectorE.at(i);
	                         hEreco50->Fill(vectorE.at(i));}
	   if ( i < percent[5] ){sample601[i]=vectorE.at(i);
	                         hEreco60->Fill(vectorE.at(i));}
	   if ( i < percent[6] ){sample701[i]=vectorE.at(i);
	                         hEreco70->Fill(vectorE.at(i));}
	   if ( i < percent[7] ){sample801[i]=vectorE.at(i);
	                         hEreco80->Fill(vectorE.at(i));}
	   if ( i < percent[8] ){sample901[i]=vectorE.at(i);
	                         hEreco90->Fill(vectorE.at(i));}
	   if ( i > percent[4] ){sample501s[i]=vectorE.at(i);
	                         hEreco50s->Fill(vectorE.at(i));}				 
	   sample1[i]=vectorE.at(i);
	   //std::cout<<sample1[i]<<std::endl;
	}
   //std::cout<<"check3";
   for (int i = 0; i < vectorE2.size(); i++) 
    {
	   if ( i < percent2[0] ){sample102[i]=vectorE2.at(i);
	                          hEreco102->Fill(vectorE2.at(i));}
	   if ( i < percent2[1] ){sample202[i]=vectorE2.at(i);
	                          hEreco202->Fill(vectorE2.at(i));}
	   if ( i < percent2[2] ){sample302[i]=vectorE2[i];
	                          hEreco302->Fill(vectorE2.at(i));}
	   if ( i < percent2[3] ){sample402[i]=vectorE2[i];
	                          hEreco402->Fill(vectorE2.at(i));}
	   if ( i < percent2[4] ){sample502[i]=vectorE2.at(i);
	                          hEreco502->Fill(vectorE2.at(i));}
	   if ( i < percent2[5] ){sample602[i]=vectorE2.at(i);
	                          hEreco602->Fill(vectorE2.at(i));}
	   if ( i < percent2[6] ){sample702[i]=vectorE2.at(i);
	                          hEreco702->Fill(vectorE2.at(i));}
	   if ( i < percent2[7] ){sample802[i]=vectorE2.at(i);
	                          hEreco802->Fill(vectorE2.at(i));}
	   if ( i < percent2[8] ){sample902[i]=vectorE2.at(i);
	                          hEreco902->Fill(vectorE2.at(i));}
	   if ( i > percent2[4] ){sample502s[i]=vectorE2.at(i);
	                         hEreco502s->Fill(vectorE2.at(i));}
	   sample2[i]=vectorE2.at(i);
    }
	
   //std::cout<<"check4";
   double pvalueAD[10];
   double pvalueKS[10];
   double pvaluechi[10];
   double pvalueADtest;
   double pvalueKStest;
   double pvaluechitest;
   
   std::cout<<"Starting tests:"<<std::endl;
   
   /*
   ROOT::Math::GoFTest* goftest10 = new ROOT::Math::GoFTest(percent[0],sample101, percent2[0], sample102);
   pvalueAD[0] = goftest10-> AndersonDarling2SamplesTest();
   pvalueKS[0] = goftest10-> KolmogorovSmirnov2SamplesTest();
   std::cout<<"AD and KS for 0.1 sample done"<<std::endl;
   
   ROOT::Math::GoFTest* goftest20 = new ROOT::Math::GoFTest(percent[1],sample201, percent2[1], sample202);
   pvalueAD[1] = goftest20-> AndersonDarling2SamplesTest();
   pvalueKS[1] = goftest20-> KolmogorovSmirnov2SamplesTest();
   std::cout<<"AD and KS for 0.2 sample done"<<std::endl;
   
   ROOT::Math::GoFTest* goftest30 = new ROOT::Math::GoFTest(percent[2],sample301, percent2[2], sample302);
   pvalueAD[2] = goftest30-> AndersonDarling2SamplesTest();
   pvalueKS[2] = goftest30-> KolmogorovSmirnov2SamplesTest();
   std::cout<<"AD and KS for 0.3 sample done"<<std::endl;
   
   ROOT::Math::GoFTest* goftest40 = new ROOT::Math::GoFTest(percent[3],sample401, percent2[3], sample402);
   pvalueAD[3] = goftest40-> AndersonDarling2SamplesTest();
   pvalueKS[3] = goftest40-> KolmogorovSmirnov2SamplesTest();
   std::cout<<"AD and KS for 0.4 sample done"<<std::endl;
   
   ROOT::Math::GoFTest* goftest50 = new ROOT::Math::GoFTest(percent[4],sample501, percent2[4], sample502);
   pvalueAD[4] = goftest50-> AndersonDarling2SamplesTest();
   pvalueKS[4] = goftest50-> KolmogorovSmirnov2SamplesTest();
   std::cout<<"AD and KS for 0.5 sample done"<<std::endl;
   
   ROOT::Math::GoFTest* goftest60 = new ROOT::Math::GoFTest(percent[5],sample601, percent2[5], sample602);
   pvalueAD[5] = goftest60-> AndersonDarling2SamplesTest();
   pvalueKS[5] = goftest60-> KolmogorovSmirnov2SamplesTest();
   std::cout<<"AD and KS for 0.6 sample done"<<std::endl;
   
   ROOT::Math::GoFTest* goftest70 = new ROOT::Math::GoFTest(percent[6],sample701, percent2[6], sample702);
   pvalueAD[6] = goftest70-> AndersonDarling2SamplesTest();
   pvalueKS[6] = goftest70-> KolmogorovSmirnov2SamplesTest();
   std::cout<<"AD and KS for 0.7 sample done"<<std::endl;
   
   ROOT::Math::GoFTest* goftest80 = new ROOT::Math::GoFTest(percent[7],sample801, percent2[7], sample802);
   pvalueAD[7] = goftest80-> AndersonDarling2SamplesTest();
   pvalueKS[7] = goftest80-> KolmogorovSmirnov2SamplesTest();
   std::cout<<"AD and KS for 0.8 sample done"<<std::endl;
   
   ROOT::Math::GoFTest* goftest90 = new ROOT::Math::GoFTest(percent[8],sample901, percent2[8], sample902);
   pvalueAD[8] = goftest90-> AndersonDarling2SamplesTest();
   pvalueKS[8] = goftest90-> KolmogorovSmirnov2SamplesTest();
   std::cout<<"AD and KS for 0.9 sample done"<<std::endl;
   
   ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(vectorE.size(),sample1, vectorE2.size(), sample2);
   pvalueAD[9] = goftest-> AndersonDarling2SamplesTest();
   pvalueKS[9] = goftest-> KolmogorovSmirnov2SamplesTest();
   std::cout<<"AD and KS for tot sample done"<<std::endl;
   
   ROOT::Math::GoFTest* goftestcontrol = new ROOT::Math::GoFTest(vectorE.size()-percent[4],sample501s, percent[4], sample501);
   pvalueADtest = goftestcontrol-> AndersonDarling2SamplesTest();
   pvalueKStest = goftestcontrol-> KolmogorovSmirnov2SamplesTest();
   std::cout<<"AD and KS for test sample done"<<std::endl;
   
   */
   hEreco50s->Scale(hEreco50->Integral()/hEreco50s->Integral());
   pvaluechi[0] = hEreco10->Chi2Test(hEreco102);
   std::cout<<"chi2 for 0.1 sample done"<<std::endl;
   pvaluechi[1] = hEreco20->Chi2Test(hEreco202);
   std::cout<<"chi2 for 0.2 sample done"<<std::endl;
   pvaluechi[2] = hEreco30->Chi2Test(hEreco302);
   std::cout<<"chi2 for 0.3 sample done"<<std::endl;
   pvaluechi[3] = hEreco40->Chi2Test(hEreco402);
   std::cout<<"chi2 for 0.4 sample done"<<std::endl;
   pvaluechi[4] = hEreco50->Chi2Test(hEreco502);
   std::cout<<"chi2 for 0.5 sample done"<<std::endl;
   pvaluechi[5] = hEreco60->Chi2Test(hEreco602);
   std::cout<<"chi2 for 0.6 sample done"<<std::endl;
   pvaluechi[6] = hEreco70->Chi2Test(hEreco702);
   std::cout<<"chi2 for 0.7 sample done"<<std::endl;
   pvaluechi[7] = hEreco80->Chi2Test(hEreco802);
   std::cout<<"chi2 for 0.8 sample done"<<std::endl;
   pvaluechi[8] = hEreco90->Chi2Test(hEreco902);
   std::cout<<"chi2 for 0.9 sample done"<<std::endl;
   pvaluechi[9] = hEreco->Chi2Test(hEreco2);
   std::cout<<"chi2 for tot sample done"<<std::endl;
   pvaluechitest = hEreco50->Chi2Test(hEreco50s);
   
   std::cout<<"chi2 for test sample done"<<std::endl;
  
   
   ////////////////////////////////////////PETTI TEST
   double chi2Petti = 0;
   double chi2Pettitest = 0;
   
   
   
   for (int i=1; i <= hEreco->GetNbinsX(); i++)
   {
	  double N=(hEreco->GetBinContent(i)-hEreco2->GetBinContent(i))*(hEreco->GetBinContent(i)-hEreco2->GetBinContent(i));
	  double sigma=hEreco->GetBinContent(i)+hEreco2->GetBinWidth(i);
	  chi2Petti += N / sigma;
   }
   for (int i=1; i <= hEreco50->GetNbinsX(); i++)
   {
	  double N=(hEreco50->GetBinContent(i)-hEreco50s->GetBinContent(i))*(hEreco50->GetBinContent(i)-hEreco50s->GetBinContent(i));
	  double sigma=hEreco50->GetBinContent(i)+hEreco50s->GetBinWidth(i);
	  chi2Pettitest += N / sigma;
   }
   
   Double_t x[10] = {1,0.1,0.01,0.001,0.0001,0.00001,0.0000001,0.00000001,0.000000001,0.0000000001};
   /*
   TCanvas *c0 = new TCanvas("AD","AD",200,10,900,800);
   TGraph* grAD = new TGraph(10,x,pvalueAD);
   grAD->SetTitle("");
   grAD->GetYaxis()->SetTitle("p-value");
   grAD->GetXaxis()->SetTitle("Sample fraction");
   grAD->SetMarkerColor(kBlue);
   grAD->Draw("AC*");

   TCanvas *c1 = new TCanvas("KS","KS",200,10,900,800);
   TGraph* grKS = new TGraph(10,x,pvalueKS);
   grKS->SetTitle("");
   grKS->GetYaxis()->SetTitle("p-value");
   grKS->GetXaxis()->SetTitle("Sample fraction");
   grKS->SetMarkerColor(kBlue);
   grKS->Draw("AC*");
   */
   TCanvas *c3 = new TCanvas("chi","chi",200,10,900,800);
   gPad->SetLogx();
   gPad->SetLogy();
   TGraph* grchi = new TGraph(10,x,pvaluechi);
   grchi->SetTitle("");
   grchi->GetYaxis()->SetTitle("p-value");
   grchi->GetXaxis()->SetTitle("Sample fraction");
   grchi->SetMarkerColor(kBlue);
   grchi->Draw("AC*");
   
   
   
   std::cout<< " pvalueAD = " <<  pvalueAD[9] << ";   pvalueKS = "<< pvalueKS[9]<<  ";   pvaluechi = "<< pvaluechi[9]<< std::endl;
   std::cout<< " pvalueADtest = " <<  pvalueADtest << ";   pvalueKStest = "<< pvalueKStest<<  ";   pvaluechitest = "<< pvaluechitest<< std::endl;
   
   std::cout<< " chi2Petti = " <<  chi2Petti << ";   chi2Pettitest = "<< chi2Pettitest<<  std::endl;
   //hEreco->Chi2Test(hEreco2,"P");
   //hEreco50->Chi2Test(hEreco50s,"P");
   
   TFile* fo= new TFile(fOut,"RECREATE");
   fo->cd();
   
   //c0->Write("",TObject::kOverwrite);
   //c1->Write("",TObject::kOverwrite);
   c2->Write("",TObject::kOverwrite);
   c2b->Write("",TObject::kOverwrite);
   //c3->Write("",TObject::kOverwrite);
   /*
   pvalueAD.Write("",TObject::kOverwrite);
   pvalueKS.Write("",TObject::kOverwrite);
   pvaluechi.Write("",TObject::kOverwrite);
   pvalueADtest.Write("",TObject::kOverwrite);
   pvalueKStest.Write("",TObject::kOverwrite);
   pvaluechitest.Write("",TObject::kOverwrite);
   chi2Petti.Write("",TObject::kOverwrite);
   chi2Pettitest.Write("",TObject::kOverwrite);
   
   c2->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/ErecoBoth.png");
   c0->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/pvalueAD.png");
   c1->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/pvalueKD.png");
   c3->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/pvaluechi.png");
   */
   
   
}