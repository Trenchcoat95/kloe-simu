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
#include <time.h>

#include "/mnt/c/Linux/Dune/edep-sim/edep-gcc-7-x86_64-linux-gnu/include/EDepSim/TG4Event.h"
#include "/mnt/c/Linux/Dune/edep-sim/edep-gcc-7-x86_64-linux-gnu/include/EDepSim/TG4HitSegment.h"

#include "/mnt/c/Linux/Dune/kloe-simu/include/struct.h"
#include "/mnt/c/Linux/Dune/kloe-simu/include/utils.h"

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>

void draw(const char* fin, const char* fout) 
{
	//numu_nom_10k_0.0.ana.root
	TChain* tout = new TChain("tout");
	tout->Add(fin);
	
	double Enu_true, Enu_reco, xv_true, yv_true, zv_true, xv_reco, pxmu_true, pymu_true, pzmu_true, pxmu_reco, pymu_reco, pzmu_reco, maxde_frontoutlayer;
	int isCC, isPmuOK;
	
	tout->SetBranchAddress("Enu_true",&Enu_true);
	tout->SetBranchAddress("Enu_reco",&Enu_reco);
	tout->SetBranchAddress("xv_true",&xv_true);
	tout->SetBranchAddress("yv_true",&yv_true);
	tout->SetBranchAddress("zv_true",&zv_true);
	tout->SetBranchAddress("xv_reco",&zv_true);
	tout->SetBranchAddress("pxmu_true",&pxmu_true);
	tout->SetBranchAddress("pymu_true",&pymu_true);
	tout->SetBranchAddress("pzmu_true",&pzmu_true);
	tout->SetBranchAddress("pxmu_reco",&pxmu_reco);
	tout->SetBranchAddress("pymu_reco",&pymu_reco);
	tout->SetBranchAddress("pzmu_reco",&pzmu_reco);
	tout->SetBranchAddress("isCC",&isCC);
	tout->SetBranchAddress("isPmuOK",&isPmuOK);
	tout->SetBranchAddress("maxde_frontoutlayer",&maxde_frontoutlayer);
	
	
	vector<double> ptrue_muonreco_In;
	vector<double> preco_muonreco_In;
	
	TH1D hptrue_muon ("hptrue_muon","",100,0,16000);
	TH1D hptrue_muonreco ("hptrue_muonreco","",100,0,16000);
	TH1D hptrue_muonreco_In ("hptrue_muonreco_In","",100,0,16000);
	TH1D hpreco_muonreco_In ("hpreco_muonreco_In","",100,0,16000);
	TH1D hEtrue_nu ("hEtrue_nu","",100,0,40000);
	TH1D hEtrue_nu_In ("hEtrue_nu_In","",100,0,40000);
	
	
	const int nev = tout->GetEntries();

	for(int i=0; i<nev; i++)
	{
		tout->GetEntry(i);
		
		hEtrue_nu.Fill(Enu_true);
		//std::cout<< maxde_frontoutlayer<<std::endl;
		
		if(maxde_frontoutlayer < 40 && xv_reco <= 1690 && xv_reco >= -1690)
		{
			hEtrue_nu_In.Fill(Enu_true);
		}
		
		if(isCC>0)
		{
			double ptot = std::sqrt(pxmu_true*pxmu_true + pymu_true*pymu_true + pzmu_true*pzmu_true);
			hptrue_muon.Fill(ptot);
			
			double ptot_reco = std::sqrt(pxmu_reco*pxmu_reco + pymu_reco*pymu_reco + pzmu_reco*pzmu_reco);
			if(isPmuOK>0)
			{
				hptrue_muonreco.Fill(ptot);
				
				
				if (maxde_frontoutlayer < 40 && xv_reco <= 1690 && xv_reco >= -1690)
				{
					
					hptrue_muonreco_In.Fill(ptot);
					hpreco_muonreco_In.Fill(ptot_reco);
					ptrue_muonreco_In.push_back(ptot);
					preco_muonreco_In.push_back(ptot_reco);
				}
				
			}
		}
		
	}
	
	
	TFile* fo= new TFile(fout,"RECREATE");
    fo->cd();
	
	hptrue_muon.Write("",TObject::kOverwrite);
	hptrue_muonreco.Write("",TObject::kOverwrite);
	hptrue_muonreco_In.Write("",TObject::kOverwrite);
	hpreco_muonreco_In.Write("",TObject::kOverwrite);
	hEtrue_nu.Write("",TObject::kOverwrite);
	hEtrue_nu_In.Write("",TObject::kOverwrite);
	
	////////////////////////////////////CUT Plots//////////////////////////////////////////////////////////////////////
	
	TCanvas *c1 = new TCanvas("cneutrino_E_cut","cneutrino_E_cut",200,10,900,600); ///////////Energy distribution of all neutrinos + neutrinos surviving the cut
	gStyle->SetOptStat(0);
	hEtrue_nu.SetTitle("");//Neutrino energy distribution 
	hEtrue_nu.GetXaxis()->SetTitle("E_{#nu}[MeV]");
	hEtrue_nu.GetYaxis()->SetTitle("events");
	hEtrue_nu.SetFillColor(kGray);
	hEtrue_nu.SetLineColor(kBlack);
	hEtrue_nu.SetFillStyle(3001);
	hEtrue_nu_In.SetFillColor(kYellow-10);
	hEtrue_nu_In.SetLineColor(kYellow);
	hEtrue_nu_In.SetFillStyle(3001);
	c1->Modified();
	c1->Update();
	hEtrue_nu.DrawClone();
	hEtrue_nu_In.DrawClone("same");
	
	TLegend* legend1 = new TLegend(0.62,0.8,0.894,0.9);
	legend1->AddEntry("hEtrue_nu", "All neutrinos", "f");
    legend1->AddEntry("hEtrue_nu_In", "Neutrinos surviving the cut", "f");
	legend1->SetFillColor(0);
	legend1->SetBorderSize(1);
	legend1->Draw("");
	
	c1->Write("",TObject::kOverwrite);
	
	TEfficiency* pEffCut = 0;                       ///////////////////////////////Efficienza Cut
	if(TEfficiency::CheckConsistency(hEtrue_nu_In,hEtrue_nu))
	{
		pEffCut = new TEfficiency(hEtrue_nu_In,hEtrue_nu);
		TCanvas *c1a = new TCanvas("ccut_efficiency","ccut_efficiency",200,10,900,600);
		pEffCut->SetTitle(";E_{#nu}[MeV];Efficiency"); //Cut efficiency as a function of the neutrino energy
		pEffCut->Draw();
		c1a->Modified();
		c1a->Update();
		c1a->Update();
		c1a->Write("",TObject::kOverwrite);
		//pEffCut->Write("",TObject::kOverwrite);
		//std::string CutEfficiencyad = "/mnt/c/Linux/Dune/kloe-simu/plots/CutEfficiency";
		//CutEfficiencyad += std::to_string(nev)+".png";
		//c8->SaveAs(CutEfficiencyad.c_str());
	}
	
	///////////////////////////////////////////////////////RECO PLOTS///////////////////////////////////////////////////////////////////////
	
	TCanvas *c2 = new TCanvas("cmuon_p","cmuon_p",200,10,900,600); ////////////// true p + reconstructed p of muons surviving the cut
	gStyle->SetOptStat(0);
	hptrue_muonreco_In.SetTitle("");
	hptrue_muonreco_In.GetXaxis()->SetTitle("p_{#mu}[MeV/c]");
	hptrue_muonreco_In.GetYaxis()->SetTitle("events");
	hptrue_muonreco_In.SetFillColor(kBlue-10);
	hptrue_muonreco_In.SetLineColor(kBlue);
	hptrue_muonreco_In.SetFillStyle(3001);
	hpreco_muonreco_In.SetFillColor(kGreen-10);
	hpreco_muonreco_In.SetLineColor(kGreen);
	hpreco_muonreco_In.SetFillStyle(3001);
	c2->Modified();
	c2->Update();
	hpreco_muonreco_In.DrawClone();
	hptrue_muonreco_In.DrawClone("same");
	
	
	TLegend* legend2 = new TLegend(0.62,0.8,0.894,0.9);
	legend2->AddEntry("hptrue_muonreco_In", "True p", "f");
    legend2->AddEntry( "hpreco_muonreco_In", "Reconstructed p", "f");
	legend2->SetFillColor(0);
	legend2->SetBorderSize(1);
	legend2->Draw("");
	
	c2->Write("",TObject::kOverwrite);
	
	
	
	
	Double_t ptrue[ptrue_muonreco_In.size()],preco[preco_muonreco_In.size()];
	for(int i=0;i<ptrue_muonreco_In.size();i++)
	{
		ptrue[i]=ptrue_muonreco_In.at(i)*0.001;
		preco[i]=preco_muonreco_In.at(i)*0.001;
	}
	TCanvas *c2a = new TCanvas("ptrueVSpreco","ptrueVSpreco",200,10,900,700); ////////////// true p VS reconstructed p of muons surviving the cut
	//gStyle->SetOptFit();
    TGraph* grreco = new TGraph(ptrue_muonreco_In.size(),ptrue,preco);
    grreco->SetTitle("");
    grreco->GetYaxis()->SetTitle("p^{reco}_{#mu} [GeV/c]");
    grreco->GetXaxis()->SetTitle("p^{true}_{#mu} [GeV/c]");
    grreco->SetMarkerColor(kBlue);
	//grreco->Fit("pol1");
    grreco->Draw("A*"); //P
	
	c2a->Write("",TObject::kOverwrite);
	
	
	TEfficiency* pEffReco = 0;                             ///////////////////////////////Efficienza ricostruzione
	if(TEfficiency::CheckConsistency(hptrue_muonreco,hptrue_muon))
	{
		pEffReco = new TEfficiency(hptrue_muonreco,hptrue_muon);
		TCanvas *c3 = new TCanvas("creconstruction_efficiency","creconstruction_efficiency",200,10,900,600);
		pEffReco->SetTitle(";E [MeV];#varepsilon_{reco}(E)");//("Reconstruction efficiency as a function of the muon momentum;Muon momentum [MeV];Efficiency");
		pEffReco->Draw();
		c3->Modified();
		c3->Update();
		c3->Write("",TObject::kOverwrite);
		pEffReco->Write("",TObject::kOverwrite);
		//std::string RecoEfficiencyad = "/mnt/c/Linux/Dune/kloe-simu/plots/RecoEfficiency";
		//RecoEfficiencyad += std::to_string(nev)+".png";
		//c11->SaveAs(RecoEfficiencyad.c_str());
	}
	
	
}