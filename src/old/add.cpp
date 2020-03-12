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

void add(const char* ftemplate, const int n, const char* fout) //"numu_nom_10k_"
{
	//numu_nom_10k_0.0.ana.root
	
	std::string filename = ftemplate + std::to_string(0) + ".0.ana.png";
	TFile f(filename.c_str());
	
	TH1D* hmuonErecoInner = (TH1D*)f.Get("hmuonErecoInner");
	TH1D* hmuonrecoEinner = (TH1D*)f.Get("hmuonrecoEInner");
	TH1D* hEneutrinoInner = (TH1D*)f.Get("hEneutrinoInner");
	TH1D* hEneutrinoInnerMC = (TH1D*)f.Get("hEneutrinoInnerMC");
	TH1D* hpmuon = (TH1D*)f.Get("hpmuon");
	TH1D* hpmuonreco = (TH1D*)f.Get("hpmuon");
	
	
	for(int i=1; i<=n; i++)
	{
		std::string filenametemp = ftemplate + std::to_string(n) + ".0.ana.png";
		TFile ftemp(filenametemp.c_str());
		
		TH1D* hmuonErecoInnertemp = (TH1D*)ftemp.Get("hmuonErecoInner");
		TH1D* hmuonrecoEinnertemp = (TH1D*)ftemp.Get("hmuonrecoEInner");
		TH1D* hEneutrinoInnertemp = (TH1D*)ftemp.Get("hEneutrinoInner");
		TH1D* hEneutrinoInnerMCtemp = (TH1D*)ftemp.Get("hEneutrinoInnerMC");
		TH1D* hpmuontemp = (TH1D*)ftemp.Get("hpmuon");
		TH1D* hpmuonrecotemp = (TH1D*)ftemp.Get("hpmuon");
		
		hmuonErecoInner->Add(hmuonErecoInnertemp);
		hmuonrecoEinner->Add(hmuonrecoEinnertemp);
		hEneutrinoInner->Add(hEneutrinoInnertemp);
		hEneutrinoInnerMC->Add(hEneutrinoInnerMCtemp);
		hpmuon->Add(hpmuontemp);
		hpmuonreco->Add(hpmuontemp);
	}
	
	
	TFile* fo= new TFile(fOut,"RECREATE");
    fo->cd();
	
	hmuonErecoInner->Write("",TObject::kOverwrite);
	hmuonrecoEinner->Write("",TObject::kOverwrite);
	hEneutrinoInner->Write("",TObject::kOverwrite);
	hEneutrinoInnerMC->Write("",TObject::kOverwrite);
	hpmuon->Write("",TObject::kOverwrite);
	hpmuonreco->Write("",TObject::kOverwrite);
	
	TCanvas *c1 = new TCanvas("cneutrino_E_cut","cneutrino_E_cut",200,10,900,600);
	gStyle->SetOptStat(0);
	hEneutrinoInner->SetTitle("");//Neutrino energy distribution of the events surviving the cut
	hEneutrinoInner->GetXaxis()->SetTitle("GeV/c");
	hEneutrinoInner->GetYaxis()->SetTitle("events");
	hEneutrinoInner->SetFillColor(kGray);
	hEneutrinoInner->SetLineColor(kBlack);
	hEneutrinoInner->SetFillStyle(3001);
	hEneutrinoInnerMC->SetFillColor(kYellow-10);
	hEneutrinoInnerMC->SetLineColor(kYellow);
	c1->Modified();
	c1->Update();
	hEneutrinoInnerMC->Draw();
	hEneutrinoInner->Draw("same");
	
	TLegend* legend1 = new TLegend(0.62,0.8,0.894,0.9);
	legend1->AddEntry("hEneutrinoInnerMC", "Neutrino having vertex in the inner layers", "f");
    legend1->AddEntry("hEneutrinoInner", "Neutrinos surviving the cut", "f");
	legend1->SetFillColor(0);
	legend1->SetBorderSize(1);
	legend1->Draw("");
	
	c1->Write("",TObject::kOverwrite);
	
	TCanvas *c2 = new TCanvas("cmuon_p","cmuon_p",200,10,900,600);
	gStyle->SetOptStat(0);
	hpmuon->SetTitle("");
	hpmuon->GetXaxis()->SetTitle("E [MeV]");
	hpmuon->GetYaxis()->SetTitle("events");
	hpmuon->SetFillColor(kBlue-10);
	hpmuon->SetLineColor(kBlue);
	hpmuonreco->SetFillColor(kGreen-10);
	hpmuonreco->SetLineColor(kGreen);
	c2->Modified();
	c2->Update();
	hpmuon->DrawClone();
	hpmuonreco->DrawClone("same");
	
	
	TLegend* legend4 = new TLegend(0.62,0.8,0.894,0.9);
	legend4->AddEntry("hpmuon", "All muons", "f");
    legend4->AddEntry( "hpmuonreco", "Correctly reconstructed muons", "f");
	legend4->SetFillColor(0);
	legend4->SetBorderSize(1);
	legend4->Draw("");
	
	c2->Write("",TObject::kOverwrite);
	
	
	
	
	TEfficiency* pEffReco = 0;                             ///////////////////////////////Efficienza ricostruzione
	if(TEfficiency::CheckConsistency(hpmuonreco,hpmuon))
	{
		pEffReco = new TEfficiency(hpmuonreco,hpmuon);
		TCanvas *c3 = new TCanvas("creconstruction_efficiency","creconstruction_efficiency",200,10,900,600);
		pEffReco->SetTitle(";E [MeV];#varepsilon_{reco}(E)");//("Reconstruction efficiency as a function of the muon momentum;Muon momentum [MeV];Efficiency");
		pEffReco->Draw();
		c3->Modified();
		c3->Update();
		c3->Write("",TObject::kOverwrite);
		pEffReco->Write("",TObject::kOverwrite)
		//std::string RecoEfficiencyad = "/mnt/c/Linux/Dune/kloe-simu/plots/RecoEfficiency";
		//RecoEfficiencyad += std::to_string(nev)+".png";
		//c11->SaveAs(RecoEfficiencyad.c_str());
	}
	
	TEfficiency* pEffCut = 0;                       ///////////////////////////////Efficienza Cut
	if(TEfficiency::CheckConsistency(hEneutrinoInner,hEneutrinoInnerMC))
	{
		pEffCut = new TEfficiency(hEneutrinoInner,hEneutrinoInnerMC);
		TCanvas *c4 = new TCanvas("ccut_efficiency","ccut_efficiency",200,10,900,600);
		pEffCut->SetTitle(";Neutrino Energy [GeV/c];Efficiency"); //Cut efficiency as a function of the neutrino energy
		pEffCut->Draw();
		c4->Modified();
		c4->Update();
		c4->Update();
		c4->Write("",TObject::kOverwrite);
		pEffCut->Write("",TObject::kOverwrite);
		//std::string CutEfficiencyad = "/mnt/c/Linux/Dune/kloe-simu/plots/CutEfficiency";
		//CutEfficiencyad += std::to_string(nev)+".png";
		//c8->SaveAs(CutEfficiencyad.c_str());
	}
	
}