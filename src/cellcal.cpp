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

void cellcal(const char* fIn) 
{
	
	gSystem->Load("/mnt/c/Linux/Dune/kloe-simu/lib/libStruct.so");
	
	TFile* f= new TFile(fIn,"");
	
	auto t = (TTree*) f->Get("tDigit");
	
	
	
	
	std::vector<cell>* vec_cell = 0;
	
	t->SetBranchAddress("cell", &vec_cell);

	
	TH1D *hEcell= new TH1D("hEcell","",200,0,400);
	
	const int nev = t->GetEntries();
	
	for(int j = 0; j < nev; j++)
    {
		
		t->GetEntry(j);
		
		for(int i = 0; i < vec_cell->size(); i++)
	    {
			cell c;
			c = vec_cell->at(i);
			hEcell->Fill((c.adc1+c.adc2));
		}
		
		
	}
	
	TCanvas *c2 = new TCanvas("c2","hEcell",200,10,900,600);
	gStyle->SetOptFit();
	hEcell->GetXaxis()->SetTitle("p.e.(adc1+adc2)");
    hEcell->GetYaxis()->SetTitle("n^{#circ} of cells");
	hEcell->SetFillColor(kBlue-10);
	//TF1 *f1 = new TF1("f1", "landau", 50, 400);
	//hEcell->Fit("f1","R");
	hEcell->Draw();
	c2->Modified();
	c2->Update();
	gSystem->ProcessEvents();
	
	c2->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/hEcell.png");
	
}