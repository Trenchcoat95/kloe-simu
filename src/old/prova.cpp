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

void prova(char* ftemplate, int n) //"numu_nom_10k_"
{
	//numu_nom_10k_0.0.ana.root
	
	std::string filename = ftemplate + std::to_string(0) + ".0.ana.png";
	TFile f(filename.c_str());
	TH1D* hEreco = (TH1D*)f.Get("hmuonErecoInner");
	
	for(int i=1; i<=n; i++)
	{
		std::string filenametemp = ftemplate + std::to_string(n) + ".0.ana.png";
		TFile ftemp(filenametemp.c_str());
		TH1D* hErecotemp = (TH1D*)ftemp.Get("hmuonErecoInner");
		hEreco->Add(hErecotemp)
	}
	
}