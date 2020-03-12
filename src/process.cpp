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

#include "/mnt/c/Linux/Dune/edep-sim/edep-gcc-7-x86_64-linux-gnu/include/EDepSim/TG4Event.h"
#include "/mnt/c/Linux/Dune/edep-sim/edep-gcc-7-x86_64-linux-gnu/include/EDepSim/TG4HitSegment.h"

#include "/mnt/c/Linux/Dune/kloe-simu/include/struct.h"
#include "/mnt/c/Linux/Dune/kloe-simu/include/utils.h"

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>

const double m_to_mm = 1000.;

void CellXYZTE(cell c, double& x, double& y, double& z, double& t, double& e)
{      
  if(c.id < 25000) //Barrel
  {
    x = 0.5 * (c.tdc1 - c.tdc2)/ns_Digit::vlfb * m_to_mm + c.x;
    y = c.y;
  }
  else
  {
    x = c.x;
    y = -0.5 * (c.tdc1 - c.tdc2)/ns_Digit::vlfb * m_to_mm + c.y;
  }
  z = c.z;
  t = 0.5 * (c.tdc1 + c.tdc2 - ns_Digit::vlfb * c.l / m_to_mm );
  e = c.adc1 + c.adc2;
}


void ana(const char* fIn, const char* fOut)    ///////(First file, List of file in directory using wildcard *, output file)
{
	int change = 1;
	
	//////////////////LOADING THE LIBRARIES AND TREES///////////////////////////////
	
	gSystem->Load("/mnt/c/Linux/Dune/kloe-simu/lib/libStruct.so");
	
	TFile* fo= new TFile(fOut,"RECREATE");
	fo->cd();
	
	TChain* tTrueMC = new TChain("EDepSimEvents");
	tTrueMC->Add(fIn);
	TChain* tDigit = new TChain("tDigit");
	tDigit->Add(fIn);
	TChain* tEvent = new TChain("tEvent");
	tEvent->Add(fIn);
	
	TTree* tout = new TTree("tout","tout");
        double Enu_true;
	double Enu_reco;
	double xv_true, yv_true, zv_true;
        double xv_reco;
        double maxde_frontoutlayer;
        int isCC;
	double pxmu_true, pymu_true, pzmu_true;
	double pxmu_reco, pymu_reco, pzmu_reco;
	int isPmuOK;

	tout->Branch("Enu_true",&Enu_true,"Enu_true/D");
	tout->Branch("Enu_reco",&Enu_reco,"Enu_reco/D");
	tout->Branch("xv_true",&xv_true,"xv_true/D");
	tout->Branch("yv_true",&yv_true,"yv_true/D");
	tout->Branch("zv_true",&zv_true,"zv_true/D");
        tout->Branch("xv_reco",&xv_reco,"xv_reco/D");
	tout->Branch("isCC",&isCC,"isCC/I");
	tout->Branch("pxmu_true",&pxmu_true,"pxmu_true/D");
	tout->Branch("pymu_true",&pymu_true,"pymu_true/D");
	tout->Branch("pzmu_true",&pzmu_true,"pzmu_true/D");
	tout->Branch("pxmu_reco",&pxmu_reco,"pxmu_reco/D");
	tout->Branch("pymu_reco",&pymu_reco,"pymu_reco/D");
	tout->Branch("pzmu_reco",&pzmu_reco,"pzmu_reco/D");
	tout->Branch("isPmuOK",&isPmuOK,"isPmuOK/I");
	tout->Branch("maxde_frontoutlayer",&maxde_frontoutlayer,"maxde_frontoutlayer/D");
	
	tDigit->AddFriend(tTrueMC);
	tDigit->AddFriend(tEvent);
  
	TTree* t = tDigit;	
	
	
	TG4Event* ev = new TG4Event;
	std::vector<cell>* vec_cell = new std::vector<cell>;
	event* e = new event;
	
	t->SetBranchAddress("Event",&ev);
	t->SetBranchAddress("cell",&vec_cell);
	t->SetBranchAddress("event",&e);

	const int nev = t->GetEntries();

	double x,y,z,tm,de,etot;

	std::cout << "Events: " << nev << " [";
	std::cout << std::setw(3) << int(0) << "%]" << std::flush;

	for (int i = 0; i < nev; i++) 
	{ 
      
	  t->GetEntry(i);

	  std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
	  
	  Enu_true = e->Enu;
	  Enu_reco = e->Enureco;

	  xv_true = e->x;
	  yv_true = e->y;
	  zv_true = e->z;

	  xv_reco = 0.;

	  isPmuOK = 0;
	  isCC = 0;

	  pxmu_true = 0;
	  pymu_true = 0;
	  pzmu_true = 0;

          pxmu_reco = 0;
          pymu_reco = 0;
          pzmu_reco = 0;

	  TString reaction = ev->Primaries.at(0).Reaction;

	  maxde_frontoutlayer = 0.;
	  etot = 0.;

	  if(reaction.Contains("CC"))
	  {
	    isCC = 1;
	  }

	  for (unsigned int j = 0; j < e->particles.size(); j++)
	  {
	    if (e->particles.at(j).primary == 1 && e->particles.at(j).pdg == 13)
	    {
	      pxmu_true = e->particles.at(j).pxtrue;
	      pymu_true = e->particles.at(j).pytrue;
	      pzmu_true = e->particles.at(j).pztrue;

	      pxmu_reco = e->particles.at(j).pxreco;
	      pymu_reco = e->particles.at(j).pyreco;
	      pzmu_reco = e->particles.at(j).pzreco;

	      if(e->particles.at(j).tr.ret_cr == 0 && e->particles.at(j).tr.ret_ln == 0)
	      {
	        isPmuOK = 1;
	      }

	      break;
	    }
	  }

	  for (int j = 0; j < vec_cell->size(); j++)   ///////////////////////Ciclo sulle celle per la flag sugli eventi esterni
	  {
		 
		  int layer = vec_cell->at(j).lay;
		  int mod = vec_cell->at(j).mod;
		  int cell = vec_cell->at(j).cel;
		  
	          CellXYZTE(vec_cell->at(j), x, y, z, tm, de);

		  if (mod >= 2 && mod <=10 && layer == 4 && de >= maxde_frontoutlayer)
		  {
			  maxde_frontoutlayer = de;
		  }
	          
	          if (mod >= 2 && mod <=10)
	          {
		  	  xv_reco += de * x;
			  etot += de;
	          }	  
		  
	  }

	  xv_reco /= etot;
	  xv_reco *= -1.;
	 
	  tout->Fill();
      
	}

	std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
	std::cout << std::endl;
	
	tout->Write();
	fo->Close();
}


int main(int argc, char* argv[])
{
	if(argc < 3)
	  return -1;

	ana(argv[1],argv[2]);
	return 0;	
}
