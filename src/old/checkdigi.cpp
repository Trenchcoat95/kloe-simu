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

void checkdigi(const char* fIn, const char* fOut)
{
	gSystem->Load("/mnt/c/Linux/Dune/kloe-simu/lib/libStruct.so");
	
	
	TFile f(fIn,"UPDATE");	
	TTree* tTrueMC = (TTree*) f.Get("EDepSimEvents");
	TTree* tDigit = (TTree*) f.Get("tDigit");
	
	
	tDigit->AddFriend(tTrueMC);
    
  
    TTree* t = tDigit;	
	
	
	
	TG4Event* ev = new TG4Event;
	std::vector<cell>* vec_cell = new std::vector<cell>;
	std::vector<digit>* vec_stt = new std::vector<digit>;
	
	
	
	t->SetBranchAddress("Event",&ev);
	t->SetBranchAddress("cell",&vec_cell);
	t->SetBranchAddress("Stt",&vec_stt);
	
	
	const int nev = t->GetEntries();
	
	
	
	
	
	
	
	TH1D hdistYZ("hdistYZ","",1000,0,400);
	TH2D hrelposYZ("hrelposYZ","",10000,-400,400,10000,-400,400);
	TH2D habsposYZ("habsposYZ","",10000,20000,30000,10000,-5000,1000);
	TH2D hcellposYZ("hcellposYZ","",10000,20000,30000,10000,-5000,1000);
	
	TH1D hdistXZ("hdistXZ","",1000,0,400);
	TH2D hrelposXZ("hrelposXZ","",10000,-400,400,10000,-400,400);
	TH2D habsposXZ("habsposXZ","",10000,20000,30000,10000,-5000,1000);
	TH2D hcellposXZ("hcellposXZ","",10000,20000,30000,10000,-5000,1000);
	
	
	
	for (int i = 0; i < nev; i++) 
	{
      t->GetEntry(i);
	  std::cout<<"Entry number: "<<i<<std::endl;
	  
	  
	  
	    for (int j = 0; j < vec_cell->size(); j++)
	    {
		  int layer = vec_cell->at(j).lay;
		  int mod = vec_cell->at(j).mod;
		  int cell = vec_cell->at(j).cel;
		  
		    for(int k = 0; k < vec_cell->at(j).hindex1.size(); k++)
			{
			  double hx = ev->SegmentDetectors["EMCalSci"].at(vec_cell->at(j).hindex1.at(k)).GetStart().X();
			  double hy = ev->SegmentDetectors["EMCalSci"].at(vec_cell->at(j).hindex1.at(k)).GetStart().Y();
			  double hz = ev->SegmentDetectors["EMCalSci"].at(vec_cell->at(j).hindex1.at(k)).GetStart().Z();
			  double distX = vec_cell->at(j).x - hx;
			  double distY = vec_cell->at(j).y - hy;
			  double distZ = vec_cell->at(j).z - hz;
			  
			  /*std::cout << "hindex1: " << mod << " " << layer << " " << cell << " " << 
			  vec_cell->at(j).x << " " << 
			  vec_cell->at(j).y << " " << 
			  vec_cell->at(j).z << " " << 
			  k << " " << vec_cell->at(j).hindex1.at(k) << " " << 
			  ev->SegmentDetectors["EMCalSci"].at(vec_cell->at(j).hindex1.at(k)).GetStart().X() << " " <<
			  ev->SegmentDetectors["EMCalSci"].at(vec_cell->at(j).hindex1.at(k)).GetStart().Y() << " " <<
			  ev->SegmentDetectors["EMCalSci"].at(vec_cell->at(j).hindex1.at(k)).GetStart().Z() << std::endl;*/
			  
			  double dist2DYZ = sqrt(distY*distY + distZ*distZ);
			  double dist2DXZ = sqrt(distX*distX + distZ*distZ);
			  
			  if(mod < 24)
			  {
				hdistYZ.Fill(dist2DYZ);
				hrelposYZ.Fill(distZ,distY);
				habsposYZ.Fill(hz,hy);
				hcellposYZ.Fill(vec_cell->at(j).z,vec_cell->at(j).y);
			  }
			  else if (mod == 30 || mod == 40)
			  {
				hdistXZ.Fill(dist2DXZ);
				hrelposXZ.Fill(distZ,distX);
				habsposXZ.Fill(hx,hy);
				hcellposXZ.Fill(vec_cell->at(j).z,vec_cell->at(j).x);
			  }
			  
			}
	    }
  
      
      
	}
   
	
	f.Close() ;
	
	TFile fout(fOut,"RECREATE");
	hdistYZ.Write();
	hrelposYZ.Write();
	habsposYZ.Write();
	hcellposYZ.Write();
	hdistXZ.Write();
	hrelposXZ.Write();
	habsposXZ.Write();
	hcellposXZ.Write();
	fout.Close();
	
    
    
}




