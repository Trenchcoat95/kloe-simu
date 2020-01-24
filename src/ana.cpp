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

void ana(const char* fIn)
{
	gSystem->Load("/mnt/c/Linux/Dune/kloe-simu/lib/libStruct.so");
	
	//TFile f(("");
	TFile f(fIn,"UPDATE");	
    TTree* tTracker = (TTree*) f.Get("gRooTracker");
	TTree* tReco = (TTree*) f.Get("tReco");
	TTree* tTrueMC = (TTree*) f.Get("EDepSimEvents");
	TTree* tDigit = (TTree*) f.Get("tDigit");
	TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry");
	
	tReco->AddFriend(tTrueMC);
    tReco->AddFriend(tTracker);
	tReco->AddFriend(tDigit);
  
    TTree* t = tReco;	
	
	const int kMaxStdHepN = 200;
	
	TG4Event* ev = new TG4Event;
	std::vector<track>* vec_tr = new std::vector<track>;
	std::vector<cluster>* vec_cl = new std::vector<cluster>;
	std::vector<cell>* vec_cell = new std::vector<cell>;
	std::vector<digit>* vec_stt = new std::vector<digit>;
	double part_mom[kMaxStdHepN][4];
	double evt_vtx[4];
	int part_pdg[kMaxStdHepN];
	
	
	t->SetBranchAddress("Event",&ev);
	t->SetBranchAddress("track",&vec_tr);
	t->SetBranchAddress("cluster",&vec_cl);
	t->SetBranchAddress("cell",&vec_cell);
	t->SetBranchAddress("Stt",&vec_stt);
	t->SetBranchAddress("StdHepP4",part_mom);
	t->SetBranchAddress("EvtVtx",evt_vtx);
	t->SetBranchAddress("StdHepPdg",part_pdg);
	
	const int nev = t->GetEntries();
	
	double centerKLOE[3];
	double dummyLoc[3];

	geo->cd("volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/volSTTFULL_PV_0");

	dummyLoc[0] = 0.;
	dummyLoc[1] = 0.;
	dummyLoc[2] = 0.;
	geo->LocalToMaster(dummyLoc, centerKLOE);
	
	TTree *tcut = t->CloneTree(0);
	
	TH2D hedeplay("hedeplay","",5,0,5,50,0,10000);
	TH2D hedepradius("hedepradius","",50,2000,2235,100,0,1000);
	
	TH1D hdistYZ("hdistYZ","",1000,0,400);
	TH2D hrelposYZ("hrelposYZ","",10000,-400,400,10000,-400,400);
	TH2D habsposYZ("habsposYZ","",10000,20000,30000,10000,-5000,1000);
	TH2D hcellposYZ("hcellposYZ","",10000,20000,30000,10000,-5000,1000);
	
	TH1D hdistXZ("hdistXZ","",1000,0,400);
	TH2D hrelposXZ("hrelposXZ","",10000,-400,400,10000,-400,400);
	TH2D habsposXZ("habsposXZ","",10000,20000,30000,10000,-5000,1000);
	TH2D hcellposXZ("hcellposXZ","",10000,20000,30000,10000,-5000,1000);
	
	int nlay = -1;
	double edepmax = 0.;
	
	for (int i = 0; i < nev; i++) 
	{
      t->GetEntry(i);
	  std::cout<<"Entry number: "<<i<<std::endl;
	  
	  bool firstlayer = false;
	  
	  //std::cout<<"Number of cells hit: "<<vec_cell->size()<<std::endl;
	  
	  double dy = evt_vtx[1]*1000 - centerKLOE[1];
	  double dz = evt_vtx[2]*1000 - centerKLOE[2];
	  double radius = sqrt(dy*dy+dz*dz);
	  double dr = radius - 2000.; 
	  
	  //std::cout << evt_vtx[1] << " " << evt_vtx[2] << " " << centerKLOE[1] << " " << centerKLOE[2] << std::endl;
	  
	  nlay = -1;
	  
	  if(dr < 44.)
	  {
		  nlay = 0;
	  }
	  else if(dr >= 44. && dr < 88.)
	  {
		  nlay = 1;
	  }
	  else if(dr >= 88. && dr < 132.)
	  {
		  nlay = 2;
	  }
	  else if(dr >= 132. && dr < 176.)
	  {
		  nlay = 3;
	  }
	  else if(dr >= 176. && dr < 235.)
	  {
		  nlay = 4;
	  }
	  else
	  {
		  std::cout << "dr = " << dr << std::endl;
	  }

	  edepmax = 0.;
	  
	  for (int j = 0; j < vec_cell->size(); j++)
	  {
		  int layer = vec_cell->at(j).lay;
		  int mod = vec_cell->at(j).mod;
		  int cell = vec_cell->at(j).cel;
		  
		  if(layer == 4 && mod >= 2 && mod <=10 && vec_cell->at(j).adc1 + vec_cell->at(j).adc2 > edepmax)
		  {
			  edepmax = vec_cell->at(j).adc1 + vec_cell->at(j).adc2;
			  //std::cout << "entrato: " << vec_cell->at(j).adc1 + vec_cell->at(j).adc2 << std::endl;
		  }
		  
		  //std::cout<<"mod id: "<<mod<<std::endl;
		  //std::cout<<"cell layer: "<<layer<<std::endl;
		  
		  double edep = vec_cell->at(j).adc1 + vec_cell->at(j).adc2;
		  
		  if (mod >= 2 && mod <=10)
		  {
			if (layer == 4 && edep >= 35)
			{
			  firstlayer = true;
			  //break;
			}
		  }
		  
		  
		  
		  if(true)
		  {
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
  
      if (firstlayer == false) {tcut->Fill();}
	  std::cout << nlay << " " << edepmax << std::endl;
	  hedeplay.Fill(nlay,edepmax);
	  hedepradius.Fill(radius,edepmax);
      
	}
   
	TCanvas *c1 = new TCanvas("c1","c1",200,10,900,600);

	tcut->Draw("StdHepP4[0][3]","","", 1000, 0);
	TH1D *h = (TH1D*)gPad->GetPrimitive("htemp");
	h->SetTitle("Neutrino energy distribution");
	h->GetXaxis()->SetTitle("E(GeV/c)");
	h->GetYaxis()->SetTitle("number of events");
	c1->Modified();
	c1->Update();
	
	TCanvas *c2 = new TCanvas("c2","c2",200,10,900,600);
	
	tcut->Draw("EvtVtx[1]:EvtVtx[0]","","", 1000, 0);
	
	
	
	TH1D *h1 = (TH1D*)gPad->GetPrimitive("htemp");
	h1->SetTitle("Vertex xy");
	h1->GetXaxis()->SetTitle("x(m)");
	h1->GetYaxis()->SetTitle("y(m)");
	
	c2->Modified();
	c2->Update();
	
	TCanvas *c3 = new TCanvas("c3","c3",200,10,900,600);
	
	
	tcut->Draw("EvtVtx[2]:EvtVtx[1]>>hyx(1,-10,10,1,0,50)","","", 1000, 0);
	TH1D *h2 = (TH1D*)gPad->GetPrimitive("hyx");
	h2->SetTitle("Vertex yz");
	h2->GetXaxis()->SetTitle("y(m)");
	h2->GetYaxis()->SetTitle("z(m)");
	c3->Modified();
	c3->Update();
	
	
	c1->SaveAs("/mnt/c/Linux/Dune/kloe-simu/E.png");
	c2->SaveAs("/mnt/c/Linux/Dune/kloe-simu/xy.png");
	c3->SaveAs("/mnt/c/Linux/Dune/kloe-simu/yz.png");
	
	
	
	TFile fout("ana.root","RECREATE");
	hedeplay.Write();
	hedepradius.Write();
	hdistYZ.Write();
	hrelposYZ.Write();
	habsposYZ.Write();
	hcellposYZ.Write();
	hdistXZ.Write();
	hrelposXZ.Write();
	habsposXZ.Write();
	hcellposXZ.Write();
	fout.Close();
	
    
    f.Close();
}


void help_ana()
{
  std::cout << "Ana <input file>" << std::endl;
} 

int main(int argc, char* argv[])
{
  gSystem->Load("libStruct.so");
  if(argc != 2)
    help_ana();
  else
    ana(argv[1]);
}

