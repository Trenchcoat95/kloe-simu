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
	
	//////////////////LOADING THE LIBRARIES AND TREES///////////////////////////////
	
	gSystem->Load("/mnt/c/Linux/Dune/kloe-simu/lib/libStruct.so");
	
	
	TFile* f= new TFile(fIn,"UPDATE");	
	
    TTree* tTracker = (TTree*) f->Get("gRooTracker");
	TTree* tTrueMC = (TTree*) f->Get("EDepSimEvents");
	TTree* tDigit = (TTree*) f->Get("tDigit");
	TTree* tReco = (TTree*) f->Get("tReco");
	TTree* tEvent = (TTree*) f->Get("tEvent");
	
	TTree* tout = new TTree("tLay","tLay");
	int firstlayer;
    tout->Branch("firstlayer",&firstlayer,"firstlayer/I");
	
	
	tDigit->AddFriend(tTrueMC);
    tDigit->AddFriend(tTracker);
	tDigit->AddFriend(tReco);
	tDigit->AddFriend(tEvent);
	tDigit->AddFriend(tout);
  
    TTree* t = tDigit;	
	
	const int kMaxStdHepN = 200;
	
	TG4Event* ev = new TG4Event;
	std::vector<cell>* vec_cell = new std::vector<cell>;
	std::vector<digit>* vec_stt = new std::vector<digit>;
	double part_mom[kMaxStdHepN][4];
	double evt_vtx[4];
	event* e = new event;
	int part_pdg[kMaxStdHepN];
	std::vector<track>* vec_tr = new std::vector<track>;
	std::vector<cluster>* vec_cl = new std::vector<cluster>;
	
	
	t->SetBranchAddress("Event",&ev);
	t->SetBranchAddress("cell",&vec_cell);
	t->SetBranchAddress("Stt",&vec_stt);
	t->SetBranchAddress("StdHepP4",&part_mom);
	t->SetBranchAddress("EvtVtx",&evt_vtx);
	t->SetBranchAddress("StdHepPdg",&part_pdg);
	t->SetBranchAddress("event",&e);
	t->SetBranchAddress("track",&vec_tr);
	t->SetBranchAddress("cluster",&vec_cl);
	
    TCut internalEvent = "firstlayer==0";
	
	
	TH1D* hmuonpe= new TH1D("hmuonpe","",100,0,16000);
	TH1D* hmuonEreco =  new TH1D("hmuonEreco","",100,0,16000);
	
	int ncut = 0;
	int muons = 0;
	int muonreco = 0;
	const int nev = t->GetEntries();
	
	
	
	///////////////////////////////////////Ciclo sugli eventi ///////////////////////////////////	
	
	for (int i = 0; i < nev; i++) 
	{ 
      
	  
	  
      t->GetEntry(i);
	  
	  
	  firstlayer = 0;
	  int pemuon = 0;
	  
	  for (int j = 0; j < vec_cell->size(); j++)   ///////////////////////Ciclo sulle celle per la flag sugli eventi esterni
	  {
		 
		  int layer = vec_cell->at(j).lay;
		  int mod = vec_cell->at(j).mod;
		  int cell = vec_cell->at(j).cel;
		  
		  //std::cout<<"mod id: "<<mod<<std::endl;
		  //std::cout<<"cell layer: "<<layer<<std::endl;
		  
		  double edep = vec_cell->at(j).adc1 + vec_cell->at(j).adc2;
		  
		  if (mod >= 2 && mod <=10 && layer == 4 && edep >= 75)
		  {
			  firstlayer = 1;
			  //std::cout<< "recognized as having vertex in the outer layer" << std::endl;
			  break;
		  }
		  
		  
	  }
	  
	  for (int j = 0; j < vec_cell->size(); j++)   //////////////////////Ciclo sulle celle per totale foto-elettroni prodotti dal muone: pemuon 
	    {
			for(int k = 0; k < vec_cell->at(j).hindex1.size(); k++)
			{
				int Trackid = ev->SegmentDetectors["EMCalSci"].at(vec_cell->at(j).hindex1.at(k)).GetPrimaryId();
				int pdg = ev->Trajectories.at(Trackid).GetPDGCode();
				if (pdg == 13) {pemuon ++;}
			}
			
			for(int k = 0; k < vec_cell->at(j).hindex2.size(); k++)
			{
				int Trackid = ev->SegmentDetectors["EMCalSci"].at(vec_cell->at(j).hindex2.at(k)).GetPrimaryId();
				int pdg = ev->Trajectories.at(Trackid).GetPDGCode();
				if (pdg == 13) {pemuon ++;}
			}
			
		}
	  //std::cout<<"size "<<vec_partcut->size()<<std::endl;
	  if (firstlayer == 0)
	  {
		  ncut++;
		  for (int k = 0; k < e->particles.size(); k++)	
		  {
			  if (e->particles.at(k).tr.pid == 13 ) 
			  {
				  muons++;
				  
				  if(e->particles.at(k).tr.ret_ln == 0 && e->particles.at(k).tr.ret_cr == 0)
					{
						muonreco++;
						hmuonEreco->Fill(e->particles.at(k).Ereco);   /////////////////////////////////Fill dell'istogramma dell'energia ricostruita
					}
				  else {std::cout<<"Problematic event number:"<<i<<std::endl;}
			  }
		  }		  
      }
	  if (firstlayer == 0 && pemuon != 0){hmuonpe->Fill(pemuon);}
	  
	  tout->Fill();
      
	}
	
	///////////////STAMPA INFO SUL CUT /////////////////////////
	
	
	std::cout<<"Number of events post cut:"<<ncut<<std::endl;
	std::cout<<"Number of muons:"<<muons<<std::endl;
	std::cout<<"Number of correctly reconstructed muons:"<<muonreco<<std::endl;
	
	
	///////////////D/////////////////////////
	
	TCanvas *c1 = new TCanvas("c1","c1",200,10,900,600);

	t->Draw("StdHepP4[0][3]","","", 1000, 0);
	TH1D *h = (TH1D*)gPad->GetPrimitive("htemp");
	h->SetTitle("Neutrino energy distribution");
	h->GetXaxis()->SetTitle("E(GeV/c)");
	h->GetYaxis()->SetTitle("number of events");
	c1->Modified();
	c1->Update();
	
	TCanvas *c2 = new TCanvas("c2","c2",200,10,900,600);
	
	t->SetMarkerColor(kBlack);
	t->Draw("EvtVtx[1]:EvtVtx[0]>>hyx(1,-2.5,2.5,1,0,-5)",internalEvent,"*", 1000, 0); //internalEvent
	TH1D *h1 = (TH1D*)gPad->GetPrimitive("hyx");
	h1->SetTitle("Vertex xy");
	h1->GetXaxis()->SetTitle("x(m)");
	h1->GetYaxis()->SetTitle("y(m)");
	t->SetMarkerColor(kRed);
	t->Draw("EvtVtx[1]:EvtVtx[0]",!internalEvent,"*same", 1000, 0); ////!internalEvent
	
	c2->Modified();
	c2->Update();
	
	TCanvas *c3 = new TCanvas("c3","c3",200,10,900,600);
	
	t->SetMarkerColor(kBlack);
	t->Draw("EvtVtx[1]:EvtVtx[2]>>hyz(1,20,25,1,0,-5)",internalEvent,"*", 1000, 0); //internalEvent
	TH1D *h2 = (TH1D*)gPad->GetPrimitive("hyz");
	h2->SetTitle("Vertex yz");
	h2->GetXaxis()->SetTitle("y(m)");
	h2->GetYaxis()->SetTitle("z(m)");
	t->SetMarkerColor(kRed);
	t->Draw("EvtVtx[1]:EvtVtx[2]",!internalEvent,"*same", 1000, 0); //!internalEvent
	c3->Modified();
	c3->Update();
	

	gStyle->SetOptStat(kTRUE);
	gStyle->SetOptFit(kTRUE);
	
	TCanvas *c4 = new TCanvas("c4","c4",200,10,900,600);
	hmuonpe->SetTitle("Muon p.e. production");
	hmuonpe->GetXaxis()->SetTitle("number of photoelectrons");
	hmuonpe->GetYaxis()->SetTitle("events");
	hmuonpe->SetFillColor(kBlue-10);
	hmuonpe->Fit("landau");
	hmuonpe->Draw();
	c4->Modified();
	c4->Update();
	c4->Update();
	
	TCanvas *c5 = new TCanvas("c5","c5",200,10,900,600);
	hmuonEreco->SetTitle("Muon reconstructed energy");
	hmuonEreco->GetXaxis()->SetTitle("MeV");
	hmuonEreco->GetYaxis()->SetTitle("events");
	hmuonEreco->SetFillColor(kBlue-10);
	hmuonEreco->Fit("landau");
	hmuonEreco->Draw();
	c5->Modified();
	c5->Update();
	c5->Update();
	
	/*
	TH1D* hist_new=(TH1D*)hmuonpe->Clone();
    hist_new->SetName("hist_new");
	
	hmuonpe -> Chi2Test(hist_new,"P");
	*/
	
	///////////////////////////////SALVA TUTTO///////////////////////////
	std::string Ead = "/mnt/c/Linux/Dune/kloe-simu/plots/E";
	std::string xyad = "/mnt/c/Linux/Dune/kloe-simu/plots/xy";
	std::string yzad = "/mnt/c/Linux/Dune/kloe-simu/plots/yz";
	std::string Muonpead = "/mnt/c/Linux/Dune/kloe-simu/plots/Muonpe";
	std::string Erecoad = "/mnt/c/Linux/Dune/kloe-simu/plots/Ereco";
	
	Ead += std::to_string(nev)+".png";
	xyad += std::to_string(nev)+".png";
	yzad += std::to_string(nev)+".png";
	Muonpead += std::to_string(nev)+".png";
	Erecoad += std::to_string(nev)+".png";
	
	c1->SaveAs(Ead.c_str());
	c2->SaveAs(xyad.c_str());
	c3->SaveAs(yzad.c_str());
	c4->SaveAs(Muonpead.c_str());
	c5->SaveAs(Erecoad.c_str());
	
	tout->Write("",TObject::kOverwrite);
	
	hmuonpe->Write("",TObject::kOverwrite);
	hmuonEreco->Write("",TObject::kOverwrite);
	
	
	
	
	//delete f;
	//delete fout;
	
    
    
}



