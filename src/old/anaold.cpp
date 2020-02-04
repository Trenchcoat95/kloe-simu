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

void ana(const char* fIn, const char* fOut)
{
	
	//////////////////LOADING THE LIBRARIES AND TREES///////////////////////////////
	
	gSystem->Load("/mnt/c/Linux/Dune/kloe-simu/lib/libStruct.so");
	
	
	TFile* f= new TFile(fIn,"");	
    TTree* tTracker = (TTree*) f->Get("gRooTracker");
	TTree* tTrueMC = (TTree*) f->Get("EDepSimEvents");
	TTree* tDigit = (TTree*) f->Get("tDigit");
	TGeoManager* geo = (TGeoManager*) f->Get("EDepSimGeometry");
	TTree* tReco = (TTree*) f->Get("tReco");
	TTree* tEvent = (TTree*) f->Get("tEvent");
	//if (!tEvent){std::cout<<"Tree brutto"<<std::endl;}
	
	tDigit->AddFriend(tTrueMC);
    tDigit->AddFriend(tTracker);
	tDigit->AddFriend(tReco);
	tDigit->AddFriend(tEvent);
  
    TTree* t = tDigit;	
	
	const int kMaxStdHepN = 200;
	
	TG4Event* ev = new TG4Event;
	std::vector<cell>* vec_cell = new std::vector<cell>;
	std::vector<digit>* vec_stt = new std::vector<digit>;
	double part_mom[kMaxStdHepN][4];
	double evt_vtx[4];
	std::vector<particle>* vec_part = new std::vector<particle>;
	event* e = new event;
	int part_pdg[kMaxStdHepN];
	std::vector<track>* vec_tr = new std::vector<track>;
	std::vector<cluster>* vec_cl = new std::vector<cluster>;
	int firstlayer;
	
	
	t->SetBranchAddress("Event",&ev);
	t->SetBranchAddress("cell",&vec_cell);
	t->SetBranchAddress("Stt",&vec_stt);
	t->SetBranchAddress("StdHepP4",&part_mom);
	t->SetBranchAddress("EvtVtx",&evt_vtx);
	t->SetBranchAddress("StdHepPdg",&part_pdg);
	t->SetBranchAddress("event",&e);
	t->SetBranchAddress("track",&vec_tr);
	t->SetBranchAddress("cluster",&vec_cl);
	
	const int nev = t->GetEntries();
	
	
	TFile* fout = new TFile(fOut,"RECREATE");
	TTree *tcutDigit = tDigit->CloneTree(0);
	TTree *tcutTracker = tTracker->CloneTree(0);
	TTree *tcutTrueMC = tTrueMC->CloneTree(0);
	TTree *tcutReco = tReco->CloneTree(0);
	//std::cout<<"Tree brutto"<<std::endl;
	TTree *tcutEvent = tEvent->CloneTree(0);
	//std::cout<<"check"<<std::endl;
	
	tcutDigit->SetDirectory(0);
	tcutTracker->SetDirectory(0);
	tcutTrueMC->SetDirectory(0);
	tcutReco->SetDirectory(0);
	tcutEvent->SetDirectory(0);
	
	tcutDigit->Branch("firstlayer",&firstlayer,"firstlayer/I");
		
	/*
	double centerKLOE[3];
	double dummyLoc[3];

	geo->cd("volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/volSTTFULL_PV_0");

	dummyLoc[0] = 0.;
	dummyLoc[1] = 0.;
	dummyLoc[2] = 0.;
	geo->LocalToMaster(dummyLoc, centerKLOE);
	
	
	
	TH2D hedeplay("hedeplay","",5,0,5,50,0,10000);
	TH2D hedepradius("hedepradius","",50,2000,2235,100,0,1000);
	
	
	int nlay = -1;
	double edepmax = 0.;
	
	int correctcut = 0;
	int wrongcut = 0;
	int nev_outerlayer = 0;
	int nev_cut = 0;
	*/
	for (int i = 0; i < nev; i++) 
	{
	  if (i==104 || i==105) {std::cout<<"Check 00 "<<i<<std::endl;}
      t->GetEntry(i);
	  //std::cout<<"Entry number: "<<i<<std::endl;
	  
	  firstlayer = 0;
	  /*
	  
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
	  else if(dr >= 176. && dr < 245.)
	  {
		  nlay = 4;
	  }
	  else
	  {
		  std::cout << "dr = " << dr << std::endl;
	  }

	  edepmax = 0.;
	  
	  //std::cout<< "nlay= "<< nlay << std::endl;
	  */
	  if (i==104 || i==105) {std::cout<<"Check 01 "<<i<<std::endl;}
	  for (int j = 0; j < vec_cell->size(); j++)
	  {
		  if (i==104 || i==105) {std::cout<<"Check 1 "<<j<<std::endl;}
		  int layer = vec_cell->at(j).lay;
		  int mod = vec_cell->at(j).mod;
		  int cell = vec_cell->at(j).cel;
		  if (i==104 || i==105) {std::cout<<"Check 2 "<<j<<std::endl;}
		  /*
		  if(layer == 4 && mod >= 2 && mod <=10 && vec_cell->at(j).adc1 + vec_cell->at(j).adc2 > edepmax)
		  {
			  edepmax = vec_cell->at(j).adc1 + vec_cell->at(j).adc2;
			  //std::cout << "entrato: " << vec_cell->at(j).adc1 + vec_cell->at(j).adc2 << std::endl;
		  }
		  */
		  
		  //std::cout<<"mod id: "<<mod<<std::endl;
		  //std::cout<<"cell layer: "<<layer<<std::endl;
		  
		  double edep = vec_cell->at(j).adc1 + vec_cell->at(j).adc2;
		  if (i==104 || i==105) {std::cout<<"Check 3 "<<j<<std::endl;}
		  if (mod >= 2 && mod <=10 && layer == 4 && edep >= 75)
		  {
			  firstlayer = 1;
			  //std::cout<< "recognized as having vertex in the outer layer" << std::endl;
			  //break;
		  }
		  if (i==104 || i==105) {std::cout<<"Check 4 "<<j<<std::endl;}
		  
		  
	  }
	  if (i==104 || i==105) {std::cout<<"Check 02 "<<i<<std::endl;}
	  
	  tcutDigit->Fill();
	  if (i==104 || i==105) {std::cout<<"Check 03 "<<i<<std::endl;}
	  //tcutTrueMC->Fill();
	  if (i==104 || i==105) {std::cout<<"Check 04 "<<i<<std::endl;}
	  tcutTracker->Fill();
	  if (i==104 || i==105) {std::cout<<"Check 05 "<<i<<std::endl;}
	  tcutReco->Fill();
	  if (i==104 || i==105) {std::cout<<"Check 06 "<<i<<std::endl;}
	  tcutEvent->Fill();
	  
	  if (i==106 || i==105) {std::cout<<"Check 07 "<<i<<std::endl;}
	  std::cout<<"Check 07 "<<i<<std::endl;
	  //std::cout<<"Tree brutto"<<std::endl;
	  /*
	  if (nlay == 4)
	  {
		  nev_outerlayer ++;
		  if (firstlayer == true) {correctcut ++;}
	  }
	  
	  if (nlay != 4 && firstlayer == true) {wrongcut ++;}
	  
	  if (firstlayer == true) {nev_cut ++ ;}
	  
	  //std::cout << nlay << " " << edepmax << std::endl;
	  hedeplay.Fill(nlay,edepmax);
	  hedepradius.Fill(radius,edepmax);
	  */
	  
      
	}
	
	
	//tcutDigit->AddFriend(tcutTrueMC);
	std::cout<<"Check 08 "<<std::endl;
    tcutDigit->AddFriend(tcutTracker);
	std::cout<<"Check 09 "<<std::endl;
	tcutDigit->AddFriend(tcutReco);
	std::cout<<"Check 10 "<<std::endl;
	tcutDigit->AddFriend(tcutEvent);
	std::cout<<"Check 10 "<<std::endl;
	TTree* tcut = tcutDigit;
	TCut internalEvent = "firstlayer==0";
	TG4Event* evcut = new TG4Event;
	std::vector<cell>* vec_cellcut = new std::vector<cell>;
	std::vector<digit>* vec_sttcut = new std::vector<digit>;
	double part_momcut[kMaxStdHepN][4];
	double evt_vtxcut[4];
	int part_pdgcut[kMaxStdHepN];
	event* ecut = new event;
	//event* evtcut;
	
	tcut->SetBranchAddress("Event",&evcut);
	tcut->SetBranchAddress("cell",&vec_cellcut);
	tcut->SetBranchAddress("Stt",&vec_sttcut);
	//tcut->SetBranchAddress("StdHepP4",part_momcut);
	//tcut->SetBranchAddress("EvtVtx",evt_vtxcut);
	//tcut->SetBranchAddress("StdHepPdg",part_pdgcut);
	tcut->SetBranchAddress("event",&ecut);
    
	
	TCanvas *c1 = new TCanvas("c1","c1",200,10,900,600);

	tcut->Draw("StdHepP4[0][3]","","", 1000, 0);
	TH1D *h = (TH1D*)gPad->GetPrimitive("htemp");
	h->SetTitle("Neutrino energy distribution");
	h->GetXaxis()->SetTitle("E(GeV/c)");
	h->GetYaxis()->SetTitle("number of events");
	c1->Modified();
	c1->Update();
	
	TCanvas *c2 = new TCanvas("c2","c2",200,10,900,600);
	
	tcut->SetMarkerColor(kBlack);
	tcut->Draw("EvtVtx[1]:EvtVtx[0]>>hyx(1,-2.5,2.5,1,0,-5)",internalEvent,"*", 1000, 0); //internalEvent
	TH1D *h1 = (TH1D*)gPad->GetPrimitive("hyx");
	h1->SetTitle("Vertex xy");
	h1->GetXaxis()->SetTitle("x(m)");
	h1->GetYaxis()->SetTitle("y(m)");
	tcut->SetMarkerColor(kRed);
	tcut->Draw("EvtVtx[1]:EvtVtx[0]",!internalEvent,"*same", 1000, 0); ////!internalEvent
	
	c2->Modified();
	c2->Update();
	
	TCanvas *c3 = new TCanvas("c3","c3",200,10,900,600);
	
	tcut->SetMarkerColor(kBlack);
	tcut->Draw("EvtVtx[1]:EvtVtx[2]>>hyz(1,20,25,1,0,-5)",internalEvent,"*", 1000, 0); //internalEvent
	TH1D *h2 = (TH1D*)gPad->GetPrimitive("hyz");
	h2->SetTitle("Vertex yz");
	h2->GetXaxis()->SetTitle("y(m)");
	h2->GetYaxis()->SetTitle("z(m)");
	tcut->SetMarkerColor(kRed);
	tcut->Draw("EvtVtx[1]:EvtVtx[2]",!internalEvent,"*same", 1000, 0); //!internalEvent
	c3->Modified();
	c3->Update();
	
	
	
	
	
	const int nevnew = tcut->GetEntries();
	TH1D* hmuonpe= new TH1D("hmuonpe","",100,0,16000);
	TH1D* hmuonEreco =  new TH1D("hmuonEreco","",100,0,16000);
	int ncut, cc, nc, muons, muonreco = 0;
	for (int i = 0; i < nevnew; i++) 
	{
		
	  //std::cout<<firstlayer<<std::endl;	
      tcut->GetEntry(i);
	  int pemuon = 0;
	  
	  for (int j = 0; j < vec_cellcut->size(); j++)
	    {
			for(int k = 0; k < vec_cellcut->at(j).hindex1.size(); k++)
			{
				//int pdg = evcut->SegmentDetectors["EMCalSci"].at(vec_cellcut->at(j).hindex1.at(0)).GetStart().X();
				//std::cout<<"number of tracks for this hit"<<evcut->SegmentDetectors["EMCalSci"].at(vec_cellcut->at(j).hindex1.at(k)).Contrib.size()<<std::endl;
				int Trackid = evcut->SegmentDetectors["EMCalSci"].at(vec_cellcut->at(j).hindex1.at(k)).GetPrimaryId();
				int pdg = evcut->Trajectories.at(Trackid).GetPDGCode();
				if (pdg == 13) {pemuon ++;}
			}
			
			for(int k = 0; k < vec_cellcut->at(j).hindex2.size(); k++)
			{
				//int pdg = evcut->SegmentDetectors["EMCalSci"].at(vec_cellcut->at(j).hindex1.at(0)).GetStart().X();
				//std::cout<<"number of tracks for this hit"<<evcut->SegmentDetectors["EMCalSci"].at(vec_cellcut->at(j).hindex1.at(k)).Contrib.size()<<std::endl;
				int Trackid = evcut->SegmentDetectors["EMCalSci"].at(vec_cellcut->at(j).hindex2.at(k)).GetPrimaryId();
				int pdg = evcut->Trajectories.at(Trackid).GetPDGCode();
				if (pdg == 13) {pemuon ++;}
			}
			
		}
	  //std::cout<<"size "<<vec_partcut->size()<<std::endl;
	  if (firstlayer == 0)
	  {
		  ncut++;
		  for (int k = 0; k < ecut->particles.size(); k++)	
		  {
			  if (ecut->particles.at(k).tr.pid == 13 ) 
			  {
				  muons++;
				  cc++;
				  if(ecut->particles.at(k).tr.ret_ln == 0 && ecut->particles.at(k).tr.ret_cr == 0)
					{
						muonreco++;
						hmuonEreco->Fill(ecut->particles.at(k).Ereco);
					}
				  else {std::cout<<"Problematic event number:"<<i<<std::endl;}
			  }
		  }		  
      }
	  if (firstlayer == 0 && pemuon != 0){hmuonpe->Fill(pemuon);}
	  
	  
	}
	
	std::cout<<"Number of events post cut:"<<ncut<<std::endl;
	//std::cout<<"Number of cc interactions:"<<cc<<std::endl;
	//std::cout<<"Number of nc interactions:"<<ncut-cc<<std::endl;
	std::cout<<"Number of muons:"<<muons<<std::endl;
	std::cout<<"Number of correctly reconstructed muons:"<<muonreco<<std::endl;
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
	
	c1->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/E.png");
	c2->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/xy.png");
	c3->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/yz.png");
	c4->SaveAs("/mnt/c/Linux/Dune/kloe-simu/plots/Muonpe.png");
	
	fout->cd();
	tcutDigit->Write();
	tcutTrueMC->Write();
	tcutTracker->Write();
	
	//hedeplay.Write();
	//hedepradius.Write();
	hmuonpe->Write();
	hmuonEreco->Write();
	
	tTracker->SetDirectory(0);
	tTrueMC->SetDirectory(0);
	tDigit->SetDirectory(0);
	
	/*
	std::cout << "number of events cut: " << nev_cut << std::endl;
	std::cout << "number of events in the outer layer: " << nev_outerlayer << std::endl;
	std::cout << "number of events in the outer layer that have been cut: " << correctcut << std::endl;
	std::cout << "number of events in the outer layer that have not been cut: " << nev_outerlayer-correctcut << std::endl;
	std::cout << "number of events not in the outer layer that have been cut: " << wrongcut	<< std::endl;
	*/
	
	delete f;
	//delete fout;
	
    
    
}



