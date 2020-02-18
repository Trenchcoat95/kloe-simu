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

void ana(const char* fgeo, const char* fIn, const char* fOut)    ///////(First file, List of file in directory using wildcard *, output file)
{
	int change = 1;
	
	//////////////////LOADING THE LIBRARIES AND TREES///////////////////////////////
	
	gSystem->Load("/mnt/c/Linux/Dune/kloe-simu/lib/libStruct.so");
	
	//TChain *f=new TChain("f");
	//f->Add(fIn);
	TFile* fg= new TFile(fgeo);
    TFile* fo= new TFile(fOut,"RECREATE");
    fo->cd();	
	//TTree* tTracker = (TTree*) f->Get("gRooTracker");
	
	TChain* tTracker = new TChain("gRooTracker");
	tTracker->Add(fIn);
	TChain* tTrueMC = new TChain("EDepSimEvents");
	tTrueMC->Add(fIn);
	TChain* tDigit = new TChain("tDigit");
	tDigit->Add(fIn);
	TChain* tReco = new TChain("tReco");
	tReco->Add(fIn);
	TChain* tEvent = new TChain("tEvent");
	tEvent->Add(fIn);
	TGeoManager* geo = (TGeoManager*) fg->Get("EDepSimGeometry");
	
	TTree* tout = new TTree("tLay","tLay");
	int outerlayer;
	int outerlayerMC;
	int MuonReco;
	double Ereco;
    tout->Branch("outerlayer",&outerlayer,"outerlayer/I");
	tout->Branch("outerlayerMC",&outerlayerMC,"outerlayerMC/I");
	tout->Branch("MuonReco",&MuonReco,"MuonReco/I");
	tout->Branch("Ereco",&Ereco,"Ereco/D");
	
	
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
	
    TCut internalEvent = "outerlayer==0";
	TCut internalEventMC = "outerlayerMC==0";
	
	
	TH1D* hmuonpeInner= new TH1D("hmuonpeInner","",100,0,16000);
	TH1D* hmuonErecoInner =  new TH1D("hmuonErecoInner","",100,0,16000);
	
	TH1D hEneutrinoInner ("hEneutrinoInner","",100,0,60);
	TH1D hEneutrinoInnerMC ("hEneutrinoInnerMC","",100,0,60);
	
	TH1D hpmuon ("hpmuon","",100,0,16000);
	TH1D hpmuonreco ("hpmuonreco","",100,0,16000);
	TH1D hmuonrecoEInner ("hmuonrecoEInner","",100,0,16000);
	
	int ncutMC = 0;
	int ncut = 0;
	int muons = 0;
	int muonreco = 0;
	int wrongcut = 0;
	int missedcut = 0;
	const int nev = t->GetEntries();
	
	double centerKLOE[3];
	double dummyLoc[3];
    geo->cd("volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/volSTTFULL_PV_0");
    dummyLoc[0] = 0.;
	dummyLoc[1] = 0.;
	dummyLoc[2] = 0.;
	geo->LocalToMaster(dummyLoc, centerKLOE);
	int nlayMC;
	
	///////////////////////////////////////Ciclo sugli eventi ///////////////////////////////////	
	
	for (int i = 0; i < nev; i++) 
	{ 
      
	  t->GetEntry(i);
	  
	  MuonReco = 0;
	  outerlayer = 0;
	  outerlayerMC = 0;
	  int pemuon = 0;
	  
	  double dy = evt_vtx[1]*1000 - centerKLOE[1];
	  double dz = evt_vtx[2]*1000 - centerKLOE[2];
	  double radius = sqrt(dy*dy+dz*dz);
	  double dr = radius - 2000.; 
	  
	  
	  nlayMC = -1;
	  
	  if(dr < 44.)
	  {
		  nlayMC = 0;
	  }
	  else if(dr >= 44. && dr < 88.)
	  {
		  nlayMC = 1;
	  }
	  else if(dr >= 88. && dr < 132.)
	  {
		  nlayMC = 2;
	  }
	  else if(dr >= 132. && dr < 176.)
	  {
		  nlayMC = 3;
	  }
	  else if(dr >= 176. && dr < 245.)
	  {
		  nlayMC = 4;
	  }
	  else
	  {
		  std::cout << "dr = " << dr << std::endl;
	  }
	  
	  if (nlayMC == 4) {outerlayerMC=1;}
	  
		
		
	
	  
	  
	  
	  
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
			  outerlayer = 1;
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
	  if (outerlayer == 0)
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
						hmuonErecoInner->Fill(e->particles.at(k).Ereco);   /////////////////////////////////Fill dell'istogramma dell'energia ricostruita
						Ereco=e->particles.at(k).Ereco;
						MuonReco=1;
					}
				  else {std::cout<<"Problematic event number:"<<i<<std::endl;}
			  }
		  }		  
      }
	  
	  
	  for(int j = 0; j < e->particles.size(); j++)
	  {
		  int pdg = e->particles.at(j).tr.pid;
		  int tid = e->particles.at(j).tr.tid;
		  
		  if (pdg == 13)
		  {
			  TLorentzVector p = ev->Trajectories.at(tid).GetInitialMomentum();
			  double pmod = std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
			  //std::cout<<"pmod"<<pmod<<std::endl;               
			  hpmuon.Fill(pmod);                          ////////////////////////////////////////////Fill istogramma momento muoni
			  
			  if(e->particles.at(j).tr.ret_ln == 0 && e->particles.at(j).tr.ret_cr == 0)
			  {
				  //std::cout<<"pmod"<<pmod<<std::endl;
				  hpmuonreco.Fill(pmod);                  ////////////////////////////////////////////Fill istogramma momento muoni ricostruiti
				  if (outerlayer == 0){hmuonrecoEInner.Fill(p[3]);}
			  }
			  
		  }
	  }
	  
	  if (outerlayer == 0 && pemuon != 0){hmuonpeInner->Fill(pemuon);}
	  
	  if (outerlayer==0){hEneutrinoInner.Fill(part_mom[0][3]);}
	  if (outerlayerMC==0){hEneutrinoInnerMC.Fill(part_mom[0][3]);}
	  if (outerlayer==1 && outerlayerMC==0){wrongcut++;}
	  if (outerlayer==0 && outerlayerMC==1){missedcut++;}
	  
	  tout->Fill();
      
	}
	
	///////////////STAMPA INFO SUL CUT /////////////////////////
	
	
	std::cout<<"Number of events post cut:"<<ncut<<std::endl;
	std::cout<<"Number of events cut not in the outer layer:"<<wrongcut<<std::endl;
	std::cout<<"Number of events not cut in the outer layer:"<<missedcut<<std::endl;
	std::cout<<"Number of events in the inner layers:"<<hEneutrinoInnerMC.GetEntries()<<std::endl;
	std::cout<<"Number of muons:"<<muons<<std::endl;
	std::cout<<"Number of correctly reconstructed muons:"<<muonreco<<std::endl;
	
	
	///////////////DRAW/////////////////////////
	
	TCanvas *c1 = new TCanvas("cneutrino_E","cneutrino_E",200,10,900,600);

	t->Draw("StdHepP4[0][3]","","", 1000, 0);
	TH1D *h = (TH1D*)gPad->GetPrimitive("htemp");
	h->SetTitle("Neutrino energy distribution");
	h->GetXaxis()->SetTitle("E(GeV/c)");
	h->GetYaxis()->SetTitle("number of events");
	c1->Modified();
	c1->Update();
	
	TCanvas *c2 = new TCanvas("cyx","cyx",200,10,900,900);
	
	t->SetMarkerColor(kBlack);
	t->Draw("EvtVtx[1]:EvtVtx[0]>>hyx",internalEvent,"*", 1000, 0); //internalEvent
	TH1D *h1 = (TH1D*)gPad->GetPrimitive("hyx");
	h1->SetTitle("");
	//h1->GetXaxis()->SetRangeUser(-2.5, 2.5);
	//h1->GetYaxis()->SetRangeUser(-5, 0);
	h1->GetXaxis()->SetTitle("x(m)");
	h1->GetYaxis()->SetTitle("y(m)");
	t->SetMarkerColor(kRed);
	t->Draw("EvtVtx[1]:EvtVtx[0]",!internalEvent,"*same", 1000, 0); ////!internalEvent
	
	c2->Modified();
	c2->Update();
	
	TCanvas *c3 = new TCanvas("cyz","cyz",200,10,900,900);
	
	t->SetMarkerColor(kBlack);
	t->Draw("EvtVtx[1]:EvtVtx[2]>>hyz",internalEvent,"*", 1000, 0); //internalEvent
	TH1D *h2 = (TH1D*)gPad->GetPrimitive("hyz");
	h2->SetTitle("");
	h2->GetXaxis()->SetTitle("z(m)");
	h2->GetYaxis()->SetTitle("y(m)");
	t->SetMarkerColor(kRed);
	t->Draw("EvtVtx[1]:EvtVtx[2]",!internalEvent,"*same", 1000, 0); //!internalEvent
	c3->Modified();
	c3->Update();
	

	gStyle->SetOptStat(kTRUE);
	gStyle->SetOptFit(kTRUE);
	
	TCanvas *c4 = new TCanvas("cmuon_pe","cmuon_pe",200,10,900,600);
	hmuonpeInner->SetTitle("Muon p.e. production");
	hmuonpeInner->GetXaxis()->SetTitle("number of photoelectrons");
	hmuonpeInner->GetYaxis()->SetTitle("events");
	hmuonpeInner->SetFillColor(kBlue-10);
	hmuonpeInner->Fit("landau");
	hmuonpeInner->Draw();
	c4->Modified();
	c4->Update();

	
	TCanvas *c5 = new TCanvas("cmuon_reconstructed_E","cmuon_reconstructed_E",200,10,900,600);
	hmuonErecoInner->SetTitle("Muon reconstructed energy");
	hmuonErecoInner->GetXaxis()->SetTitle("MeV");
	hmuonErecoInner->GetYaxis()->SetTitle("events");
	hmuonErecoInner->SetFillColor(kBlue-10);
	//hmuonErecoInner->Fit("landau");
	hmuonErecoInner->Draw();
	c5->Modified();
	c5->Update();
	
	
	
	TCanvas *c6 = new TCanvas("cneutrino_E_cut","cneutrino_E_cut",200,10,900,600);
	gStyle->SetOptStat(0);
	hEneutrinoInner.SetTitle("");//Neutrino energy distribution of the events surviving the cut
	hEneutrinoInner.GetXaxis()->SetTitle("GeV/c");
	hEneutrinoInner.GetYaxis()->SetTitle("events");
	hEneutrinoInner.SetFillColor(kGray);
	hEneutrinoInner.SetLineColor(kBlack);
	hEneutrinoInner.SetFillStyle(3001);
	hEneutrinoInnerMC.SetFillColor(kYellow-10);
	hEneutrinoInnerMC.SetLineColor(kYellow);
		//hEneutrinoInner.Fit("landau");
	c6->Modified();
	c6->Update();
	hEneutrinoInnerMC.DrawClone();
	hEneutrinoInner.DrawClone("same");
	
	TLegend* legend6 = new TLegend(0.62,0.8,0.894,0.9);
	legend6->AddEntry("hEneutrinoInnerMC", "Neutrino having vertex in the inner layers", "f");
    legend6->AddEntry( "hEneutrinoInner", "Neutrinos surviving the cut", "f");
	legend6->SetFillColor(0);
	legend6->SetBorderSize(1);
	legend6->Draw("");
	
	//TCanvas *c7 = new TCanvas("c7","c7",200,10,900,600);
	
	//hEneutrinoInnerMC.SetTitle("");//Neutrino energy distribution of the events with vertex in the inner layers
	//hEneutrinoInnerMC.GetXaxis()->SetTitle("GeV/c");
	//hEneutrinoInnerMC.GetYaxis()->SetTitle("events");
	//hEneutrinoInnerMC.SetFillColor(kRed-10);
	//hEneutrinoInnerMC.SetLineColor(kRed);
	//hEneutrinoInnerMC.Fit("landau");
	//c7->Modified();
	//c7->Update();
	//hEneutrinoInnerMC.DrawClone("same");
	
	
	
	TEfficiency* pEffCut = 0;                       ///////////////////////////////Efficienza Cut
	if(TEfficiency::CheckConsistency(hEneutrinoInner,hEneutrinoInnerMC))
	{
		pEffCut = new TEfficiency(hEneutrinoInner,hEneutrinoInnerMC);
		TCanvas *c8 = new TCanvas("ccut_efficiency","ccut_efficiency",200,10,900,600);
		pEffCut->SetTitle(";Neutrino Energy [GeV/c];Efficiency"); //Cut efficiency as a function of the neutrino energy
		//pEffCut->GetXaxis()->SetTitle("Neutrino Energy (GeV/c)");
		//pEffCut->GetYaxis()->SetTitle("Cut Efficiency");
		//pEffCut->SetFillColor(kBlue-10);
		//pEffCut->Fit("landau");
		pEffCut->Draw();
		c8->Modified();
		c8->Update();
		c8->Update();
		
		std::string CutEfficiencyad = "/mnt/c/Linux/Dune/kloe-simu/plots/CutEfficiency";
		CutEfficiencyad += std::to_string(nev)+".png";
		//c8->SaveAs(CutEfficiencyad.c_str());
		c8->Write("",TObject::kOverwrite);
		pEffCut->Write("",TObject::kOverwrite);
	}
	
	TCanvas *c9 = new TCanvas("cmuon_p","cmuon_p",200,10,900,600);
	gStyle->SetOptStat(0);
	hpmuon.SetTitle("");
	hpmuon.GetXaxis()->SetTitle("E [MeV]");
	hpmuon.GetYaxis()->SetTitle("events");
	hpmuon.SetFillColor(kBlue-10);
	hpmuon.SetLineColor(kBlue);
	//hpmuon.Fit("landau");
	c9->Modified();
	c9->Update();
	hpmuon.DrawClone();
	
	
	//TCanvas *c10 = new TCanvas("c10","c10",200,10,900,600);
	hpmuonreco.SetTitle("");
	hpmuonreco.GetXaxis()->SetTitle("MeV");
	hpmuonreco.GetYaxis()->SetTitle("events");
	hpmuonreco.SetFillColor(kGreen-10);
	hpmuonreco.SetLineColor(kGreen);
	//hpmuonreco.Fit("landau");
	//c10->Modified();
	//c10->Update();
	hpmuonreco.DrawClone("same");
	
	TLegend* legend4 = new TLegend(0.62,0.8,0.894,0.9);
	legend4->AddEntry("hpmuon", "All muons", "f");
    legend4->AddEntry( "hpmuonreco", "Correctly reconstructed muons", "f");
	legend4->SetFillColor(0);
	legend4->SetBorderSize(1);
	legend4->Draw("");
	
	TCanvas *c10a = new TCanvas("cmuon_E","cmuon_E",200,10,900,600);
	hmuonrecoEInner.SetTitle("Muon Real Energy distribution (Correctly Recunstructed and survived the cut)");
	hmuonrecoEInner.GetXaxis()->SetTitle("MeV");
	hmuonrecoEInner.GetYaxis()->SetTitle("events");
	hmuonrecoEInner.SetFillColor(kBlue-10);
	//hpmuonreco.Fit("landau");
	c10a->Modified();
	c10a->Update();
	hmuonrecoEInner.DrawClone();
	
	
	
	//hmuonrecoEInner.Chi2Test(hmuonErecoInner,"P");
	
	
	TEfficiency* pEffReco = 0;                             ///////////////////////////////Efficienza ricostruzione
	if(TEfficiency::CheckConsistency(hpmuonreco,hpmuon))
	{
		pEffReco = new TEfficiency(hpmuonreco,hpmuon);
		TCanvas *c11 = new TCanvas("creconstruction_efficiency","creconstruction_efficiency",200,10,900,600);
		pEffReco->SetTitle(";E [MeV];#varepsilon_{reco}(E)");//("Reconstruction efficiency as a function of the muon momentum;Muon momentum [MeV];Efficiency");
		//pEffCut->GetXaxis()->SetTitle("Neutrino Energy (GeV/c)");
		//pEffCut->GetYaxis()->SetTitle("Cut Efficiency");
		//pEffCut->SetFillColor(kBlue-10);
		//pEffCut->Fit("landau");
		pEffReco->Draw();
		c11->Modified();
		c11->Update();
		c11->Update();
		
		std::string RecoEfficiencyad = "/mnt/c/Linux/Dune/kloe-simu/plots/RecoEfficiency";
		RecoEfficiencyad += std::to_string(nev)+".png";
		//c11->SaveAs(RecoEfficiencyad.c_str());
		c11->Write("",TObject::kOverwrite);
		pEffReco->Write("",TObject::kOverwrite);
	}
	
	
	/*
	TH1D* hist_new=(TH1D*)hmuonpe->Clone();
    hist_new->SetName("hist_new");
	
	hmuonpe -> Chi2Test(hist_new,"P");
	*/
	
	///////////////////////////////SALVA TUTTO////////////////////////////////////
	/*
	if (change==0)
	{
		std::string Eneutrinoad = "/mnt/c/Linux/Dune/kloe-simu/plots/E";
		std::string xyad = "/mnt/c/Linux/Dune/kloe-simu/plots/xy";
		std::string yzad = "/mnt/c/Linux/Dune/kloe-simu/plots/yz";
		std::string Muonpead = "/mnt/c/Linux/Dune/kloe-simu/plots/Muonpe";
		std::string Erecoad = "/mnt/c/Linux/Dune/kloe-simu/plots/Ereco";
		std::string EneutrinoInnerad = "/mnt/c/Linux/Dune/kloe-simu/plots/EneutrinoInner";
		std::string EneutrinoInnerMCad = "/mnt/c/Linux/Dune/kloe-simu/plots/EneutrinoInnerMC";
		std::string PMuonad = "/mnt/c/Linux/Dune/kloe-simu/plots/PMuon";
		std::string PMuonRecoad = "/mnt/c/Linux/Dune/kloe-simu/plots/PMuonReco";
		std::string PMuonRecoInnerad = "/mnt/c/Linux/Dune/kloe-simu/plots/PMuonRecoInner";
	}
	else if (change==1)
	*/
		std::string Eneutrinoad = "/mnt/c/Linux/Dune/kloe-simu/plots/CHE";
		std::string xyad = "/mnt/c/Linux/Dune/kloe-simu/plots/CHxy";
		std::string yzad = "/mnt/c/Linux/Dune/kloe-simu/plots/CHyz";
		std::string Muonpead = "/mnt/c/Linux/Dune/kloe-simu/plots/CHMuonpe";
		std::string Erecoad = "/mnt/c/Linux/Dune/kloe-simu/plots/CHEreco";
		std::string EneutrinoInnerad = "/mnt/c/Linux/Dune/kloe-simu/plots/CHEneutrinoInner";
		std::string EneutrinoInnerMCad = "/mnt/c/Linux/Dune/kloe-simu/plots/CHEneutrinoInnerMC";
		std::string PMuonad = "/mnt/c/Linux/Dune/kloe-simu/plots/CHPMuon";
		std::string PMuonRecoad = "/mnt/c/Linux/Dune/kloe-simu/plots/CHPMuonReco";
		std::string PMuonRecoInnerad = "/mnt/c/Linux/Dune/kloe-simu/plots/CHPMuonRecoInner";
	
	
	Eneutrinoad += std::to_string(nev)+".png";
	xyad += std::to_string(nev)+".png";
	yzad += std::to_string(nev)+".png";
	Muonpead += std::to_string(nev)+".png";
	Erecoad += std::to_string(nev)+".png";
	EneutrinoInnerad += std::to_string(nev)+".png";
	EneutrinoInnerMCad += std::to_string(nev)+".png";
	PMuonad += std::to_string(nev)+".png";
	PMuonRecoad += std::to_string(nev)+".png";
	PMuonRecoInnerad += std::to_string(nev)+".png";
	
	/*
	c1->SaveAs(Eneutrinoad.c_str());
	c2->SaveAs(xyad.c_str());
	c3->SaveAs(yzad.c_str());
	c4->SaveAs(Muonpead.c_str());
	c5->SaveAs(Erecoad.c_str());
	c6->SaveAs(EneutrinoInnerad.c_str());
	//c7->SaveAs(EneutrinoInnerMCad.c_str());
	c9->SaveAs(PMuonad.c_str());
	//c10->SaveAs(PMuonRecoad.c_str());
	c10a->SaveAs(PMuonRecoad.c_str());
	*/
	
	c1->Write("",TObject::kOverwrite);
	c2->Write("",TObject::kOverwrite);
	c3->Write("",TObject::kOverwrite);
	c4->Write("",TObject::kOverwrite);
	c5->Write("",TObject::kOverwrite);
	c6->Write("",TObject::kOverwrite);
	c9->Write("",TObject::kOverwrite);
	c10a->Write("",TObject::kOverwrite);
	
	tout->Write("",TObject::kOverwrite);
	
	hmuonpeInner->Write("",TObject::kOverwrite);
	hmuonErecoInner->Write("",TObject::kOverwrite);
	hEneutrinoInner.Write("",TObject::kOverwrite);
	hEneutrinoInnerMC.Write("",TObject::kOverwrite);
	hpmuon.Write("",TObject::kOverwrite);
	hpmuonreco.Write("",TObject::kOverwrite);
	hmuonrecoEInner.Write("",TObject::kOverwrite);
	
	
	//delete f;
	//delete fout;
	
    
    
}



