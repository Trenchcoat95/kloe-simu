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

#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <vector>
#include <map>
#include <iostream>

#include "struct.h"
#include "utils.h"

// Energy MeV
// Distance mm
// Time ns


void LayerTime(TH1D* hEcell,std::vector<cell>* vec_cell, double tlay [5] ,double en [5], double pe [5], double & tcell )
{
	//CellTime(vec_cell, vec_celltime);
	int ncell = vec_cell->size(); 
	//std::cout << "ncell=" << ncell << std::endl;
	int laycell [5] ={0,0,0,0,0};
	
	for(int i = 0; i < ncell; i++)
	{
		cell c;
		
		c = vec_cell->at(i);
		
		double t0 = 0.5*( c.tdc1 + c.tdc2 - ns_Digit::lCalBarrel * ns_Digit::vlfb);
		int planeID = c.id / 100;
		int cellID = c.id - planeID*100;
		//std::cout<<"cellid="<<cellID<<std::endl;
		std::vector<double> pe1= c.pe_time1;
		std::vector<double> pe2= c.pe_time2;
		//std::cout<<"planeID:"<<planeID<<std::endl;
		//std::cout<<"t0:"<<c.t0<<std::endl;
		hEcell->Fill(c.adc1+c.adc2);
		
		if (planeID<=4 and t0>0)
		{		
			tlay[planeID] += t0 * (c.adc1+c.adc2);
			en[planeID] += c.adc1 + c.adc2;
			if(cellID==6)
			{
				pe[planeID] += 1.*(pe1.size()+pe2.size());
				laycell[planeID] ++;
				if (planeID==0)
				{tcell=t0;}
			}
		}	
		
		//if (planeID==0) {std::cout<<"npe on layer 0: "<<pe[planeID]<<std::endl;}
		
		if (t0<0)
		{
			//std::cout<<"t0="<<t0<<std::endl;
		}
	}

	for(int j = 0; j < 5; j++)
	{
		tlay[j] = tlay[j] / en[j];
	}
	
	//std::cout<< "laycell [" <<laycell[0]<<","<<laycell[1]<<","<<laycell[2]<<","<<laycell[3]<<","<<laycell[4]<<"]"<<std::endl;
}	

double RelativeTime(std::vector<cell>* vec_cell,
                    TH1D* ht01,TH1D* ht02,TH1D* ht03,TH1D* ht04,TH1D* ht05,
					TH1D* hE1,TH1D* hE2,TH1D* hE3,TH1D* hE4,TH1D* hE5, TH1D* hEcell, TH1D* hEsum,
					TH1D* hpe1,TH1D* hpe2,TH1D* hpe3,TH1D* hpe4,TH1D* hpe5,
					TH1D* ht01rel,TH1D* ht02rel,TH1D* ht03rel,TH1D* ht04rel,TH1D* ht05rel,
					TH1D* ht01singlecell)
{
	double tlay [5] = {0,0,0,0,0};
	double tcell = 0;
	double en [5] = {0,0,0,0,0};
	double pe [5] = {0,0,0,0,0};
	double dz [5] = {0,0.044,0.088,0.132,0.18095};
	
	LayerTime(hEcell,vec_cell,tlay,en,pe,tcell);
	//std::cout<<"tcell="<<tcell<<std::endl;
	ht01->Fill(tlay[0]);
	ht02->Fill(tlay[1]);
	ht03->Fill(tlay[2]);
	ht04->Fill(tlay[3]);
	ht05->Fill(tlay[4]);
	hE1->Fill(en[0]);
	hE2->Fill(en[1]);
	hE3->Fill(en[2]);
	hE4->Fill(en[3]);
	hE5->Fill(en[4]);
	hpe1->Fill(pe[0]);
	hpe2->Fill(pe[1]);
	hpe3->Fill(pe[2]);
	hpe4->Fill(pe[3]);
	hpe5->Fill(pe[4]);
	ht01singlecell->Fill(tcell);
	hEsum->Fill(en[0]+en[1]+en[2]+en[3]+en[4]);
	
	for(int j = 0; j < 5; j++)
	{
		tlay[j] = tlay[j]-dz[j]/0.2998; //t0-(dz_layer[m]/c[m/ns]) [ns]
	}
	
	ht01rel->Fill(tlay[0]);
	ht02rel->Fill(tlay[1]);
	ht03rel->Fill(tlay[2]);
	ht04rel->Fill(tlay[3]);
	ht05rel->Fill(tlay[4]);
	
	double t = 0;
	double e = 0;
	for (int j=0; j < 5; j++) 
	{
		t +=  tlay[j] * en[j];
		e +=  en[j];
	}
	
	//std::cout << "2- " << hpe1 << std::endl;
	
	return t/e;
	
	
}

void Resolution(const char* finname, const char* foutname)
{
	  
    auto f = TFile::Open(finname);
    auto t = (TTree*) f->Get("tDigit");
      
	std::vector<cell>* vec_cell = 0;
	
	t->SetBranchAddress("cell", &vec_cell);
	
	//std::vector<celltime> vec_celltime;
	
	auto fout = TFile::Open(foutname, "RECREATE");
	
    //TTree *tout= new TTree("tmu","Resolution");
	auto htmu = new TH1D("htmu", "Time resolution distribution", 200, 0, 4);
	TH1D *ht01= new TH1D("ht01","T0 first layer",200,0,4);
	TH1D *ht02= new TH1D("ht02","T0 second layer",200,0,4);
	TH1D *ht03= new TH1D("ht03","T0 third layer",200,0,4);
	TH1D *ht04= new TH1D("ht04","T0 fourth layer",200,0,4);
	TH1D *ht05= new TH1D("ht05","T0 fifth layer",200,0,4);
	
	TH1D *ht01singlecell= new TH1D("ht01singlecell","T0 first layer seventh cell",200,0,4);
	
	TH1D *hE1= new TH1D("hE1","Energy first layer",200,0,1000);
	TH1D *hE2= new TH1D("hE2","Energy second layer",200,0,1000);
	TH1D *hE3= new TH1D("hE3","Energy third layer",200,0,1000);
	TH1D *hE4= new TH1D("hE4","Energy fourth layer",200,0,1000);
	TH1D *hE5= new TH1D("hE5","Energy fifth layer",200,0,1000);
	
	TH1D *hEcell= new TH1D("hEcell","Energy per cell",200,0,1000);
	
	TH1D *hpe1= new TH1D("hpe1","Photo-electrons collected in the first layer",300,0,300);
	TH1D *hpe2= new TH1D("hpe2","Photo-electrons collected in the second layer",300,0,300);
	TH1D *hpe3= new TH1D("hpe3","Photo-electrons collected in the third layer",300,0,300);
	TH1D *hpe4= new TH1D("hpe4","Photo-electrons collected in the fourth layer",300,0,300);
	TH1D *hpe5= new TH1D("hpe5","Photo-electrons collected in the fifth layer",300,0,300);
	//std::cout << "1- " << hpe1 << std::endl;
	
	TH1D *ht01rel= new TH1D("ht01rel","T0 first layer relative",200,0,4);
	TH1D *ht02rel= new TH1D("ht02rel","T0 second layer relative",200,0,4);
	TH1D *ht03rel= new TH1D("ht03rel","T0 third layer relative",200,0,4);
	TH1D *ht04rel= new TH1D("ht04rel","T0 fourth layer relative",200,0,4);
	TH1D *ht05rel= new TH1D("ht05rel","T0 fifth layer relative",200,0,4);
	
	TH1D *hEsum= new TH1D("hEsum","Total Energy Deposit",200,0,2000);
    //tout->Branch("time","TH1D",&htmu,32000,0);
    double tmu = 0;
    
    const int nev = t->GetEntries();
    
    std::cout << "Events: " << nev << " [";
    std::cout << std::setw(3) << int(0) << "%]" << std::flush ;
    
	
	
    for(int i = 0; i < nev; i++)
    {
		
		t->GetEntry(i);
       
		std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush ;
      
	//std::cout << "*- " << hpe1 << std::endl;
		tmu=RelativeTime(vec_cell,ht01,ht02,ht03,ht04,ht05,hE1,hE2,hE3,hE4,hE5,hEcell,hEsum,hpe1,hpe2,hpe3,hpe4,hpe5,ht01rel,ht02rel,ht03rel,ht04rel,ht05rel,ht01singlecell);
      
	    //std::cout<<"tmu="<<tmu<<std::endl;
        
		htmu->Fill(tmu);
		
    }
	//tout->Fill();
    std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
    std::cout << std::endl;
	
	
	
	
	fout->cd();
	//tout->Print();
    fout->Write();
	
	//hpe1->Write();
	//hpe2->Write();
	
	
	
	TApplication theApp("App", 0, 0);
	gStyle->SetOptFit();
	TCanvas *c1 = new TCanvas("c1","htmu",200,10,900,600);
	htmu->GetXaxis()->SetTitle("t_{#mu}[ns]");
    htmu->GetYaxis()->SetTitle("n^{#circ} events");
	htmu->SetFillColor(kBlue-10);
	htmu->Fit("gaus");
	htmu->Draw();
	c1->Modified();
	c1->Update();
	gSystem->ProcessEvents();
	
	TCanvas *c2 = new TCanvas("c2","hEcell",200,10,900,600);
	hEcell->GetXaxis()->SetTitle("adc1+adc2");
    hEcell->GetYaxis()->SetTitle("n^{#circ} of cells");
	hEcell->SetFillColor(kBlue-10);
	hEcell->Fit("landau");
	hEcell->Draw();
	c2->Modified();
	c2->Update();
	gSystem->ProcessEvents();
	
	
	c1->SaveAs("/mnt/c/Linux/Dune/kloe-simu/files/htmu.png");
	c2->SaveAs("/mnt/c/Linux/Dune/kloe-simu/files/hEcell.png");
	
    fout->Close();
    f->Close();
	
	theApp.Run(true); 

}

void help_time()
{
  std::cout << "Time <input file> <output file>" << std::endl;
  std::cout << "input file name could contain wild card" << std::endl;
} 

int main(int argc, char* argv[])
{
  gSystem->Load("libStruct.so");
  if(argc != 3)
    help_time();
  else
    Resolution(argv[1], argv[2]);
}
