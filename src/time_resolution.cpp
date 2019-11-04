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


void LayerTime(std::vector<cell>* vec_cell, double tlay [5] ,double en [5] )
{
	//CellTime(vec_cell, vec_celltime);
	int ncell = vec_cell->size(); 
	//std::cout << "ncell=" << ncell << std::endl;
	
	
	for(int i = 0; i < ncell; i++)
	{
		cell c;
		
		c = vec_cell->at(i);
		
		double t0 = 0.5*( c.tdc1 + c.tdc2 - ns_Digit::lCalBarrel * ns_Digit::vlfb);
		int planeID = c.id / 100;
		//std::cout<<"planeID:"<<planeID<<std::endl;
		//std::cout<<"t0:"<<c.t0<<std::endl;
		
		if (planeID<=4 and t0>0)
		{		
			tlay[planeID] += t0 * (c.adc1+c.adc2);
			en[planeID] += c.adc1 + c.adc2;
		}		
		if (t0<0)
		{
			//std::cout<<"t0="<<t0<<std::endl;
		}
	}

	for(int j = 0; j < 5; j++)
	{
		tlay[j] = tlay[j] / en[j];
	}
	
	//std::cout<< "tlay [" <<tlay[0]<<","<<tlay[1]<<","<<tlay[2]<<","<<tlay[3]<<","<<tlay[4]<<"]"<<std::endl;
}	

double RelativeTime(std::vector<cell>* vec_cell)
{
	double tlay [5] = {0,0,0,0,0};
	double en [5] = {0,0,0,0,0};
	double dz [5] = {0,0.044,0.088,0.132,0.18095};
	LayerTime(vec_cell,tlay,en);
	
	for(int j = 0; j < 5; j++)
	{
		tlay[j] = tlay[j]-dz[j]/0.2998; //t0-(dz_layer[m]/c[m/ns]) [ns]
	}
	
	double t = 0;
	double e = 0;
	for (int j=0; j < 5; j++) 
	{
		t +=  tlay[j] * en[j];
		e +=  en[j];
	}
	
	
	
	return t/e;
	
	
}
void Resolution(const char* finname, const char* foutname)
{
	  
    TFile f(finname,"READ");
    TTree* t = (TTree*) f.Get("tDigit");
      
    std::vector<cell>* vec_cell= new std::vector<cell>;
	t->SetBranchAddress("cell",&vec_cell);
	
	//std::vector<celltime> vec_celltime;
	
	TFile fout(foutname,"RECREATE");
    //TTree *tout= new TTree("tmu","Resolution");
	TH1D *htmu= new TH1D("htmu","Time resolution distribution",200,0,4);
    //tout->Branch("time","TH1D",&htmu,32000,0);
    double tmu = 0;
    
    const int nev = t->GetEntries();
    
    std::cout << "Events: " << nev << " [";
    std::cout << std::setw(3) << int(0) << "%]" << std::flush ;
    
	
	
    for(int i = 0; i < nev; i++)
    {
		
		t->GetEntry(i);
       
		std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush ;
      
		tmu=RelativeTime(vec_cell);
      
	    //std::cout<<"tmu="<<tmu<<std::endl;
        
		htmu->Fill(tmu);
		
    }
	//tout->Fill();
    std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
    std::cout << std::endl;
	
	
	fout.cd();
	//tout->Print();
    fout.Write();
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
	
	
	c1->SaveAs("/mnt/c/Linux/Dune/kloe-simu/files/htmu.png");
	
    fout.Close();
    f.Close();
	
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
