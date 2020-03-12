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

void control(const char* fIn)
{
	gSystem->Load("/mnt/c/Linux/Dune/kloe-simu/lib/libStruct.so");
	
	TChain* tEvent = new TChain("tEvent");
	tEvent->Add(fIn);
	event* e = new event;
	tEvent->SetBranchAddress("event",&e);
	
	vector<double> ptrue_muonreco;
	vector<double> preco_muonreco;
	
	for(int i=0; i< tEvent->GetEntries(); i++)
	{
		tEvent->GetEntry(i);
		
		for (int k = 0; k < e->particles.size(); k++)	
		  {
			  //if (e->particles.at(k).tr.pid == 13 && e->particles.at(k).pxreco == 0 && e->particles.at(k).pyreco == 0 && e->particles.at(k).pzreco == 0  &&e->particles.at(k).tr.ret_ln == 0 && e->particles.at(k).tr.ret_cr == 0)
			  //{
				//  std::cout<<"N of digits = "<<e->particles.at(k).tr.digits.size()<<std::endl;
			  //}
			  if (e->particles.at(k).tr.pid == 13 && e->particles.at(k).tr.digits.size() == 0 )
			  {
				  
				  //std::cout<<"(px,py,pz) = "<<e->particles.at(k).pxreco<< " " << e->particles.at(k).pyreco << " " << e->particles.at(k).pzreco << std::endl;
				  //std::cout<<"controls = "<<e->particles.at(k).tr.ret_ln  << " " << e->particles.at(k).tr.ret_cr  << std::endl;
				  
			  }
			  if(e->particles.at(k).tr.pid == 13&&e->particles.at(k).tr.ret_ln == 0 && e->particles.at(k).tr.ret_cr == 0 && e->particles.at(k).tr.chi2_ln  <= 30)
			  {
				  double ptot_reco=std::sqrt(e->particles.at(k).pxreco*e->particles.at(k).pxreco + e->particles.at(k).pyreco*e->particles.at(k).pyreco+ e->particles.at(k).pzreco*e->particles.at(k).pzreco); 
				  double ptot=std::sqrt(e->particles.at(k).pxtrue*e->particles.at(k).pxtrue + e->particles.at(k).pytrue*e->particles.at(k).pytrue+ e->particles.at(k).pztrue*e->particles.at(k).pztrue); 
				  ptrue_muonreco.push_back(ptot);
				  preco_muonreco.push_back(ptot_reco);
				  std::cout<<"N of digits nice = "<<e->particles.at(k).tr.digits.size()<<std::endl;
			  }
			  else if(e->particles.at(k).tr.pid == 13&&e->particles.at(k).tr.ret_ln == 0 && e->particles.at(k).tr.ret_cr == 0 && e->particles.at(k).tr.chi2_ln  > 30)
			  {
				 std::cout<<"N of digits = "<<e->particles.at(k).tr.digits.size()<<std::endl;
			  }
		  }
		
	}
	
	Double_t ptrue[ptrue_muonreco.size()],preco[preco_muonreco.size()];
	for(int i=0;i<ptrue_muonreco.size();i++)
	{
		ptrue[i]=ptrue_muonreco.at(i)*0.001;
		preco[i]=preco_muonreco.at(i)*0.001;
	}
	TCanvas *c2a = new TCanvas("ptrueVSpreco","ptrueVSpreco",200,10,900,700); ////////////// true p VS reconstructed p of muons surviving the cut
	//gStyle->SetOptFit();
    TGraph* grreco = new TGraph(ptrue_muonreco.size(),ptrue,preco);
    grreco->SetTitle("");
    grreco->GetYaxis()->SetTitle("p^{reco}_{#mu} [GeV/c]");
    grreco->GetXaxis()->SetTitle("p^{true}_{#mu} [GeV/c]");
    grreco->SetMarkerColor(kBlue);
	//grreco->Fit("pol1");
    grreco->Draw("A*"); //P
	
	std::cout<<"controls = "<<ptrue_muonreco.size() << std::endl;
}