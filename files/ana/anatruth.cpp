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

TString path = "/mnt/c/Linux/Dune/kloe-simu/files/ana/";
TFile fout(path + "result.root","RECREATE");


TCut fiducial_truth = "dr < 176 && abs(xv_true) < 1500";
TCut CC = "isCC == 1";
TCut lowE = "Enu_true/1E3<20";

TCanvas c1("c1","c1",1000,1000,1000,1000);

void EvalChi2(TH1D& h1, TH1D& h2, double& chi2, int& ndof)
{
  double n1, n2;
  double vchi2 = 0.;
  
  chi2 = 0.;
  
  ndof =  h1.GetNbinsX();

  for(int i = 0; i < ndof; i++)
  {
    n1 = h1.GetBinContent(i+1);
    n2 = h2.GetBinContent(i+1);
    
    if(n1 != 0 || n2 != 0)
      vchi2 = (n1 - n2) * (n1 - n2) / (n1 + n2);
    else
      vchi2 = 0.;
    
    chi2 += vchi2;
  }
}



void anatruth()
{
  gROOT->SetBatch(1);
  TChain tNom("tout","tNom");
  TChain tH1Y("tout","tH1Y");
  tNom.Add(path + "nominal/*");
  tH1Y.Add(path + "shift/*");
  std::cout << "H1Y: " << tH1Y.GetEntries() << "; " <<std::endl;
  
  //////////////////////////////////////////////SetALIASES
  
  tNom.SetAlias("pmu_true","sqrt(pxmu_true*pxmu_true+pymu_true*pymu_true+pzmu_true*pzmu_true)");
  tH1Y.SetAlias("pmu_true","sqrt(pxmu_true*pxmu_true+pymu_true*pymu_true+pzmu_true*pzmu_true)");
  
  
  tNom.SetAlias("ptmu_true","sqrt(pymu_true*pymu_true+pzmu_true*pzmu_true)");
  tH1Y.SetAlias("ptmu_true","sqrt(pymu_true*pymu_true+pzmu_true*pzmu_true)");

  tNom.SetAlias("dr","sqrt( (yv_true+2384.73)*(yv_true+2384.73) + (zv_true-23910)*(zv_true-23910) ) - 2000");
  tH1Y.SetAlias("dr","sqrt( (yv_true+2384.73)*(yv_true+2384.73) + (zv_true-23910)*(zv_true-23910) ) - 2000");


  ////////////////////////////////////////////////////////////////////


  int long NlowEok = tNom.GetEntries(lowE);
  
  

  /////////////////////////////////////////////////////////////////MOMENTUM HISTOGRAMS FROM TWO FLUXES
  
  
  TH1D hNom_true_all("hNom_true_all","",40,0,16);
  TH1D hH1Y_true_all("hH1Y_true_all","",40,0,16);
  

  
  ////////////Two fluxes all muons with fiducial cut (ptrue)
  
  tNom.Draw("pmu_true/1E3>>hNom_true_all", fiducial_truth + CC + lowE,"E0");
  tH1Y.Draw("pmu_true/1E3>>hH1Y_true_all", fiducial_truth + CC + lowE,"E0",NlowEok,0);
  
  hNom_true_all.GetXaxis()->SetTitle("p_{#mu}[GeV/c]");
  hNom_true_all.GetYaxis()->SetTitle("events");
  //hNom_true_fid.SetFillColor(kBlue-10);
  hNom_true_all.SetLineColor(kBlue);
  hNom_true_all.SetFillStyle(3001);
  //hH1Y_true_fid.SetFillColor(kRed-10);
  hH1Y_true_all.SetLineColor(kRed);
  hH1Y_true_all.SetFillStyle(3001);
  
  c1.SetLogy(true);
  hNom_true_all.Draw();
  hH1Y_true_all.Draw("same");
  TLegend* legend5all = new TLegend(0.6,0.75,0.9,0.9);
  legend5all->AddEntry("hNom_true_all", "p_{#mu}^{true} nominal", "f");
  legend5all->AddEntry((TObject*)0, TString::Format("entries: %.0f",hNom_true_all.GetEntries()).Data(), ""); 
  legend5all->AddEntry("hH1Y_true_all", "p_{#mu}^{true} shifted", "f");
  legend5all->AddEntry((TObject*)0, TString::Format("entries: %.0f",hH1Y_true_all.GetEntries()).Data(), ""); 
  legend5all->SetFillColor(0);
  legend5all->SetBorderSize(1);
  legend5all->Draw("");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "ptrueFluxeseverything.png");
  c1.SetLogy(false);

  
    
  /////////////////////////////////////////////////////CHI2TEST//////////////////////////////////////////////////////////////////
  
  
  int ndof;
  double chi2_true_all = 0.;

  
  
  
  EvalChi2(hNom_true_all, hH1Y_true_all, chi2_true_all, ndof);
  
	
  
  
  std::cout << "chi2: " << chi2_true_all << std::endl;
  std::cout << "p-value: " << TMath::Prob(chi2_true_all,hNom_true_all.GetNbinsX()) << std::endl;
  std::cout << "Confidence level: " <<TMath::ErfInverse(1-TMath::Prob(chi2_true_all,hNom_true_all.GetNbinsX()))*sqrt(2)<< std::endl;

  
  
 
  
}
