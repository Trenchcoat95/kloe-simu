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



#include <vector>
#include <map>
#include <iostream>

void test(const char* fIn, const char* fIn2)
{
   TFile* f= new TFile(fIn,"");
   TFile* f2= new TFile(fIn2,"");
   
   TH1D* hEreco =  new TH1D("hEreco","",100,0,16000);
   hEreco = (TH1D*)f->Get("hmuonErecoInner");
   TH1D* hEreco2 =  new TH1D("hEreco2","",100,0,16000);
   hEreco2 = (TH1D*)f2->Get("hmuonErecoInner");
   
   hEreco->Chi2Test(hEreco2,"P");
   hEreco->KolmogorovTest(hEreco2,"D");
   hEreco->AndersonDarlingTest(hEreco2,"D");
}