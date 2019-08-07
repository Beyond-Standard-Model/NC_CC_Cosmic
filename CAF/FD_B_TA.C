#ifdef __CINT__
void  FD_B_TA() // For Binary Classification
{
  std::cout << "Sorry, you must run in compiled mode" << std::endl;
}
#else

#include "CAFAna/Core/SpectrumLoaderBase.h"
#include "CAFAna/Core/Progress.h"
#include "CAFAna/Core/SAMQuerySource.h"
#include "CAFAna/Core/SAMProjectSource.h"
#include "CAFAna/Core/WildcardSource.h"
#include "CAFAna/Core/IFileSource.h"
#include "CAFAna/Core/FileListSource.h"
#include "CAFAna/Core/Utilities.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Cuts/Cuts.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/OscillatableSpectrum.h"
#include "CAFAna/Vars/Vars.h"
#include "CAFAna/Core/OscCurve.h"
#include "CAFAna/Core/ReweightableSpectrum.h"
#include "CAFAna/Core/EventList.h"
#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Cuts/Cuts.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/nus/NuSLoadMacroDefs.h"
#include "TLorentzVector.h"
#include <cmath>
#include <iostream>
#include <functional>
#include <list>
#include <memory>
#include <set>
#include <string>
#include <vector>
//ROOT includes
#include "TStyle.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TVector3.h"
#include "StandardRecord/StandardRecord.h"
#include "StandardRecord/Proxy/SRProxy.h" 
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"

using namespace ana;

void FD_B_TA(){
  //Training and Testing 
  const std::string fnamecos = "prod_decaf_R17-03-01-prod3reco.h_fd_cosmic_full_nue_or_numu_or_nus_contain_v1_goodruns";
  //const std::string fnamecos = "defname:prod_decaf_R17-03-01-prod3reco.h_fd_cosmic_full_nue_or_numu_or_nus_contain_v1_goodruns with limit 10"; 
  SpectrumLoader loadercos(fnamecos);
  const std::string fnamenc = "prod_decaf_R17-03-01-prod3reco.l_fd_genie_nonswap_fhc_nova_v08_full_nue_or_numu_or_nus_contain_v1";
  //const std::string fnamenc = "defname:prod_decaf_R17-03-01-prod3reco.l_fd_genie_nonswap_fhc_nova_v08_full_nue_or_numu_or_nus_contain_v1 with limit 10";
  SpectrumLoader loadernc(fnamenc);
  //Output file
  TFile tMVALoad("TRA_FD_B.root","RECREATE");
  //Two trees for four interactions
  TTree *STree = new TTree("STree", "NC Signal tree");
  TTree *BTree = new TTree("BTree", "Backgrounds tree");
  ////////////////////////////////// Definding the useful variables to transform into Output root File ///////////////////////////////////
  //True Info
  int   Run;  
  int   SubRun;
  int   Evt;
  int   SubEvt;
  int   IsCC;
  ////////////////////////////////The below are the selected vars for training//////////////////////////////////////////////////////
  //Slice Info
  float nhit;      
  float ncalhit;   
  float nmiphit;   
  float ncontplanes;
  float ncellsfromedge;
  float slcCalE;
  float boxminy;
  float boxmaxy;
  float meanposX;   
  float meanposY;   
  float meanposZ;   
  float starttime;  
  float endtime;    
  float meantime;   
  float closestslicetime;
  float closestslicenhit;
  //NuS Sandbox 
  float sumtx;     
  float ewsumtx;   
  float angsumtx;  
  float angewsumtx;
  float cossumtx;  
  float cosewsumtx;
  float sumty;     
  float ewsumty;   
  float angsumty;  
  float angewsumty;
  float cossumty;  
  float cosewsumty;
  //Vertex Info
  float vtxTime;          
  float vtxX;             
  float vtxY;             
  float vtxZ;             
  float vtxpngNhit;       
  float vtxpngNhitX;      
  float vtxpngNhitY;      
  //FuzzyK Prongs
  float fuzzykntot;
  float fuzzyknpng;
  float fuzzyknpng2d;
  float fuzzykpngdirx;
  float fuzzykpngdiry;
  float fuzzykpngdirz;
  float fuzzykpngCalE;
  float fuzzykpnglen;
  float fuzzykpngnhit;
  float fuzzykpngnhitx;
  float fuzzykpngnhity;
  //CVN  
  float ncid; 
  float nueid;
  float numuid;
  float nutauid;
  float cosmicid;
  /////////////////////////////////////////////////////// Background Tree //////////////////////////////////////////////////////
  //True Info
  BTree->Branch("Run",&Run,"Run/I");          
  BTree->Branch("SubRun",&SubRun,"SubRun/I"); 
  BTree->Branch("Evt",&Evt,"Evt/I");          
  BTree->Branch("SubEvt",&SubEvt,"SubEvt/I"); 
  BTree->Branch("IsCC",&IsCC,"IsCC/I");       
  //Slice Info
  BTree->Branch("nhit",&nhit,"nhit/F");              
  BTree->Branch("ncalhit",&ncalhit,"ncalhit/F");              
  BTree->Branch("nmiphit",&nmiphit,"nmiphit/F");          
  BTree->Branch("ncontplanes",&ncontplanes,"ncontplanes/F");  
  BTree->Branch("ncellsfromedge",&ncellsfromedge,"ncellsfromedge/F");  
  BTree->Branch("slcCalE",&slcCalE,"slcCalE/F");        
  BTree->Branch("boxminy",&boxminy,"boxminy/F");              
  BTree->Branch("boxmaxy",&boxmaxy,"boxmaxy/F");              
  BTree->Branch("starttime",&starttime,"starttime/F");           
  BTree->Branch("meantime",&meantime,"meantime/F");           
  BTree->Branch("endtime",&endtime,"endtime/F");           
  BTree->Branch("meanposX",&meanposX,"meanposX/F");           
  BTree->Branch("meanposY",&meanposY,"meanposY/F");           
  BTree->Branch("meanposZ",&meanposZ,"meanposZ/F");           
  BTree->Branch("closestslicetime",&closestslicetime,"closestslicetime/F");           
  BTree->Branch("closestslicenhit",&closestslicenhit,"closestslicenhit/F");           
  //Nus Sandbox
  BTree->Branch("sumtx",&sumtx,"sumtx/F");                                           
  BTree->Branch("sumty",&sumty,"sumty/F");                                           
  BTree->Branch("ewsumtx",&ewsumtx,"ewsumtx/F");                                     
  BTree->Branch("ewsumty",&ewsumty,"ewsumty/F");                                     
  BTree->Branch("cossumtx",&cossumtx,"cossumtx/F");                                  
  BTree->Branch("cossumty",&cossumty,"cossumty/F");                                  
  BTree->Branch("cosewsumtx",&cosewsumtx,"cosewsumtx/F");                            
  BTree->Branch("cosewsumty",&cosewsumty,"cosewsumty/F");                            
  BTree->Branch("angsumtx",&angsumtx,"angsumtx/F");                                  
  BTree->Branch("angsumty",&angsumty,"angsumty/F");                                  
  BTree->Branch("angewsumtx",&angewsumtx,"angewsumtx/F");                            
  BTree->Branch("angewsumty",&angewsumty,"angewsumty/F");                            
  //Vertex Info
  BTree->Branch("vtxTime",&vtxTime,"vtxTime/F");  
  BTree->Branch("vtxX",&vtxX,"vtxX/F");           
  BTree->Branch("vtxY",&vtxY,"vtxY/F");           
  BTree->Branch("vtxZ",&vtxZ,"vtxZ/F");           
  BTree->Branch("vtxpngNhit",&vtxpngNhit,"vtxpngNhit/F");  
  BTree->Branch("vtxpngNhitX",&vtxpngNhitX,"vtxpngNhitX/F"); 
  BTree->Branch("vtxpngNhitY",&vtxpngNhitY,"vtxpngNhitY/F"); 
  //FuzzyK Prongs                                                                     
  BTree->Branch("fuzzykntot",&fuzzykntot,"fuzzykntot/F");                                      
  BTree->Branch("fuzzyknpng",&fuzzyknpng,"fuzzyknpng/F");                                
  BTree->Branch("fuzzyknpng2d",&fuzzyknpng2d,"fuzzyknpng2d/F");                       
  BTree->Branch("fuzzykpngdirx",&fuzzykpngdirx,"fuzzykpngdirx/F");                          
  BTree->Branch("fuzzykpngdiry",&fuzzykpngdiry,"fuzzykpngdiry/F");                          
  BTree->Branch("fuzzykpngdirz",&fuzzykpngdirz,"fuzzykpngdirz/F");                          
  BTree->Branch("fuzzykpngCalE",&fuzzykpngCalE,"fuzzykpngCalE/F");                             
  BTree->Branch("fuzzykpnglen",&fuzzykpnglen,"fuzzykpnglen/F");                             
  BTree->Branch("fuzzykpngnhit",&fuzzykpngnhit,"fuzzykpngnhit/F");                             
  BTree->Branch("fuzzykpngnhitx",&fuzzykpngnhitx,"fuzzykpngnhitx/F");                             
  BTree->Branch("fuzzykpngnhity",&fuzzykpngnhity,"fuzzykpngnhity/F");                             
  //CVN
  BTree->Branch("ncid",&ncid,"ncid/F");                
  BTree->Branch("nueid",&nueid,"nueid/F");      
  BTree->Branch("numuid",&numuid,"numuid/F");   
  BTree->Branch("nutauid",&nutauid,"nutauid/F");
  BTree->Branch("cosmicid",&cosmicid,"cosmicid/F");    
  /////////////////////////////////////////////////////// Signal Tree ////////////////////////////////////////////////////////////
  //True Info   
  STree->Branch("Run",&Run,"Run/I");                                          
  STree->Branch("SubRun",&SubRun,"SubRun/I");                                 
  STree->Branch("Evt",&Evt,"Evt/I");                                          
  STree->Branch("SubEvt",&SubEvt,"SubEvt/I");                                 
  STree->Branch("IsCC",&IsCC,"IsCC/I");                                       
  //Slice Info
  STree->Branch("nhit",&nhit,"nhit/F");              
  STree->Branch("ncalhit",&ncalhit,"ncalhit/F");              
  STree->Branch("nmiphit",&nmiphit,"nmiphit/F");          
  STree->Branch("ncontplanes",&ncontplanes,"ncontplanes/F");  
  STree->Branch("ncellsfromedge",&ncellsfromedge,"ncellsfromedge/F");  
  STree->Branch("slcCalE",&slcCalE,"slcCalE/F");        
  STree->Branch("boxminy",&boxminy,"boxminy/F");              
  STree->Branch("boxmaxy",&boxmaxy,"boxmaxy/F");              
  STree->Branch("starttime",&starttime,"starttime/F");           
  STree->Branch("meantime",&meantime,"meantime/F");           
  STree->Branch("endtime",&endtime,"endtime/F");           
  STree->Branch("meanposX",&meanposX,"meanposX/F");           
  STree->Branch("meanposY",&meanposY,"meanposY/F");           
  STree->Branch("meanposZ",&meanposZ,"meanposZ/F");           
  STree->Branch("closestslicetime",&closestslicetime,"closestslicetime/F");           
  STree->Branch("closestslicenhit",&closestslicenhit,"closestslicenhit/F");           
  //Nus Sand Box
  STree->Branch("sumtx",&sumtx,"sumtx/F");                                    
  STree->Branch("sumty",&sumty,"sumty/F");                                    
  STree->Branch("ewsumtx",&ewsumtx,"ewsumtx/F");                              
  STree->Branch("ewsumty",&ewsumty,"ewsumty/F");                              
  STree->Branch("cossumtx",&cossumtx,"cossumtx/F");                           
  STree->Branch("cossumty",&cossumty,"cossumty/F");                           
  STree->Branch("cosewsumtx",&cosewsumtx,"cosewsumtx/F");                     
  STree->Branch("cosewsumty",&cosewsumty,"cosewsumty/F");                     
  STree->Branch("angsumtx",&angsumtx,"angsumtx/F");                           
  STree->Branch("angsumty",&angsumty,"angsumty/F");                           
  STree->Branch("angewsumtx",&angewsumtx,"angewsumtx/F");                     
  STree->Branch("angewsumty",&angewsumty,"angewsumty/F");                     
  //Vertex Info
  STree->Branch("vtxTime",&vtxTime,"vtxTime/F");                        
  STree->Branch("vtxX",&vtxX,"vtxX/F");                                 
  STree->Branch("vtxY",&vtxY,"vtxY/F");                                 
  STree->Branch("vtxZ",&vtxZ,"vtxZ/F");                                 
  STree->Branch("vtxpngNhit",&vtxpngNhit,"vtxpngNhit/F");               
  STree->Branch("vtxpngNhitX",&vtxpngNhitX,"vtxpngNhitX/F");            
  STree->Branch("vtxpngNhitY",&vtxpngNhitY,"vtxpngNhitY/F");            
  //FuzzyK Prongs                                                                     
  STree->Branch("fuzzykntot",&fuzzykntot,"fuzzykntot/F");                                      
  STree->Branch("fuzzyknpng",&fuzzyknpng,"fuzzyknpng/F");                                
  STree->Branch("fuzzyknpng2d",&fuzzyknpng2d,"fuzzyknpng2d/F");                       
  STree->Branch("fuzzykpngdirx",&fuzzykpngdirx,"fuzzykpngdirx/F");                          
  STree->Branch("fuzzykpngdiry",&fuzzykpngdiry,"fuzzykpngdiry/F");                          
  STree->Branch("fuzzykpngdirz",&fuzzykpngdirz,"fuzzykpngdirz/F");                          
  STree->Branch("fuzzykpngCalE",&fuzzykpngCalE,"fuzzykpngCalE/F");                             
  STree->Branch("fuzzykpnglen",&fuzzykpnglen,"fuzzykpnglen/F");                             
  STree->Branch("fuzzykpngnhit",&fuzzykpngnhit,"fuzzykpngnhit/F");                             
  STree->Branch("fuzzykpngnhitx",&fuzzykpngnhitx,"fuzzykpngnhitx/F");                             
  STree->Branch("fuzzykpngnhity",&fuzzykpngnhity,"fuzzykpngnhity/F");                             
  //CVN
  STree->Branch("ncid",&ncid,"ncid/F");                                 
  STree->Branch("nueid",&nueid,"nueid/F");                                          
  STree->Branch("numuid",&numuid,"numuid/F");                                       
  STree->Branch("nutauid",&nutauid,"nutauid/F");                                    
  STree->Branch("cosmicid",&cosmicid,"cosmicid/F");                     
  ///////////////////////////////////////////////////////////////////////////////////////////////////  
  IFileSource* nonswap=loadernc.WildcardOrSAMQuery(fnamenc);
  int AllNon = nonswap->NFiles();
  Progress* prog = 0;
  TFile* nonfile;
  int NumNon = -1;
  while(TFile* nonfile =(TFile*)nonswap->GetNextFile()){               //Straring to loop over files
    ++NumNon;
    std::cout<<"******************************NonSwap File NUMBER: "<<NumNon<<std::endl;
    if(AllNon >= 0 && !prog) prog = new Progress(TString::Format("Filling tree from %d files", AllNon).Data());
    loadernc.HandleFile(nonfile, AllNon == 1 ? prog : 0);
    if(AllNon > 1 && prog) prog->SetProgress((NumNon+1.)/AllNon);
    TTree *recTreenon = (TTree*)nonfile->Get("recTree");
    caf::StandardRecord* recTreeObject = 0;
    recTreenon->SetBranchAddress("rec", &recTreeObject);
    // Turn off all the branches first
    recTreenon->SetBranchStatus("*",0);                                      
    //Turn on the useful branches from the beam file
    //Ture Info
    recTreenon->SetBranchStatus("hdr.run",1);                                      
    recTreenon->SetBranchStatus("hdr.subrun",1);                                   
    recTreenon->SetBranchStatus("hdr.evt",1);                                      
    recTreenon->SetBranchStatus("hdr.subevt",1);                                   
    recTreenon->SetBranchStatus("mc.nnu",1);                                       
    recTreenon->SetBranchStatus("mc.nu",1);                                        
    recTreenon->SetBranchStatus("mc.nu.iscc",1);                                   
    recTreenon->SetBranchStatus("mc.nu.pdg",1);                                    
    //Slice Info
    recTreenon->SetBranchStatus("slc.nhit",1);          
    recTreenon->SetBranchStatus("slc.ncalhit",1);        
    recTreenon->SetBranchStatus("slc.nmiphit",1);          
    recTreenon->SetBranchStatus("slc.ncontplanes",1);   
    recTreenon->SetBranchStatus("slc.ncellsfromedge",1);
    recTreenon->SetBranchStatus("slc.starttime",1);     
    recTreenon->SetBranchStatus("slc.endtime",1);       
    recTreenon->SetBranchStatus("slc.meantime",1);      
    recTreenon->SetBranchStatus("slc.meanpos.x",1);    
    recTreenon->SetBranchStatus("slc.meanpos.y",1);    
    recTreenon->SetBranchStatus("slc.meanpos.z",1);    
    recTreenon->SetBranchStatus("slc.calE",1);    
    recTreenon->SetBranchStatus("slc.boxmin.y",1);    
    recTreenon->SetBranchStatus("slc.boxmax.y",1);    
    recTreenon->SetBranchStatus("slc.closestslicetime",1);    
    recTreenon->SetBranchStatus("slc.closestslicenhit",1);    
    //Nus Sandbox 
    recTreenon->SetBranchStatus("sand.nus.sumtx",1);       
    recTreenon->SetBranchStatus("sand.nus.sumty",1);       
    recTreenon->SetBranchStatus("sand.nus.ewsumtx",1);     
    recTreenon->SetBranchStatus("sand.nus.ewsumty",1);     
    recTreenon->SetBranchStatus("sand.nus.cossumtx",1);    
    recTreenon->SetBranchStatus("sand.nus.cossumty",1);    
    recTreenon->SetBranchStatus("sand.nus.cosewsumtx",1);  
    recTreenon->SetBranchStatus("sand.nus.cosewsumty",1);  
    recTreenon->SetBranchStatus("sand.nus.angsumtx",1);    
    recTreenon->SetBranchStatus("sand.nus.angsumty",1);    
    recTreenon->SetBranchStatus("sand.nus.angewsumtx",1);  
    recTreenon->SetBranchStatus("sand.nus.angewsumty",1);  
    //Vertex Info 
    recTreenon->SetBranchStatus("vtx.elastic.time",1);                             
    recTreenon->SetBranchStatus("vtx.elastic.vtx.x",1);                           
    recTreenon->SetBranchStatus("vtx.elastic.vtx.y",1);                           
    recTreenon->SetBranchStatus("vtx.elastic.vtx.z",1);                           
    recTreenon->SetBranchStatus("vtx.nelastic",1);
    //FuzzyK Prongs                                                                                                       
    recTreenon->SetBranchStatus("vtx.elastic.fuzzyk.ntot",1);                             
    recTreenon->SetBranchStatus("vtx.elastic.fuzzyk.npng",1);                           
    recTreenon->SetBranchStatus("vtx.elastic.fuzzyk.npng2d",1);                           
    recTreenon->SetBranchStatus("vtx.elastic.fuzzyk.png",1);                                   
    //CVN
    recTreenon->SetBranchStatus("sel.cvn.ncid",1);                               
    recTreenon->SetBranchStatus("sel.cvn.nueid",1);                               
    recTreenon->SetBranchStatus("sel.cvn.numuid",1);                              
    recTreenon->SetBranchStatus("sel.cvn.nutauid",1);                             
    recTreenon->SetBranchStatus("sel.cvn.output",1);  
    //NueCosRej         
    recTreenon->SetBranchStatus("sel.nuecosrej.distallpngtop"   ,1);                         
    recTreenon->SetBranchStatus("sel.nuecosrej.distallpngbottom",1);                         
    recTreenon->SetBranchStatus("sel.nuecosrej.distallpngeast"  ,1);                        
    recTreenon->SetBranchStatus("sel.nuecosrej.distallpngwest"  ,1);                         
    recTreenon->SetBranchStatus("sel.nuecosrej.distallpngfront" ,1);                         
    recTreenon->SetBranchStatus("sel.nuecosrej.distallpngback"  ,1);                         
    Int_t nbeam = recTreenon->GetEntries(); 
    for (Int_t j = 0; j < nbeam; ++j) {                                        //Starting to loop over cosmic events  
      recTreenon->GetEntry(j);      
      cerr << "\r-- Processing event " << j << " of " << nbeam;
      Short_t nnu   = recTreeObject->mc.nnu;
      if(nnu != 1) continue;      
      Int_t   nccc  = recTreeObject->mc.nu[0].iscc;
      Short_t pdg   = recTreeObject->mc.nu[0].pdg;
      // NC Offical Quality Cuts                               
      if(recTreeObject->vtx.nelastic == 0) continue;
      if(recTreeObject->sel.nuecosrej.hitsperplane >= 8) continue;
      if(recTreeObject->vtx.elastic[0].fuzzyk.npng == 0) continue;
      if(recTreeObject->slc.ncontplanes <= 2) continue;
      // NC Official Containment Cuts 
      if(recTreeObject->sel.nuecosrej.distallpngtop    < 100) continue;
      if(recTreeObject->sel.nuecosrej.distallpngbottom < 10 ) continue;
      if(recTreeObject->sel.nuecosrej.distallpngeast   < 50 ) continue;
      if(recTreeObject->sel.nuecosrej.distallpngwest   < 50 ) continue;
      if(recTreeObject->sel.nuecosrej.distallpngfront  < 50 ) continue;
      if(recTreeObject->sel.nuecosrej.distallpngback   < 50 ) continue; 
      // NC Official NC/CC rejection             
      if(recTreeObject->slc.nhit < 25)  continue;               
      // Energy Cut
      if(recTreeObject->slc.calE < 0.25) continue;
      //if(recTreeObject->slc.calE > 10.0) continue;
      ///////////////////For Nus SandBox////////////////
      if(std::isnan(recTreeObject->sand.nus.angsumty))   continue;      
      if(std::isnan(recTreeObject->sand.nus.angewsumty)) continue;      
      if(std::isnan(recTreeObject->sand.nus.cossumty))   continue;      
      if(std::isnan(recTreeObject->sand.nus.cosewsumty)) continue;            
      if(std::isnan(recTreeObject->sand.nus.angsumtx))   continue;
      if(std::isnan(recTreeObject->sand.nus.angewsumtx)) continue;
      if(std::isnan(recTreeObject->sand.nus.cossumtx))   continue;
      if(std::isnan(recTreeObject->sand.nus.cosewsumtx)) continue;
      //True Info    
      Run               = recTreeObject->hdr.run;      
      SubRun            = recTreeObject->hdr.subrun;   
      Evt               = recTreeObject->hdr.evt;      
      SubEvt            = recTreeObject->hdr.subevt;   
      //Slice Info   
      nhit              = recTreeObject->slc.nhit;       
      ncalhit           = recTreeObject->slc.ncalhit;    
      nmiphit           = recTreeObject->slc.nmiphit;    
      ncontplanes       = recTreeObject->slc.ncontplanes;      
      ncellsfromedge    = recTreeObject->slc.ncellsfromedge;  
      slcCalE           = recTreeObject->slc.calE;
      boxminy           = recTreeObject->slc.boxmin.y;
      boxmaxy           = recTreeObject->slc.boxmax.y;
      meanposX          = recTreeObject->slc.meanpos.x;
      meanposY          = recTreeObject->slc.meanpos.y;
      meanposZ          = recTreeObject->slc.meanpos.z;
      starttime         = recTreeObject->slc.starttime; 
      endtime           = recTreeObject->slc.endtime;   
      meantime          = recTreeObject->slc.meantime;  
      closestslicetime  = recTreeObject->slc.closestslicetime;
      closestslicenhit  = recTreeObject->slc.closestslicenhit;
      //NuS Sandbox
      sumtx             = recTreeObject->sand.nus.sumtx;                             
      sumty             = recTreeObject->sand.nus.sumty;                             
      ewsumtx           = recTreeObject->sand.nus.ewsumtx;                           
      ewsumty           = recTreeObject->sand.nus.ewsumty;                           
      cossumtx          = recTreeObject->sand.nus.cossumtx;                          
      cossumty          = recTreeObject->sand.nus.cossumty;                          
      cosewsumtx        = recTreeObject->sand.nus.cosewsumtx;                        
      cosewsumty        = recTreeObject->sand.nus.cosewsumty;                        
      angsumtx          = recTreeObject->sand.nus.angsumtx;                          
      angsumty          = recTreeObject->sand.nus.angsumty;                          
      angewsumtx        = recTreeObject->sand.nus.angewsumtx;                        
      angewsumty        = recTreeObject->sand.nus.angewsumty;                        
      //Vertex Info 
      vtxTime           = recTreeObject->vtx.elastic[0].time;
      vtxX              = recTreeObject->vtx.elastic[0].vtx.x;
      vtxY              = recTreeObject->vtx.elastic[0].vtx.y; 
      vtxZ              = recTreeObject->vtx.elastic[0].vtx.z; 
      //FuzzyK Prongs
      fuzzykntot        = recTreeObject->vtx.elastic[0].fuzzyk.ntot;
      fuzzyknpng        = recTreeObject->vtx.elastic[0].fuzzyk.npng;
      fuzzyknpng2d      = recTreeObject->vtx.elastic[0].fuzzyk.npng2d;
      fuzzykpngdirx     = recTreeObject->vtx.elastic[0].fuzzyk.png[0].dir.x;
      fuzzykpngdiry     = recTreeObject->vtx.elastic[0].fuzzyk.png[0].dir.y;
      fuzzykpngdirz     = recTreeObject->vtx.elastic[0].fuzzyk.png[0].dir.z;
      fuzzykpngCalE     = recTreeObject->vtx.elastic[0].fuzzyk.png[0].calE;
      fuzzykpnglen      = recTreeObject->vtx.elastic[0].fuzzyk.png[0].len;
      fuzzykpngnhit     = recTreeObject->vtx.elastic[0].fuzzyk.png[0].nhit;
      fuzzykpngnhitx    = recTreeObject->vtx.elastic[0].fuzzyk.png[0].nhitx;
      fuzzykpngnhity    = recTreeObject->vtx.elastic[0].fuzzyk.png[0].nhity;
      //CVN  
      ncid              = recTreeObject->sel.cvn.ncid;         
      nueid             = recTreeObject->sel.cvn.nueid; 
      numuid            = recTreeObject->sel.cvn.numuid; 
      nutauid           = recTreeObject->sel.cvn.nutauid;
      cosmicid          = recTreeObject->sel.cvn.output[14];         
      if(nccc == 0){ STree->Fill(); }
      else {BTree->Fill();}
    }                                                                         // Ending to loop over beam events  
  }                                                                           // Ending to loop over files
  
  IFileSource* cosmicdata=loadercos.WildcardOrSAMQuery(fnamecos);
  int AllCos = cosmicdata->NFiles();
  //std::cout<<"*****************************************NUMBER OF FILES: "<<AllCos<<std::endl;           
  //Progress* prog = 0;
  TFile* cosfile;
  int NumCos = -1;
  while(TFile* cosfile =(TFile*)cosmicdata->GetNextFile()){               //Straring to loop over files                                              
    ++NumCos;
    std::cout<<"*******************************************************Cosmic File NUMBER: "<<NumCos<<std::endl;                             
    if(AllCos >= 0 && !prog) prog = new Progress(TString::Format("Filling tree from %d files", AllCos).Data());
    loadercos.HandleFile(cosfile, AllCos == 1 ? prog : 0);
    if(AllCos > 1 && prog) prog->SetProgress((NumCos+1.)/AllCos);
    TTree *recTreecos = (TTree*)cosfile->Get("recTree");
    //caf::SRProxy* recTreeObject = 0;
    caf::StandardRecord* recTreeObject = 0;
    recTreecos->SetBranchAddress("rec", &recTreeObject);    
    // Turn off all the branches first
    recTreecos->SetBranchStatus("*",0);                                      
    // Turn on the useful branches from the cosmic file
    //Ture Info 
    recTreecos->SetBranchStatus("hdr.run",1);
    recTreecos->SetBranchStatus("hdr.subrun",1);
    recTreecos->SetBranchStatus("hdr.evt",1);
    recTreecos->SetBranchStatus("hdr.subevt",1);
    recTreecos->SetBranchStatus("mc.nnu",1);
    recTreecos->SetBranchStatus("mc.nu",1);
    recTreecos->SetBranchStatus("mc.nu.iscc",1);
    recTreecos->SetBranchStatus("mc.nu.pdg",1);
    //Slice Info
    recTreecos->SetBranchStatus("slc.nhit",1);            
    recTreecos->SetBranchStatus("slc.ncalhit",1);          
    recTreecos->SetBranchStatus("slc.nmiphit",1);         
    recTreecos->SetBranchStatus("slc.ncontplanes",1);     
    recTreecos->SetBranchStatus("slc.ncellsfromedge",1);  
    recTreecos->SetBranchStatus("slc.starttime",1);       
    recTreecos->SetBranchStatus("slc.endtime",1);         
    recTreecos->SetBranchStatus("slc.meantime",1);        
    recTreecos->SetBranchStatus("slc.meanpos.x",1);    
    recTreecos->SetBranchStatus("slc.meanpos.y",1);    
    recTreecos->SetBranchStatus("slc.meanpos.z",1);    
    recTreecos->SetBranchStatus("slc.calE",1);
    recTreecos->SetBranchStatus("slc.boxmin.y",1);
    recTreecos->SetBranchStatus("slc.boxmax.y",1);
    recTreecos->SetBranchStatus("slc.closestslicetime",1);
    recTreecos->SetBranchStatus("slc.closestslicenhit",1);
    //Nus Sandbox 
    recTreecos->SetBranchStatus("sand.nus.sumtx",1);                               
    recTreecos->SetBranchStatus("sand.nus.sumty",1);                               
    recTreecos->SetBranchStatus("sand.nus.ewsumtx",1);                             
    recTreecos->SetBranchStatus("sand.nus.ewsumty",1);                             
    recTreecos->SetBranchStatus("sand.nus.cossumtx",1);                            
    recTreecos->SetBranchStatus("sand.nus.cossumty",1);                            
    recTreecos->SetBranchStatus("sand.nus.cosewsumtx",1);                          
    recTreecos->SetBranchStatus("sand.nus.cosewsumty",1);                          
    recTreecos->SetBranchStatus("sand.nus.angsumtx",1);                            
    recTreecos->SetBranchStatus("sand.nus.angsumty",1);                            
    recTreecos->SetBranchStatus("sand.nus.angewsumtx",1);                          
    recTreecos->SetBranchStatus("sand.nus.angewsumty",1);                          
    //Vertex Info              
    recTreecos->SetBranchStatus("vtx.elastic.time",1);                        
    recTreecos->SetBranchStatus("vtx.elastic.vtx.x",1);                      
    recTreecos->SetBranchStatus("vtx.elastic.vtx.y",1);                      
    recTreecos->SetBranchStatus("vtx.elastic.vtx.z",1);                      
    recTreecos->SetBranchStatus("vtx.nelastic",1);                            
    //FuzzyK Prongs                                                                                                                                    
    recTreecos->SetBranchStatus("vtx.elastic.fuzzyk.ntot",1);
    recTreecos->SetBranchStatus("vtx.elastic.fuzzyk.npng",1);
    recTreecos->SetBranchStatus("vtx.elastic.fuzzyk.npng2d",1);
    recTreecos->SetBranchStatus("vtx.elastic.fuzzyk.png",1);
    //CVN   
    recTreecos->SetBranchStatus("sel.cvn.ncid",1);                                 
    recTreecos->SetBranchStatus("sel.cvn.nueid",1);                                
    recTreecos->SetBranchStatus("sel.cvn.numuid",1);                               
    recTreecos->SetBranchStatus("sel.cvn.nutauid",1);                              
    recTreecos->SetBranchStatus("sel.cvn.output",1);             
    //NueCosRej                                                                                                                             
    recTreecos->SetBranchStatus("sel.nuecosrej.distallpngtop"   ,1);
    recTreecos->SetBranchStatus("sel.nuecosrej.distallpngbottom",1);
    recTreecos->SetBranchStatus("sel.nuecosrej.distallpngeast"  ,1);
    recTreecos->SetBranchStatus("sel.nuecosrej.distallpngwest"  ,1);
    recTreecos->SetBranchStatus("sel.nuecosrej.distallpngfront" ,1);
    recTreecos->SetBranchStatus("sel.nuecosrej.distallpngback"  ,1);
    Int_t ncosmic = recTreecos->GetEntriesFast(); 
    for (Int_t i = 0; i < ncosmic; ++i) {                                       //Starting to loop over cosmic events  
      recTreecos->GetEntry(i);      
      cerr << "\r-- Processing event " << i << " of " << ncosmic;
      // NC Offical Quality Cuts                                                                                                       
      if(recTreeObject->vtx.nelastic == 0) continue;
      if(recTreeObject->sel.nuecosrej.hitsperplane >= 8) continue;
      if(recTreeObject->vtx.elastic[0].fuzzyk.npng == 0) continue;
      if(recTreeObject->slc.ncontplanes <= 2) continue;
      // NC Official Containment Cuts                                                                                            
      if(recTreeObject->sel.nuecosrej.distallpngtop    < 100) continue;
      if(recTreeObject->sel.nuecosrej.distallpngbottom < 10 ) continue;
      if(recTreeObject->sel.nuecosrej.distallpngeast   < 50 ) continue;
      if(recTreeObject->sel.nuecosrej.distallpngwest   < 50 ) continue;
      if(recTreeObject->sel.nuecosrej.distallpngfront  < 50 ) continue;
      if(recTreeObject->sel.nuecosrej.distallpngback   < 50 ) continue;
      // NC Official NC/CC rejection                                                                                                                           
      if (recTreeObject->slc.nhit < 25) continue;
      // Energy Cut                                                                                                                             
      if(recTreeObject->slc.calE < 0.25) continue;
      //if(recTreeObject->slc.calE > 10.0) continue;
      ///////////////////Nus SandBox////////////////                                                                                               
      if(std::isnan(recTreeObject->sand.nus.angsumty))   continue;
      if(std::isnan(recTreeObject->sand.nus.angewsumty)) continue;
      if(std::isnan(recTreeObject->sand.nus.cossumty))   continue;
      if(std::isnan(recTreeObject->sand.nus.cosewsumty)) continue;
      if(std::isnan(recTreeObject->sand.nus.angsumtx))   continue;
      if(std::isnan(recTreeObject->sand.nus.angewsumtx)) continue;
      if(std::isnan(recTreeObject->sand.nus.cossumtx))   continue;
      if(std::isnan(recTreeObject->sand.nus.cosewsumtx)) continue;
      //True Info 
      Run    = recTreeObject->hdr.run;
      SubRun = recTreeObject->hdr.subrun;
      Evt    = recTreeObject->hdr.evt;
      SubEvt = recTreeObject->hdr.subevt;
      //Slice Info
      nhit        = recTreeObject->slc.nhit;       
      ncalhit     = recTreeObject->slc.ncalhit;    
      nmiphit     = recTreeObject->slc.nmiphit;    
      ncontplanes = recTreeObject->slc.ncontplanes;
      ncellsfromedge = recTreeObject->slc.ncellsfromedge;
      slcCalE     = recTreeObject->slc.calE;
      boxminy     = recTreeObject->slc.boxmin.y;
      boxmaxy     = recTreeObject->slc.boxmax.y;        
      starttime   = recTreeObject->slc.starttime;  
      endtime     = recTreeObject->slc.endtime;    
      meantime    = recTreeObject->slc.meantime;   
      meanposX    = recTreeObject->slc.meanpos.x;  
      meanposY    = recTreeObject->slc.meanpos.y;  
      meanposZ    = recTreeObject->slc.meanpos.z;
      closestslicetime  = recTreeObject->slc.closestslicetime;
      closestslicenhit  = recTreeObject->slc.closestslicenhit;  
      //NuS Sandbox
      sumtx       = recTreeObject->sand.nus.sumtx;                           
      sumty       = recTreeObject->sand.nus.sumty;                           
      ewsumtx     = recTreeObject->sand.nus.ewsumtx;                         
      ewsumty     = recTreeObject->sand.nus.ewsumty;                         
      cossumtx    = recTreeObject->sand.nus.cossumtx;                        
      cossumty    = recTreeObject->sand.nus.cossumty;                        
      cosewsumtx  = recTreeObject->sand.nus.cosewsumtx;                      
      cosewsumty  = recTreeObject->sand.nus.cosewsumty;                      
      angsumtx    = recTreeObject->sand.nus.angsumtx;                        
      angsumty    = recTreeObject->sand.nus.angsumty;                        
      angewsumtx  = recTreeObject->sand.nus.angewsumtx;                      
      angewsumty  = recTreeObject->sand.nus.angewsumty;                      
      //Vertex Info 
      vtxTime     = recTreeObject->vtx.elastic[0].time;     
      vtxX        = recTreeObject->vtx.elastic[0].vtx.x;   
      vtxY        = recTreeObject->vtx.elastic[0].vtx.y;   
      vtxZ        = recTreeObject->vtx.elastic[0].vtx.z;   
      //FuzzyK Prongs                                                                    
      fuzzykntot      = recTreeObject->vtx.elastic[0].fuzzyk.ntot;
      fuzzyknpng      = recTreeObject->vtx.elastic[0].fuzzyk.npng;
      fuzzyknpng2d    = recTreeObject->vtx.elastic[0].fuzzyk.npng2d;
      fuzzykpngdirx   = recTreeObject->vtx.elastic[0].fuzzyk.png[0].dir.x;
      fuzzykpngdiry   = recTreeObject->vtx.elastic[0].fuzzyk.png[0].dir.y;
      fuzzykpngdirz   = recTreeObject->vtx.elastic[0].fuzzyk.png[0].dir.z;
      fuzzykpngCalE   = recTreeObject->vtx.elastic[0].fuzzyk.png[0].calE;
      fuzzykpnglen    = recTreeObject->vtx.elastic[0].fuzzyk.png[0].len;
      fuzzykpngnhit   = recTreeObject->vtx.elastic[0].fuzzyk.png[0].nhit;
      fuzzykpngnhitx  = recTreeObject->vtx.elastic[0].fuzzyk.png[0].nhitx;
      fuzzykpngnhity  = recTreeObject->vtx.elastic[0].fuzzyk.png[0].nhity;
      //CVN  
      ncid        = recTreeObject->sel.cvn.ncid;     
      nueid       = recTreeObject->sel.cvn.nueid;    
      numuid      = recTreeObject->sel.cvn.numuid;   
      nutauid     = recTreeObject->sel.cvn.nutauid;  
      cosmicid    = recTreeObject->sel.cvn.output[14]; 
      BTree->Fill();
    }                                                                         // Ending to loop over cosmic events 
  }                                                                           // Ending to loop over files
  tMVALoad.Write();
  tMVALoad.Close();
}
#endif
