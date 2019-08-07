#ifdef __CINT__
void  FD_M_TA() // For MultiClass Classification
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

#include "CAFAna/Core/EventList.h"

using namespace ana;

void FD_M_TA(){
  //Training and Testing 
  const std::string fnamecos = "prod_decaf_R17-03-01-prod3reco.h_fd_cosmic_full_nue_or_numu_or_nus_contain_v1_goodruns";
  //const std::string fnamecos = "defname:prod_decaf_R17-03-01-prod3reco.h_fd_cosmic_full_nue_or_numu_or_nus_contain_v1_goodruns with limit 10";
  SpectrumLoader loadercos(fnamecos);
  const std::string fnamenc  = "prod_decaf_R17-03-01-prod3reco.l_fd_genie_nonswap_fhc_nova_v08_full_nue_or_numu_or_nus_contain_v1";
  //const std::string fnamenc = "defname:prod_decaf_R17-03-01-prod3reco.l_fd_genie_nonswap_fhc_nova_v08_full_nue_or_numu_or_nus_contain_v1 with limit 10";
  SpectrumLoader loadernc(fnamenc);
  //Output file
  TFile tMVALoad("TRA_FD_M.root","RECREATE");
  //Four trees for four interactions
  TTree *ncTree = new TTree("ncTree", "NC event tree");
  TTree *muTree = new TTree("muTree", "Numu event tree");
  TTree *neTree = new TTree("neTree", "Nue event tree");
  TTree *csTree = new TTree("csTree", "Cosmic event tree");
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
  ///////////////////////There are 50 variables, about 30 input vars will be selected from them.////////////////////////////
  /////////////////////////////////////////////////////// Cosmic Tree //////////////////////////////////////////////////////
  //True Info
  csTree->Branch("Run",&Run,"Run/I");          
  csTree->Branch("SubRun",&SubRun,"SubRun/I"); 
  csTree->Branch("Evt",&Evt,"Evt/I");          
  csTree->Branch("SubEvt",&SubEvt,"SubEvt/I"); 
  csTree->Branch("IsCC",&IsCC,"IsCC/I");       
  //Slice Info
  csTree->Branch("nhit",&nhit,"nhit/F");              
  csTree->Branch("ncalhit",&ncalhit,"ncalhit/F");              
  csTree->Branch("nmiphit",&nmiphit,"nmiphit/F");          
  csTree->Branch("ncontplanes",&ncontplanes,"ncontplanes/F");  
  csTree->Branch("ncellsfromedge",&ncellsfromedge,"ncellsfromedge/F");  
  csTree->Branch("slcCalE",&slcCalE,"slcCalE/F");        
  csTree->Branch("boxminy",&boxminy,"boxminy/F");              
  csTree->Branch("boxmaxy",&boxmaxy,"boxmaxy/F");              
  csTree->Branch("starttime",&starttime,"starttime/F");           
  csTree->Branch("meantime",&meantime,"meantime/F");           
  csTree->Branch("endtime",&endtime,"endtime/F");           
  csTree->Branch("meanposX",&meanposX,"meanposX/F");           
  csTree->Branch("meanposY",&meanposY,"meanposY/F");           
  csTree->Branch("meanposZ",&meanposZ,"meanposZ/F");           
  csTree->Branch("closestslicetime",&closestslicetime,"closestslicetime/F");           
  csTree->Branch("closestslicenhit",&closestslicenhit,"closestslicenhit/F");           
  //Nus Sandbox
  csTree->Branch("sumtx",&sumtx,"sumtx/F");                                           
  csTree->Branch("sumty",&sumty,"sumty/F");                                           
  csTree->Branch("ewsumtx",&ewsumtx,"ewsumtx/F");                                     
  csTree->Branch("ewsumty",&ewsumty,"ewsumty/F");                                     
  csTree->Branch("cossumtx",&cossumtx,"cossumtx/F");                                  
  csTree->Branch("cossumty",&cossumty,"cossumty/F");                                  
  csTree->Branch("cosewsumtx",&cosewsumtx,"cosewsumtx/F");                            
  csTree->Branch("cosewsumty",&cosewsumty,"cosewsumty/F");                            
  csTree->Branch("angsumtx",&angsumtx,"angsumtx/F");                                  
  csTree->Branch("angsumty",&angsumty,"angsumty/F");                                  
  csTree->Branch("angewsumtx",&angewsumtx,"angewsumtx/F");                            
  csTree->Branch("angewsumty",&angewsumty,"angewsumty/F");                            
  //Vertex Info
  csTree->Branch("vtxTime",&vtxTime,"vtxTime/F");  
  csTree->Branch("vtxX",&vtxX,"vtxX/F");           
  csTree->Branch("vtxY",&vtxY,"vtxY/F");           
  csTree->Branch("vtxZ",&vtxZ,"vtxZ/F");           
  csTree->Branch("vtxpngNhit",&vtxpngNhit,"vtxpngNhit/F");  
  csTree->Branch("vtxpngNhitX",&vtxpngNhitX,"vtxpngNhitX/F"); 
  csTree->Branch("vtxpngNhitY",&vtxpngNhitY,"vtxpngNhitY/F"); 
  //FuzzyK Prongs                                                                     
  csTree->Branch("fuzzykntot",    &fuzzykntot,     "fuzzykntot/F");                                      
  csTree->Branch("fuzzyknpng",    &fuzzyknpng,     "fuzzyknpng/F");                                
  csTree->Branch("fuzzyknpng2d",  &fuzzyknpng2d,   "fuzzyknpng2d/F");                       
  csTree->Branch("fuzzykpngdirx", &fuzzykpngdirx,  "fuzzykpngdirx/F");                          
  csTree->Branch("fuzzykpngdiry", &fuzzykpngdiry,  "fuzzykpngdiry/F");                          
  csTree->Branch("fuzzykpngdirz", &fuzzykpngdirz,  "fuzzykpngdirz/F");                          
  csTree->Branch("fuzzykpngCalE", &fuzzykpngCalE,  "fuzzykpngCalE/F");                             
  csTree->Branch("fuzzykpnglen",  &fuzzykpnglen,   "fuzzykpnglen/F");                             
  csTree->Branch("fuzzykpngnhit", &fuzzykpngnhit,  "fuzzykpngnhit/F");                             
  csTree->Branch("fuzzykpngnhitx",&fuzzykpngnhitx, "fuzzykpngnhitx/F");                             
  csTree->Branch("fuzzykpngnhity",&fuzzykpngnhity, "fuzzykpngnhity/F");                             
  //CVN
  csTree->Branch("ncid",&ncid,"ncid/F");                
  csTree->Branch("nueid",&nueid,"nueid/F");      
  csTree->Branch("numuid",&numuid,"numuid/F");   
  csTree->Branch("nutauid",&nutauid,"nutauid/F");
  csTree->Branch("cosmicid",&cosmicid,"cosmicid/F");    
  /////////////////////////////////////////////////////// NC Tree ////////////////////////////////////////////////////////////
  //True Info   
  ncTree->Branch("Run",&Run,"Run/I");                                          
  ncTree->Branch("SubRun",&SubRun,"SubRun/I");                                 
  ncTree->Branch("Evt",&Evt,"Evt/I");                                          
  ncTree->Branch("SubEvt",&SubEvt,"SubEvt/I");                                 
  ncTree->Branch("IsCC",&IsCC,"IsCC/I");                                       
  //Slice Info
  ncTree->Branch("nhit",&nhit,"nhit/F");              
  ncTree->Branch("ncalhit",&ncalhit,"ncalhit/F");              
  ncTree->Branch("nmiphit",&nmiphit,"nmiphit/F");          
  ncTree->Branch("ncontplanes",&ncontplanes,"ncontplanes/F");  
  ncTree->Branch("ncellsfromedge",&ncellsfromedge,"ncellsfromedge/F");  
  ncTree->Branch("slcCalE",&slcCalE,"slcCalE/F");        
  ncTree->Branch("boxminy",&boxminy,"boxminy/F");              
  ncTree->Branch("boxmaxy",&boxmaxy,"boxmaxy/F");              
  ncTree->Branch("starttime",&starttime,"starttime/F");           
  ncTree->Branch("meantime",&meantime,"meantime/F");           
  ncTree->Branch("endtime",&endtime,"endtime/F");           
  ncTree->Branch("meanposX",&meanposX,"meanposX/F");           
  ncTree->Branch("meanposY",&meanposY,"meanposY/F");           
  ncTree->Branch("meanposZ",&meanposZ,"meanposZ/F");           
  ncTree->Branch("closestslicetime",&closestslicetime,"closestslicetime/F");           
  ncTree->Branch("closestslicenhit",&closestslicenhit,"closestslicenhit/F");           
  //Nus Sand Box
  ncTree->Branch("sumtx",&sumtx,"sumtx/F");                                    
  ncTree->Branch("sumty",&sumty,"sumty/F");                                    
  ncTree->Branch("ewsumtx",&ewsumtx,"ewsumtx/F");                              
  ncTree->Branch("ewsumty",&ewsumty,"ewsumty/F");                              
  ncTree->Branch("cossumtx",&cossumtx,"cossumtx/F");                           
  ncTree->Branch("cossumty",&cossumty,"cossumty/F");                           
  ncTree->Branch("cosewsumtx",&cosewsumtx,"cosewsumtx/F");                     
  ncTree->Branch("cosewsumty",&cosewsumty,"cosewsumty/F");                     
  ncTree->Branch("angsumtx",&angsumtx,"angsumtx/F");                           
  ncTree->Branch("angsumty",&angsumty,"angsumty/F");                           
  ncTree->Branch("angewsumtx",&angewsumtx,"angewsumtx/F");                     
  ncTree->Branch("angewsumty",&angewsumty,"angewsumty/F");                     
  //Vertex Info
  ncTree->Branch("vtxTime",&vtxTime,"vtxTime/F");                        
  ncTree->Branch("vtxX",&vtxX,"vtxX/F");                                 
  ncTree->Branch("vtxY",&vtxY,"vtxY/F");                                 
  ncTree->Branch("vtxZ",&vtxZ,"vtxZ/F");                                 
  ncTree->Branch("vtxpngNhit",&vtxpngNhit,"vtxpngNhit/F");               
  ncTree->Branch("vtxpngNhitX",&vtxpngNhitX,"vtxpngNhitX/F");            
  ncTree->Branch("vtxpngNhitY",&vtxpngNhitY,"vtxpngNhitY/F");            
  //FuzzyK Prongs                                                                     
  ncTree->Branch("fuzzykntot",    &fuzzykntot,     "fuzzykntot/F");                                      
  ncTree->Branch("fuzzyknpng",    &fuzzyknpng,     "fuzzyknpng/F");                                
  ncTree->Branch("fuzzyknpng2d",  &fuzzyknpng2d,   "fuzzyknpng2d/F");                       
  ncTree->Branch("fuzzykpngdirx", &fuzzykpngdirx,  "fuzzykpngdirx/F");                          
  ncTree->Branch("fuzzykpngdiry", &fuzzykpngdiry,  "fuzzykpngdiry/F");                          
  ncTree->Branch("fuzzykpngdirz", &fuzzykpngdirz,  "fuzzykpngdirz/F");                          
  ncTree->Branch("fuzzykpngCalE", &fuzzykpngCalE,  "fuzzykpngCalE/F");                             
  ncTree->Branch("fuzzykpnglen",  &fuzzykpnglen,   "fuzzykpnglen/F");                             
  ncTree->Branch("fuzzykpngnhit", &fuzzykpngnhit,  "fuzzykpngnhit/F");                             
  ncTree->Branch("fuzzykpngnhitx",&fuzzykpngnhitx, "fuzzykpngnhitx/F");                             
  ncTree->Branch("fuzzykpngnhity",&fuzzykpngnhity, "fuzzykpngnhity/F");                             
  //CVN
  ncTree->Branch("ncid",&ncid,"ncid/F");                                 
  ncTree->Branch("nueid",&nueid,"nueid/F");                                          
  ncTree->Branch("numuid",&numuid,"numuid/F");                                       
  ncTree->Branch("nutauid",&nutauid,"nutauid/F");                                    
  ncTree->Branch("cosmicid",&cosmicid,"cosmicid/F");                     
  ///////////////////////////////////////////////// Numu Tree ///////////////////////////////////////////////
  //True Info
  muTree->Branch("Run",&Run,"Run/I");             
  muTree->Branch("SubRun",&SubRun,"SubRun/I");    
  muTree->Branch("Evt",&Evt,"Evt/I");             
  muTree->Branch("SubEvt",&SubEvt,"SubEvt/I");    
  muTree->Branch("IsCC",&IsCC,"IsCC/I");          
  //Slice Info
  muTree->Branch("nhit",&nhit,"nhit/F");              
  muTree->Branch("ncalhit",&ncalhit,"ncalhit/F");              
  muTree->Branch("nmiphit",&nmiphit,"nmiphit/F");          
  muTree->Branch("ncontplanes",&ncontplanes,"ncontplanes/F");  
  muTree->Branch("ncellsfromedge",&ncellsfromedge,"ncellsfromedge/F");  
  muTree->Branch("slcCalE",&slcCalE,"slcCalE/F");        
  muTree->Branch("boxminy",&boxminy,"boxminy/F");              
  muTree->Branch("boxmaxy",&boxmaxy,"boxmaxy/F");              
  muTree->Branch("starttime",&starttime,"starttime/F");           
  muTree->Branch("meantime",&meantime,"meantime/F");           
  muTree->Branch("endtime",&endtime,"endtime/F");           
  muTree->Branch("meanposX",&meanposX,"meanposX/F");           
  muTree->Branch("meanposY",&meanposY,"meanposY/F");           
  muTree->Branch("meanposZ",&meanposZ,"meanposZ/F");           
  muTree->Branch("closestslicetime",&closestslicetime,"closestslicetime/F");           
  muTree->Branch("closestslicenhit",&closestslicenhit,"closestslicenhit/F");           
  //Nus Sand Box
  muTree->Branch("sumtx",&sumtx,"sumtx/F");       
  muTree->Branch("sumty",&sumty,"sumty/F");       
  muTree->Branch("ewsumtx",&ewsumtx,"ewsumtx/F"); 
  muTree->Branch("ewsumty",&ewsumty,"ewsumty/F"); 
  muTree->Branch("cossumtx",&cossumtx,"cossumtx/F");
  muTree->Branch("cossumty",&cossumty,"cossumty/F");
  muTree->Branch("cosewsumtx",&cosewsumtx,"cosewsumtx/F");
  muTree->Branch("cosewsumty",&cosewsumty,"cosewsumty/F");
  muTree->Branch("angsumtx",&angsumtx,"angsumtx/F");      
  muTree->Branch("angsumty",&angsumty,"angsumty/F");      
  muTree->Branch("angewsumtx",&angewsumtx,"angewsumtx/F");
  muTree->Branch("angewsumty",&angewsumty,"angewsumty/F");
  //Vertex Info
  muTree->Branch("vtxTime",&vtxTime,"vtxTime/F");
  muTree->Branch("vtxX",&vtxX,"vtxX/F");         
  muTree->Branch("vtxY",&vtxY,"vtxY/F");         
  muTree->Branch("vtxZ",&vtxZ,"vtxZ/F");         
  muTree->Branch("vtxpngNhit",&vtxpngNhit,"vtxpngNhit/F"); 
  muTree->Branch("vtxpngNhitX",&vtxpngNhitX,"vtxpngNhitX/F");
  muTree->Branch("vtxpngNhitY",&vtxpngNhitY,"vtxpngNhitY/F");
  //FuzzyK Prongs                                                                     
  muTree->Branch("fuzzykntot",     &fuzzykntot,     "fuzzykntot/F");                                      
  muTree->Branch("fuzzyknpng",     &fuzzyknpng,     "fuzzyknpng/F");                                
  muTree->Branch("fuzzyknpng2d",   &fuzzyknpng2d,   "fuzzyknpng2d/F");                       
  muTree->Branch("fuzzykpngdirx",  &fuzzykpngdirx,  "fuzzykpngdirx/F");                          
  muTree->Branch("fuzzykpngdiry",  &fuzzykpngdiry,  "fuzzykpngdiry/F");                          
  muTree->Branch("fuzzykpngdirz",  &fuzzykpngdirz,  "fuzzykpngdirz/F");                          
  muTree->Branch("fuzzykpngCalE",  &fuzzykpngCalE,  "fuzzykpngCalE/F");                             
  muTree->Branch("fuzzykpnglen",   &fuzzykpnglen,   "fuzzykpnglen/F");                             
  muTree->Branch("fuzzykpngnhit",  &fuzzykpngnhit,  "fuzzykpngnhit/F");                             
  muTree->Branch("fuzzykpngnhitx", &fuzzykpngnhitx, "fuzzykpngnhitx/F");                             
  muTree->Branch("fuzzykpngnhity", &fuzzykpngnhity, "fuzzykpngnhity/F");                             
  //CVN
  muTree->Branch("ncid",&ncid,"ncid/F");                                             
  muTree->Branch("nueid",&nueid,"nueid/F");               
  muTree->Branch("numuid",&numuid,"numuid/F");            
  muTree->Branch("nutauid",&nutauid,"nutauid/F");         
  muTree->Branch("cosmicid",&cosmicid,"cosmicid/F");                                 
  //////////////////////////////////////////////////// Nue Tree ///////////////////////////////////////////
  //True Info 
  neTree->Branch("Run",&Run,"Run/I");           
  neTree->Branch("SubRun",&SubRun,"SubRun/I");  
  neTree->Branch("Evt",&Evt,"Evt/I");           
  neTree->Branch("SubEvt",&SubEvt,"SubEvt/I");  
  neTree->Branch("IsCC",&IsCC,"IsCC/I");        
  //Slice Info
  neTree->Branch("nhit",&nhit,"nhit/F");              
  neTree->Branch("ncalhit",&ncalhit,"ncalhit/F");              
  neTree->Branch("nmiphit",&nmiphit,"nmiphit/F");          
  neTree->Branch("ncontplanes",&ncontplanes,"ncontplanes/F");  
  neTree->Branch("ncellsfromedge",&ncellsfromedge,"ncellsfromedge/F");  
  neTree->Branch("slcCalE",&slcCalE,"slcCalE/F");        
  neTree->Branch("boxminy",&boxminy,"boxminy/F");              
  neTree->Branch("boxmaxy",&boxmaxy,"boxmaxy/F");              
  neTree->Branch("starttime",&starttime,"starttime/F");           
  neTree->Branch("meantime",&meantime,"meantime/F");           
  neTree->Branch("endtime",&endtime,"endtime/F");           
  neTree->Branch("meanposX",&meanposX,"meanposX/F");           
  neTree->Branch("meanposY",&meanposY,"meanposY/F");           
  neTree->Branch("meanposZ",&meanposZ,"meanposZ/F");           
  neTree->Branch("closestslicetime",&closestslicetime,"closestslicetime/F");           
  neTree->Branch("closestslicenhit",&closestslicenhit,"closestslicenhit/F");           
  //Nus Sand Box
  neTree->Branch("sumtx",&sumtx,"sumtx/F");             
  neTree->Branch("sumty",&sumty,"sumty/F");             
  neTree->Branch("ewsumtx",&ewsumtx,"ewsumtx/F");       
  neTree->Branch("ewsumty",&ewsumty,"ewsumty/F");       
  neTree->Branch("cossumtx",&cossumtx,"cossumtx/F");    
  neTree->Branch("cossumty",&cossumty,"cossumty/F");    
  neTree->Branch("cosewsumtx",&cosewsumtx,"cosewsumtx/F"); 
  neTree->Branch("cosewsumty",&cosewsumty,"cosewsumty/F"); 
  neTree->Branch("angsumtx",&angsumtx,"angsumtx/F");       
  neTree->Branch("angsumty",&angsumty,"angsumty/F");       
  neTree->Branch("angewsumtx",&angewsumtx,"angewsumtx/F"); 
  neTree->Branch("angewsumty",&angewsumty,"angewsumty/F"); 
  //Vertex Info
  neTree->Branch("vtxTime",&vtxTime,"vtxTime/F");                                 
  neTree->Branch("vtxX",&vtxX,"vtxX/F");                                          
  neTree->Branch("vtxY",&vtxY,"vtxY/F");                                          
  neTree->Branch("vtxZ",&vtxZ,"vtxZ/F");                                          
  neTree->Branch("vtxpngNhit",&vtxpngNhit,"vtxpngNhit/F");                        
  neTree->Branch("vtxpngNhitX",&vtxpngNhitX,"vtxpngNhitX/F");                     
  neTree->Branch("vtxpngNhitY",&vtxpngNhitY,"vtxpngNhitY/F");                     
  //FuzzyK Prongs                                                                     
  neTree->Branch("fuzzykntot",    &fuzzykntot,     "fuzzykntot/F");                                      
  neTree->Branch("fuzzyknpng",    &fuzzyknpng,     "fuzzyknpng/F");                                
  neTree->Branch("fuzzyknpng2d",  &fuzzyknpng2d,   "fuzzyknpng2d/F");                       
  neTree->Branch("fuzzykpngdirx", &fuzzykpngdirx,  "fuzzykpngdirx/F");                          
  neTree->Branch("fuzzykpngdiry", &fuzzykpngdiry,  "fuzzykpngdiry/F");                          
  neTree->Branch("fuzzykpngdirz", &fuzzykpngdirz,  "fuzzykpngdirz/F");                          
  neTree->Branch("fuzzykpngCalE", &fuzzykpngCalE,  "fuzzykpngCalE/F");                             
  neTree->Branch("fuzzykpnglen",  &fuzzykpnglen,   "fuzzykpnglen/F");                             
  neTree->Branch("fuzzykpngnhit", &fuzzykpngnhit,  "fuzzykpngnhit/F");                             
  neTree->Branch("fuzzykpngnhitx",&fuzzykpngnhitx, "fuzzykpngnhitx/F");                             
  neTree->Branch("fuzzykpngnhity",&fuzzykpngnhity, "fuzzykpngnhity/F");                             
  //CVN
  neTree->Branch("ncid",&ncid,"ncid/F");                                          
  neTree->Branch("nueid",&nueid,"nueid/F");                                          
  neTree->Branch("numuid",&numuid,"numuid/F");                                       
  neTree->Branch("nutauid",&nutauid,"nutauid/F");                                     
  neTree->Branch("cosmicid",&cosmicid,"cosmicid/F");                                
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
      ////////////////For Nus SandBox////////////////
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
      
      if(nccc == 0){ ncTree->Fill(); }
      else if(pdg == 14){muTree->Fill();}
      else {neTree->Fill();}
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
      csTree->Fill();
    }                                                                         // Ending to loop over cosmic events 
  }                                                                           // Ending to loop over files
  tMVALoad.Write();
  tMVALoad.Close();
}
#endif
