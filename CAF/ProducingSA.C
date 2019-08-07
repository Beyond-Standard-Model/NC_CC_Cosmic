#ifdef __CINT__
void  ProducingSA()
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
// ROOT includes
#include "TStyle.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TVector3.h"
#include "StandardRecord/StandardRecord.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"

using namespace ana;

void ProducingSA(){

  //Training and Testing                                                                                                                                                       
  const std::string fnamecos = "prod_decaf_R17-03-01-prod3reco.h_fd_cosmic_full_nue_or_numu_or_nus_contain_v1_goodruns";
  //const std::string fnamecos = "defname:prod_decaf_R17-03-01-prod3reco.h_fd_cosmic_full_nue_or_numu_or_nus_contain_v1_goodruns with limit 10";          
  SpectrumLoader loadercos(fnamecos);
  const std::string fnamenc = "prod_decaf_R17-03-01-prod3reco.l_fd_genie_nonswap_fhc_nova_v08_full_nue_or_numu_or_nus_contain_v1";
  //const std::string fnamenc = "defname:prod_decaf_R17-03-01-prod3reco.l_fd_genie_nonswap_fhc_nova_v08_full_nue_or_numu_or_nus_contain_v1 with limit 10";                    
  SpectrumLoader loadernc(fnamenc);

  
  ///Training and Testing  
  //const std::string fnamecos = "defname:prod_limitedcaf_R16-03-03-prod2reco.a_fd_cosmic_full_nueveto_v1_goodruns with limit 100";
  //SpectrumLoader loadercos(fnamecos);
  //const std::string fnamenc = "defname:prod_caf_R16-03-03-prod2reco.f_fd_genie_nonswap_fhc_nova_v08_period2_v1_prod2-snapshot with limit 100 ";
  //SpectrumLoader loadernc(fnamenc);
 
  // Output file
  TFile tMVALoad("TrainingSA.root","RECREATE");
  // Four trees for four interactions
  TTree *ncTree = new TTree("ncTree", "NC event tree");
  TTree *muTree = new TTree("muTree", "Numu event tree");
  TTree *neTree = new TTree("neTree", "Nue event tree");
  TTree *csTree = new TTree("csTree", "Cosmic event tree");
  // Definding the useful variables to transform into Input.root
  int   Run;  
  int   SubRun;
  int   Evt;
  int   SubEvt;
  int   IsCC;

  float ncid;                                                                         // CVN                          
  float nueid;                                                                        // CVN
  float numuid;                                                                       // CVN 
  float nutauid;                                                                      // CVN 
  float cosmicid;                                                                     // CVN             SY 1

  float partptp;                                                                      // NueCosRej       SY 2 We need to choose one from them    
  float pngptp;                                                                       // NueCosRej       SY 2                   
  float starttop;                                                                     // NueCosRej
  float startbottom;                                                                  // NueCosRej
  float startfront;                                                                   // NueCosRej
  float startback;                                                                    // NueCosRej
  float startwest;                                                                    // NueCosRej
  float starteast;                                                                    // NueCosRej
  float stoptop;                                                                      // NueCosRej
  float stopbottom;                                                                   // NueCosRej
  float stopfront;                                                                    // NueCosRej
  float stopback;                                                                     // NueCosRej
  float stopwest;                                                                     // NueCosRej
  float stopeast;                                                                     // NueCosRej
  float cosdang;                                                                      // NueCosRej
  float vtxdoca;                                                                      // NueCosRej
  float prongmaxx;                                                                    // NueCosRej
  float prongmaxy;                                                                    // NueCosRej
  float prongmaxz;                                                                    // NueCosRej
  float prongminx;                                                                    // NueCosRej
  float prongminy;                                                                    // NueCosRej
  float prongminz;                                                                    // NueCosRej

  float sparsenessasymm;                                                              // NueCosRej      I wanna check this one 

  float hitsperplaneasymm;                                                            // NueCosRej      I wanna check this one  
  float hitsperplane;                                                                 // NueCosRej
  float musliceidx;                                                                   // NueCosRej
  float muanglediff;                                                                  // NueCosRej
  float mutimediff;                                                                   // NueCosRej
  float muclosestapproach;                                                            // NueCosRej

  /// Do we want to use Prong to replace Showe at all!!!
  float shwnhit;                                                                      // Shower Branch                SY 3
  float shwnhitx;                                                                     // Shower Branch                SY 4
  float shwnhity;                                                                     // Shower Branch                SY 5
  float shwxminusy;                                                                   // shwnhitx - shwnhity          SY 6
  float shwxplusy;                                                                    // shwnhitx + shwnhity          SY 7
  float shwxovery;                                                                    // shwnhitx/shwnhity            SY 8
  float shwnplane;                                                                    // Shower Branch                
  float shwmaxplanecont;                                                              // Shower Branch                
  float shwmaxplanegap;                                                               // Shower Branch                
  float shwcalE;                                                                      // Shower Branch                SY 9               
  float shwstartX;                                                                    // Shower Branch
  float shwstartY;                                                                    // Shower Branch                
  float shwstartZ;                                                                    // Shower Branch                
  float shwdirX;                                                                      // Shower Branch
  float shwdirY;                                                                      // Shower Branch                SY 10
  float shwdirZ;                                                                      // Shower Branch
  float shwlen;                                                                       // Shower Branch                SY 11
  float shwwwidth;                                                                    // Shower Branch                SY 12
  float shwNplaneX;                                                                   // Shower Branch               
  float shwNplaneY;                                                                   // Shower Branch
  float shwGap;                                                                       // Shower Branch                SY 13
  float shwstopX;                                                                     // Shower Branch
  float shwstopY;                                                                     // Shower Branch
  float shwstopZ;                                                                     // Shower Branch
  float nshwlid;                                                                      // Shower Branch                SY 14

  //float nhit;                                                                       // Slice Branch
  float ncalhit;                                                                      // Slice Branch                 SY 15
  /// What is the different between the ncalhit and nmiphit
  float nmiphit;                                                                      // Slice Branch
  float ncontplanes;                                                                  // Slice Branch
  float starttime;                                                                    // Slice Branch
  float endtime;                                                                      // Slice Branch
  float meantime;                                                                     // Slice Branch
  float meanposX;                                                                     // Slice Branch
  float meanposY;                                                                     // Slice Branch
  float meanposZ;                                                                     // Slice Branch

  /// I do want to use our NUS Sandbox Vars
  float sumtx;                                                                        // NuS Sandbox
  float sumty;                                                                        // NuS Sandbox
  float ewsumtx;                                                                      // NuS Sandbox
  float ewsumty;                                                                      // NuS Sandbox
  float cossumtx;                                                                     // NuS Sandbox
  float cossumty;                                                                     // NuS Sandbox
  float cosewsumtx;                                                                   // NuS Sandbox
  float cosewsumty;                                                                   // NuS Sandbox
  float angsumtx;                                                                     // NuS Sandbox
  float angsumty;                                                                     // NuS Sandbox
  float angewsumtx;                                                                   // NuS Sandbox
  float angewsumty;                                                                   // NuS Sandbox

  float numupid;                                                                      // CosBDT Score by Kirk
  
  float rempid;                                                                       // RemID    /// We will try it as a selection CUT 
  float remlen;                                                                       // RemID

  float vtxTime;                                                                      // Vertex Branch
  float vtxX;                                                                         // Vertex Branch
  float vtxY;                                                                         // Vertex Branch
  float vtxZ;                                                                         // Vertex Branch
  float vtxpngNhit;                                                                   // Vertex Branch Prong Hits
  float vtxpngNhitX;                                                                  // Vertex Branch Prong Hits X View
  float vtxpngNhitY;                                                                  // Vertex Branch Prong Hits Y View

  float anglekal;                                                                     // Kirk 1
  float dirFY;                                                                        // Kirk 2
  float boxmaxFY;                                                                     // Kirk 3 
  float kalnhit;                                                                      // Kirk 4                           
  float kallen;                                                                       // Kirk 5
  float kalfwdcell;                                                                   // Kirk 6
  float kalbakcell;                                                                   // Kirk 6
  float scatt;                                                                        // Kirk 7
  float nhit;                                                                         // Kirk 8
  float energy;                                                                       // Kirk 9
  float boxminFY;                                                                     // Kirk 10                  
  float nkal;                                                                         // Kirk 11  

  // cosmic tree
  csTree->Branch("Run",&Run,"Run/I");                                                 // True Info
  csTree->Branch("SubRun",&SubRun,"SubRun/I");                                        // True Info
  csTree->Branch("Evt",&Evt,"Evt/I");                                                 // True Info
  csTree->Branch("SubEvt",&SubEvt,"SubEvt/I");                                        // True Info
  csTree->Branch("IsCC",&IsCC,"IsCC/I");                                              // True Info

  csTree->Branch("nueid",&nueid,"nueid/F");                                           // CVN
  csTree->Branch("numuid",&numuid,"numuid/F");                                        // CVN
  csTree->Branch("nutauid",&nutauid,"nutauid/F");                                     // CVN
  
  csTree->Branch("pngptp",&pngptp,"pngptp/F");                                        // NueCosRej
  csTree->Branch("starttop",&starttop,"starttop/F");                                  // NueCosRej
  csTree->Branch("startbottom",&startbottom,"startbottom/F");                         // NueCosRej
  csTree->Branch("startfront",&startfront,"startfront/F");                            // NueCosRej
  csTree->Branch("startback",&startback,"startback/F");                               // NueCosRej
  csTree->Branch("startwest",&startwest,"startwest/F");                               // NueCosRej
  csTree->Branch("starteast",&starteast,"starteast/F");                               // NueCosRej
  csTree->Branch("stoptop",&stoptop,"stoptop/F");                                     // NueCosRej
  csTree->Branch("stopbottom",&stopbottom,"stopbottom/F");                            // NueCosRej
  csTree->Branch("stopfront",&stopfront,"stopfront/F");                               // NueCosRej
  csTree->Branch("stopback",&stopback,"stopback/F");                                  // NueCosRej
  csTree->Branch("stopwest",&stopwest,"stopwest/F");                                  // NueCosRej
  csTree->Branch("stopeast",&stopeast,"stopeast/F");                                  // NueCosRej
  csTree->Branch("cosdang",&cosdang,"cosdang/F");                                     // NueCosRej
  csTree->Branch("vtxdoca",&vtxdoca,"vtxdoca/F");                                     // NueCosRej
  csTree->Branch("prongmaxx",&prongmaxx,"prongmaxx/F");                               // NueCosRej
  csTree->Branch("prongmaxy",&prongmaxy,"prongmaxy/F");                               // NueCosRej
  csTree->Branch("prongmaxz",&prongmaxz,"prongmaxz/F");                               // NueCosRej
  csTree->Branch("prongminx",&prongminx,"prongminx/F");                               // NueCosRej
  csTree->Branch("prongminy",&prongminy,"prongminy/F");                               // NueCosRej
  csTree->Branch("prongminz",&prongminz,"prongminz/F");                               // NueCosRej
  csTree->Branch("sparsenessasymm",&sparsenessasymm,"sparsenessasymm/F");             // NueCosRej
  csTree->Branch("hitsperplaneasymm",&hitsperplaneasymm,"hitsperplaneasymm/F");       // NueCosRej
  csTree->Branch("hitsperplane",&hitsperplane,"hitsperplane/F");                      // NueCosRej
  csTree->Branch("musliceidx",&musliceidx,"musliceidx/F");                            // NueCosRej
  csTree->Branch("muanglediff",&muanglediff,"muanglediff/F");                         // NueCosRej
  csTree->Branch("mutimediff",&mutimediff,"mutimediff/F");                            // NueCosRej
  csTree->Branch("muclosestapproach",&muclosestapproach,"muclosestapproach/F");       // NueCosRej

  csTree->Branch("shwnplane",&shwnplane,"shwnplane/F");                               // Shower Branch 
  csTree->Branch("shwmaxplanecont",&shwmaxplanecont,"shwmaxplanecont/F");             // Shower Branch
  csTree->Branch("shwmaxplanegap",&shwmaxplanegap,"shwmaxplanegap/F");                // Shower Branch
  csTree->Branch("shwstartX",&shwstartX,"shwstartX/F");                               // Shower Branch
  csTree->Branch("shwstartY",&shwstartY,"shwstartY/F");                               // Shower Branch
  csTree->Branch("shwstartZ",&shwstartZ,"shwstartZ/F");                               // Shower Branch
  csTree->Branch("shwdirX",&shwdirX,"shwdirX/F");                                     // Shower Branch
  csTree->Branch("shwdirZ",&shwdirZ,"shwdirZ/F");                                     // Shower Branch    
  csTree->Branch("shwNplaneX",&shwNplaneX,"shwNplaneX/F");                            // Shower Branch
  csTree->Branch("shwNplaneY",&shwNplaneY,"shwNplaneY/F");                            // Shower Branch    
  csTree->Branch("shwstopX",&shwstopX,"shwstopX/F");                                  // Shower Branch
  csTree->Branch("shwstopY",&shwstopY,"shwstopY/F");                                  // Shower Branch
  csTree->Branch("shwstopZ",&shwstopZ,"shwstopZ/F");                                  // Shower Branch

  csTree->Branch("sumtx",&sumtx,"sumtx/F");                                           // NuS Sandbox
  csTree->Branch("sumty",&sumty,"sumty/F");                                           // NuS Sandbox
  csTree->Branch("ewsumtx",&ewsumtx,"ewsumtx/F");                                     // NuS Sandbox
  csTree->Branch("ewsumty",&ewsumty,"ewsumty/F");                                     // NuS Sandbox
  csTree->Branch("cossumtx",&cossumtx,"cossumtx/F");                                  // NuS Sandbox
  csTree->Branch("cossumty",&cossumty,"cossumty/F");                                  // NuS Sandbox
  csTree->Branch("cosewsumtx",&cosewsumtx,"cosewsumtx/F");                            // NuS Sandbox
  csTree->Branch("cosewsumty",&cosewsumty,"cosewsumty/F");                            // NuS Sandbox
  csTree->Branch("angsumtx",&angsumtx,"angsumtx/F");                                  // NuS Sandbox
  csTree->Branch("angsumty",&angsumty,"angsumty/F");                                  // NuS Sandbox
  csTree->Branch("angewsumtx",&angewsumtx,"angewsumtx/F");                            // NuS Sandbox
  csTree->Branch("angewsumty",&angewsumty,"angewsumty/F");                            // NuS Sandbox  
  
  csTree->Branch("numupid",&numupid,"numupid/F");                                     // CosBDT score by Kirk

  csTree->Branch("rempid",&rempid,"rempid/F");                                        // RemID
  csTree->Branch("remlen",&remlen,"remlen/F");                                        // RemID 

  csTree->Branch("ncalhit",&ncalhit,"ncalhit/F");                                     // Slice Branch

  csTree->Branch("ncontplanes",&ncontplanes,"ncontplanes/F");                         // Slice Branch
  csTree->Branch("starttime",&starttime,"starttime/F");                               // Slice Branch
  csTree->Branch("endtime",&endtime,"endtime/F");                                     // Slice Branch
  csTree->Branch("meantime",&meantime,"meantime/F");                                  // Slice Branch
  csTree->Branch("meanposX",&meanposX,"meanposX/F");                                  // Slice Branch
  csTree->Branch("meanposY",&meanposY,"meanposY/F");                                  // Slice Branch
  csTree->Branch("meanposZ",&meanposZ,"meanposZ/F");                                  // Slice Branch

  csTree->Branch("vtxTime",&vtxTime,"vtxTime/F");                                     // Vertex Branch
  csTree->Branch("vtxX",&vtxX,"vtxX/F");                                              // Vertex Branch
  csTree->Branch("vtxY",&vtxY,"vtxY/F");                                              // Vertex Branch
  csTree->Branch("vtxZ",&vtxZ,"vtxZ/F");                                              // Vertex Branch
  csTree->Branch("vtxpngNhit",&vtxpngNhit,"vtxpngNhit/F");                            // Vertex Branch Prong Hits
  csTree->Branch("vtxpngNhitX",&vtxpngNhitX,"vtxpngNhitX/F");                         // Vertex Branch Prong Hits X View
  csTree->Branch("vtxpngNhitY",&vtxpngNhitY,"vtxpngNhitY/F");                         // Vertex Branch Prong Hits Y View

  csTree->Branch("ncid",&ncid,"ncid/F");                                              // CVN              SY 1
  csTree->Branch("cosmicid",&cosmicid,"cosmicid/F");                                  // CVN              SY 1
  csTree->Branch("partptp",&partptp,"partptp/F");                                     // NueCosRej        SY 2   LIDShower PtP 
  csTree->Branch("shwnhit",&shwnhit,"shwnhit/F");                                     // Shower Branch    SY 3
  csTree->Branch("shwnhitx",&shwnhitx,"shwnhitx/F");                                  // Shower Branch    SY 4
  csTree->Branch("shwnhity",&shwnhity,"shwnhity/F");                                  // Shower Branch    SY 5
  csTree->Branch("shwxminusy",&shwxminusy,"shwxminusy/F");                            // Shower Branch    SY 6
  csTree->Branch("shwxplusy",&shwxplusy,"shwxplusy/F");                               // Shower Branch    SY 7
  csTree->Branch("shwxovery",&shwxovery,"shwxovery/F");                               // Shower Branch    SY 8  
  csTree->Branch("shwcalE",&shwcalE,"shwcalE/F");                                     // Shower Branch    SY 9
  csTree->Branch("shwdirY",&shwdirY,"shwdirY/F");                                     // Shower Branch    SY 10
  csTree->Branch("shwlen",&shwlen,"shwlen/F");                                        // Shower Branch    SY 11
  csTree->Branch("shwwwidth",&shwwwidth,"shwwwidth/F");                               // Shower Branch    SY 12
  csTree->Branch("shwGap",&shwGap,"shwGap/F");                                        // Shower Branch    SY 13
  csTree->Branch("nshwlid",&nshwlid,"nshwlid/F");                                     // Shower Branch    SY 14
  csTree->Branch("nmiphit",&nmiphit,"nmiphit/F");                                     // Slice Branch     SY 15
  
  csTree->Branch("anglekal",&anglekal,"anglekal/F");                                  // Kirk 1
  csTree->Branch("dirFY",&dirFY,"dirFY/F");                                           // Kirk 2
  csTree->Branch("boxmaxFY",&boxmaxFY,"boxmaxFY/F");                                  // Kirk 3
  csTree->Branch("kalnhit",&kalnhit,"kalnhit/F");                                     // Kirk 4
  csTree->Branch("kallen",&kallen,"kallen/F");                                        // Kirk 5
  csTree->Branch("kalfwdcell",&kalfwdcell,"kalfwdcell/F");                            // Kirk 6
  csTree->Branch("kalbakcell",&kalbakcell,"kalbakcell/F");                            // Kirk 6
  csTree->Branch("scatt",&scatt,"scatt/F");                                           // Kirk 7
  csTree->Branch("nhit",&nhit,"nhit/F");                                              // Kirk 8
  csTree->Branch("energy",&energy,"energy/F");                                        // kirk 9
  csTree->Branch("boxminFY",&boxminFY,"boxminFY/F");                                  // Kirk 10
  csTree->Branch("nkal",&nkal,"nkal/F");                                              // Kirk 11

  // nc tree
  ncTree->Branch("Run",&Run,"Run/I");                                                 // Ture Info                    
  ncTree->Branch("SubRun",&SubRun,"SubRun/I");                                        // Ture Info                    
  ncTree->Branch("Evt",&Evt,"Evt/I");                                                 // Ture Info                    
  ncTree->Branch("SubEvt",&SubEvt,"SubEvt/I");                                        // Ture Info                    
  ncTree->Branch("IsCC",&IsCC,"IsCC/I");                                              // Ture Info                    

  ncTree->Branch("sumtx",&sumtx,"sumtx/F");                                           // NuS Sandbox
  ncTree->Branch("sumty",&sumty,"sumty/F");                                           // NuS Sandbox
  ncTree->Branch("ewsumtx",&ewsumtx,"ewsumtx/F");                                     // NuS Sandbox
  ncTree->Branch("ewsumty",&ewsumty,"ewsumty/F");                                     // NuS Sandbox
  ncTree->Branch("cossumtx",&cossumtx,"cossumtx/F");                                  // NuS Sandbox
  ncTree->Branch("cossumty",&cossumty,"cossumty/F");                                  // NuS Sandbox
  ncTree->Branch("cosewsumtx",&cosewsumtx,"cosewsumtx/F");                            // NuS Sandbox
  ncTree->Branch("cosewsumty",&cosewsumty,"cosewsumty/F");                            // NuS Sandbox
  ncTree->Branch("angsumtx",&angsumtx,"angsumtx/F");                                  // NuS Sandbox
  ncTree->Branch("angsumty",&angsumty,"angsumty/F");                                  // NuS Sandbox
  ncTree->Branch("angewsumtx",&angewsumtx,"angewsumtx/F");                            // NuS Sandbox
  ncTree->Branch("angewsumty",&angewsumty,"angewsumty/F");                            // NuS Sandbox  
  
  ncTree->Branch("starttop",&starttop,"starttop/F");                                  // NueCosRej
  ncTree->Branch("startbottom",&startbottom,"startbottom/F");                         // NueCosRej
  ncTree->Branch("startfront",&startfront,"startfront/F");                            // NueCosRej
  ncTree->Branch("startback",&startback,"startback/F");                               // NueCosRej
  ncTree->Branch("startwest",&startwest,"startwest/F");                               // NueCosRej
  ncTree->Branch("starteast",&starteast,"starteast/F");                               // NueCosRej
  ncTree->Branch("stoptop",&stoptop,"stoptop/F");                                     // NueCosRej
  ncTree->Branch("stopbottom",&stopbottom,"stopbottom/F");                            // NueCosRej
  ncTree->Branch("stopfront",&stopfront,"stopfront/F");                               // NueCosRej
  ncTree->Branch("stopback",&stopback,"stopback/F");                                  // NueCosRej
  ncTree->Branch("stopwest",&stopwest,"stopwest/F");                                  // NueCosRej
  ncTree->Branch("stopeast",&stopeast,"stopeast/F");                                  // NueCosRej
  ncTree->Branch("cosdang",&cosdang,"cosdang/F");                                     // NueCosRej
  ncTree->Branch("vtxdoca",&vtxdoca,"vtxdoca/F");                                     // NueCosRej
  ncTree->Branch("prongmaxx",&prongmaxx,"prongmaxx/F");                               // NueCosRej
  ncTree->Branch("prongmaxy",&prongmaxy,"prongmaxy/F");                               // NueCosRej
  ncTree->Branch("prongmaxz",&prongmaxz,"prongmaxz/F");                               // NueCosRej
  ncTree->Branch("prongminx",&prongminx,"prongminx/F");                               // NueCosRej
  ncTree->Branch("prongminy",&prongminy,"prongminy/F");                               // NueCosRej
  ncTree->Branch("prongminz",&prongminz,"prongminz/F");                               // NueCosRej
  ncTree->Branch("sparsenessasymm",&sparsenessasymm,"sparsenessasymm/F");             // NueCosRej
  ncTree->Branch("hitsperplaneasymm",&hitsperplaneasymm,"hitsperplaneasymm/F");       // NueCosRej
  ncTree->Branch("musliceidx",&musliceidx,"musliceidx/F");                            // NueCosRej
  ncTree->Branch("muanglediff",&muanglediff,"muanglediff/F");                         // NueCosRej
  ncTree->Branch("mutimediff",&mutimediff,"mutimediff/F");                            // NueCosRej
  ncTree->Branch("muclosestapproach",&muclosestapproach,"muclosestapproach/F");       // NueCosRej
  //ncTree->Branch("partptp",&partptp,"partptp/F");                                   // NueCosRej                         
  ncTree->Branch("pngptp",&pngptp,"pngptp/F");                                        // NueCosRej  
  
  ncTree->Branch("numupid",&numupid,"numupid/F");

  ncTree->Branch("nueid",&nueid,"nueid/F");                                           // CVN
  ncTree->Branch("numuid",&numuid,"numuid/F");                                        // CVN
  ncTree->Branch("nutauid",&nutauid,"nutauid/F");                                     // CVN

  ncTree->Branch("rempid",&rempid,"rempid/F");                                        // RemID
  ncTree->Branch("remlen",&remlen,"remlen/F");                                        // RemID 

  ncTree->Branch("ncalhit",&ncalhit,"ncalhit/F");                                     // Slice Branch
  //ncTree->Branch("nmiphit",&nmiphit,"nmiphit/F");                                   // Slice Branch
  ncTree->Branch("ncontplanes",&ncontplanes,"ncontplanes/F");                         // Slice Branch
  ncTree->Branch("starttime",&starttime,"starttime/F");                               // Slice Branch
  ncTree->Branch("endtime",&endtime,"endtime/F");                                     // Slice Branch
  ncTree->Branch("meantime",&meantime,"meantime/F");                                  // Slice Branch
  ncTree->Branch("meanposX",&meanposX,"meanposX/F");                                  // Slice Branch
  ncTree->Branch("meanposY",&meanposY,"meanposY/F");                                  // Slice Branch
  ncTree->Branch("meanposZ",&meanposZ,"meanposZ/F");                                  // Slice Branch

  //ncTree->Branch("shwnhit",&shwnhit,"shwnhit/F");                                   // Shower Branch
  //ncTree->Branch("shwnhitx",&shwnhitx,"shwnhitx/F");                                // Shower Branch
  //ncTree->Branch("shwnhity",&shwnhity,"shwnhity/F");                                // Shower Branch
  ncTree->Branch("shwnplane",&shwnplane,"shwnplane/F");                               // Shower Branch
  ncTree->Branch("shwmaxplanecont",&shwmaxplanecont,"shwmaxplanecont/F");             // Shower Branch
  ncTree->Branch("shwmaxplanegap",&shwmaxplanegap,"shwmaxplanegap/F");                // Shower Branch
  //ncTree->Branch("shwcalE",&shwcalE,"shwcalE/F");                                   // Shower Branch
  ncTree->Branch("shwstartX",&shwstartX,"shwstartX/F");                               // Shower Branch
  ncTree->Branch("shwstartY",&shwstartY,"shwstartY/F");                               // Shower Branch
  ncTree->Branch("shwstartZ",&shwstartZ,"shwstartZ/F");                               // Shower Branch
  ncTree->Branch("shwdirX",&shwdirX,"shwdirX/F");                                     // Shower Branch
  //ncTree->Branch("shwdirY",&shwdirY,"shwdirY/F");                                   // Shower Branch
  ncTree->Branch("shwdirZ",&shwdirZ,"shwdirZ/F");                                     // Shower Branch
  //ncTree->Branch("shwlen",&shwlen,"shwlen/F");                                      // Shower Branch
  //ncTree->Branch("shwwwidth",&shwwwidth,"shwwwidth/F");                             // Shower Branch
  ncTree->Branch("shwNplaneX",&shwNplaneX,"shwNplaneX/F");                            // Shower Branch
  ncTree->Branch("shwNplaneY",&shwNplaneY,"shwNplaneY/F");                            // Shower Branch
  //ncTree->Branch("shwGap",&shwGap,"shwGap/F");                                      // Shower Branch
  ncTree->Branch("shwstopX",&shwstopX,"shwstopX/F");                                  // Shower Branch
  ncTree->Branch("shwstopY",&shwstopY,"shwstopY/F");                                  // Shower Branch
  ncTree->Branch("shwstopZ",&shwstopZ,"shwstopZ/F");                                  // Shower Branch
  //ncTree->Branch("nshwlid",&nshwlid,"nshwlid/F");                                   // Shower Branch

  ncTree->Branch("vtxTime",&vtxTime,"vtxTime/F");                                     // Vertex Branch
  ncTree->Branch("vtxX",&vtxX,"vtxX/F");                                              // Vertex Branch
  ncTree->Branch("vtxY",&vtxY,"vtxY/F");                                              // Vertex Branch
  ncTree->Branch("vtxZ",&vtxZ,"vtxZ/F");                                              // Vertex Branch
  ncTree->Branch("vtxpngNhit",&vtxpngNhit,"vtxpngNhit/F");                            // Vertex Branch Prong Hits
  ncTree->Branch("vtxpngNhitX",&vtxpngNhitX,"vtxpngNhitX/F");                         // Vertex Branch Prong Hits X View
  ncTree->Branch("vtxpngNhitY",&vtxpngNhitY,"vtxpngNhitY/F");                         // Vertex Branch Prong Hits Y View

  ncTree->Branch("ncid",&ncid,"ncid/F");                                              // CVN              SY 1
  ncTree->Branch("cosmicid",&cosmicid,"cosmicid/F");                                  // CVN              SY 1
  ncTree->Branch("partptp",&partptp,"partptp/F");                                     // NueCosRej        SY 2  
  ncTree->Branch("shwnhit",&shwnhit,"shwnhit/F");                                     // Shower Branch    SY 3
  ncTree->Branch("shwnhitx",&shwnhitx,"shwnhitx/F");                                  // Shower Branch    SY 4
  ncTree->Branch("shwnhity",&shwnhity,"shwnhity/F");                                  // Shower Branch    SY 5
  ncTree->Branch("shwxminusy",&shwxminusy,"shwxminusy/F");                            // Shower Branch    SY 6
  ncTree->Branch("shwxplusy",&shwxplusy,"shwxplusy/F");                               // Shower Branch    SY 7
  ncTree->Branch("shwxovery",&shwxovery,"shwxovery/F");                               // Shower Branch    SY 8  
  ncTree->Branch("shwcalE",&shwcalE,"shwcalE/F");                                     // Shower Branch    SY 9
  ncTree->Branch("shwdirY",&shwdirY,"shwdirY/F");                                     // Shower Branch    SY 10
  ncTree->Branch("shwlen",&shwlen,"shwlen/F");                                        // Shower Branch    SY 11
  ncTree->Branch("shwwwidth",&shwwwidth,"shwwwidth/F");                               // Shower Branch    SY 12
  ncTree->Branch("shwGap",&shwGap,"shwGap/F");                                        // Shower Branch    SY 13
  ncTree->Branch("nshwlid",&nshwlid,"nshwlid/F");                                     // Shower Branch    SY 14
  ncTree->Branch("nmiphit",&nmiphit,"nmiphit/F");                                     // Slice Branch     SY 15

  ncTree->Branch("anglekal",&anglekal,"anglekal/F");                                  // Kirk 1
  ncTree->Branch("dirFY",&dirFY,"dirFY/F");                                           // Kirk 2
  ncTree->Branch("boxmaxFY",&boxmaxFY,"boxmaxFY/F");                                  // Kirk 3
  ncTree->Branch("kalnhit",&kalnhit,"kalnhit/F");                                     // Kirk 4
  ncTree->Branch("kallen",&kallen,"kallen/F");                                        // Kirk 5
  ncTree->Branch("kalfwdcell",&kalfwdcell,"kalfwdcell/F");                            // Kirk 6
  ncTree->Branch("kalbakcell",&kalbakcell,"kalbakcell/F");                            // Kirk 6
  ncTree->Branch("scatt",&scatt,"scatt/F");                                           // Kirk 7
  ncTree->Branch("nhit",&nhit,"nhit/F");                                              // Kirk 8
  ncTree->Branch("energy",&energy,"energy/F");                                        // Kirk 9
  ncTree->Branch("boxminFY",&boxminFY,"boxminFY/F");                                  // Kirk 10
  ncTree->Branch("nkal",&nkal,"nkal/F");                                              // Kirk 11

  // numu tree
  muTree->Branch("Run",&Run,"Run/I");                                                 // Ture Info                    
  muTree->Branch("SubRun",&SubRun,"SubRun/I");                                        // Ture Info                    
  muTree->Branch("Evt",&Evt,"Evt/I");                                                 // Ture Info                    
  muTree->Branch("SubEvt",&SubEvt,"SubEvt/I");                                        // Ture Info                    
  muTree->Branch("IsCC",&IsCC,"IsCC/I");                                              // Ture Info                    

  muTree->Branch("sumtx",&sumtx,"sumtx/F");                                           // NuS Sandbox
  muTree->Branch("sumty",&sumty,"sumty/F");                                           // NuS Sandbox
  muTree->Branch("ewsumtx",&ewsumtx,"ewsumtx/F");                                     // NuS Sandbox
  muTree->Branch("ewsumty",&ewsumty,"ewsumty/F");                                     // NuS Sandbox
  muTree->Branch("cossumtx",&cossumtx,"cossumtx/F");                                  // NuS Sandbox
  muTree->Branch("cossumty",&cossumty,"cossumty/F");                                  // NuS Sandbox
  muTree->Branch("cosewsumtx",&cosewsumtx,"cosewsumtx/F");                            // NuS Sandbox
  muTree->Branch("cosewsumty",&cosewsumty,"cosewsumty/F");                            // NuS Sandbox
  muTree->Branch("angsumtx",&angsumtx,"angsumtx/F");                                  // NuS Sandbox
  muTree->Branch("angsumty",&angsumty,"angsumty/F");                                  // NuS Sandbox
  muTree->Branch("angewsumtx",&angewsumtx,"angewsumtx/F");                            // NuS Sandbox
  muTree->Branch("angewsumty",&angewsumty,"angewsumty/F");                            // NuS Sandbox  
  
  muTree->Branch("starttop",&starttop,"starttop/F");                                  // NueCosRej
  muTree->Branch("startbottom",&startbottom,"startbottom/F");                         // NueCosRej
  muTree->Branch("startfront",&startfront,"startfront/F");                            // NueCosRej
  muTree->Branch("startback",&startback,"startback/F");                               // NueCosRej
  muTree->Branch("startwest",&startwest,"startwest/F");                               // NueCosRej
  muTree->Branch("starteast",&starteast,"starteast/F");                               // NueCosRej
  muTree->Branch("stoptop",&stoptop,"stoptop/F");                                     // NueCosRej
  muTree->Branch("stopbottom",&stopbottom,"stopbottom/F");                            // NueCosRej
  muTree->Branch("stopfront",&stopfront,"stopfront/F");                               // NueCosRej
  muTree->Branch("stopback",&stopback,"stopback/F");                                  // NueCosRej
  muTree->Branch("stopwest",&stopwest,"stopwest/F");                                  // NueCosRej
  muTree->Branch("stopeast",&stopeast,"stopeast/F");                                  // NueCosRej
  muTree->Branch("cosdang",&cosdang,"cosdang/F");                                     // NueCosRej
  muTree->Branch("vtxdoca",&vtxdoca,"vtxdoca/F");                                     // NueCosRej
  muTree->Branch("prongmaxx",&prongmaxx,"prongmaxx/F");                               // NueCosRej
  muTree->Branch("prongmaxy",&prongmaxy,"prongmaxy/F");                               // NueCosRej
  muTree->Branch("prongmaxz",&prongmaxz,"prongmaxz/F");                               // NueCosRej
  muTree->Branch("prongminx",&prongminx,"prongminx/F");                               // NueCosRej
  muTree->Branch("prongminy",&prongminy,"prongminy/F");                               // NueCosRej
  muTree->Branch("prongminz",&prongminz,"prongminz/F");                               // NueCosRej
  muTree->Branch("sparsenessasymm",&sparsenessasymm,"sparsenessasymm/F");             // NueCosRej
  muTree->Branch("hitsperplaneasymm",&hitsperplaneasymm,"hitsperplaneasymm/F");       // NueCosRej
  muTree->Branch("hitsperplane",&hitsperplane,"hitsperplane/F");                      // NueCosRej
  muTree->Branch("musliceidx",&musliceidx,"musliceidx/F");                            // NueCosRej
  muTree->Branch("muanglediff",&muanglediff,"muanglediff/F");                         // NueCosRej
  muTree->Branch("mutimediff",&mutimediff,"mutimediff/F");                            // NueCosRej
  muTree->Branch("muclosestapproach",&muclosestapproach,"muclosestapproach/F");       // NueCosRej
  //muTree->Branch("partptp",&partptp,"partptp/F");                                   // NueCosRej        
  muTree->Branch("pngptp",&pngptp,"pngptp/F");                                        // NueCosRej
                                                    

  muTree->Branch("nueid",&nueid,"nueid/F");                                           // CVN
  muTree->Branch("numuid",&numuid,"numuid/F");                                        // CVN
  muTree->Branch("nutauid",&nutauid,"nutauid/F");                                     // CVN

  muTree->Branch("numupid",&numupid,"numupid/F");

  muTree->Branch("rempid",&rempid,"rempid/F");                                        // RemID
  muTree->Branch("remlen",&remlen,"remlen/F");                                        // RemID 

  muTree->Branch("ncalhit",&ncalhit,"ncalhit/F");                                     // Slice Branch
  //muTree->Branch("nmiphit",&nmiphit,"nmiphit/F");                                   // Slice Branch
  muTree->Branch("ncontplanes",&ncontplanes,"ncontplanes/F");                         // Slice Branch
  muTree->Branch("starttime",&starttime,"starttime/F");                               // Slice Branch
  muTree->Branch("endtime",&endtime,"endtime/F");                                     // Slice Branch
  muTree->Branch("meantime",&meantime,"meantime/F");                                  // Slice Branch
  muTree->Branch("meanposX",&meanposX,"meanposX/F");                                  // Slice Branch
  muTree->Branch("meanposY",&meanposY,"meanposY/F");                                  // Slice Branch
  muTree->Branch("meanposZ",&meanposZ,"meanposZ/F");                                  // Slice Branch

  //muTree->Branch("shwnhit",&shwnhit,"shwnhit/F");                                   // Shower Branch
  //muTree->Branch("shwnhitx",&shwnhitx,"shwnhitx/F");                                // Shower Branch
  //muTree->Branch("shwnhity",&shwnhity,"shwnhity/F");                                // Shower Branch
  muTree->Branch("shwnplane",&shwnplane,"shwnplane/F");                               // Shower Branch
  muTree->Branch("shwmaxplanecont",&shwmaxplanecont,"shwmaxplanecont/F");             // Shower Branch
  muTree->Branch("shwmaxplanegap",&shwmaxplanegap,"shwmaxplanegap/F");                // Shower Branch
  //muTree->Branch("shwcalE",&shwcalE,"shwcalE/F");                                   // Shower Branch
  muTree->Branch("shwstartX",&shwstartX,"shwstartX/F");                               // Shower Branch
  muTree->Branch("shwstartY",&shwstartY,"shwstartY/F");                               // Shower Branch
  muTree->Branch("shwstartZ",&shwstartZ,"shwstartZ/F");                               // Shower Branch
  muTree->Branch("shwdirX",&shwdirX,"shwdirX/F");                                     // Shower Branch
  //muTree->Branch("shwdirY",&shwdirY,"shwdirY/F");                                   // Shower Branch
  muTree->Branch("shwdirZ",&shwdirZ,"shwdirZ/F");                                     // Shower Branch
  //muTree->Branch("shwlen",&shwlen,"shwlen/F");                                      // Shower Branch
  //muTree->Branch("shwwwidth",&shwwwidth,"shwwwidth/F");                             // Shower Branch
  muTree->Branch("shwNplaneX",&shwNplaneX,"shwNplaneX/F");                            // Shower Branch
  muTree->Branch("shwNplaneY",&shwNplaneY,"shwNplaneY/F");                            // Shower Branch
  //muTree->Branch("shwGap",&shwGap,"shwGap/F");                                      // Shower Branch
  muTree->Branch("shwstopX",&shwstopX,"shwstopX/F");                                  // Shower Branch
  muTree->Branch("shwstopY",&shwstopY,"shwstopY/F");                                  // Shower Branch
  muTree->Branch("shwstopZ",&shwstopZ,"shwstopZ/F");                                  // Shower Branch
  //muTree->Branch("nshwlid",&nshwlid,"nshwlid/F");                                   // Shower Branch

  muTree->Branch("vtxTime",&vtxTime,"vtxTime/F");                                     // Vertex Branch
  muTree->Branch("vtxX",&vtxX,"vtxX/F");                                              // Vertex Branch
  muTree->Branch("vtxY",&vtxY,"vtxY/F");                                              // Vertex Branch
  muTree->Branch("vtxZ",&vtxZ,"vtxZ/F");                                              // Vertex Branch
  muTree->Branch("vtxpngNhit",&vtxpngNhit,"vtxpngNhit/F");                            // Vertex Branch Prong Hits
  muTree->Branch("vtxpngNhitX",&vtxpngNhitX,"vtxpngNhitX/F");                         // Vertex Branch Prong Hits X View
  muTree->Branch("vtxpngNhitY",&vtxpngNhitY,"vtxpngNhitY/F");                         // Vertex Branch Prong Hits Y View

  muTree->Branch("ncid",&ncid,"ncid/F");                                              // CVN              SY 1
  muTree->Branch("cosmicid",&cosmicid,"cosmicid/F");                                  // CVN              SY 1
  muTree->Branch("partptp",&partptp,"partptp/F");                                     // NueCosRej        SY 2  
  muTree->Branch("shwnhit",&shwnhit,"shwnhit/F");                                     // Shower Branch    SY 3
  muTree->Branch("shwnhitx",&shwnhitx,"shwnhitx/F");                                  // Shower Branch    SY 4
  muTree->Branch("shwnhity",&shwnhity,"shwnhity/F");                                  // Shower Branch    SY 5
  muTree->Branch("shwxminusy",&shwxminusy,"shwxminusy/F");                            // Shower Branch    SY 6
  muTree->Branch("shwxplusy",&shwxplusy,"shwxplusy/F");                               // Shower Branch    SY 7
  muTree->Branch("shwxovery",&shwxovery,"shwxovery/F");                               // Shower Branch    SY 8  
  muTree->Branch("shwcalE",&shwcalE,"shwcalE/F");                                     // Shower Branch    SY 9
  muTree->Branch("shwdirY",&shwdirY,"shwdirY/F");                                     // Shower Branch    SY 10
  muTree->Branch("shwlen",&shwlen,"shwlen/F");                                        // Shower Branch    SY 11
  muTree->Branch("shwwwidth",&shwwwidth,"shwwwidth/F");                               // Shower Branch    SY 12
  ncTree->Branch("shwGap",&shwGap,"shwGap/F");                                        // Shower Branch    SY 13
  muTree->Branch("nshwlid",&nshwlid,"nshwlid/F");                                     // Shower Branch    SY 14
  muTree->Branch("nmiphit",&nmiphit,"nmiphit/F");                                     // Slice Branch     SY 15
  
  muTree->Branch("anglekal",&anglekal,"anglekal/F");                                  // Kirk 1
  muTree->Branch("dirFY",&dirFY,"dirFY/F");                                           // Kirk 2
  muTree->Branch("boxmaxFY",&boxmaxFY,"boxmaxFY/F");                                  // Krik 3
  muTree->Branch("kalnhit",&kalnhit,"kalnhit/F");                                     // Kirk 4
  muTree->Branch("kallen",&kallen,"kallen/F");                                        // Kirk 5
  muTree->Branch("kalfwdcell",&kalfwdcell,"kalfwdcell/F");                            // Kirk 6
  muTree->Branch("kalbakcell",&kalbakcell,"kalbakcell/F");                            // Kirk 6
  muTree->Branch("scatt",&scatt,"scatt/F");                                           // Kirk 7
  muTree->Branch("nhit",&nhit,"nhit/F");                                              // Kirk 8
  muTree->Branch("energy",&energy,"energy/F");                                        // Kirk 9
  muTree->Branch("boxminFY",&boxminFY,"boxminFY/F");                                  // Kirk 10
  muTree->Branch("nkal",&nkal,"nkal/F");                                              // Krik 11

  // nue tree
  neTree->Branch("Run",&Run,"Run/I");                                                 // Ture Info                    
  neTree->Branch("SubRun",&SubRun,"SubRun/I");                                        // Ture Info                    
  neTree->Branch("Evt",&Evt,"Evt/I");                                                 // Ture Info                    
  neTree->Branch("SubEvt",&SubEvt,"SubEvt/I");                                        // Ture Info                    
  neTree->Branch("IsCC",&IsCC,"IsCC/I");                                              // Ture Info                    

  neTree->Branch("sumtx",&sumtx,"sumtx/F");                                           // NuS Sandbox
  neTree->Branch("sumty",&sumty,"sumty/F");                                           // NuS Sandbox
  neTree->Branch("ewsumtx",&ewsumtx,"ewsumtx/F");                                     // NuS Sandbox
  neTree->Branch("ewsumty",&ewsumty,"ewsumty/F");                                     // NuS Sandbox
  neTree->Branch("cossumtx",&cossumtx,"cossumtx/F");                                  // NuS Sandbox
  neTree->Branch("cossumty",&cossumty,"cossumty/F");                                  // NuS Sandbox
  neTree->Branch("cosewsumtx",&cosewsumtx,"cosewsumtx/F");                            // NuS Sandbox
  neTree->Branch("cosewsumty",&cosewsumty,"cosewsumty/F");                            // NuS Sandbox
  neTree->Branch("angsumtx",&angsumtx,"angsumtx/F");                                  // NuS Sandbox
  neTree->Branch("angsumty",&angsumty,"angsumty/F");                                  // NuS Sandbox
  neTree->Branch("angewsumtx",&angewsumtx,"angewsumtx/F");                            // NuS Sandbox
  neTree->Branch("angewsumty",&angewsumty,"angewsumty/F");                            // NuS Sandbox  
  
  neTree->Branch("starttop",&starttop,"starttop/F");                                  // NueCosRej
  neTree->Branch("startbottom",&startbottom,"startbottom/F");                         // NueCosRej
  neTree->Branch("startfront",&startfront,"startfront/F");                            // NueCosRej
  neTree->Branch("startback",&startback,"startback/F");                               // NueCosRej
  neTree->Branch("startwest",&startwest,"startwest/F");                               // NueCosRej
  neTree->Branch("starteast",&starteast,"starteast/F");                               // NueCosRej
  neTree->Branch("stoptop",&stoptop,"stoptop/F");                                     // NueCosRej
  neTree->Branch("stopbottom",&stopbottom,"stopbottom/F");                            // NueCosRej
  neTree->Branch("stopfront",&stopfront,"stopfront/F");                               // NueCosRej
  neTree->Branch("stopback",&stopback,"stopback/F");                                  // NueCosRej
  neTree->Branch("stopwest",&stopwest,"stopwest/F");                                  // NueCosRej
  neTree->Branch("stopeast",&stopeast,"stopeast/F");                                  // NueCosRej
  neTree->Branch("cosdang",&cosdang,"cosdang/F");                                     // NueCosRej
  neTree->Branch("vtxdoca",&vtxdoca,"vtxdoca/F");                                     // NueCosRej
  neTree->Branch("prongmaxx",&prongmaxx,"prongmaxx/F");                               // NueCosRej
  neTree->Branch("prongmaxy",&prongmaxy,"prongmaxy/F");                               // NueCosRej
  neTree->Branch("prongmaxz",&prongmaxz,"prongmaxz/F");                               // NueCosRej
  neTree->Branch("prongminx",&prongminx,"prongminx/F");                               // NueCosRej
  neTree->Branch("prongminy",&prongminy,"prongminy/F");                               // NueCosRej
  neTree->Branch("prongminz",&prongminz,"prongminz/F");                               // NueCosRej
  neTree->Branch("sparsenessasymm",&sparsenessasymm,"sparsenessasymm/F");             // NueCosRej
  neTree->Branch("hitsperplaneasymm",&hitsperplaneasymm,"hitsperplaneasymm/F");       // NueCosRej
  neTree->Branch("musliceidx",&musliceidx,"musliceidx/F");                            // NueCosRej
  neTree->Branch("muanglediff",&muanglediff,"muanglediff/F");                         // NueCosRej
  neTree->Branch("mutimediff",&mutimediff,"mutimediff/F");                            // NueCosRej
  neTree->Branch("muclosestapproach",&muclosestapproach,"muclosestapproach/F");       // NueCosRej
  //neTree->Branch("partptp",&partptp,"partptp/F");                                   // NueCosRej         
  neTree->Branch("pngptp",&pngptp,"pngptp/F");                                        // NueCosRej                                                    


  neTree->Branch("nueid",&nueid,"nueid/F");                                           // CVN
  neTree->Branch("numuid",&numuid,"numuid/F");                                        // CVN
  neTree->Branch("nutauid",&nutauid,"nutauid/F");                                     // CVN

  neTree->Branch("numupid",&numupid,"numupid/F");

  neTree->Branch("rempid",&rempid,"rempid/F");                                        // RemID
  neTree->Branch("remlen",&remlen,"remlen/F");                                        // RemID 

  neTree->Branch("ncalhit",&ncalhit,"ncalhit/F");                                     // Slice Branch
  //neTree->Branch("nmiphit",&nmiphit,"nmiphit/F");                                     // Slice Branch
  neTree->Branch("ncontplanes",&ncontplanes,"ncontplanes/F");                         // Slice Branch
  neTree->Branch("starttime",&starttime,"starttime/F");                               // Slice Branch
  neTree->Branch("endtime",&endtime,"endtime/F");                                     // Slice Branch
  neTree->Branch("meantime",&meantime,"meantime/F");                                  // Slice Branch
  neTree->Branch("meanposX",&meanposX,"meanposX/F");                                  // Slice Branch
  neTree->Branch("meanposY",&meanposY,"meanposY/F");                                  // Slice Branch
  neTree->Branch("meanposZ",&meanposZ,"meanposZ/F");                                  // Slice Branch

  //neTree->Branch("shwnhit",&shwnhit,"shwnhit/F");                                   // Shower Branch
  //neTree->Branch("shwnhitx",&shwnhitx,"shwnhitx/F");                                // Shower Branch
  //neTree->Branch("shwnhity",&shwnhity,"shwnhity/F");                                // Shower Branch
  neTree->Branch("shwnplane",&shwnplane,"shwnplane/F");                               // Shower Branch
  neTree->Branch("shwmaxplanecont",&shwmaxplanecont,"shwmaxplanecont/F");             // Shower Branch
  neTree->Branch("shwmaxplanegap",&shwmaxplanegap,"shwmaxplanegap/F");                // Shower Branch
  //neTree->Branch("shwcalE",&shwcalE,"shwcalE/F");                                   // Shower Branch
  neTree->Branch("shwstartX",&shwstartX,"shwstartX/F");                               // Shower Branch
  neTree->Branch("shwstartY",&shwstartY,"shwstartY/F");                               // Shower Branch
  neTree->Branch("shwstartZ",&shwstartZ,"shwstartZ/F");                               // Shower Branch
  neTree->Branch("shwdirX",&shwdirX,"shwdirX/F");                                     // Shower Branch
  //neTree->Branch("shwdirY",&shwdirY,"shwdirY/F");                                   // Shower Branch
  neTree->Branch("shwdirZ",&shwdirZ,"shwdirZ/F");                                     // Shower Branch
  //neTree->Branch("shwlen",&shwlen,"shwlen/F");                                      // Shower Branch
  //neTree->Branch("shwwwidth",&shwwwidth,"shwwwidth/F");                             // Shower Branch
  neTree->Branch("shwNplaneX",&shwNplaneX,"shwNplaneX/F");                            // Shower Branch
  neTree->Branch("shwNplaneY",&shwNplaneY,"shwNplaneY/F");                            // Shower Branch
  //neTree->Branch("shwGap",&shwGap,"shwGap/F");                                      // Shower Branch
  neTree->Branch("shwstopX",&shwstopX,"shwstopX/F");                                  // Shower Branch
  neTree->Branch("shwstopY",&shwstopY,"shwstopY/F");                                  // Shower Branch
  neTree->Branch("shwstopZ",&shwstopZ,"shwstopZ/F");                                  // Shower Branch
  //neTree->Branch("nshwlid",&nshwlid,"nshwlid/F");                                   // Shower Branch

  neTree->Branch("vtxTime",&vtxTime,"vtxTime/F");                                     // Vertex Branch
  neTree->Branch("vtxX",&vtxX,"vtxX/F");                                              // Vertex Branch
  neTree->Branch("vtxY",&vtxY,"vtxY/F");                                              // Vertex Branch
  neTree->Branch("vtxZ",&vtxZ,"vtxZ/F");                                              // Vertex Branch
  neTree->Branch("vtxpngNhit",&vtxpngNhit,"vtxpngNhit/F");                            // Vertex Branch Prong Hits
  neTree->Branch("vtxpngNhitX",&vtxpngNhitX,"vtxpngNhitX/F");                         // Vertex Branch Prong Hits X View
  neTree->Branch("vtxpngNhitY",&vtxpngNhitY,"vtxpngNhitY/F");                         // Vertex Branch Prong Hits Y View

  neTree->Branch("ncid",&ncid,"ncid/F");                                              // CVN              SY 1
  neTree->Branch("cosmicid",&cosmicid,"cosmicid/F");                                  // CVN              SY 1
  neTree->Branch("partptp",&partptp,"partptp/F");                                     // NueCosRej        SY 2  
  neTree->Branch("shwnhit",&shwnhit,"shwnhit/F");                                     // Shower Branch    SY 3
  neTree->Branch("shwnhitx",&shwnhitx,"shwnhitx/F");                                  // Shower Branch    SY 4
  neTree->Branch("shwnhity",&shwnhity,"shwnhity/F");                                  // Shower Branch    SY 5
  neTree->Branch("shwxminusy",&shwxminusy,"shwxminusy/F");                            // Shower Branch    SY 6
  neTree->Branch("shwxplusy",&shwxplusy,"shwxplusy/F");                               // Shower Branch    SY 7
  neTree->Branch("shwxovery",&shwxovery,"shwxovery/F");                               // Shower Branch    SY 8  
  neTree->Branch("shwcalE",&shwcalE,"shwcalE/F");                                     // Shower Branch    SY 9
  neTree->Branch("shwdirY",&shwdirY,"shwdirY/F");                                     // Shower Branch    SY 10
  neTree->Branch("shwlen",&shwlen,"shwlen/F");                                        // Shower Branch    SY 11
  neTree->Branch("shwwwidth",&shwwwidth,"shwwwidth/F");                               // Shower Branch    SY 12
  neTree->Branch("shwGap",&shwGap,"shwGap/F");                                        // Shower Branch    SY 13
  neTree->Branch("nshwlid",&nshwlid,"nshwlid/F");                                     // Shower Branch    SY 14
  neTree->Branch("nmiphit",&nmiphit,"nmiphit/F");                                     // Slice Branch     SY 15
  
  neTree->Branch("anglekal",&anglekal,"anglekal/F");                                  // Kirk 1
  neTree->Branch("dirFY",&dirFY,"dirFY/F");                                           // Kirk 2
  neTree->Branch("boxmaxFY",&boxmaxFY,"boxmaxFY/F");                                  // Kirk 3
  neTree->Branch("kalnhit",&kalnhit,"kalnhit/F");                                     // Kirk 4
  neTree->Branch("kallen",&kallen,"kallen/F");                                        // Kirk 5
  neTree->Branch("kalfwdcell",&kalfwdcell,"kalfwdcell/F");                            // Kirk 6
  neTree->Branch("kalbakcell",&kalbakcell,"kalbakcell/F");                            // Kirk 6
  neTree->Branch("scatt",&scatt,"scatt/F");                                           // Kirk 7
  neTree->Branch("nhit",&nhit,"nhit/F");                                              // Kirk 8
  neTree->Branch("energy",&energy,"energy/F");                                        // Kirk 9
  neTree->Branch("boxminFY",&boxminFY,"boxminFY/F");                                  // Kirk 10
  neTree->Branch("nkal",&nkal,"nkal/F");                                              // Kirk 11

  //recTreecos->SetBranchAddress("rec", &recTreeObject);
  //recTreenon->SetBranchAddress("rec", &recTreeObject);

  IFileSource* nonswap=loadernc.WildcardOrSAMQuery(fnamenc);
  int AllNon = nonswap->NFiles();
  //std::cout<<"NUMBER OF FILES: "<<Nfiles<<std::endl;
  
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
    // Turn on the useful branches from the beam file
    recTreenon->SetBranchStatus("hdr.run",1);                                      // Ture Info 
    recTreenon->SetBranchStatus("hdr.subrun",1);                                   // Ture Info 
    recTreenon->SetBranchStatus("hdr.evt",1);                                      // Ture Info 
    recTreenon->SetBranchStatus("hdr.subevt",1);                                   // Ture Info 

    recTreenon->SetBranchStatus("mc.nnu",1);                                       // Ture Info 
    recTreenon->SetBranchStatus("mc.nu",1);                                        // Ture Info 
    recTreenon->SetBranchStatus("mc.nu.iscc",1);                                   // Ture Info 
    recTreenon->SetBranchStatus("mc.nu.pdg",1);                                    // Ture Info 
    recTreenon->SetBranchStatus("trk.ncosmic",1);                                  // Ture Info 
    
    recTreenon->SetBranchStatus("sel.nuecosrej.starttop",1);                        // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.startbottom",1);                     // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.startfront",1);                      // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.startback",1);                       // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.startwest",1);                       // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.starteast",1);                       // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.stoptop",1);                         // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.stopbottom",1);                      // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.stopfront",1);                       // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.stopback",1);                        // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.stopwest",1);                        // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.stopeast",1);                        // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.cosdang",1);                         // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.vtxdoca",1);                         // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.prongmaxx",1);                       // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.prongmaxy",1);                       // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.prongmaxz",1);                       // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.prongminx",1);                       // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.prongminy",1);                       // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.prongminz",1);                       // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.sparsenessasymm",1);                 // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.hitsperplane",1);                    // NueCosRej  /////// 
    recTreenon->SetBranchStatus("sel.nuecosrej.hitsperplaneasymm",1);               // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.musliceidx",1);                      // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.muanglediff",1);                     // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.mutimediff",1);                      // NueCosRej 
    recTreenon->SetBranchStatus("sel.nuecosrej.muclosestapproach",1);               // NueCosRej

    //recTreenon->SetBranchStatus("sel.nuecosrej.pngptp",1);                        // NueCosRej              Shaokai 1
    //recTreenon->SetBranchStatus("shw.shwlid.nhit",1);                             // Shower Branch          Shaokai 2
    //recTreenon->SetBranchStatus("shw.shwlid.nhitx",1);                            // Shower Branch          Shaokai 3
    //recTreenon->SetBranchStatus("shw.shwlid.nhity",1);                            // Shower Branch          Shaokai 3
    recTreenon->SetBranchStatus("shw.shwlid.nplane",1);                             // Shower Branch          Shaokai 4
    recTreenon->SetBranchStatus("shw.shwlid.maxplanecont",1);                       // Shower Branch          Shaokai 5
    recTreenon->SetBranchStatus("shw.shwlid.maxplanegap",1);                        // Shower Branch          Shaokai 6
    
    recTreenon->SetBranchStatus("shw.shwlid.start.fX",1);                           // Shower Branch
    recTreenon->SetBranchStatus("shw.shwlid.start.fY",1);                           // Shower Branch  
    recTreenon->SetBranchStatus("shw.shwlid.start.fZ",1);                           // Shower Branch
    recTreenon->SetBranchStatus("shw.shwlid.dir.fX",1);                             // Shower Branch
    recTreenon->SetBranchStatus("shw.shwlid.dir.fZ",1);                             // Shower Branch
    recTreenon->SetBranchStatus("shw.shwlid.nplanex",1);                            // Shower Branch
    recTreenon->SetBranchStatus("shw.shwlid.nplaney",1);                            // Shower Branch

    recTreenon->SetBranchStatus("shw.shwlid.stop.fX",1);                            // Shower Branch
    recTreenon->SetBranchStatus("shw.shwlid.stop.fY",1);                            // Shower Branch
    recTreenon->SetBranchStatus("shw.shwlid.stop.fZ",1);                            // Shower Branch

    //recTreenon->SetBranchStatus("shw.shwlid.gap",1);                              // Shower Branch

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    recTreenon->SetBranchStatus("sel.elecid.shwlid.ismuon",1);                        //sel.elecid.shwlid[0].ismuon
    recTreenon->SetBranchStatus("sel.cosrej.numucontpid",1);
    //sr->sel.cosrej.numucontpid    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    recTreenon->SetBranchStatus("sel.cvn.nueid",1);                                 // CVN
    recTreenon->SetBranchStatus("sel.cvn.numuid",1);                                // CVN
    recTreenon->SetBranchStatus("sel.cvn.nutauid",1);                               // CVN

    recTreenon->SetBranchStatus("slc.nhit",1);                                      // Slice Branch
    //recTreenon->SetBranchStatus("slc.calhit",1);                                  // Slice Branch

    recTreenon->SetBranchStatus("slc.ncontplanes",1);                               // Slice Branch
    recTreenon->SetBranchStatus("slc.ncellsfromedge",1);                            // Slice Branch
    recTreenon->SetBranchStatus("slc.starttime",1);                                 // Slice Branch
    recTreenon->SetBranchStatus("slc.endtime",1);                                   // Slice Branch
    recTreenon->SetBranchStatus("slc.meantime",1);                                  // Slice Branch
    recTreenon->SetBranchStatus("slc.meanpos.fX",1);                                // Slice Branch
    recTreenon->SetBranchStatus("slc.meanpos.fY",1);                                // Slice Branch
    recTreenon->SetBranchStatus("slc.meanpos.fZ",1);                                // Slice Branch

    recTreenon->SetBranchStatus("sand.nus.sumtx",1);                                // NuS Sandbox
    recTreenon->SetBranchStatus("sand.nus.sumty",1);                                // NuS Sandbox
    recTreenon->SetBranchStatus("sand.nus.ewsumtx",1);                              // NuS Sandbox
    recTreenon->SetBranchStatus("sand.nus.ewsumty",1);                              // NuS Sandbox
    recTreenon->SetBranchStatus("sand.nus.cossumtx",1);                             // NuS Sandbox
    recTreenon->SetBranchStatus("sand.nus.cossumty",1);                             // NuS Sandbox
    recTreenon->SetBranchStatus("sand.nus.cosewsumtx",1);                           // NuS Sandbox
    recTreenon->SetBranchStatus("sand.nus.cosewsumty",1);                           // NuS Sandbox
    recTreenon->SetBranchStatus("sand.nus.angsumtx",1);                             // NuS Sandbox
    recTreenon->SetBranchStatus("sand.nus.angsumty",1);                             // NuS Sandbox
    recTreenon->SetBranchStatus("sand.nus.angewsumtx",1);                           // NuS Sandbox
    recTreenon->SetBranchStatus("sand.nus.angewsumty",1);                           // NuS Sandbox    
    
    //float shwGap;                                                                 // Shower Branch                SY 13                           
    recTreenon->SetBranchStatus("sel.remid.pid",1);                                 // RemID
    recTreenon->SetBranchStatus("sel.remid.len",1);                                 // RemID

    recTreenon->SetBranchStatus("vtx.elastic.time",1);                              // Vertex Branch
    recTreenon->SetBranchStatus("vtx.elastic.vtx.fX",1);                            // Vertex Branch
    recTreenon->SetBranchStatus("vtx.elastic.vtx.fY",1);                            // Vertex Branch
    recTreenon->SetBranchStatus("vtx.elastic.vtx.fZ",1);                            // Vertex Branch
    //recTreenon->SetBranchStatus("vtx.elastic.png3d.nhit",1);                      // Vertex Branch
    //recTreenon->SetBranchStatus("vtx.elastic.png3d.nhitx",1);                     // Vertex Branch
    //recTreenon->SetBranchStatus("vtx.elastic.png3d.nhity",1);                     // Vertex Branch
    recTreenon->SetBranchStatus("vtx.nelastic",1);                                  // Vertex Branch

    recTreenon->SetBranchStatus("sel.cvn.ncid",1);                                  // CVN                    SY 1
    //recTreenon->SetBranchStatus("sel.cvn.output[14]",1);                          // CVN                    SY 1 
    recTreenon->SetBranchStatus("sel.cvn.output",1);                                // CVN                    SY 1 
    recTreenon->SetBranchStatus("sel.nuecosrej.partptp",1);                         // NueCosRej              SY 2            
    recTreenon->SetBranchStatus("shw.shwlid.nhit",1);                               // Shower Branch          SY 3
    recTreenon->SetBranchStatus("shw.shwlid.nhitx",1);                              // Shower Branch          SY 4
    recTreenon->SetBranchStatus("shw.shwlid.nhity",1);                              // Shower Branch          SY 5
    recTreenon->SetBranchStatus("shw.shwlid.calE",1);                               // Shower Branch          SY 6            
    recTreenon->SetBranchStatus("shw.shwlid.dir.fY",1);                             // Shower Branch          SY 7
    recTreenon->SetBranchStatus("shw.shwlid.len",1);                                // Shower Branch          SY 8
    recTreenon->SetBranchStatus("shw.shwlid.width",1);                              // Shower Branch          SY 9
    recTreenon->SetBranchStatus("shw.shwlid.gap",1);                                // Shower Branch          SY 11
    recTreenon->SetBranchStatus("shw.nshwlid",1);                                   // Shower Branch          SY 12  
    recTreenon->SetBranchStatus("slc.nmiphit",1);                                   // Slice Branch           SY 13

    recTreenon->SetBranchStatus("sel.cosrej.anglekal",1);                           // Kirk 1  
    //recTreenon->SetBranchStatus("trk.kalman.dir.*",1);                            // Kirk 2   
    recTreenon->SetBranchStatus("trk.kalman.dir.fY",1);                             // Kirk 2   
    recTreenon->SetBranchStatus("slc.boxmax.*",1);                                  // Kirk 3
    recTreenon->SetBranchStatus("trk.kalman.nhit",1);                               // Kirk 4
    recTreenon->SetBranchStatus("trk.kalman.len",1);                                // Kirk 5
    recTreenon->SetBranchStatus("sel.contain.kalfwdcell",1);                        // Kirk 6 
    recTreenon->SetBranchStatus("sel.contain.kalbakcell",1);                        // Kirk 6
    recTreenon->SetBranchStatus("sel.cosrej.scatt",1);                              // Kirk 7
    recTreenon->SetBranchStatus("slc.nhit",1);                                      // Kirk 8
    recTreenon->SetBranchStatus("slc.calE",1);                                      // Kirk 9
    recTreenon->SetBranchStatus("slc.boxmin.*",1);                                  // Kirk 10
    recTreenon->SetBranchStatus("trk.nkalman",1 );                                  // Kirk 11
    recTreenon->SetBranchStatus("trk.kalman.calE",1);
    recTreenon->SetBranchStatus("sel.cosrej.numucontpid",1);
    recTreenon->SetBranchStatus("sel.contain.cosfwdcell",1);
    recTreenon->SetBranchStatus("sel.contain.cosbakcell",1);
    recTreenon->SetBranchStatus("sel.contain.nplanestofront",1);
    recTreenon->SetBranchStatus("sel.contain.nplanestoback",1);
    
    Int_t nbeam = recTreenon->GetEntries(); 
    //Int_t nbeam = recTreenon->GetEntriesFast(); 
    for (Int_t j = 0; j < nbeam; ++j) {                                        //Starting to loop over cosmic events  
      recTreenon->GetEntry(j);      
      cerr << "\r-- Processing event " << j << " of " << nbeam;
      
      //if(recTreeObject->trk.nkalman == 0) continue;
      Short_t nnu   = recTreeObject->mc.nnu;
      //std::cout<<"------------nnu is "<<nnu<<std::endl;    
      if(nnu != 1) continue;      
      Int_t   nccc  = recTreeObject->mc.nu[0].iscc;
      Short_t pdg   = recTreeObject->mc.nu[0].pdg;
      //Short_t CalE = recTreeObject->slc.calE;
      ///(recTreeObject->sel.cvn.ncid)
      // NC Offical Quality Cuts                                                                                                       
      if(recTreeObject->sel.nuecosrej.hitsperplane >= 8) continue;
      if(recTreeObject->shw.nshwlid == 0 ) continue;
      if(recTreeObject->shw.shwlid[0].gap >= 100.) continue;
      if(recTreeObject->slc.ncontplanes <= 2) continue;
      if(recTreeObject->vtx.nelastic == 0) continue;
      // NC Official Fiducial Cuts
      if (recTreeObject->vtx.elastic[0].vtx.fX < -680.0) continue;
      if (recTreeObject->vtx.elastic[0].vtx.fX >  650.0) continue;
      if (recTreeObject->vtx.elastic[0].vtx.fY < -720.0) continue;
      if (recTreeObject->vtx.elastic[0].vtx.fY >  500.0) continue; ///////Loosen Cuts
      if (recTreeObject->vtx.elastic[0].vtx.fZ <   50.0) continue;
      if (recTreeObject->vtx.elastic[0].vtx.fZ >  5450.0) continue;
      // NC Official Containment Cuts                                                                                                
      if(recTreeObject->sel.nuecosrej.startfront < 10) continue;
      if(recTreeObject->sel.nuecosrej.stopfront < 10) continue;  
      if(recTreeObject->sel.nuecosrej.startback < 10) continue;
      if(recTreeObject->sel.nuecosrej.stopback < 10) continue;
      if(recTreeObject->sel.nuecosrej.starteast < 10) continue;
      if(recTreeObject->sel.nuecosrej.stopeast < 10) continue;
      if(recTreeObject->sel.nuecosrej.startwest < 10) continue;
      if(recTreeObject->sel.nuecosrej.stopwest < 10) continue;  
      if(recTreeObject->sel.nuecosrej.starttop < 10) continue;
      if(recTreeObject->sel.nuecosrej.stoptop < 10) continue;
      if(recTreeObject->sel.nuecosrej.startbottom < 10) continue;
      if(recTreeObject->sel.nuecosrej.stopbottom < 10) continue;

      // NC Official NC/CC rejection                                                                                                                
      //if (recTreeObject->slc.nhit >= 200) continue;                                                                                                               
      //if (recTreeObject->slc.nhit < 20) continue;                                                                                                                 
      //if (recTreeObject->sel.cvn.ncid < 0.2)  continue;

      // Muon Removed Cut
      if(recTreeObject->sel.remid.pid >= 0.9 ) continue;
      // Energy Cut
      if(recTreeObject->slc.calE > 4) continue;
      if(recTreeObject->slc.calE < 0.5) continue;

      Run         = recTreeObject->hdr.run;                                  // True Info
      SubRun      = recTreeObject->hdr.subrun;                               // True Info
      Evt         = recTreeObject->hdr.evt;                                  // True Info
      SubEvt      = recTreeObject->hdr.subevt;                               // True Info

      starttop    = recTreeObject->sel.nuecosrej.starttop;                   // NueCosRej
      startbottom = recTreeObject->sel.nuecosrej.startbottom;                // NueCosRe
      startfront  = recTreeObject->sel.nuecosrej.startfront;                 // NueCosRej
      startback   = recTreeObject->sel.nuecosrej.startback;                  // NueCosRej
      startwest   = recTreeObject->sel.nuecosrej.startwest;                  // NueCosRej
      starteast   = recTreeObject->sel.nuecosrej.starteast;                  // NueCosRej
      stoptop     = recTreeObject->sel.nuecosrej.stoptop;                    // NueCosRej
      stopbottom  = recTreeObject->sel.nuecosrej.stopbottom;                 // NueCosRej
      stopfront   = recTreeObject->sel.nuecosrej.stopfront;                  // NueCosRej
      stopback    = recTreeObject->sel.nuecosrej.stopback;                   // NueCosRej
      stopwest    = recTreeObject->sel.nuecosrej.stopwest;                   // NueCosRej
      stopeast    = recTreeObject->sel.nuecosrej.stopeast;                   // NueCosRej
      cosdang     = recTreeObject->sel.nuecosrej.cosdang;                    // NueCosRej
      vtxdoca     = recTreeObject->sel.nuecosrej.vtxdoca;                    // NueCosRej
      prongmaxx   = recTreeObject->sel.nuecosrej.prongmaxx;                  // NueCosRej
      prongmaxy   = recTreeObject->sel.nuecosrej.prongmaxy;                  // NueCosRej
      prongmaxz   = recTreeObject->sel.nuecosrej.prongmaxz;                  // NueCosRej
      prongminx   = recTreeObject->sel.nuecosrej.prongminx;                  // NueCosRej
      prongminy   = recTreeObject->sel.nuecosrej.prongminy;                  // NueCosRej
      prongminz   = recTreeObject->sel.nuecosrej.prongminz;                  // NueCosRej
      sparsenessasymm  = recTreeObject->sel.nuecosrej.sparsenessasymm;       // NueCosRej
      hitsperplane = recTreeObject->sel.nuecosrej.hitsperplane;              // NueCosRej
      hitsperplaneasymm = recTreeObject->sel.nuecosrej.hitsperplaneasymm;    // NueCosRej
      musliceidx  = recTreeObject->sel.nuecosrej.musliceidx;                 // NueCosRej
      muanglediff = recTreeObject->sel.nuecosrej.muanglediff;                // NueCosRej
      mutimediff  = recTreeObject->sel.nuecosrej.mutimediff;                 // NueCosRej
      muclosestapproach  = recTreeObject->sel.nuecosrej.muclosestapproach;   // NueCosRej
      numupid     = recTreeObject->sel.cosrej.numucontpid;  

      pngptp      = recTreeObject->sel.nuecosrej.pngptp;                          // NueCosRej                                                             
      //sumtx       = recTreeObject->sand.nus.sumtx;                              // NuS Sandbox
      //sumty       = recTreeObject->sand.nus.sumty;                              // NuS Sandbox
      //ewsumtx     = recTreeObject->sand.nus.ewsumtx;                            // NuS Sandbox
      //ewsumty     = recTreeObject->sand.nus.ewsumty;                            // NuS Sandbox
      //cossumtx    = recTreeObject->sand.nus.cossumtx;                           // NuS Sandbox
      //cossumty    = recTreeObject->sand.nus.cossumty;                           // NuS Sandbox
      //cosewsumtx  = recTreeObject->sand.nus.cosewsumtx;                         // NuS Sandbox
      //cosewsumty  = recTreeObject->sand.nus.cosewsumty;                         // NuS Sandbox
      //angsumtx    = recTreeObject->sand.nus.angsumtx;                           // NuS Sandbox
      //angsumty    = recTreeObject->sand.nus.angsumty;                           // NuS Sandbox
      //angewsumtx  = recTreeObject->sand.nus.angewsumtx;                         // NuS Sandbox
      //angewsumty  = recTreeObject->sand.nus.angewsumty;                         // NuS Sandbox

      shwnplane   = recTreeObject->shw.shwlid[0].nplane;                           // Shower Branch
      shwmaxplanecont = recTreeObject->shw.shwlid[0].maxplanecont;                 // Shower Branch
      shwmaxplanegap  = recTreeObject->shw.shwlid[0].maxplanegap;                  // Shower Branch
      shwstartX   = recTreeObject->shw.shwlid[0].start.fX;                         // Shower Branch
      shwstartY   = recTreeObject->shw.shwlid[0].start.fY;                         // Shower Branch
      shwstartZ   = recTreeObject->shw.shwlid[0].start.fZ;                         // Shower Branch
      shwdirX     = recTreeObject->shw.shwlid[0].dir.fX;                           // Shower Branch
      shwdirZ     = recTreeObject->shw.shwlid[0].dir.fZ;                           // Shower Branch
      //shwNplaneX  = recTreeObject->shw.shwlid[0].nplanex;                        // Shower Branch
      //shwNplaneY  = recTreeObject->shw.shwlid[0].nplaney;                        // Shower Branch
      //shwstopX    = recTreeObject->shw.shwlid[0].stop.fX;                        // Shower Branch
      //shwstopY    = recTreeObject->shw.shwlid[0].stop.fY;                        // Shower Branch
      //shwstopZ    = recTreeObject->shw.shwlid[0].stop.fZ;                        // Shower Branch

      //nueid       = recTreeObject->sel.cvn.nueid;                                // CVN
      //numuid      = recTreeObject->sel.cvn.numuid;                               // CVN
      //nutauid     = recTreeObject->sel.cvn.nutauid;                              // CVN

      nhit        = recTreeObject->slc.nhit;                                       // Slice Branch
      //ncalhit     = recTreeObject->slc.ncalhit;                                  // Slice Branch
      //ncontplanes = recTreeObject->slc.ncontplanes;                              // Slice Branch
      //ncellsfromedge = recTreeObject->slc.;                                      // Slice Branch
      //starttime   = recTreeObject->slc.starttime;                                // Slice Branch
      //endtime     = recTreeObject->slc.endtime;                                  // Slice Branch
      //meantime    = recTreeObject->slc.meantime;                                 // Slice Branch
      //meanposX    = recTreeObject->slc.meanpos.fX;                               // Slice Branch   
      //meanposY    = recTreeObject->slc.meanpos.fY;                               // Slice Branch     
      //meanposZ    = recTreeObject->slc.meanpos.fZ;                               // Slice Branch  
      
      rempid      = recTreeObject->sel.remid.pid;                                  // RemID
      remlen      = recTreeObject->sel.remid.len;                                  // RemID

      //vtxTime     = recTreeObject->vtx.elastic[0].time;                          // Vertex Branch
      //vtxX        = recTreeObject->vtx.elastic[0].vtx.fX;                        // Vertex Branch
      //vtxY        = recTreeObject->vtx.elastic[0].vtx.fY;                        // Vertex Branch
      //vtxZ        = recTreeObject->vtx.elastic[0].vtx.fZ;                        // Vertex Branch
      //vtxX        = recTreeObject->vtx.elastic;                                  // Vertex Branch

      ncid        = recTreeObject->sel.cvn.ncid;                                   // CVN                 SY 1
      //cosmicid    = recTreeObject->sel.cvn.kCosmic; 
      cosmicid    = recTreeObject->sel.cvn.output[14]; 
      partptp     = recTreeObject->sel.nuecosrej.partptp;                          // NueCosRej           SY 2                                     
      shwnhit     = recTreeObject->shw.shwlid[0].nhit;                             // Shower Branch       SY 3
      shwnhitx    = recTreeObject->shw.shwlid[0].nhitx;                            // Shower Branch       SY 4
      shwnhity    = recTreeObject->shw.shwlid[0].nhity;                            // Shower Branch       SY 5
      shwxminusy  = shwnhitx - shwnhity;                                           //                     SY 6
      shwxplusy   = shwnhitx + shwnhity;                                           //                     SY 7
      shwxovery   = shwxminusy/shwxplusy;                                          //                     SY 8
      shwcalE     = recTreeObject->shw.shwlid[0].calE;                             // Shower Branch       SY 9
      shwdirY     = recTreeObject->shw.shwlid[0].dir.fY;                           // Shower Branch       SY 10
      shwlen      = recTreeObject->shw.shwlid[0].len;                              // Shower Branch       SY 11
      shwwwidth   = recTreeObject->shw.shwlid[0].width;                            // Shower Branch       SY 12
      shwGap      = recTreeObject->shw.shwlid[0].gap;                              // Shower Branch       SY 13
      nshwlid     = recTreeObject->shw.nshwlid;                                    // Shower Branch       SY 14
      nmiphit     = recTreeObject->slc.nmiphit;                                    // Slice Branch        SY 15

      anglekal    = recTreeObject->sel.cosrej.anglekal;                            // Kirk 1            
      dirFY       = recTreeObject->trk.kalman[0].dir.fY;                           // Kirk 2                                                         
      boxmaxFY    = recTreeObject->slc.boxmax.fY;                                  // Kirk 3                                                         
      kalnhit     = recTreeObject->trk.kalman[0].nhit;                             // Kirk 4     
      kallen      = recTreeObject->trk.kalman[0].len;                              // Kirk 5                 
      kalfwdcell  = recTreeObject->sel.contain.kalfwdcell;                         // Kirk 6                     
      kalbakcell  = recTreeObject->sel.contain.kalbakcell;                         // Kirk 6                     
      scatt       = recTreeObject->sel.cosrej.scatt;                               // Kirk 7                                                  
      nhit       = recTreeObject->slc.nhit;                                        // Kirk 8                              
      energy     = recTreeObject->slc.calE;                                        // Kirk 9                                          
      boxminFY   = recTreeObject->slc.boxmin.fY;                                   // Kirk 10                                    
      nkal       = recTreeObject->trk.nkalman;                                     // Kirk 11         
      
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
    recTreecos->SetBranchStatus("hdr.run",1);
    recTreecos->SetBranchStatus("hdr.subrun",1);
    recTreecos->SetBranchStatus("hdr.evt",1);
    recTreecos->SetBranchStatus("hdr.subevt",1);
    recTreecos->SetBranchStatus("mc.nnu",1);
    recTreecos->SetBranchStatus("mc.nu",1);
    recTreecos->SetBranchStatus("mc.nu.iscc",1);
    recTreecos->SetBranchStatus("mc.nu.pdg",1);
    recTreecos->SetBranchStatus("trk.ncosmic",1);
    
    recTreecos->SetBranchStatus("sand.nus.sumtx",1);                               // NuS Sandbox
    recTreecos->SetBranchStatus("sand.nus.sumty",1);                               // NuS Sandbox
    recTreecos->SetBranchStatus("sand.nus.ewsumtx",1);                             // NuS Sandbox
    recTreecos->SetBranchStatus("sand.nus.ewsumty",1);                             // NuS Sandbox
    recTreecos->SetBranchStatus("sand.nus.cossumtx",1);                            // NuS Sandbox
    recTreecos->SetBranchStatus("sand.nus.cossumty",1);                            // NuS Sandbox
    recTreecos->SetBranchStatus("sand.nus.cosewsumtx",1);                          // NuS Sandbox
    recTreecos->SetBranchStatus("sand.nus.cosewsumty",1);                          // NuS Sandbox
    recTreecos->SetBranchStatus("sand.nus.angsumtx",1);                            // NuS Sandbox
    recTreecos->SetBranchStatus("sand.nus.angsumty",1);                            // NuS Sandbox
    recTreecos->SetBranchStatus("sand.nus.angewsumtx",1);                          // NuS Sandbox
    recTreecos->SetBranchStatus("sand.nus.angewsumty",1);                          // NuS Sandbox    
    
    recTreecos->SetBranchStatus("sel.nuecosrej.starttop",1);                        // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.startbottom",1);                     // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.startfront",1);                      // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.startback",1);                       // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.startwest",1);                       // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.starteast",1);                       // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.stoptop",1);                         // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.stopbottom",1);                      // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.stopfront",1);                       // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.stopback",1);                        // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.stopwest",1);                        // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.stopeast",1);                        // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.cosdang",1);                         // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.vtxdoca",1);                         // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.prongmaxx",1);                       // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.prongmaxy",1);                       // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.prongmaxz",1);                       // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.prongminx",1);                       // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.prongminy",1);                       // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.prongminz",1);                       // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.sparsenessasymm",1);                 // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.hitsperplane",1);                    // NueCosRej  /////// 
    recTreecos->SetBranchStatus("sel.nuecosrej.hitsperplaneasymm",1);               // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.musliceidx",1);                      // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.muanglediff",1);                     // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.mutimediff",1);                      // NueCosRej 
    recTreecos->SetBranchStatus("sel.nuecosrej.muclosestapproach",1);               // NueCosRej
    recTreecos->SetBranchStatus("sel.nuecosrej.partptp",1);                         // NueCosRej                                                     
    recTreecos->SetBranchStatus("sel.nuecosrej.pngptp",1);                          // NueCosRej          

    recTreecos->SetBranchStatus("sel.cvn.ncid",1);                                  // CVN
    recTreecos->SetBranchStatus("sel.cvn.nueid",1);                                 // CVN
    recTreecos->SetBranchStatus("sel.cvn.numuid",1);                                // CVN
    recTreecos->SetBranchStatus("sel.cvn.nutauid",1);                               // CVN

    recTreecos->SetBranchStatus("sel.remid.pid",1);                                 // RemID
    recTreecos->SetBranchStatus("sel.remid.len",1);                                 // RemID

    recTreecos->SetBranchStatus("slc.nhit",1);                                      // Slice Branch
    //recTreecos->SetBranchStatus("slc.calhit",1);                                  // Slice Branch
    recTreecos->SetBranchStatus("slc.nmiphit",1);                                   // Slice Branch
    recTreecos->SetBranchStatus("slc.ncontplanes",1);                               // Slice Branch
    recTreecos->SetBranchStatus("slc.ncellsfromedge",1);                            // Slice Branch
    recTreecos->SetBranchStatus("slc.starttime",1);                                 // Slice Branch
    recTreecos->SetBranchStatus("slc.endtime",1);                                   // Slice Branch
    recTreecos->SetBranchStatus("slc.meantime",1);                                  // Slice Branch
    recTreecos->SetBranchStatus("slc.meanpos.fX",1);                                // Slice Branch
    recTreecos->SetBranchStatus("slc.meanpos.fY",1);                                // Slice Branch
    recTreecos->SetBranchStatus("slc.meanpos.fZ",1);                                // Slice Branch
   
    recTreecos->SetBranchStatus("shw.shwlid.nhit",1);                               // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.nhitx",1);                              // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.nhity",1);                              // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.nplane",1);                             // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.maxplanecont",1);                       // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.maxplanegap",1);                        // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.calE",1);                               // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.start.fX",1);                           // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.start.fY",1);                           // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.start.fZ",1);                           // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.dir.fX",1);                             // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.dir.fY",1);                             // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.dir.fZ",1);                             // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.len",1);                                // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.width",1);                              // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.nplanex",1);                            // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.nplaney",1);                            // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.gap",1);                                // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.stop.fX",1);                            // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.stop.fY",1);                            // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.stop.fZ",1);                            // Shower Branch
    recTreecos->SetBranchStatus("shw.nshwlid",1);                                   // Shower Branch
    recTreecos->SetBranchStatus("shw.shwlid.gap",1);                                // Shower Branch

    recTreecos->SetBranchStatus("vtx.elastic.time",1);                              // Vertex Branch
    recTreecos->SetBranchStatus("vtx.elastic.vtx.fX",1);                            // Vertex Branch
    recTreecos->SetBranchStatus("vtx.elastic.vtx.fY",1);                            // Vertex Branch
    recTreecos->SetBranchStatus("vtx.elastic.vtx.fZ",1);                            // Vertex Branch
    //recTreecos->SetBranchStatus("vtx.elastic.png3d.nhit",1);                      // Vertex Branch
    //recTreecos->SetBranchStatus("vtx.elastic.png3d.nhitx",1);                     // Vertex Branch
    //recTreecos->SetBranchStatus("vtx.elastic.png3d.nhity",1);                     // Vertex Branch
    recTreecos->SetBranchStatus("vtx.nelastic",1);                                  // Vertex Branch
    
    recTreecos->SetBranchStatus("sel.cvn.output",1);                                // CVN                    SY 1                       
    //recTreecos->SetBranchStatus("sel.cvn.output[14]",1);                          // CVN                    SY 1                       
    recTreecos->SetBranchStatus("sel.nuecosrej.partptp",1);                         // NueCosRej              SY 2            
    recTreecos->SetBranchStatus("shw.shwlid.nhit",1);                               // Shower Branch          SY 3
    recTreecos->SetBranchStatus("shw.shwlid.nhitx",1);                              // Shower Branch          SY 4
    recTreecos->SetBranchStatus("shw.shwlid.nhity",1);                              // Shower Branch          SY 5
    recTreecos->SetBranchStatus("shw.shwlid.calE",1);                               // Shower Branch          SY 6            
    recTreecos->SetBranchStatus("shw.shwlid.dir.fY",1);                             // Shower Branch          SY 7
    recTreecos->SetBranchStatus("shw.shwlid.len",1);                                // Shower Branch          SY 8
    recTreecos->SetBranchStatus("shw.shwlid.width",1);                              // Shower Branch          SY 9
    recTreecos->SetBranchStatus("shw.shwlid.gap",1);                                // Shower Branch          SY 11
    recTreecos->SetBranchStatus("shw.nshwlid",1);                                   // Shower Branch          SY 12  
    recTreecos->SetBranchStatus("slc.nmiphit",1);                                   // Slice Branch           SY 13

    recTreecos->SetBranchStatus("sel.cosrej.anglekal",1);             // Kirk 1
    recTreecos->SetBranchStatus("trk.kalman.dir.*",1);                // Kirk 2   
    recTreecos->SetBranchStatus("slc.boxmax.*",1);                    // Kirk 3
    recTreecos->SetBranchStatus("trk.kalman.nhit",1);                 // Kirk 4
    recTreecos->SetBranchStatus("trk.kalman.len",1);                  // Kirk 5
    recTreecos->SetBranchStatus("sel.contain.kalfwdcell",1);          // Kirk 6
    recTreecos->SetBranchStatus("sel.contain.kalbakcell",1);          // Kirk 6
    recTreecos->SetBranchStatus("sel.cosrej.scatt",1);                // Kirk 7
    recTreecos->SetBranchStatus("slc.nhit",1);                        // Kirk 8
    recTreecos->SetBranchStatus("slc.calE",1);                        // Kirk 9
    recTreecos->SetBranchStatus("slc.boxmin.*",1);                    // Kirk 10
    recTreecos->SetBranchStatus("trk.nkalman",1 );                    // Kirk 11

    recTreecos->SetBranchStatus("trk.kalman.calE",1);

    recTreecos->SetBranchStatus("sel.cosrej.numucontpid",1);
    recTreecos->SetBranchStatus("sel.contain.cosfwdcell",1);
    recTreecos->SetBranchStatus("sel.contain.cosbakcell",1);
    recTreecos->SetBranchStatus("sel.contain.nplanestofront",1);
    recTreecos->SetBranchStatus("sel.contain.nplanestoback",1);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    recTreecos->SetBranchStatus("sel.elecid.shwlid.ismuon",1);                        //sel.elecid.shwlid[0].ismuon
    recTreecos->SetBranchStatus("sel.cosrej.numucontpid",1);
    //sr->sel.cosrej.numucontpid    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    Int_t ncosmic = recTreecos->GetEntriesFast(); 
    for (Int_t i = 0; i < ncosmic; ++i) {                                       //Starting to loop over cosmic events  
      recTreecos->GetEntry(i);      
      cerr << "\r-- Processing event " << i << " of " << ncosmic;
      //if(recTreeObject->mc.nnu == 0) continue;
      //if(recTreeObject->trk.nkalman == 0) continue;

      //Short_t nnu   = recTreeObject->mc.nnu;
      //std::cout<<"------------nnu is "<<nnu<<std::endl;    
      //if(nnu != 1) continue;      
      //Int_t   nccc  = recTreeObject->mc.nu[0].iscc;
      //Short_t pdg   = recTreeObject->mc.nu[0].pdg;

      // NC Offical Quality Cuts                                                                                                       
      if(recTreeObject->sel.nuecosrej.hitsperplane >= 8) continue;
      if(recTreeObject->shw.nshwlid == 0 ) continue;
      if(recTreeObject->shw.shwlid[0].gap >= 100.) continue;
      if(recTreeObject->slc.ncontplanes <= 2) continue;
      if(recTreeObject->vtx.nelastic == 0) continue;

      // NC Official Fiducial Cuts
      if (recTreeObject->vtx.elastic[0].vtx.fX < -680.0) continue;
      if (recTreeObject->vtx.elastic[0].vtx.fX >  650.0) continue;
      if (recTreeObject->vtx.elastic[0].vtx.fY < -720.0) continue;
      if (recTreeObject->vtx.elastic[0].vtx.fY >  500.0) continue; ///////Loosen Cuts                                                                          
      if (recTreeObject->vtx.elastic[0].vtx.fZ <   50.0) continue;
      if (recTreeObject->vtx.elastic[0].vtx.fZ >  5450.0) continue;

      // NC Official Containment Cuts                                                                                                                           
      if(recTreeObject->sel.nuecosrej.startfront < 10) continue;
      if(recTreeObject->sel.nuecosrej.stopfront < 10) continue;
      if(recTreeObject->sel.nuecosrej.startback < 10) continue;
      if(recTreeObject->sel.nuecosrej.stopback < 10) continue;
      if(recTreeObject->sel.nuecosrej.starteast < 10) continue;
      if(recTreeObject->sel.nuecosrej.stopeast < 10) continue;
      if(recTreeObject->sel.nuecosrej.startwest < 10) continue;
      if(recTreeObject->sel.nuecosrej.stopwest < 10) continue;
      if(recTreeObject->sel.nuecosrej.starttop < 10) continue;
      if(recTreeObject->sel.nuecosrej.stoptop < 10) continue;
      if(recTreeObject->sel.nuecosrej.startbottom < 10) continue;
      if(recTreeObject->sel.nuecosrej.stopbottom < 10) continue;

      // NC Official NC/CC rejection                                                                                                                              
      //if (recTreeObject->slc.nhit >= 200) continue;
      //if (recTreeObject->slc.nhit < 20) continue;
      //if (recTreeObject->sel.cvn.ncid < 0.2)  continue;

      // Muon Removed Cut
      if(recTreeObject->sel.remid.pid >= 0.9 ) continue;

      // Energy Cut                                                                                                                             
      if(recTreeObject->slc.calE > 4) continue;
      if(recTreeObject->slc.calE < 0.5) continue;

      Run    = recTreeObject->hdr.run;
      SubRun = recTreeObject->hdr.subrun;
      Evt    = recTreeObject->hdr.evt;
      SubEvt = recTreeObject->hdr.subevt;

      //sumtx       = recTreeObject->sand.nus.sumtx;                           // NuS Sandbox
      //sumty       = recTreeObject->sand.nus.sumty;                           // NuS Sandbox
      //ewsumtx     = recTreeObject->sand.nus.ewsumtx;                         // NuS Sandbox
      //ewsumty     = recTreeObject->sand.nus.ewsumty;                         // NuS Sandbox
      //cossumtx    = recTreeObject->sand.nus.cossumtx;                        // NuS Sandbox
      //cossumty    = recTreeObject->sand.nus.cossumty;                        // NuS Sandbox
      //cosewsumtx  = recTreeObject->sand.nus.cosewsumtx;                      // NuS Sandbox
      //cosewsumty  = recTreeObject->sand.nus.cosewsumty;                      // NuS Sandbox
      //angsumtx    = recTreeObject->sand.nus.angsumtx;                        // NuS Sandbox
      //angsumty    = recTreeObject->sand.nus.angsumty;                        // NuS Sandbox
      //angewsumtx  = recTreeObject->sand.nus.angewsumtx;                      // NuS Sandbox
      //angewsumty  = recTreeObject->sand.nus.angewsumty;                      // NuS Sandbox
      
      starttop    = recTreeObject->sel.nuecosrej.starttop;                   // NueCosRej
      startbottom = recTreeObject->sel.nuecosrej.startbottom;                // NueCosRe
      startfront  = recTreeObject->sel.nuecosrej.startfront;                 // NueCosRej
      startback   = recTreeObject->sel.nuecosrej.startback;                  // NueCosRej
      startwest   = recTreeObject->sel.nuecosrej.startwest;                  // NueCosRej
      starteast   = recTreeObject->sel.nuecosrej.starteast;                  // NueCosRej
      stoptop     = recTreeObject->sel.nuecosrej.stoptop;                    // NueCosRej
      stopbottom  = recTreeObject->sel.nuecosrej.stopbottom;                 // NueCosRej
      stopfront   = recTreeObject->sel.nuecosrej.stopfront;                  // NueCosRej
      stopback    = recTreeObject->sel.nuecosrej.stopback;                   // NueCosRej
      stopwest    = recTreeObject->sel.nuecosrej.stopwest;                   // NueCosRej
      stopeast    = recTreeObject->sel.nuecosrej.stopeast;                   // NueCosRej
      cosdang     = recTreeObject->sel.nuecosrej.cosdang;                    // NueCosRej
      vtxdoca     = recTreeObject->sel.nuecosrej.vtxdoca;                    // NueCosRej
      prongmaxx   = recTreeObject->sel.nuecosrej.prongmaxx;                  // NueCosRej
      prongmaxy   = recTreeObject->sel.nuecosrej.prongmaxy;                  // NueCosRej
      prongmaxz   = recTreeObject->sel.nuecosrej.prongmaxz;                  // NueCosRej
      prongminx   = recTreeObject->sel.nuecosrej.prongminx;                  // NueCosRej
      prongminy   = recTreeObject->sel.nuecosrej.prongminy;                  // NueCosRej
      prongminz   = recTreeObject->sel.nuecosrej.prongminz;                  // NueCosRej
      sparsenessasymm  = recTreeObject->sel.nuecosrej.sparsenessasymm;       // NueCosRej
      hitsperplane = recTreeObject->sel.nuecosrej.hitsperplane;              // NueCosRej
      hitsperplaneasymm = recTreeObject->sel.nuecosrej.hitsperplaneasymm;    // NueCosRej
      musliceidx  = recTreeObject->sel.nuecosrej.musliceidx;                 // NueCosRej
      muanglediff = recTreeObject->sel.nuecosrej.muanglediff;                // NueCosRej
      mutimediff  = recTreeObject->sel.nuecosrej.mutimediff;                 // NueCosRej
      muclosestapproach  = recTreeObject->sel.nuecosrej.muclosestapproach;   // NueCosRej
      
      //partptp     = recTreeObject->sel.nuecosrej.partptp;                    // NueCosRej                                                            
      pngptp      = recTreeObject->sel.nuecosrej.pngptp;                     // NueCosRej                                                             
      
      numupid     = recTreeObject->sel.cosrej.numucontpid;  

      //ncid        = recTreeObject->sel.cvn.ncid;                             // CVN
      //nueid       = recTreeObject->sel.cvn.nueid;                            // CVN
      //numuid      = recTreeObject->sel.cvn.numuid;                           // CVN
      //nutauid     = recTreeObject->sel.cvn.nutauid;                          // CVN

      rempid      = recTreeObject->sel.remid.pid;                              // RemID
      remlen      = recTreeObject->sel.remid.len;                              // RemID

      // nhit        = recTreeObject->slc.nhit;                                // Slice Branch
      //ncalhit     = recTreeObject->slc.ncalhit;                              // Slice Branch
      //nmiphit     = recTreeObject->slc.nmiphit;                              // Slice Branch
      //ncontplanes = recTreeObject->slc.ncontplanes;                          // Slice Branch
      //ncellsfromedge = recTreeObject->slc.;                                  // Slice Branch
      //starttime   = recTreeObject->slc.starttime;                            // Slice Branch
      //endtime     = recTreeObject->slc.endtime;                              // Slice Branch
      //meantime    = recTreeObject->slc.meantime;                             // Slice Branch
      //meanposX    = recTreeObject->slc.meanpos.fX;                           // Slice Branch   
      //meanposY    = recTreeObject->slc.meanpos.fY;                           // Slice Branch     
      //meanposZ    = recTreeObject->slc.meanpos.fZ;                           // Slice Branch  
      
      shwnhit     = recTreeObject->shw.shwlid[0].nhit;                           // Shower Branch
      shwnhitx    = recTreeObject->shw.shwlid[0].nhitx;                          // Shower Branch
      //shwnhity    = recTreeObject->shw.shwlid[0].nhity;                          // Shower Branch
      shwnplane   = recTreeObject->shw.shwlid[0].nplane;                         // Shower Branch
      shwmaxplanecont = recTreeObject->shw.shwlid[0].maxplanecont;               // Shower Branch
      shwmaxplanegap  = recTreeObject->shw.shwlid[0].maxplanegap;                // Shower Branch
      //shwcalE     = recTreeObject->shw.shwlid[0].calE;                           // Shower Branch
      shwstartX   = recTreeObject->shw.shwlid[0].start.fX;                       // Shower Branch
      shwstartY   = recTreeObject->shw.shwlid[0].start.fY;                       // Shower Branch
      shwstartZ   = recTreeObject->shw.shwlid[0].start.fZ;                       // Shower Branch
      shwdirX     = recTreeObject->shw.shwlid[0].dir.fX;                         // Shower Branch
      //shwdirY     = recTreeObject->shw.shwlid[0].dir.fY;                         // Shower Branch
      shwdirZ     = recTreeObject->shw.shwlid[0].dir.fZ;                         // Shower Branch
      //shwlen      = recTreeObject->shw.shwlid[0].len;                            // Shower Branch
      //shwwwidth   = recTreeObject->shw.shwlid[0].width;                          // Shower Branch
      shwNplaneX  = recTreeObject->shw.shwlid[0].nplanex;                        // Shower Branch
      shwNplaneY  = recTreeObject->shw.shwlid[0].nplaney;                        // Shower Branch
      //shwGap      = recTreeObject->shw.shwlid[0].gap;                            // Shower Branch
      //shwstopX    = recTreeObject->shw.shwlid[0].stop.fX;                        // Shower Branch
      //shwstopY    = recTreeObject->shw.shwlid[0].stop.fY;                        // Shower Branch
      //shwstopZ    = recTreeObject->shw.shwlid[0].stop.fZ;                        // Shower Branch
      //nshwlid     = recTreeObject->shw.nshwlid;                               // Shower Branch

      //vtxTime     = recTreeObject->vtx.elastic[0].time;                          // Vertex Branch
      //vtxX        = recTreeObject->vtx.elastic[0].vtx.fX;                        // Vertex Branch
      //vtxY        = recTreeObject->vtx.elastic[0].vtx.fY;                        // Vertex Branch
      //vtxZ        = recTreeObject->vtx.elastic[0].vtx.fZ;                        // Vertex Branch

      ncid        = recTreeObject->sel.cvn.ncid;                                   // CVN                 SY 1
      cosmicid    = recTreeObject->sel.cvn.output[14]; 
      partptp     = recTreeObject->sel.nuecosrej.partptp;                          // NueCosRej           SY 2                                     
      shwnhit     = recTreeObject->shw.shwlid[0].nhit;                             // Shower Branch       SY 3
      shwnhitx    = recTreeObject->shw.shwlid[0].nhitx;                            // Shower Branch       SY 4
      shwnhity    = recTreeObject->shw.shwlid[0].nhity;                            // Shower Branch       SY 5
      shwxminusy  = shwnhitx - shwnhity;                                           //                     SY 6
      shwxplusy   = shwnhitx + shwnhity;                                           //                     SY 7
      shwxovery   = shwxminusy/shwxplusy;                                          //                     SY 8
      shwcalE     = recTreeObject->shw.shwlid[0].calE;                             // Shower Branch       SY 9
      shwdirY     = recTreeObject->shw.shwlid[0].dir.fY;                           // Shower Branch       SY 10
      shwlen      = recTreeObject->shw.shwlid[0].len;                              // Shower Branch       SY 11
      shwwwidth   = recTreeObject->shw.shwlid[0].width;                            // Shower Branch       SY 12
      shwGap      = recTreeObject->shw.shwlid[0].gap;                              // Shower Branch       SY 13
      nshwlid     = recTreeObject->shw.nshwlid;                                    // Shower Branch       SY 14
      nmiphit     = recTreeObject->slc.nmiphit;                                    // Slice Branch        SY 15

      anglekal   = recTreeObject->sel.cosrej.anglekal;                            // Kirk 1        
      dirFY      = recTreeObject->trk.kalman[0].dir.fY;                           // Kirk 2        
      boxmaxFY   = recTreeObject->slc.boxmax.fY;                                  // Kirk 3 
      kalnhit    = recTreeObject->trk.kalman[0].nhit;                             // Kirk 4   
      kallen     = recTreeObject->trk.kalman[0].len;                              // Kirk 5    
      kalfwdcell = recTreeObject->sel.contain.kalfwdcell;                         // Kirk 6
      kalbakcell = recTreeObject->sel.contain.kalbakcell;                         // Kirk 6     
      scatt      = recTreeObject->sel.cosrej.scatt;                               // Kirk 7                       
      nhit       = recTreeObject->slc.nhit;                                       // Kirk 8     
      energy     = recTreeObject->slc.calE;                                       // Kirk 9 
      boxminFY   = recTreeObject->slc.boxmin.fY;                                  // Kirk 10              
      nkal       = recTreeObject->trk.nkalman;                                    // Kirk 11            
      csTree->Fill();
    }                                                                             // Ending to loop over cosmic events 
  }                                                                              //Ending to loop over files

  tMVALoad.Write();
  tMVALoad.Close();
}
#endif
