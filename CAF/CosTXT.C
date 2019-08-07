#ifdef __CINT__
void CosTXT()
{
  std::cout << "Sorry, you must run in compiled mode" << std::endl;
}
#else

#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/EventList.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Cuts/NumuCuts.h"
#include "CAFAna/Cuts/NusCuts.h"
#include "CAFAna/Cuts/NusCuts17.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Cuts/TimingCuts.h"
#include "CAFAna/Decomp/ProportionalDecomp.h"
#include "CAFAna/Extrap/ExtrapSterile.h"
#include "CAFAna/Prediction/PredictionExtrap.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionSterile.h"
#include "CAFAna/Vars/GenieWeights.h"
#include "CAFAna/Vars/HistAxes.h"
#include "CAFAna/nus/Nus17/NusLoadProd3.h"
#include "CAFAna/Vars/GenieWeights.h"
#include "CAFAna/Analysis/Style.h"
#include "OscLib/func/OscCalculatorSterile.h"
#include "CAFAna/Analysis/Calcs.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"

#include "CAFAna/nus/NuSPlotFunctions.h"

#include "StandardRecord/Proxy/SRProxy.h"
#include "StandardRecord/StandardRecord.h"

using namespace ana;

void CosTXT()
{  
  const std::string fFDcosmicdata_test = "defname: prod_caf_R17-03-01-prod3reco.h_fd_cosmic_full_v1_goodruns with limit 1";

  //Input Variables
  const Var nhit            ([](const caf::SRProxy* sr){return (float)sr->slc.nhit;});
  const Var ncalhit         ([](const caf::SRProxy* sr){return (float)sr->slc.ncalhit;});
  const Var nmiphit         ([](const caf::SRProxy* sr){return (float)sr->slc.nmiphit;});
  const Var ncontplanes     ([](const caf::SRProxy* sr){return (float)sr->slc.ncontplanes;});
  const Var ncellsfromedge  ([](const caf::SRProxy* sr){return (float)sr->slc.ncellsfromedge;});
  const Var slcCalE         ([](const caf::SRProxy* sr){return (float)sr->slc.calE;});
  const Var boxminy         ([](const caf::SRProxy* sr){return (float)sr->slc.boxmin.y;});
  const Var boxmaxy         ([](const caf::SRProxy* sr){return (float)sr->slc.boxmax.y;});
  const Var closestslicetime([](const caf::SRProxy* sr){return (float)sr->slc.closestslicetime;});
  const Var closestslicenhit([](const caf::SRProxy* sr){return (float)sr->slc.closestslicenhit;});
  const Var sumtx           ([](const caf::SRProxy* sr){if(std::isnan(sr->sand.nus.sumtx))     {return -5.f;}return (float)sr->sand.nus.sumtx;});
  const Var sumty           ([](const caf::SRProxy* sr){if(std::isnan(sr->sand.nus.sumty))     {return -5.f;}return (float)sr->sand.nus.sumty;});
  const Var ewsumtx         ([](const caf::SRProxy* sr){if(std::isnan(sr->sand.nus.ewsumtx))   {return -5.f;}return (float)sr->sand.nus.ewsumtx;});
  const Var ewsumty         ([](const caf::SRProxy* sr){if(std::isnan(sr->sand.nus.ewsumty))   {return -5.f;}return (float)sr->sand.nus.ewsumty;});
  const Var cossumtx        ([](const caf::SRProxy* sr){if(std::isnan(sr->sand.nus.cossumtx))  {return -5.f;}return (float)sr->sand.nus.cossumtx;});
  const Var cossumty        ([](const caf::SRProxy* sr){if(std::isnan(sr->sand.nus.cossumty))  {return -5.f;}return (float)sr->sand.nus.cossumty;});
  const Var cosewsumtx      ([](const caf::SRProxy* sr){if(std::isnan(sr->sand.nus.cosewsumtx)){return -5.f;}return (float)sr->sand.nus.cosewsumtx;});
  const Var cosewsumty      ([](const caf::SRProxy* sr){if(std::isnan(sr->sand.nus.cosewsumty)){return -5.f;}return (float)sr->sand.nus.cosewsumty;});
  const Var angsumtx        ([](const caf::SRProxy* sr){if(std::isnan(sr->sand.nus.angsumtx))  {return -5.f;}return (float)sr->sand.nus.angsumtx;});
  const Var angsumty        ([](const caf::SRProxy* sr){if(std::isnan(sr->sand.nus.angsumty))  {return -5.f;}return (float)sr->sand.nus.angsumty;});
  const Var angewsumtx      ([](const caf::SRProxy* sr){if(std::isnan(sr->sand.nus.angewsumtx)){return -5.f;}return (float)sr->sand.nus.angewsumtx;});
  const Var angewsumty      ([](const caf::SRProxy* sr){if(std::isnan(sr->sand.nus.angewsumty)){return -5.f;}return (float)sr->sand.nus.angewsumty;});
  const Var fuzzykntot      ([](const caf::SRProxy* sr){if(sr->vtx.nelastic <1){return 900.f;}return (float)sr->vtx.elastic[0].fuzzyk.ntot;});
  const Var fuzzyknpng      ([](const caf::SRProxy* sr){if(sr->vtx.nelastic <1){return 900.f;}return (float)sr->vtx.elastic[0].fuzzyk.npng;});
  const Var fuzzyknpng2d    ([](const caf::SRProxy* sr){if(sr->vtx.nelastic <1){return 900.f;}return (float)sr->vtx.elastic[0].fuzzyk.npng2d;});
  const Var fuzzykpngdirx   ([](const caf::SRProxy* sr){if(sr->vtx.nelastic <1){return 900.f;}if(sr->vtx.elastic[0].fuzzyk.npng < 1)return 900.0f;return(float)sr->vtx.elastic[0].fuzzyk.png[0].dir.x;});
  const Var fuzzykpngdiry   ([](const caf::SRProxy* sr){if(sr->vtx.nelastic <1){return 900.f;}if(sr->vtx.elastic[0].fuzzyk.npng < 1)return 900.0f;return(float)sr->vtx.elastic[0].fuzzyk.png[0].dir.y;});
  const Var fuzzykpngdirz   ([](const caf::SRProxy* sr){if(sr->vtx.nelastic <1){return 900.f;}if(sr->vtx.elastic[0].fuzzyk.npng < 1)return 900.0f;return(float)sr->vtx.elastic[0].fuzzyk.png[0].dir.z;});
  const Var fuzzykpngCalE   ([](const caf::SRProxy* sr){if(sr->vtx.nelastic <1){return 900.f;}if(sr->vtx.elastic[0].fuzzyk.npng < 1)return 900.0f;return(float)sr->vtx.elastic[0].fuzzyk.png[0].calE;});
  const Var fuzzykpnglen    ([](const caf::SRProxy* sr){if(sr->vtx.nelastic <1){return 900.f;}if(sr->vtx.elastic[0].fuzzyk.npng < 1)return 900.0f;return(float)sr->vtx.elastic[0].fuzzyk.png[0].len;});
  const Var fuzyykpngnhit   ([](const caf::SRProxy* sr){if(sr->vtx.nelastic <1){return 900.f;}if(sr->vtx.elastic[0].fuzzyk.npng < 1)return 900.0f;return(float)sr->vtx.elastic[0].fuzzyk.png[0].nhit;});
  const Var fuzyykpngnhitx  ([](const caf::SRProxy* sr){if(sr->vtx.nelastic <1){return 900.f;}if(sr->vtx.elastic[0].fuzzyk.npng < 1)return 900.0f;return(float)sr->vtx.elastic[0].fuzzyk.png[0].nhitx;});
  const Var fuzyykpngnhity  ([](const caf::SRProxy* sr){if(sr->vtx.nelastic <1){return 900.f;}if(sr->vtx.elastic[0].fuzzyk.npng < 1)return 900.0f;return(float)sr->vtx.elastic[0].fuzzyk.png[0].nhity;});
  const Var cvnncid         ([](const caf::SRProxy* sr){if(sr->sel.cvn.noutput == 0) {return 900.f;}return (float)sr->sel.cvn.ncid;});
  
  const Cut kNus17HitCut(
			[](const caf::SRProxy* sr){
			  if(sr->slc.nhit     < 25)  return false;
			  return true;
			});
  
  std::vector<const Var*> vars = {&ncalhit,          &nmiphit,       &ncellsfromedge,   &boxminy,        &boxmaxy,      
				  &closestslicenhit, &sumtx,         &ewsumtx,          &angsumtx,       &angewsumtx, 
				  &cossumtx,         &cosewsumtx,    &sumty,            &ewsumty,        &angsumty,     
				  &angewsumty,       &cossumty,      &cosewsumty,       &fuzzykntot,     &fuzzyknpng, 
				  &fuzzyknpng2d,     &fuzzykpngdirx, &fuzzykpngdiry,    &fuzzykpngdirz,  &fuzzykpngCalE, 
				  &fuzzykpnglen,     &fuzyykpngnhit, &fuzyykpngnhitx,   &fuzyykpngnhity, &cvnncid};
  
  MakeTextListFile(fFDcosmicdata_test, {kNus17EnergyCut && kNus17EventQuality && kNus17FDContain && kNus17HitCut}, {"Cos.txt"}, vars, &kStandardSpillCuts);
}
#endif

