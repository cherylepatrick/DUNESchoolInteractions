// To run this, type: cafe Exercise2aSolution.C

// These are standard header files from the CAFAna analysis tool
// They allow you to load and plot variables for each interaction event
// in your simulation file (and later, in data files)
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"

// As we are working with simulation, we have access to information about the
// TRUE event - what GENIE, the simulation program, simulated for the event
// this is separate to the RECONSTRUCTED information - what the simulation 
// program thinks the DUNE detector would have seen
#include "CAFAna/Vars/Vars.h" // Variables
#include "CAFAna/Cuts/TruthCuts.h" // Cuts

#include "StandardRecord/SRProxy.h" // A wrapper for the CAF format

// These files come from the ROOT data analysis package
// This is used in many particle-physics experiments to make plots
// and do some basic statistics, cuts etc. The CAF simulation and data files
// are a special DUNE-specific format or a ROOT data file
// ROOT classes all start with a T, for some reason. It makes them easy to search
// for, except for things like the unfortunate TAxis...
#include "TCanvas.h" // Plots are drawn on a "canvas"
#include "TH1.h" // 1-dimensional histogram
#include "TPad.h" // Canvases are divided into pads. That could let you draw more than one plot on a canvas, if you wanted to, by using multiple pads. We will not bother with that today.
#include "TLegend.h" // Lets us draw a legend
#include "TMath.h" // I'll use some basic math functions
#include "THStack.h" // We'll have a stacked histogram this time

// Standard C++ library for input and output
#include <iostream>


/* *****************
 Define some GENIE interaction modes.
 See full list at https://wiki.dunescience.org/wiki/Scattering_mode
 Use these to make your TRUTH cuts - this the interaction type GENIE simulated
 */

const int MODE_QE = 1;
const int MODE_RES = 4;
const int MODE_DIS = 3;
const int MODE_MEC = 10;

/* ********
 Define some PDG codes (particle identifiers from https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf)
 You can also find the ID's for things like protons, neutrons, pions and even whole nuclei in the list!
 */
const int PDG_MU=13;
const int PDG_E=11;
const int PDG_NUMU=14;
const int PDG_NUE=12;


using namespace ana;


// This is the main function. To use ROOT's interpreted interface, you need to define a function
// with the same name as your file (minus the .C file extension)
void Exercise2aSolution()
{
  // Define four options for our input CAF samples.
  // Environment variables and wildcards work, as do SAM datasets
  // (a metadata database Fermilab uses to organise large volumes of data and simulation files).
  
  const std::string NDGAR_FHC = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5gas/CAF_FHC_90*.root"; //ND-GAr FHC
  const std::string NDGAR_RHC = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5gas/CAF_RHC_90*.root"; //ND-GAr RHC
  const std::string NDLAR_FHC = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5/00/CAF_FHC_90*.root"; //ND-LAr FHC
  const std::string NDLAR_RHC = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5/00/CAF_RHC_90*.root"; //ND-LAr RHC
  
  // Source of events - load them from the one of the sets of files
  SpectrumLoader loader(NDGAR_FHC); // ***** Change this to use a different sample ***

  // We want to plot a histogram with 40 bins, covering the range 0 to 10 GeV
  const Binning binsEnergy = Binning::Simple(40, 0, 10);

  // Define the label, binning, and contents that we want for our first histogram
  // The axis label can be whatever you like.
  // The binning needs to be a Binning object like the one we just made.
  // The Variable that we are plotting can be a single variable or function of variables from
  // the CAFs. See https://wiki.dunescience.org/wiki/CAFAna_Variables
  // We want to plot true energy
  const HistAxis axTrue("True neutrino energy (GeV)", binsEnergy, kTrueEnergy);

  // Now we have defined an axis and a variable we want to plot on it, let's decide which
  // events to plot. So here we are telling the loader to load a spectrum with our defined axis
  // Here, we are defining the selections or "cuts" using the variable names in CAFAna https://wiki.dunescience.org/wiki/CAFAna_Cuts
  
  // This is the cut we used before to select CC nu_mu interactions
  // You need to update it to select only CCQE.
  //Spectrum sTrueENumu(loader, axTrue, kIsNumuCC && !kIsAntiNu);
  
  // Select the true interaction types
  const Cut kIsQE = SIMPLEVAR(mode) == MODE_QE; //The modes are defined at the top
  const Cut kIsRES = SIMPLEVAR(mode) == MODE_RES;
  const Cut kIsDIS = SIMPLEVAR(mode) == MODE_DIS;
  const Cut kIsMEC = SIMPLEVAR(mode) == MODE_MEC;

  // This time, we are looking for CC0pi - one negative muon, at least one proton, and no pions
  // Define the cut...
  const Cut kHasCC0PiFinalState([](const caf::SRProxy* sr)
                                {
                                  const int totPi = sr->nipip + sr->nipim + sr->nipi0;
                                  return sr->LepPDG == PDG_MU && sr->nP >= 1 && totPi == 0;
                                });
  
  // 4 Spectrum objects for the 4 true cuts
  Spectrum sCC0piQE(loader, axTrue, kIsQE && kHasCC0PiFinalState);
  Spectrum sCC0piRES(loader, axTrue, kIsRES && kHasCC0PiFinalState);
  Spectrum sCC0piMEC(loader, axTrue, kIsMEC && kHasCC0PiFinalState);
  Spectrum sCC0piDIS(loader, axTrue, kIsDIS && kHasCC0PiFinalState);
  
  // Fill all the Spectrum objects
  loader.Go();

  /* 
     Set to the same exposure as before
  */  
  const double pot = 1e20;

  // Convert and draw
  TCanvas *canvas = new TCanvas; // Make a canvas
  
  // Make them all into histograms
  // ROOT colors are defined at https://root.cern.ch/doc/master/classTColor.
  TH1D *hCC0piQE = sCC0piQE.ToTH1(pot, kAzure-7);
  TH1D *hCC0piRES =sCC0piRES.ToTH1(pot, kOrange-2);
  TH1D *hCC0piMEC =sCC0piMEC.ToTH1(pot, kOrange+7);
  TH1D *hCC0piDIS =sCC0piDIS.ToTH1(pot, kAzure-9);
  
  // This makes unfilled histograms, so let's fill 'em up!
  hCC0piQE->SetFillColor(kAzure-7);
  hCC0piRES->SetFillColor(kOrange-2);
  hCC0piMEC->SetFillColor(kOrange+7);
  hCC0piDIS->SetFillColor(kAzure-9);
  
  // Make a stacked histogram
  THStack *stack = new THStack("stack","");
  stack->Add(hCC0piDIS);
  stack->Add(hCC0piRES);
  stack->Add(hCC0piMEC);
  stack->Add(hCC0piQE);
  stack->Draw("hist");
  
  gPad->SetLogy(false);
  
  auto legend = new TLegend(0.65,0.65,0.9,0.9); // x and y coordinates of corners
  legend->SetHeader("Legend","C"); // option "C" to center the header
  legend->AddEntry(hCC0piQE,"QE","f");
  legend->AddEntry(hCC0piMEC,"MEC","f");
  legend->AddEntry(hCC0piRES,"RES","f");
  legend->AddEntry(hCC0piDIS,"DIS","f");
  legend->Draw();
  
  canvas->SaveAs("Exercise2a.png"); // Save the result
}
