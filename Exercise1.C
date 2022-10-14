// To run this, type: cafe  Exercise1.C
// To use only a fraction of the input files (e.g. 1 in 20), type: cafe --stride 20 Exercise1.C


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

// Standard C++ library for input and output
#include <iostream>

using namespace ana;


// This is the main function. To use ROOT's interpreted interface, you need to define a function
// with the same name as your file (minus the .C file extension)
void Exercise1()
{
  // Define four options for our input CAF samples.
  // Environment variables and wildcards work. As do SAM datasets
  // (a metadata database Fermilab uses to organise large volumes of data and simulation files).
  
  const std::string NDGAR_FHC = "/Disk/ds-sopa-group/PPE/dune/DuneSchool/CAFs/NDGAr/CAF_FHC_90*.root"; //ND-GAr FHC
  const std::string NDGAR_RHC = "/Disk/ds-sopa-group/PPE/dune/DuneSchool/CAFs/NDGAr/CAF_RHC_90*.root"; //ND-GAr RHC
  const std::string NDLAR_FHC = "//Disk/ds-sopa-group/PPE/dune/DuneSchool/CAFs/NDLAr/CAF_FHC_90*.root"; //ND-LAr FHC
  const std::string NDLAR_RHC = "/Disk/ds-sopa-group/PPE/dune/DuneSchool/CAFs/NDLAr/AF_RHC_90*.root"; //ND-LAr RHC
  
  // Source of events - load them from the one of the sets of files
  SpectrumLoader loader(NDGAR_FHC);
  // ******** Change this to use a different sample ***

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
  
  // Cut for True muon neutrino charged current interactions
  Spectrum sTrueENumu(loader, axTrue, kIsNumuCC && !kIsAntiNu);  
  // kIsNumuCC = Muon neutrino charged-current interactions.
  // kIsAntiNu = Interaction initiated by an antineutrino. !kIsAntiNu means it is NOT an antineutrino.
  // && means we want both the first AND the second condition to be true. If we wanted the first OR the second condition, we would use ||

  // *********  ADD MORE Spectrum OBJECTS HERE **********
  //Spectrum sTrueENumubar(??, ??, ??); // Muon antineutrino charged-current interactions


  // This is the call that actually fills in those spectra
  loader.Go();

  /* 
     The amount of simulation we have depends on how many files we're using.
     To get an idea of what DUNE detectors might see, we want to scale it to a DUNE exposure.
     We normally think about that in terms of "protons on target" (POT) - the total number of protons
     delivered by the Fermilab accelerator to generate our neutrino beam. Let's use the same
     number of POT that MINERvA used in its initial 5-month run - 10^20! That's a lot of protons...
  */  
  const double pot = 1e20;

  // Convert each spectrum to a histogram and scale it to our POT value - then draw it!
  TCanvas *canvas = new TCanvas; // Make a canvas
  
  TH1D *hTrueENumu = sTrueENumu.ToTH1(pot, kBlue);// Draw our spectrum in blue. ROOT colors are defined at https://root.cern.ch/doc/master/classTColor.html
  hTrueENumu->Draw("HIST"); // This time we turn our spectrum into a ROOT histogram, and draw that. It means we can use the histogram for other things - like a legend.
  TH1D *hTrueENumubar=sTrueENumubar.ToTH1(pot, kBlue, 7);//Antineutrinos are getting a dashed line.
  hTrueENumubar->Draw("HIST SAME"); // SAME canvas as the previous spectrum
  
  // ******** ADD THE OTHER SPECTRA HERE *********
  // How about red for the electron component...? Make the antineutrinos dashed again
  
  /*    ********* RHC MODE *********
   What happens when you draw in RHC mode? The maximum value of the y axis is set by the first histogram you draw. You might want to change the order you draw them in... (or if you're a ROOT pro, you can do some clever checks on the histograms and set the maximum manually... but that's a skill to learn another day)
   */
  
  
  gPad->SetLogy(); // To turn off, use gPad->SetLogy(false);
  auto legend = new TLegend(0.7,0.7,0.9,0.9); // x and y coordinates of corners
  legend->SetHeader("Legend","C"); // option "C" to center the header
  legend->AddEntry(hTrueENumu,"#nu_{#mu}","l");
  // ****** FILL IN THE LEGEND *******
  legend->Draw();
  canvas->SaveAs("Exercise1.png"); // Save the result
}
