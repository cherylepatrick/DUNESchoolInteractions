// To run this, type: cafe --stride 20 Example1.C

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

// Standard C++ library for input and output
#include <iostream>

using namespace ana;


// This is the main function. To use ROOT's interpreted interface, you need to define a function
// with the same name as your file (minus the .C file extension)
void Example1()
{
  // This points to the source files we will use. These ones are for the ND-GAr.
  // Environment variables and wildcards work, as do SAM datasets
  // (a metadata database Fermilab uses to organise large volumes of data and simulation files).
  const std::string fname = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5gas/CAF_FHC_90*.root";

  // Source of events - load them from the files
  SpectrumLoader loader(fname);

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
  sTrueENumu.ToTH1(pot, kBlue)->Draw("hist"); // Draw our spectrum in blue. ROOT colors are defined at https://root.cern.ch/doc/master/classTColor.html
  
  canvas->SaveAs("Example1.png"); // Save the result
}
