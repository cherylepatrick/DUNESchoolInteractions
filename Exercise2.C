// To run this, type: cafe Exercise2.C

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


/* *****************
 Define some GENIE interaction modes.
 See full list at https://wiki.dunescience.org/wiki/Scattering_mode
 Use these to make your TRUTH cuts - this the interaction type GENIE simulated
 */

const int MODE_QE = 1;
const int MODE_RES = 4;
const int MODE_DIS = 3;
const int MODE_MEC = 10;

using namespace ana;


// This is the main function. To use ROOT's interpreted interface, you need to define a function
// with the same name as your file (minus the .C file extension)
void Exercise2Solution()
{
  // Define four options for our input CAF samples.
  // Environment variables and wildcards work, as do SAM datasets
  // (a metadata database Fermilab uses to organise large volumes of data and simulation files).
  
  const std::string NDGAR_FHC = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5gas/CAF_FHC_90*.root"; //ND-GAr FHC
  const std::string NDGAR_RHC = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5gas/CAF_RHC_90*.root"; //ND-GAr RHC
  const std::string NDLAR_FHC = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5/00/CAF_FHC_90*.root"; //ND-LAr FHC
  const std::string NDLAR_RHC = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5/00/CAF_RHC_90*.root"; //ND-LAr RHC
  
  // Source of events - load them from the one of the sets of files
  SpectrumLoader loader(NDGAR_FHC); // ***** Change this to use a different sample *****

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
  
  // This cut selects true QE interactions
  const Cut kIsQE = SIMPLEVAR(mode) == MODE_QE; //The modes are defined at the top
  /* ***************
  But for the sample we want:
   - Muon neutrinos
   - Not antineutrinos
   - QE interactions
   Define a cut that requires all of those on the line below
    */
  //const Cut kIsCCQE = ************;
  
  // Make your spectrum
  //Spectrum sTrueEQE(****************);


  /* ******* THIS CUT DEFINITION IS FOR THE SECOND PART OF EXERCISE 2 ***
   The CCQE final state is 1 proton and 1 muon. The code below uses the CAF
   identify this state.  As it is more complicated, instead of just defining it as
   a true/false boolean, we are going to have multiple lines of code.
   The little chunk of code here needs to return true if the event will pass
   the cut, or false if it will fail.
   Pass: 1 proton and 1 muon, no other particles
   The input for this function is a CAF "Standard record"
   */
   const Cut kHasQEFinalState([](const caf::SRProxy* sr)
                             {
             // This counts all the particles that aren't protons or muons: neutron, pi plus, pi minus, pi 0, positive kaon, negative kaon, neutral kaon, electromagnetic (gammas, electrons) and nuclear fragments. We want NONE of those!
                               const int totOthers = sr->nN + sr->nipip + sr->nipim + sr->nipi0 + sr->nikp + sr->nikm + sr->nik0 + sr->niem + sr->nNucleus;
             // pass if the lepton is a mu- (PDG code 13), number of protons is 1, and total other particles is zero
                               return sr->LepPDG == PDG_MU && sr->nP == 1 && totOthers == 0;
                             });
  
  
  
  // Fill all the Spectrum objects
  loader.Go();

  /* 
     Set to the same exposure as before
  */  
  const double pot = 1e20;

  // Convert and draw
  TCanvas *canvas = new TCanvas; // Make a canvas
  
  TH1D *hTrueEQE = sTrueEQE.ToTH1(pot, kBlue);// Turn spectrum to histogram and draw it
  hTrueEQE->Draw("HIST");
  //ROOT colors are defined at https://root.cern.ch/doc/master/classTColor.html


  canvas->SaveAs("Exercise2.png"); // Save the result
}
