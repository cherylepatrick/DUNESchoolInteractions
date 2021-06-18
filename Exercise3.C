// To run this, type: cafe Exercise3.C

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

// Standard C++ library for input and output
#include <iostream>

/* *********
 Define some physical constants
 */
const double M_P = .938; // Proton mass in GeV
const double M_N = .939; // Neutron mass in GeV
const double M_MU = .106; // Muon mass in GeV
const double E_B = .028; // Binding energy for nucleons in argon-40 in GeV


/*
 Define some GENIE interaction modes.
 See full list at https://wiki.dunescience.org/wiki/Scattering_mode
 Use these to make your TRUTH cuts - this the interaction type GENIE simulated
 */

const int MODE_QE = 1;
const int MODE_RES = 4;
const int MODE_DIS = 3;
const int MODE_MEC = 10;

/*
 Define some PDG codes (particle identifiers from https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf)
 You can also find the ID's for things like protons, neutrons, pions and even whole nuclei in the list!
 */
const int PDG_MU=13;
const int PDG_E=11;
const int PDG_NUMU=14;
const int PDG_NUE=12;

using namespace ana;
using util::sqr; // Square

// Define the quasi-elastic formula for neutrino energy. Feel free to use it or code it yourself
double QEFormula(double Emu, double cosmu) // Muon energy and cosine of muon angle
{
  //Muon momentum
  const double pmu = sqrt(sqr(Emu) - sqr(M_MU)); // Use the relativity formula E^2 = p^2 + m^2
  // This is the neutrino-mode version of the formula. For antineutrino mode, you'd swap neutron and proton masses.
  const double num = sqr(M_P) - sqr(M_N - E_B) - sqr(M_MU) + 2 * (M_N - E_B) * Emu;
  const double denom = 2 * (M_N - E_B - Emu + pmu * cosmu);
  return num/denom;
}


// This is the main function. To use ROOT's interpreted interface, you need to define a function
// with the same name as your file (minus the .C file extension)
void Exercise3()
{
  // Define four options for our input CAF samples.
  // Environment variables and wildcards work, as do SAM datasets
  // (a metadata database Fermilab uses to organise large volumes of data and simulation files).
  // Change 90* to 9* for 10 times as many files!
  const std::string NDGAR_FHC = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5gas/CAF_FHC_90*.root"; //ND-GAr FHC
  const std::string NDGAR_RHC = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5gas/CAF_RHC_90*.root"; //ND-GAr RHC
  const std::string NDLAR_FHC = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5/00/CAF_FHC_90*.root"; //ND-LAr FHC
  const std::string NDLAR_RHC = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5/00/CAF_RHC_90*.root"; //ND-LAr RHC
  
  // Source of events - load them from the one of the sets of files
  SpectrumLoader loader(NDGAR_FHC); // ***** Change this to use a different sample ***

  // We want to plot a histogram with 40 bins, covering the range 0 to 10 GeV
  const Binning binsEnergy = Binning::Simple(40, 0, 10);

  // ***************   Define our own Variables
  // Conservation of energy from the true final-state particle energies
  // You can find the standard record variable names at
  // https://wiki.dunescience.org/wiki/CAF_ntuple_format
  const Var kConservedETrue([](const caf::SRProxy* sr)
                                {
                              const double Emu = sr->*****; // Enter the variable name for (final-state) true lepton energy
                              const double protonKE=sr->*******; // Enter the variable name for the true proton kinetic energy. We only have 1 proton for this example
                              double nuEnergy=*******;
    /* ***** Put in the formula for the neutrino energy, using conservation of energy
      You'll need the two things you just calculated plus some particle masses
      and the binding energy (defined at the top of the file */
                              return nuEnergy;
  });
  
  
 /*
  // Reconstructed neutrino energy reported by the CAF. Replace the **** with the right variable name. Note the trick to remove any that couldn't reconstruct it!
  const Var kRecoE([](const caf::SRProxy* sr)
                   {
                     // As you know - if we can't understand the final state - we can't reconstruct neutrino energy!
                     if(std::isnan(sr->******)) return 0.; // This line just deals with records where the energy reconstruction didn't work
                     return sr->******;
                   });
  */

  // Define our axes: title, Binning, Variable
  const HistAxis axConservedETrue("E_#nu (conserve true energies) (GeV)", binsEnergy, kConservedETrue);
  // ***** You'll be adding more of these!


  

  /* For exercise 3, we use the cut for a CCQE final state of 1 proton and 1 muon.
   */
   const Cut kHasQEFinalState([](const caf::SRProxy* sr)
                             {
                 // I'm cutting events where the reconstructed neutrino energy or muon energy is zero or not a number, to make this easier to interpret!
              if(std::isnan(sr->Ev_reco)) return false;
              if(sr->Ev_reco<=0.) return false;
              if(sr->Elep_reco<=0.) return false;
     
             // This counts all the particles that aren't protons or muons: neutron, pi plus, pi minus, pi 0, positive kaon, negative kaon, neutral kaon, electromagnetic (gammas, electrons) and nuclear fragments. We want NONE of those!
                               const int totOthers = sr->nN + sr->nipip + sr->nipim + sr->nipi0 + sr->nikp + sr->nikm + sr->nik0 + sr->niem + sr->nNucleus;
             // pass if the lepton is a mu- (PDG code 13), number of protons is 1, and total other particles is zero

                               return sr->LepPDG == PDG_MU && sr->nP == 1 && totOthers == 0 ;
                             });
  
  // Now our cut's defined, we can make all of our Spectrum objects
  
  // ***** You'll be adding more of these!
  Spectrum sConservedETrue (loader, axConservedETrue, kHasQEFinalState);

  
  // Fill all the Spectrum objects
  loader.Go();

  /* 
     Set to the same exposure as before
  */  
  const double pot = 1e20;

  // Convert and draw
  TCanvas *canvas = new TCanvas; // Make a canvas
  
  // Convert Spectrum objects to histograms
  // ROOT colors are defined at https://root.cern.ch/doc/master/classTColor.html
  TH1D *hConservedETrue = sConservedETrue.ToTH1(pot, kAzure-7);
  
  // ***** You'll be adding more of these!
  
  // These next lines will set the scale so nothing falls off the top
  // You might want to uncomment and adapt them!
  double height=hConservedETrue->GetMaximum();
  
  /*height= TMath::Max(height,hConservedEReco->GetMaximum());
  height=max(height,hEQE->GetMaximum());
  height=max(height,hEReco->GetMaximum());
  height=max(height,hETrue->GetMaximum());*/
  
  hConservedETrue->GetYaxis()->SetRangeUser(0,height * 1.1); // set the y axis range to 1.1 times the height
  // Reformatting this cos the first one you draw sets the title for everything!
  hConservedETrue->GetXaxis()->SetTitle("Energy calculated various ways (GeV)");

  
  hConservedETrue->Draw("HIST");  // ***** You'll be adding more of these!
  
  gPad->SetLogy(false);
  
  auto legend = new TLegend(0.65,0.65,0.9,0.9); // x and y coordinates of corners
  legend->SetHeader("Reconstruction method","C"); // option "C" to center the header
  legend->AddEntry(hConservedETrue,"Energy cons. (true fs)","l");
  // ***** You'll be adding more of these!
  legend->Draw();
  
  canvas->SaveAs("Exercise3.png"); // Save the result
}
