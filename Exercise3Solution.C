// To run this, type: cafe --stride 20 Example3.C

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

// Define the quasi-elastic formula for neutrino energy
double QEFormula(double Emu, double cosmu) // Muon energy and cosine of muon angle
{
  //Muon momentum
  const double pmu = sqrt(sqr(Emu) - sqr(M_MU)); // Use the relativity formula E^2 = p^2 + m^2
  // This is the neutrino-mode version of the formula. For antineutrino mode, swap neutron and proton masses.
  const double num = sqr(M_P) - sqr(M_N - E_B) - sqr(M_MU) + 2 * (M_N - E_B) * Emu;
  const double denom = 2 * (M_N - E_B - Emu + pmu * cosmu);
  return num/denom;
}


// This is the main function. To use ROOT's interpreted interface, you need to define a function
// with the same name as your file (minus the .C file extension)
void Exercise3Solution()
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

  // Define our own Variables


  // Conservation of energy from the true final-state particle energies
  const Var kConservedETrue([](const caf::SRProxy* sr)
                                {
                              const double Emu = sr->LepE; // LepE is (final-state) true lepton energy
                              const double protonKE=sr->eP; // kinetic energy
    // Final energy (proton and muon) - other initial particles' energy (bound stationary neutron)
                              return Emu + (protonKE + M_P) - (M_N - E_B);
  });
  
  // Conservation of energy for the reconstructed final-state particle energies
  const Var kConservedEReco([](const caf::SRProxy* sr)
                                {
                              const double Emu = sr->Elep_reco; // Final-state reconstructed lepton (muon) energy
                              const double protonKE=sr->eRecoP; // kinetic energy
                              // Final energy (proton and muon) - other initial particles' energy (bound stationary neutron)
                              return Emu + (protonKE + M_P) - (M_N - E_B);
  });
  
  
  // Reconstructed energy reported by the CAF.
  const Var kRecoE([](const caf::SRProxy* sr)
                   {
                     // As you know - if we can't understand the final state - we can't reconstruct neutrino energy!
                     if(std::isnan(sr->Ev_reco)) return 0.; // This line just deals with records where the energy reconstruction didn't work
                     return sr->Ev_reco;
                   });
  
  // The quasi-elastic formula was defined earlier.
  // The inputs this time are the true muon energy and cosine of the muon angle
  const Var kQEFormulaEnergy([](const caf::SRProxy* sr)
  {
    const double Emu = sr->Elep_reco;
    const double cosmu = cos(sr->theta_reco);
    // Sometimes it can't reconstruct the muon at all!
    // That's simulating a real detector where we sometimes won't be able to detect/identify a particle.
    //In that case, we'll just return 0.
    if(Emu == 0) return 0.;
    return QEFormula(Emu, cosmu);
  });

  // Define our axes: title, Binning, Variable
  const HistAxis axConservedETrue("E_#nu (conserve true energies) (GeV)", binsEnergy, kConservedETrue);
  const HistAxis axConservedEReco("E_#nu (conserve reco energies) (GeV)", binsEnergy, kConservedEReco);
  const HistAxis axEQE("E_#nu (QE formula) (GeV)", binsEnergy, kQEFormulaEnergy);
  const HistAxis axEReco("E_#nu reco (GeV)", binsEnergy, kRecoE);
  const HistAxis axETrue("True neutrino energy (GeV)", binsEnergy, kTrueEnergy);


  

  /* For exercise 3, we use the cut for a CCQE final state of 1 proton and 1 muon.
   */
   const Cut kHasQEFinalState([](const caf::SRProxy* sr)
                             {
             // This counts all the particles that aren't protons or muons: neutron, pi plus, pi minus, pi 0, positive kaon, negative kaon, neutral kaon, electromagnetic (gammas, electrons) and nuclear fragments. We want NONE of those!
                               const int totOthers = sr->nN + sr->nipip + sr->nipim + sr->nipi0 + sr->nikp + sr->nikm + sr->nik0 + sr->niem + sr->nNucleus;
             // pass if the lepton is a mu- (PDG code 13), number of protons is 1, and total other particles is zero
                               return abs(sr->LepPDG) == PDG_MU && sr->nP == 1 && totOthers == 0;
                             });
  // Now our cut's defined, we can make all of our Spectrum objects
  Spectrum sConservedETrue (loader, axConservedETrue, kHasQEFinalState);
  Spectrum sConservedEReco (loader, axConservedEReco, kHasQEFinalState);
  Spectrum sEQE (loader, axEQE, kHasQEFinalState);
  Spectrum sEReco (loader, axEReco, kHasQEFinalState);
  Spectrum sETrue (loader, axETrue, kHasQEFinalState);
  
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
  TH1D *hConservedEReco  =sConservedEReco .ToTH1(pot, kOrange-2);
  TH1D *hEQE =sEQE.ToTH1(pot, kOrange+7);
  TH1D *hEReco =sEReco.ToTH1(pot, kAzure-9);
  TH1D *hETrue =sETrue.ToTH1(pot, kGray+1);
  
  // These next lines will set the scale so nothing falls off the top
  double height= TMath::Max(hConservedETrue->GetMaximum(),hConservedEReco->GetMaximum());
  height=max(height,hEQE->GetMaximum());
  height=max(height,hEReco->GetMaximum());
  height=max(height,hETrue->GetMaximum());
  
  hConservedETrue->GetYaxis()->SetRangeUser(0,height * 1.1); // set the y axis range to 1.1 times the height
  hConservedETrue->GetXaxis()->SetTitle("Energy calculated various ways (GeV)");
  
  hConservedETrue->Draw("HIST");
  hConservedEReco->Draw("HIST SAME");
  hEQE->Draw("HIST SAME");
  hEReco->Draw("HIST SAME");
  hETrue->Draw("HIST SAME");
  
  gPad->SetLogy(false);
  
  auto legend = new TLegend(0.65,0.65,0.9,0.9); // x and y coordinates of corners
  legend->SetHeader("Legend","C"); // option "C" to center the header
  legend->AddEntry(hConservedETrue,"Energy cons. (true fs)","l");
  legend->AddEntry(hConservedEReco,"Energy cons. (reco fs)","l");
  legend->AddEntry(hEQE,"QE formula","l");
  legend->AddEntry(hEReco,"Reco from CAF","l");
  legend->AddEntry(hETrue,"True","l");
  legend->Draw();
  
  canvas->SaveAs("Exercise3.png"); // Save the result
}
