// To run this, type: cafe --stride 100 firstAnalysisSolution.C

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
#include "CAFAna/Vars/Vars.h" // for kTrueEnergy
#include "CAFAna/Cuts/TruthCuts.h" // kNumuCC and friends

#include "StandardRecord/SRProxy.h"

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
using util::sqr; // Square


// Define some physical constants
const double M_P = .938; // Proton mass in GeV
const double M_N = .939; // Neutron mass in GeV
const double M_MU = .106; // Muon mass in GeV
const double E_B = .028; // Binding energy for nucleons in argon-40 in GeV

// Define some GENIE interaction modes.
// See full list at https://wiki.dunescience.org/wiki/Scattering_mode
const int MODE_QE = 1;
const int MODE_RES = 4;
const int MODE_DIS = 3;
const int MODE_MEC = 10;

// Define some PDG codes (particle identifiers from https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf)
// You can also find the ID's for things like protons, neutrons, pions and even whole nuclei in the list!
const int PDG_MU=13;
const int PDG_E=11;
const int PDG_NUMU=14;
const int PDG_NUE=12;

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
void firstAnalysisSolution()
{
  // This is points to the source files we will use. These ones are for the ND-GAr.
  // Environment variables and wildcards work. As do SAM datasets 
  // (a metadata database Fermilab uses to organise large volumes of data and simulation files).
  const std::string fname = "/pnfs/dune/persistent/users/marshalc/CAF/CAFv5gas/CAF_FHC_*.root";

  // Source of events - load them from the files
  SpectrumLoader loader(fname);

  // We want to plot a histogram with 40 bins, covering the range 0 to 10 GeV
  const Binning binsEnergy = Binning::Simple(40, 0, 10);

  // Define the label, binning, and contents that we want for our first histogram
  // the axis label can be whatever you like. 
  // The binning needs to be a Binning object like the one we just made.
  // The Variable that we are plotting can be a single variable or function of variables from
  // the CAFs. The choices are in the file xxxx
  const HistAxis axTrue("True neutrino energy (GeV)", binsEnergy, kTrueEnergy);

  // Cheryl: These cut names aren't very intuitive, you just have to read them
  // out of Vars.h. Would it be better to have them define the cuts themselves?
  // But then in a real macro I'd want them to use the standard ones...

  // Now we have defined an axis and a variable we want to plot on it, let's decide which
  // events to plot. So here we are telling the loader to load a spectrum with our defined
  // axis, for each of a few different selections

  // Here, we are defining the selections or "cuts" using the variable names in CAFAna

  Spectrum sTrueENumu(loader, axTrue, kIsNumuCC && !kIsAntiNu);  
  // kIsNumuCC = Muon neutrino charged-current interactions.
  // kIsAntiNu = Interaction initiated by an antineutrino. !kIsAntiNu means it is NOT an antineutrino.
  // && means we want both the first AND the second condition to be true. If we wanted the first OR the second condition, we would use ||

  Spectrum sTrueENumubar(loader, axTrue, kIsNumuCC && kIsAntiNu); // Muon antineutrino charged-current interactions

  Spectrum sTrueENue(loader, axTrue, kIsBeamNue && !kIsAntiNu); 
  // For some reason the electron neutrino naming convention is different 
  // Unfortunately that kind of thing is something you need to get used to in physics - many of the file formats have been developed in
  // non-linear ways, by many different people, so not everything is as consistent as you might like.
  Spectrum sTrueENuebar(loader, axTrue, kIsBeamNue && kIsAntiNu);

  // Here, we are defining our own cuts from combinations of CAFAna variables

  // First, let's select a true interaction type
  // True quasi-elastic events

  const Cut kIsQE = SIMPLEVAR(mode) == 1; // mode 1 is CCQE
  // TODO is there an enum we should be using?
  // We can combine the cuts!
  const Cut kIsCCQE = kIsNumuCC && kIsQE; // Muon charge-current events that are QE

  // Make another spectrum as before: loader, axis to plot (that's the bins and variable), cut - but this time we're using the cut we defined ourselves
  Spectrum sTrueEQE(loader, axTrue, kIsCCQE);

  // This cut is more complicated, so instead of just defining it as a true/false boolean, 
  // we are going to have multiple lines of code. The little chunk of code here needs to
  // return true if the event will pass the cut, or false if it will fail
  // Pass: 1 proton and 1 muon, no other particles
  // The input for this function is a CAF "Standard record"
  const Cut kHasQEFinalState([](const caf::SRProxy* sr)
                             {
			       // This counts all the particles that aren't protons or muons: neutron, pi plus, pi minus, pi 0, positive kaon, negative kaon, neutral kaon, electromagnetic (gammas, electrons) and nuclear fragments. We want NONE of those!
                               const int totOthers = sr->nN + sr->nipip + sr->nipim + sr->nipi0 + sr->nikp + sr->nikm + sr->nik0 + sr->niem + sr->nNucleus;
			       // pass if the lepton is a mu- (PDG code 13), number of protons is 1, and total other particles is zero
                               return abs(sr->LepPDG) == PDG_MU && sr->nP == 1 && totOthers == 0;
                             });

  // Now our cut's defined, we can make a spectrum as before:
  Spectrum sTrueEQEfs(loader, axTrue, kHasQEFinalState);

  // This time, we are looking for CC0pi - one negative muon, at least one proton, and no pions
  // Define the cut...
  const Cut kHasCC0PiFinalState([](const caf::SRProxy* sr)
                                {
                                  const int totPi = sr->nipip + sr->nipim + sr->nipi0;
                                  return abs(sr->LepPDG) == PDG_MU && sr->nP >= 1 && totPi == 0;
                                });
  // ...and make the spectrum.
  Spectrum sTrueE0pifs(loader, axTrue, kHasCC0PiFinalState);

  // Define cuts for the other interaction modes, based on how the CCQE one was defined earlier
  const Cut kIsCCMEC = SIMPLEVAR(mode) == 10 && kIsNumuCC;
  const Cut kIsCCRES = SIMPLEVAR(mode) == 4  && kIsNumuCC;
  const Cut kIsCCDIS = SIMPLEVAR(mode) == 3  && kIsNumuCC;

  // And a spectrum for each one... remember we can COMBINE two cuts here, so we are looking
  // in each case for CC0pi final states AND a particular interaction mode
  Spectrum sTrueQETrueE0pifs (loader, axTrue, kHasCC0PiFinalState && kIsCCQE);
  Spectrum sTrueMECTrueE0pifs(loader, axTrue, kHasCC0PiFinalState && kIsCCMEC);
  Spectrum sTrueRESTrueE0pifs(loader, axTrue, kHasCC0PiFinalState && kIsCCRES);
  Spectrum sTrueDISTrueE0pifs(loader, axTrue, kHasCC0PiFinalState && kIsCCDIS);


  // Now we will see what happens if, instead of the TRUE energy, we plot the energy
  // we would get if we used the quasi-elastic formula on the muon's true kinematics

  // The quasi-elastic formula was defined earlier. 
  // The inputs this time are the true muon energy and cosine of the muon angle
  const Var kQEFormulaEnergy([](const caf::SRProxy* sr)
                             {
                               const double Emu = sr->LepE; // LepE is (final-state) true lepton energy
                               const double cosmu = cos(sr->LepNuAngle); // LepNuAngle is the angle between neutrino and final-state lepton

                               return QEFormula(Emu, cosmu);
                             });

  // Define an new axis: this time we will plot this QE formula energy rather than the true one 
  // (but keep using the bins we defined earlier)
  const HistAxis axQEFormula("Quasielastic formula energy (GeV)", binsEnergy, kQEFormulaEnergy);

  // Set up a spectrum to see how this QE formula looks for our CCQE events with a CC0pi final state
  Spectrum sTrueQEQEFormulaE0pifs(loader, axQEFormula, kHasCC0PiFinalState && kIsCCQE);

  // We used the true muon kinematics in the QE formula last time
  //  what if we used the reconstructed ones instead? Instead of being the info about
  // the real (simulated) muon, this is what the simulation thinks the detector might
  // have seen. Sometimes the detector will make a mistake, and it will not be able to 
  // reproduce the true energy and angle precisely.

  const Var kRecoQEFormulaEnergy([](const caf::SRProxy* sr)
                                 {
                                   const double Emu = sr->Elep_reco;
                                   const double cosmu = cos(sr->theta_reco);

                                   // Sometimes it can't reconstruct the muon at all!
				   // That's simulating a real detector where we sometimes won't be able to detect/identify a particle.
				   //In that case, we'll just return 0.
                                   if(Emu == 0) return 0.;
                                   
                                   return QEFormula(Emu, cosmu);
                                 });
  // We're used to this now...
  const HistAxis axRecoQEFormula("Reconstructed quasielastic formula energy (GeV)", binsEnergy, kRecoQEFormulaEnergy);
  Spectrum sTrueQERecoQEFormulaE0pifs(loader, axRecoQEFormula, kHasCC0PiFinalState && kIsCCQE);

  // Now compare this to the reconstructed energy reported by the CAF. Is it the same?
  const Var kRecoE([](const caf::SRProxy* sr)
                   {
                     // As you know - if we can't understand the final state - we can't reconstruct neutrino energy!
                     if(std::isnan(sr->Ev_reco)) return 0.; // This line just deals with records where the energy reconstruction didn't work
                     return sr->Ev_reco;
                   });

  const HistAxis axRecoE("Reconstructed energy (GeV)", binsEnergy, kRecoE);
  Spectrum sTrueQERecoE0pifs(loader, axRecoE, kHasCC0PiFinalState && kIsCCQE);

  // Removing the true QE requirement
  // Because we already defined all the axes and cuts, it's easy to mix and match them

  //  Spectrum sTrueE0pifs (loader, axTrue, kHasCC0PiFinalState); // this one is already defined above
  Spectrum sQEFormulaE0pifs(loader, axQEFormula, kHasCC0PiFinalState);
  Spectrum sRecoQEFormulaE0pifs(loader, axRecoQEFormula, kHasCC0PiFinalState);
  Spectrum sRecoE0pifs(loader, axRecoE, kHasCC0PiFinalState);

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
  sTrueENumu.ToTH1(pot, kBlue)->Draw("hist"); // Muon neutrinos in blue. ROOT colors are defined at https://root.cern.ch/doc/master/classTColor.html
  sTrueENumubar.ToTH1(pot, kBlue, 7)->Draw("hist same"); //Antineutrinos are getting a dashed line. 
  //"SAME" means it will be drawn on the same axis as the previous spectrum
  sTrueENue.ToTH1(pot, kRed)->Draw("hist same"); //Electron neutrino
  sTrueENuebar.ToTH1(pot, kRed, 7)->Draw("hist same"); //Electron antineutrino
  gPad->SetLogy();

  new TCanvas;
  sTrueEQE.ToTH1(pot)->Draw("hist");
  sTrueEQEfs.ToTH1(pot, kRed)->Draw("hist same");
  sTrueE0pifs.ToTH1(pot, kBlue)->Draw("hist same");

  new TCanvas;
  // Could potentially use THStack here, but I don't understand/trust it, so I
  // am just going to sum manually. Possibly we can put a stacked plots helper
  // into Analysis/Plots.h
  Spectrum sQE = sTrueQETrueE0pifs;
  Spectrum sQEMEC = sQE + sTrueMECTrueE0pifs;
  Spectrum sQEMECRES = sQEMEC + sTrueRESTrueE0pifs;
  Spectrum sQEMECRESDIS = sQEMECRES + sTrueDISTrueE0pifs;

  TH1* hQE  = sQE.ToTH1(pot);
  TH1* hQEMEC = sQEMEC.ToTH1(pot);
  TH1* hQEMECRES = sQEMECRES.ToTH1(pot);
  TH1* hQEMECRESDIS = sQEMECRESDIS.ToTH1(pot);

  hQE->SetFillColor(kRed);
  hQEMEC->SetFillColor(kMagenta);
  hQEMECRES->SetFillColor(kBlue);
  hQEMECRESDIS->SetFillColor(kGreen+2);

  sTrueE0pifs.ToTH1(pot)->Draw("hist"); // total
  hQEMECRESDIS->Draw("hist same");
  hQEMECRES->Draw("hist same");
  hQEMEC->Draw("hist same");
  hQE->Draw("hist same");

  new TCanvas;
  // Cheryl: this plot winds up saying "true energy" on the x-axis. Do we want
  // to make a new spectrum to say "various energy estimators", or get the TH1
  // and relabel the axis, or ignore it?
  sTrueQETrueE0pifs.ToTH1(pot)->Draw("hist");
  sTrueQEQEFormulaE0pifs.ToTH1(pot, kRed)->Draw("hist same");
  sTrueQERecoQEFormulaE0pifs.ToTH1(pot, kBlue)->Draw("hist same");
  sTrueQERecoE0pifs.ToTH1(pot, kMagenta)->Draw("hist same");

  new TCanvas;
  sTrueE0pifs.ToTH1(pot)->Draw("hist");
  sQEFormulaE0pifs.ToTH1(pot, kRed)->Draw("hist same");
  sRecoQEFormulaE0pifs.ToTH1(pot, kBlue)->Draw("hist same");
  sRecoE0pifs.ToTH1(pot, kMagenta)->Draw("hist same");
}
