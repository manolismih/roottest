#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooGenericPdf.h"
#include "RooGaussian.h"
#include "RooCategory.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooJohnson.h"
#include "RooProduct.h"
#include "RooFormulaVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooSimultaneous.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

using namespace RooFit;

void fit(){


  
  double x1 = 1, x2 = 16.5;
  RooRealVar mass( "mass", "mass", x1, x2 );

  RooCategory cat( "cat", "cat" );
  cat.defineType( "pos", 1 );
  cat.defineType( "neg", -1 );

  RooArgSet *obs = new RooArgSet();
  obs->add(mass); obs->add(cat);

  //MODEL
  //constants
  RooConstVar M_thr( "M_thr", "M_thr", 0);

  //Johnson & Gaussians means
  RooRealVar mean_j( "mean_j", "mean_j", 6.6, 5, 8 );
  RooRealVar mean_g1( "mean_g1", "mean_g1", 6.6, 5, 8 );
  RooRealVar mean_g2( "mean_g2", "mean_g2", 6.55, 5, 8 );
  RooRealVar mean_g3( "mean_g3", "mean_g3", 6.65, 5, 8 );
  RooRealVar mean_g1b( "mean_g1b", "mean_g1b", 6.6, 5, 8 );
  RooRealVar mean_g2b( "mean_g2b", "mean_g2b", 6.55, 5, 8 );
  RooRealVar mean_g3b( "mean_g3b", "mean_g3b", 6.65, 5, 8 );

  //shape parameters
  //Johnson
  RooRealVar gamma_j( "gamma_j", "gamma_j", -0.4, -1, 1 );
  RooRealVar delta_j( "delta_j", "delta_j", 0.8, 0, 10 );
  //empirical background
  RooRealVar b( "b", "b", 0.5, 0, 1 );
  RooRealVar c( "c", "c", 0.03, -0.2, 0.2 );

  //sigma parameters
  //Gaussians
  RooRealVar sigma_1( "sigma_1", "sigma_1", 0.35, 0.2, 0.4 );
  RooRealVar sigma_2( "sigma_2", "sigma_2", 0.8, 0.4, 1 );
  RooRealVar sigma_3( "sigma_3", "sigma_3", 0.1, 0, 0.2 );
  //Johnson
  RooRealVar sigma_j( "sigma_j", "sigma_j", 0.7, 0, 5 );

  //fractions
  //Johnson
  RooRealVar frac_j( "frac_j", "frac_j", 0.25, 0, 1 );
  //Gaussians
  RooRealVar frac_g1( "frac_g1", "frac_g1", 0.4, 0, 1 );
  RooRealVar frac_g2( "frac_g2", "frac_g2", 0.3, 0, 1 );

  RooRealVar asymm_sig( "asymm_sig", "asymm_sig", -0.02, -1, 1 );
  RooRealVar asymm_bkg( "asymm_bkg", "asymm_bkg", 0.05, -1, 1 );

  //signal model
  RooJohnson johnson("johnson", "Johnson pdf", mass, mean_j, sigma_j, gamma_j, delta_j, 0.);
  RooGaussian gauss_1("gauss_1", "gauss_1", mass, mean_g1, sigma_1);
  RooGaussian gauss_2("gauss_2", "gauss_2", mass, mean_g2, sigma_2);
  RooGaussian gauss_3("gauss_3", "gauss_3", mass, mean_g3, sigma_3);

//  RooAddPdf sig( "sig", "sig", RooArgList( johnson, gauss_1, gauss_2, gauss_3), RooArgList( frac_j, frac_g1, frac_g2 ), kTRUE );
  RooAddPdf sig( "sig", "sig", RooArgList( johnson), RooArgList(), kTRUE );

  RooFormulaVar f( "f", "f", "0.5*(1+(cat*asymm_sig))", RooArgSet(asymm_sig, cat) );

  //background model
  //empirical
  RooGenericPdf bkg( "bkg", "TMath::Power(mass,b)*TMath::Exp(-c*mass)", RooArgList(mass, b, c) );

  RooFormulaVar bkgtag( "bkgtag", "bkgtag", "0.5*(1+@0*@1)", RooArgSet( asymm_bkg, cat ) );

  RooRealVar N_sig( "N_sig", "N_sig", 30E6, 0, 1e9 );
  RooRealVar N_bkg( "N_bkg", "N_bkg", 12E6, 0, 1e9 );
  
  RooProduct N_sig_f("N_sig_f", "N_sig * f", RooArgSet(N_sig, f));
  RooProduct N_bkg_bkgtag("N_bkg_f", "N_bkg * bkgtag", RooArgSet(N_sig, bkgtag));

  RooAddPdf pdf_tot( "pdf_tot", "pdf_tot", RooArgSet( sig, bkg), RooArgSet( N_sig_f, N_bkg_bkgtag ) );

  RooSimultaneous sim("simPdf", "simPdf", RooArgList(pdf_tot, pdf_tot), cat);
//  sim.addPdf(pdf_tot, "pos");
//  sim.addPdf(pdf_tot, "neg");
  sim.Print("T");

  TFile * in_file = new TFile("LHCb002/data.root");
  TH2D * h_mass_plus;
  TH2D * h_mass_minus;
  h_mass_plus = (TH2D * ) in_file->Get( TString("h_pos") );
  h_mass_minus = (TH2D * ) in_file->Get( TString("h_neg") );
  RooDataHist *data_h = new RooDataHist( "data_h", "data_h", mass, Index(cat),
      Import( "pos", *h_mass_plus ),
      Import( "neg", *h_mass_minus ));

  data_h->Print();

//  RooChi2Var chi2( "chi2", "chi2", pdf_tot, *data_h, Extended(kTRUE) );
//  RooMinuit m1(chi2) ;
////  m1.setVerbose(kTRUE);
////  m1.setPrintLevel(3);
//  m1.setEps(1e-16);
//  m1.setStrategy(2);
//
//  m1.migrad();
//  m1.hesse();

  mass.setRange("sideband", 8., 16.);

//  auto fitResult = pdf_tot.chi2FitTo(*data_h, Range("sideband"), Extended(true), Strategy(2), Save());
  auto fitResult = sim.fitTo(*data_h, Extended(true), /*Range("sideband"),*/ Strategy(2), Save());

  RooArgSet* params = sim.getParameters(*obs) ;
  params->writeToFile( "params.txt" );
  

  TCanvas *c_pos = new TCanvas("c_pos", "c_pos", 900, 900);
  c_pos->cd();
  RooPlot * plot_pos = mass.frame( x1, x2);
  data_h->plotOn( plot_pos, Cut("cat==cat::pos"));
  const double ndata = data_h->sumEntries("cat==cat::pos");
  sim.plotOn( plot_pos, Slice( cat,"pos" ), ProjWData( cat, *data_h ),
//      Normalization(ndata, RooAbsReal::NumEvent),
      LineColor(kBlue) );
  plot_pos->Draw();

//  TCanvas *c_neg = new TCanvas("c_neg", "c_neg", 900, 900);
//  c_neg->cd();
//  RooPlot * plot_neg = mass.frame( x1, x2);
//  data_h->plotOn( plot_neg, Cut("cat==cat::neg"));
//  const double ndataNeg = data_h->sumEntries("cat==cat::neg");
//  pdf_tot.plotOn( plot_neg, Slice( cat,"neg" ), ProjWData( cat, *data_h ),
////      Normalization(ndataNeg, RooAbsReal::NumEvent),
//      LineColor(kBlue) );
//  plot_neg->Draw();
}
