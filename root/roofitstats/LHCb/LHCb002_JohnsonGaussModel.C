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
#include "TH1D.h"
#include "TH2D.h"

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
  RooRealVar c( "c", "c", 0.03, -0.1, 0.1 );

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
  RooRealVar frac_g1( "frac_g1", "frac_g1", 0.7, 0, 1 );
  RooRealVar frac_g2( "frac_g2", "frac_g2", 0.4, 0, 1 );

  RooRealVar asymm_sig( "asymm_sig", "asymm_sig", -0.02, -1, 1 );
  RooRealVar asymm_bkg( "asymm_bkg", "asymm_bkg", 0.05, -1, 1 );

  //signal model
  //Johnson
  RooAbsPdf *johnson = new RooGenericPdf( "johnson", "(mass>M_thr)*delta_j/(sigma_j*TMath::Sqrt(TMath::Pi()))*TMath::Exp(-0.5*(gamma_j+delta_j*TMath::ASinH((mass-mean_j)/sigma_j))*(gamma_j+delta_j*TMath::ASinH((mass-mean_j)/sigma_j)))/TMath::Sqrt(1+(mass-mean_j)*(mass-mean_j)/(sigma_j*sigma_j))", RooArgSet( mass, mean_j, sigma_j, gamma_j, delta_j, M_thr ) );
  RooAbsPdf* gauss_1_pos = new RooGenericPdf( "gauss_1_pos", "(mass>M_thr)*1./(TMath::Sqrt(2*sigma_1*sigma_1*TMath::Pi()))*(TMath::Exp(-1*(mass-mean_g1)*(mass-mean_g1)/(2*sigma_1*sigma_1)))", RooArgSet( mass, mean_g1, sigma_1, M_thr ) );
  RooAbsPdf* gauss_2_pos = new RooGenericPdf( "gauss_2_pos", "(mass>M_thr)*1./(TMath::Sqrt(2*sigma_2*sigma_2*TMath::Pi()))*(TMath::Exp(-1*(mass-mean_g2)*(mass-mean_g2)/(2*sigma_2*sigma_2)))", RooArgSet( mass, mean_g2, sigma_2, M_thr ) );
  RooAbsPdf* gauss_3_pos = new RooGenericPdf( "gauss_3_pos", "(mass>M_thr)*1./(TMath::Sqrt(2*sigma_3*sigma_3*TMath::Pi()))*(TMath::Exp(-1*(mass-mean_g3)*(mass-mean_g3)/(2*sigma_3*sigma_3)))", RooArgSet( mass, mean_g3, sigma_3, M_thr ) );
  RooAbsPdf* gauss_1_neg = new RooGenericPdf( "gauss_1_neg", "(mass>M_thr)*1./(TMath::Sqrt(2*sigma_1*sigma_1*TMath::Pi()))*(TMath::Exp(-1*(mass-mean_g1b)*(mass-mean_g1b)/(2*sigma_1*sigma_1)))", RooArgSet( mass, mean_g1b, sigma_1, M_thr ) );
  RooAbsPdf* gauss_2_neg = new RooGenericPdf( "gauss_2_neg", "(mass>M_thr)*1./(TMath::Sqrt(2*sigma_2*sigma_2*TMath::Pi()))*(TMath::Exp(-1*(mass-mean_g2b)*(mass-mean_g2b)/(2*sigma_2*sigma_2)))", RooArgSet( mass, mean_g2b, sigma_2, M_thr ) );
  RooAbsPdf* gauss_3_neg = new RooGenericPdf( "gauss_3_neg", "(mass>M_thr)*1./(TMath::Sqrt(2*sigma_3*sigma_3*TMath::Pi()))*(TMath::Exp(-1*(mass-mean_g3b)*(mass-mean_g3b)/(2*sigma_3*sigma_3)))", RooArgSet( mass, mean_g3b, sigma_3, M_thr ) );

  RooAddPdf sig_pos( "sig_pos", "sig_pos", RooArgList( *johnson, *gauss_1_pos, *gauss_2_pos, *gauss_3_pos), RooArgList( frac_j, frac_g1, frac_g2 ), kTRUE );
  RooAddPdf sig_neg( "sig_neg", "sig_neg", RooArgList( *johnson, *gauss_1_neg, *gauss_2_neg, *gauss_3_neg), RooArgList( frac_j, frac_g1, frac_g2 ), kTRUE );

  RooGenericPdf tag_pos( "tag_pos", "tag_pos", "@0==1", RooArgSet(cat) );
  RooGenericPdf tag_neg( "tag_neg", "tag_neg", "@0==-1", RooArgSet(cat) );

  RooProdPdf pdf_pos( "pdf_pos", "pdf_pos", RooArgSet( tag_pos, sig_pos ) );
  RooProdPdf pdf_neg("pdf_neg", "pdf_neg", RooArgSet( tag_neg, sig_neg ) );

  RooFormulaVar f_pos( "f_pos", "f_pos", "0.5*(1+@0)", RooArgSet(asymm_sig) );
  RooFormulaVar f_neg( "f_neg", "f_neg", "0.5*(1-@0)", RooArgSet(asymm_sig) );

  RooAddPdf pdf_sig( "pdf_sig", "pdf_sig", RooArgSet( pdf_pos, pdf_neg ), RooArgSet( f_pos, f_neg ) );

  //background model
  //empirical
  RooAbsPdf* bkg = new RooGenericPdf( "bkg", "(mass>M_thr)*(TMath::Power((mass-M_thr),b)*TMath::Exp(-c*(mass-M_thr)))", RooArgList( mass, b, c, M_thr ) );

  RooGenericPdf bkgtagpdf( "bkgtagpdf", "bkgtagpdf", "0.5*(1+@0*@1)", RooArgSet( asymm_bkg, cat ) );
  RooProdPdf pdf_bkg( "pdf_bkg", "pdf_bkg", RooArgSet( bkgtagpdf, *bkg) );

  RooRealVar N_sig( "N_sig", "N_sig", 30000000, 0, 1e9 );
  RooRealVar N_bkg( "N_bkg", "N_bkg", 12000000, 0, 1e9 );  
  
  RooAddPdf pdf_tot( "pdf_tot", "pdf_tot", RooArgSet( pdf_sig, pdf_bkg), RooArgSet( N_sig, N_bkg ) );

  TFile * in_file = new TFile("LHCb002/data.root");
  TH2D * h_mass_plus;
  TH2D * h_mass_minus;
  h_mass_plus = (TH2D * ) in_file->Get( TString("h_pos") );
  h_mass_minus = (TH2D * ) in_file->Get( TString("h_neg") );
  RooDataHist *data_h = new RooDataHist( "data_h", "data_h", mass, Index(cat), Import( "pos", *h_mass_plus ),Import( "neg", *h_mass_minus ) );

  data_h->Print();

  RooChi2Var chi2( "chi2", "chi2", pdf_tot, *data_h, Extended(kTRUE) );
  RooMinuit m1(chi2) ;
//  m1.setVerbose(kTRUE);
//  m1.setPrintLevel(3);
  m1.setEps(1e-16);
  m1.setStrategy(2);

  m1.migrad();
  m1.hesse();

  RooArgSet* params = pdf_tot.getParameters(*obs) ;
  params->writeToFile( "params.txt" );
  

  TCanvas *c_pos = new TCanvas("c_pos", "c_pos", 900, 900);
  c_pos->cd();
  RooPlot * plot_pos = mass.frame( x1, x2);
  data_h->plotOn( plot_pos, Cut("cat==cat::pos"));
  pdf_tot.plotOn( plot_pos, Slice( cat,"pos" ), ProjWData( cat, *data_h ), LineColor(kBlue) );
  plot_pos->Draw();

  TCanvas *c_neg = new TCanvas("c_neg", "c_neg", 900, 900);
  c_neg->cd();
  RooPlot * plot_neg = mass.frame( x1, x2);
  data_h->plotOn( plot_neg, Cut("cat==cat::neg"));
  pdf_tot.plotOn( plot_neg, Slice( cat,"neg" ), ProjWData( cat, *data_h ), LineColor(kBlue) );
  plot_neg->Draw();
}
