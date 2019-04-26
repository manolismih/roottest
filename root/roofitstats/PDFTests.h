// Author: Stephan Hageboeck, CERN  24 Apr 2019

/*****************************************************************************
 * RooFit
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2019, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOTTEST_ROOT_ROOFITSTATS_PDFTESTS_H_
#define ROOTTEST_ROOT_ROOFITSTATS_PDFTESTS_H_

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooGlobalFunc.h"

namespace TestPDFs {

struct TestConfig
{
  TestConfig(std::string name,
      std::shared_ptr<RooAbsPdf> pdf,
      std::initializer_list<RooAbsReal*> variables,
      std::initializer_list<RooAbsReal*> parameters,
      std::initializer_list<RooAbsReal*> liveObjects = std::initializer_list<RooAbsReal*>()) :
  _name(name),
  _pdf(pdf),
  _variables(std::make_shared<RooArgSet>()),
  _parameters(std::make_shared<RooArgSet>()),
  _origParameters(std::make_shared<RooArgSet>()),
  _otherObjects(std::make_shared<RooArgSet>())
  {
    for (auto var : variables) {
      _variables->addOwned(*var);
    }
    for (auto par : parameters) {
      _parameters->addOwned(*par);
    }
    for (auto con : liveObjects) {
      _otherObjects->addOwned(con);
    }

    _origParameters->addClone(*_parameters);
  }

  friend std::ostream& operator<<(std::ostream& os, const TestConfig& data) {
    os << std::setw(20) << std::left << data._name << "\n\tPDF = ";
    data._pdf->printStream(os, data._pdf->defaultPrintContents("I"),
        RooPrintable::kInline);
    os << "\tparameters=" << *data._parameters;
    return os;
  }


  std::string _name;
  std::shared_ptr<RooAbsPdf> _pdf;
  std::shared_ptr<RooArgSet> _variables;
  std::shared_ptr<RooArgSet> _parameters;
  std::shared_ptr<RooArgSet> _origParameters;
  std::shared_ptr<RooArgSet> _otherObjects;
};


void makeSimpleGaussian(std::vector<TestConfig>& configs) {
  // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
  auto x = new RooRealVar("x", "x", -10, 10);
  auto mean = new RooRealVar("mean", "mean of gaussian", 1, -10, 10);
  auto sigma = new RooRealVar("sigma", "width of gaussian", 1, 0.1, 10);

  // Build gaussian p.d.f in terms of x,mean and sigma
  auto gauss = std::make_shared<RooGaussian>("gauss", "gaussian PDF", *x, *mean, *sigma);

  configs.push_back(TestConfig("Simple Gaussian",
      gauss,
      {x},
      {mean, sigma}));
}

void makeGaussInMeanAndX(std::vector<TestConfig>& configs)
{
   // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
   auto x = new RooRealVar("x", "x", -100, 100);
   auto mean = new RooRealVar("mean", "mean of gaussian", 1, -10, 10);
   auto sigma = new RooRealVar("sigma", "width of gaussian", 1.337, 0.1, 10);

   // Build gaussian p.d.f in terms of x,mean and sigma
   auto gauss = std::make_shared<RooGaussian>("gauss", "gaussian PDF", *x, *mean, *sigma);

   configs.push_back(TestConfig("Gauss in mean and x",
       gauss,
       {x, mean},
       {sigma}));
}

void makeAddGaussPois(std::vector<TestConfig>& configs)
{
   // S e t u p   m o d e l
   // ---------------------

   // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
   auto x = new RooRealVar("x", "x", -100, 10000);
   auto mean = new RooRealVar("mean", "mean of gaussian", 1., -10, 10);
   auto sigma = new RooRealVar("sigma", "width of gaussian", 1.337, 0.1, 10);

   // Build gaussian p.d.f in terms of x,mean and sigma
   auto gauss = new RooGaussian("gauss", "gaussian PDF", *x, *mean, *sigma);

   auto meanPois = new RooRealVar("meanPois", "Mean of Poisson", 2., 0, 10);
   auto pois = new RooPoisson("Pois", "Poisson PDF", *x, *meanPois);

   auto fractionGaus = new RooRealVar("fractionGaus", "Fraction of Gauss component", 0.5, 0., 1.);
   auto sum = std::make_shared<RooAddPdf>("SumGausPois", "Sum of Gaus and Poisson", RooArgSet(*gauss, *pois), *fractionGaus);

   configs.push_back(TestConfig("Gauss + Poisson",
       sum,
       {x},
       {mean, sigma, meanPois, fractionGaus},
       {gauss, pois}));
}


void makeAdd2GaussExp(std::vector<TestConfig>& configs)
{
   auto x = new RooRealVar("x", "x", 0., 100.);

   auto c = new RooRealVar("c", "c", -2.3, -100., 0.);
   auto expo = new RooExponential("expo", "expo", *x, *c);


   auto mean = new RooRealVar("mean", "mean of gaussian", 5., -10, 100);
   auto sigma = new RooRealVar("sigma", "width of gaussian", 1.337, 0.1, 10);
   auto gauss = new RooGaussian("gauss", "gaussian PDF", *x, *mean, *sigma);

   auto mean2 = new RooRealVar("mean2", "mean of gaussian", 10., -10, 100);
   auto sigma2 = new RooRealVar("sigma2", "width of gaussian", 2.345, 0.1, 10);
   auto gauss2 = new RooGaussian("gauss2", "gaussian PDF", *x, *mean2, *sigma2);

   auto nGauss = new RooRealVar("nGauss", "Fraction of Gauss component", 100.);
   auto nGauss2 = new RooRealVar("nGauss2", "Fraction of Gauss component", 200.);
   auto nExp = new RooRealVar("nExp", "Number of events in exp", 1000);
   auto sum = std::make_shared<RooAddPdf>("Sum2GausExp", "Sum of Gaus and Exponentials",
       RooArgSet(std::initializer_list<RooAbsArg*>{gauss, gauss2, expo}),
       RooArgSet(std::initializer_list<RooAbsArg*>{nGauss, nGauss2, nExp}));

   configs.push_back(TestConfig("Gaus+Gaus+Exp",
       sum,
       {x},
       {c, mean, sigma, mean2, sigma2, nGauss, nGauss2, nExp},
       {expo, gauss, gauss2}));
}


void makeNestedExpGauss(std::vector<TestConfig>& configs)
{
   // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
   auto x = new RooRealVar("x", "x", 0.001, 20.);

   auto y1 = new RooRealVar("y1", "y1", -10, 0.);
   auto mean1 = new RooRealVar("mean1", "mean1", -3, -5, -0.1);
   auto sigma1 = new RooRealVar("sigma1", "sigma1", 0.4, 0.001, 5);
   auto c1 = new RooGaussian("c1", "c1", *y1, *mean1, *sigma1);

   auto y2 = new RooRealVar("y2", "y2", -10, 0.);
   auto mean2 = new RooRealVar("mean2", "mean2", -1.5, -5, -0.1);
   auto sigma2 = new RooRealVar("sigma2", "sigma2", 0.5, 0.001, 5);
   auto c2 = new RooGaussian("c2", "c2", *y2, *mean2, *sigma2);

   auto expo1 = new RooExponential("expo1", "expo1", *x, *c1);
   auto expo2 = new RooExponential("expo2", "expo2", *x, *c2);


   auto mean = new RooRealVar("mean", "mean of gaussian", 5., -10, 100);
   auto sigma = new RooRealVar("sigma", "width of gaussian", 1.337, 0.1, 10);
   auto gauss = new RooGaussian("gauss", "gaussian PDF", *x, *mean, *sigma);

   auto nGauss = new RooRealVar("nGauss", "Fraction of Gauss component", 100.);
   auto nExp1 = new RooRealVar("nExp1", "Number of events in exp1", 1000);
   auto nExp2 = new RooRealVar("nExp2", "Number of events in exp2", 2000);
   auto sum = std::make_shared<RooAddPdf>("SumGaus2Exp", "Sum of Gaus and 2 Exponentials",
       RooArgSet(std::initializer_list<RooAbsArg*>{gauss, expo1, expo2}),
       RooArgSet(std::initializer_list<RooAbsArg*>{nGauss, nExp1, nExp2}));

   configs.push_back(TestConfig("Gauss + Exp(x, c1(y1)) + Exp(x, c2(y2))",
       sum,
       {x, y1, y2},
       {mean1, sigma1, mean2, mean, sigma, nGauss, nExp1},
       {c1, c2, sigma2, expo1, expo2, gauss, nExp2}));
}



// POISSON
// --------------------------------------------------------
void makePoisson(std::vector<TestConfig>& configs) {
   auto x = new RooRealVar("x", "x", -10, 10);
   auto mean = new RooRealVar("mean", "Mean of Poisson", 2., 0., 10);
   auto pois = std::make_shared<RooPoisson>("Pois", "Poisson PDF", *x, *mean);

   configs.push_back(TestConfig("Poisson",
       pois,
       {x},
       {mean}));
}


// 2D Gauss
// --------------------------------------------------------
void make2DGauss(std::vector<TestConfig>& configs) {
   auto x = new RooRealVar("x", "x", 0, -5, 5);
   auto m1 = new RooRealVar("m1", "m1", 0 , -5, 5);
   auto s1 = new RooRealVar("s1", "s1", 1, 0.7, 1.3);
   auto y = new RooRealVar("y", "y", 0, -5., 5.);
   auto m2 = new RooRealVar("m2", "m2", 1.5 , -5., 5.);
   auto s2 = new RooRealVar("s2", "s2", 2, 0.1, 5);

   //Make a 2D PDF
   auto g1 = new RooGaussian("gaus1", "gaus1", *x, *m1, *s1);
   auto g2 = new RooGaussian("gaus2", "gaus2", *y, *m2, *s2);
   auto prod = std::make_shared<RooProdPdf>("prod", "prod", RooArgSet(*g1, *g2));

   configs.push_back(TestConfig("gaus1(x)*gaus2(y)",
       prod,
       {x, y},
       {m1, s1, m2, s2},
       {g1, g2}));
}

std::vector<TestConfig> makeTestPDFs() {
  std::vector<TestConfig> configs;

  makeSimpleGaussian(configs);
  makePoisson(configs);
  makeGaussInMeanAndX(configs);
  makeAddGaussPois(configs);
  makeAdd2GaussExp(configs);
  makeNestedExpGauss(configs);
  make2DGauss(configs);

  return configs;
}

}

#endif /* ROOTTEST_ROOT_ROOFITSTATS_PDFTESTS_H_ */
