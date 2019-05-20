// Author: Stephan Hageboeck, CERN  26 Apr 2019

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

#include "VectorisedPDFTests.h"

#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooExponential.h"

class TestGaussPlusPoisson : public PDFTest
{
  protected:
    TestGaussPlusPoisson() :
      PDFTest("Gauss + Poisson", 200000)
    {
      // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
      auto x = new RooRealVar("x", "x", -1.5, 40.5);
      x->setBins(42);//Prettier plots for Poisson

      auto mean = new RooRealVar("mean", "mean of gaussian", 20., -10, 30);
      auto sigma = new RooRealVar("sigma", "width of gaussian", 4., 0.1, 10);

      // Build gaussian p.d.f in terms of x,mean and sigma
      auto gauss = new RooGaussian("gauss", "gaussian PDF", *x, *mean, *sigma);

      auto meanPois = new RooRealVar("meanPois", "Mean of Poisson", 7, 0, 30);
      auto pois = new RooPoisson("Pois", "Poisson PDF", *x, *meanPois, true);

      auto fractionGaus = new RooRealVar("fractionGaus", "Fraction of Gauss component", 0.5, 0., 1.);
      _pdf = std::make_unique<RooAddPdf>("SumGausPois", "Sum of Gaus and Poisson",
          RooArgSet(*gauss, *pois), *fractionGaus);

      _variables.addOwned(x);

//      _variablesToPlot.add(x);

      for (auto par : {mean, sigma, meanPois, fractionGaus}) {
        _parameters.addOwned(par);
      }

      for (auto obj : std::initializer_list<RooAbsPdf*>{gauss, pois}) {
        _otherObjects.addOwned(obj);
      }
    }
};

RUN_SCALAR(TestGaussPlusPoisson, Scalar)
RUN_BATCH(TestGaussPlusPoisson, Batch)
RUN_BATCH_VS_SCALAR(TestGaussPlusPoisson, CompareBatchScalar)


class TestGaussPlusGaussPlusExp : public PDFTest
{
  protected:
    TestGaussPlusGaussPlusExp() :
      PDFTest("Gauss + Gauss + Exp")
    {
      auto x = new RooRealVar("x", "x", 0., 100.);

      auto c = new RooRealVar("c", "c", -0.2, -100., 0.);
      auto expo = new RooExponential("expo", "expo", *x, *c);


      auto mean = new RooRealVar("mean1", "mean of gaussian", 20., -10, 100);
      auto sigma = new RooRealVar("sigma1", "width of gaussian", 4., 0.1, 20);
      auto gauss = new RooGaussian("gauss1", "gaussian PDF", *x, *mean, *sigma);

      auto mean2 = new RooRealVar("mean2", "mean of gaussian", 50., -10, 100);
      auto sigma2 = new RooRealVar("sigma2", "width of gaussian", 10., 0.1, 20);
      auto gauss2 = new RooGaussian("gauss2", "gaussian PDF", *x, *mean2, *sigma2);

      auto nGauss = new RooRealVar("nGauss", "Fraction of Gauss component", 100., 1., 1.E15);
      auto nGauss2 = new RooRealVar("nGauss2", "Fraction of Gauss component", 200., 1., 1.E15);
      auto nExp = new RooRealVar("nExp", "Number of events in exp", 1000, 1, 1.E15);
      _pdf = std::make_unique<RooAddPdf>("Sum2GausExp", "Sum of Gaus and Exponentials",
          RooArgSet(std::initializer_list<RooAbsArg*>{gauss, gauss2, expo}),
          RooArgSet(std::initializer_list<RooAbsArg*>{nGauss, nGauss2, nExp}));


      _variables.addOwned(x);

      _variablesToPlot.add(x);

      for (auto par : {c, mean, sigma, mean2, sigma2}) {
        _parameters.addOwned(par);
      }

      for (auto par : {nGauss, nGauss2, nExp}) {
        _yields.addOwned(par);
      }

      for (auto obj : std::initializer_list<RooAbsPdf*>{expo, gauss, gauss2}) {
        _otherObjects.addOwned(obj);
      }
    }
};

RUN_SCALAR(TestGaussPlusGaussPlusExp, Scalar)
RUN_BATCH(TestGaussPlusGaussPlusExp, Batch)
RUN_BATCH_VS_SCALAR(TestGaussPlusGaussPlusExp, CompareBatchScalar)

