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

#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExponential.h"


class TestNestedPDFs : public PDFTest
{
  protected:
    TestNestedPDFs() :
      PDFTest("Gauss + Exp(x, Gauss(y1)) + Exp(x, Gauss(y2))")
  {
      auto x = new RooRealVar("x", "x", 0.001, 20.);
      x->setBins(10); //Speed tweak for plotting

      auto y1 = new RooRealVar("y1", "y1", -10, 0.);
      y1->setBins(10);
      auto mean1 = new RooRealVar("mean1", "mean1", -3, -5, -0.1);
      auto sigma1 = new RooRealVar("sigma1", "sigma1", 0.4, 0.001, 5);
      auto c1 = new RooGaussian("gauss1", "c1", *y1, *mean1, *sigma1);

      auto y2 = new RooRealVar("y2", "y2", -10, 0.);
      y2->setBins(10);
      auto mean2 = new RooRealVar("mean2", "mean2", -1.5, -5, -0.1);
      auto sigma2 = new RooRealVar("sigma2", "sigma2", 0.5, 0.001, 5);
      auto c2 = new RooGaussian("gauss2", "c2", *y2, *mean2, *sigma2);

      auto expo1 = new RooExponential("expo1", "expo1", *x, *mean1);
      auto expo2 = new RooExponential("expo2", "expo2", *x, *mean2);


      auto mean = new RooRealVar("mean", "mean of gaussian", 5., 0., 20.);
      auto sigma = new RooRealVar("sigma", "width of gaussian", 1.337, 0.1, 10);
      auto gauss = new RooGaussian("gauss", "gaussian PDF", *x, *mean, *sigma);

      auto nGauss = new RooRealVar("nGauss", "Fraction of Gauss component", 500.);
      auto nExp1 = new RooRealVar("nExp1", "Number of events in exp1", 1000);
      auto nExp2 = new RooRealVar("nExp2", "Number of events in exp2", 2000);
      _pdf = std::make_unique<RooAddPdf>("SumGaus2Exp", "Sum of Gaus and 2 Exponentials",
          RooArgSet(*gauss, *expo1, *expo2),
          RooArgSet(*nGauss, *nExp1, *nExp2));

      for (auto var : {x, y1, y2}) {
        _variables.addOwned(*var);
      }

      for (auto var : {x, y1, y2}) {
        _variablesToPlot.add(*var);
      }

      for (auto par : {mean1, sigma1, mean2, mean, sigma}) {
        _parameters.addOwned(*par);
      }

      for (auto par : { nGauss, nExp1, nExp2}) {
        _yields.addOwned(*par);
      }

      for (auto obj : std::initializer_list<RooAbsArg*>{c1, c2, sigma2, expo1, expo2, gauss}) {
        _otherObjects.addOwned(*obj);
      }
  }
};

FIT_TEST_SCALAR(TestNestedPDFs, RunScalar)
FIT_TEST_BATCH(TestNestedPDFs, DISABLED_RunBatch)
FIT_TEST_BATCH_VS_SCALAR(TestNestedPDFs, CompareBatchScalar)





