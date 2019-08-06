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

#include "RooJohnson.h"
#include "RooAddition.h"

class TestJohnson : public PDFTest
{
  protected:
    TestJohnson() :
      PDFTest("Johnson", 200000)
  {
      auto mass = new RooRealVar("mass", "mass", 0., -100., 500.);
      auto mu = new RooRealVar("mu", "Location parameter of normal distribution", 300., 0., 500.);
      auto lambda = new RooRealVar ("lambda", "Two sigma of normal distribution", 100., 1.E-6, 200.);
      auto gamma = new RooRealVar ("gamma", "gamma", 0.5, -10., 10.);
      auto delta = new RooRealVar ("delta", "delta", 2., 1.E-6, 10.);

      _pdf = std::make_unique<RooJohnson>("johnson", "johnson", *mass, *mu, *lambda, *gamma, *delta, -1.E300);


      _variables.addOwned(*mass);

//      _variablesToPlot.add(*mass);

      for (auto par : {mu, lambda, gamma, delta}) {
        _parameters.addOwned(*par);
      }

      _toleranceParameter = 1.E-5;
  }
};

COMPARE_FIXED_VALUES_UNNORM(TestJohnson, CompareFixedUnnorm)
COMPARE_FIXED_VALUES_NORM(TestJohnson, CompareFixedNorm)

FIT_TEST_SCALAR(TestJohnson, FitScalar)
FIT_TEST_BATCH(TestJohnson, FitBatch)
FIT_TEST_BATCH_VS_SCALAR(TestJohnson, FitBatchVsScalar)



class TestJohnsonInMassAndGamma : public PDFTest
{
  protected:
    TestJohnsonInMassAndGamma() :
      PDFTest("Johnson in mass and gamma", 200000)
  {
      auto mass = new RooRealVar("mass", "mass", 0., -100., 500.);
      auto mu = new RooRealVar("mu", "Location parameter of normal distribution", 300., -100., 500.);
      auto lambda = new RooRealVar ("lambda", "Two sigma of normal distribution", 50., 1.E-6, 100.);
      auto gamma = new RooRealVar ("gamma", "gamma", -0.7, -10., 10.);
      auto delta = new RooRealVar ("delta", "delta", 1.5, 1.E-6, 10.);

      _pdf = std::make_unique<RooJohnson>("johnson", "johnson", *mass, *mu, *lambda, *gamma, *delta, -1.E300);


      _variables.addOwned(*mass);
      _variables.addOwned(*gamma);

      //      _variablesToPlot.add(x);

      for (auto par : {mu, lambda, delta}) {
        _parameters.addOwned(*par);
      }
  }
};

COMPARE_FIXED_VALUES_UNNORM(TestJohnsonInMassAndGamma, CompareFixedUnnorm)
COMPARE_FIXED_VALUES_NORM(TestJohnsonInMassAndGamma, CompareFixedNorm)

FIT_TEST_BATCH(TestJohnsonInMassAndGamma, RunBatch)
FIT_TEST_BATCH_VS_SCALAR(TestJohnsonInMassAndGamma, CompareBatchScalar)


class TestJohnsonWithFormulaParameters : public PDFTest
{
  protected:
    TestJohnsonWithFormulaParameters() :
      PDFTest("Johnson with formula")
  {
      // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
      auto mass = new RooRealVar("mass", "mass", 0., -5000., 5000.);
      auto mu = new RooRealVar("mu", "Location parameter of normal distribution", -100, -150., 200.);
      auto lambda = new RooRealVar ("lambda", "Two sigma of normal distribution", 90., 0.1, 150.);
      auto gamma = new RooRealVar ("gamma", "gamma", -0.4, -10., 10.);
      auto delta = new RooAddition("delta", "delta",
          RooArgList(RooFit::RooConst(1.337), RooFit::RooConst(0.0002)),
          RooArgList(RooFit::RooConst(1.), *mass));
//      auto delta = new RooFormulaVar("delta", "delta", "1.337 + 0.0002 * mass", RooArgList(*mass));

      // Build gaussian p.d.f in terms of x,mean and sigma
      _pdf = std::make_unique<RooJohnson>("johnson", "johnson", *mass, *mu, *lambda, *gamma, *delta, -1.E300);

      _variablesToPlot.add(*mass);

      for (auto var : {mass}) {
        _variables.addOwned(*var);
      }

      for (auto par : {mu, lambda, gamma}) {
        _parameters.addOwned(*par);
      }

      _otherObjects.addOwned(*delta);
  }
};

COMPARE_FIXED_VALUES_UNNORM(TestJohnsonWithFormulaParameters, CompareFixedUnnorm)
COMPARE_FIXED_VALUES_NORM(TestJohnsonWithFormulaParameters, DISABLED_CompareFixedNorm)

FIT_TEST_SCALAR(TestJohnsonWithFormulaParameters, DISABLED_RunScalar)
FIT_TEST_BATCH(TestJohnsonWithFormulaParameters, DISABLED_RunBatch)
FIT_TEST_BATCH_VS_SCALAR(TestJohnsonWithFormulaParameters, DISABLED_CompareBatchScalar)
