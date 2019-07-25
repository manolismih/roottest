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

class TestJohnson : public PDFTest
{
  protected:
    TestJohnson() :
      PDFTest("Johnson", 200000)
  {
      auto mass = new RooRealVar("mass", "mass", 0., -5000., 5000.);
      auto mu = new RooRealVar("mu", "Location parameter of normal distribution", -200., -1000., 200.);
      auto sigma = new RooRealVar ("sigma", "Two sigma of normal distribution", 50., 0., 100.);
      auto gamma = new RooRealVar ("gamma", "gamma", -10., -100., 100.);
      auto delta = new RooRealVar ("delta", "delta", 3., 0., 100.);

      _pdf = std::make_unique<RooJohnson>("johnson", "johnson", *mass, *mu, *sigma, *gamma, *delta, -1.E300);


      _variables.addOwned(*mass);

      _variablesToPlot.add(*mass);

      for (auto par : {mu, sigma, gamma, delta}) {
        _parameters.addOwned(*par);
      }
  }
};

FIT_TEST_BATCH(TestJohnson, RunBatch)
FIT_TEST_BATCH_VS_SCALAR(TestJohnson, CompareBatchScalar)



class TestJohnsonInMassAndGamma : public PDFTest
{
  protected:
    TestJohnsonInMassAndGamma() :
      PDFTest("Johnson in mass and gamma", 200000)
  {
      auto mass = new RooRealVar("mass", "mass", 0., -5000., 5000.);
      auto mu = new RooRealVar("mu", "Location parameter of normal distribution", -200., -1000., 200.);
      auto sigma = new RooRealVar ("sigma", "Two sigma of normal distribution", 20., 0., 100.);
      auto gamma = new RooRealVar ("gamma", "gamma", -10., -100., 100.);
      auto delta = new RooRealVar ("delta", "delta", 3., 0., 100.);

      _pdf = std::make_unique<RooJohnson>("johnson", "johnson", *mass, *mu, *sigma, *gamma, *delta, -1.E300);


      _variables.addOwned(*mass);
      _variables.addOwned(*gamma);

      //      _variablesToPlot.add(x);

      for (auto par : {mu, sigma, delta}) {
        _parameters.addOwned(*par);
      }
  }
};

FIT_TEST_BATCH(TestJohnsonInMassAndGamma, RunBatch)
FIT_TEST_BATCH_VS_SCALAR(TestJohnsonInMassAndGamma, CompareBatchScalar)


class TestJohsonWithFormulaParameters : public PDFTest
{
  protected:
    TestJohsonWithFormulaParameters() :
      PDFTest("Gauss(x, mean)")
  {
      // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
      auto mass = new RooRealVar("mass", "mass", 0., -5000., 5000.);
      auto mu = new RooRealVar("mu", "Location parameter of normal distribution", -200, -1000., 200.);
      auto sigma = new RooRealVar ("sigma", "Two sigma of normal distribution", 2., 0., 100.);
      auto gamma = new RooRealVar ("gamma", "gamma", -10., -100., 100.);
      auto delta = new RooFormulaVar("delta", "delta", "1.7*mu", RooArgList(*mu));

      // Build gaussian p.d.f in terms of x,mean and sigma
      _pdf = std::make_unique<RooJohnson>("johnson", "johnson", *mass, *mu, *sigma, *gamma, *delta, -1.E300);


      for (auto var : {mass}) {
        _variables.addOwned(*var);
//        _variablesToPlot.add(var);
      }

      for (auto par : {mu, sigma, gamma}) {
        _parameters.addOwned(*par);
      }

      _otherObjects.addOwned(*delta);
  }
};

FIT_TEST_SCALAR(TestJohsonWithFormulaParameters, RunScalar)
FIT_TEST_BATCH(TestJohsonWithFormulaParameters, RunBatch)
FIT_TEST_BATCH_VS_SCALAR(TestJohsonWithFormulaParameters, CompareBatchScalar)
