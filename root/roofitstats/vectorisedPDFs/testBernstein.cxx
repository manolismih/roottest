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

#include "RooBernstein.h"

class TestBernstein2 : public PDFTest
{
  protected:
    TestBernstein2() :
      PDFTest("Bernstein2", 30000000)
  {
      // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
        auto x = new RooRealVar("x", "x", -10, 10);
        auto a1 = new RooRealVar("a1", "a1", 0.3, 0.01, 0.5);
        auto a2 = new RooRealVar("a2", "a2", 0.2, 0.01, 0.5);

        // Build gaussian p.d.f in terms of x,mean and sigma
        _pdf = std::make_unique<RooBernstein>("polynomial2", "polynomial PDF 2 coefficients", *x, RooArgSet(*a1,*a2));


      _variables.addOwned(*x);

      _variablesToPlot.add(*x);

      for (auto par : {a1, a2}) {
        _parameters.addOwned(*par);
      }
  }
};

COMPARE_FIXED_VALUES_UNNORM(TestBernstein2, CompareFixedValuesUnnorm)
COMPARE_FIXED_VALUES_NORM(TestBernstein2, CompareFixedValuesNorm)
COMPARE_FIXED_VALUES_NORM_LOG(TestBernstein2, CompareFixedNormLog)
//FIT_TEST_SCALAR(TestBernstein2, RunScalar)
//FIT_TEST_BATCH(TestBernstein2, RunBatch)
//FIT_TEST_BATCH_VS_SCALAR(TestBernstein2, CompareBatchScalar)

class TestBernstein5 : public PDFTest
{
  protected:
    TestBernstein5() :
      PDFTest("Bernstein5", 300000)
  {
      // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
        auto x = new RooRealVar("x", "x", -100, 50);
        auto a1 = new RooRealVar("a1", "a1", 1.0);
        auto a2 = new RooRealVar("a2", "a2", 10.0, 9.0, 12.0);
        auto a3 = new RooRealVar("a3", "a3", 0.09, 0.05, 0.1);
        auto a4 = new RooRealVar("a4", "a4", 0.00009, 0.00005, 0.0005);
        auto a5 = new RooRealVar("a5", "a5", 0.00009, 0.00005, 0.0005);

        // Build gaussian p.d.f in terms of x,mean and sigma
        _pdf = std::make_unique<RooBernstein>("bernstein5", "bernstein PDF 5 coefficients", *x, RooArgSet(*a1, *a2, *a3, *a4));


      _variables.addOwned(*x);

      _variablesToPlot.add(*x);

      for (auto par : { a2, a3, a4, a5}) {
        _parameters.addOwned(*par);
      }
  }
};

COMPARE_FIXED_VALUES_UNNORM(TestBernstein5, CompareFixedValuesUnnorm)
COMPARE_FIXED_VALUES_NORM(TestBernstein5, CompareFixedValuesNorm)
COMPARE_FIXED_VALUES_NORM_LOG(TestBernstein5, CompareFixedNormLog)
FIT_TEST_SCALAR(TestBernstein5, DISABLED_RunScalar)
FIT_TEST_BATCH(TestBernstein5, DISABLED_RunBatch)
FIT_TEST_BATCH_VS_SCALAR(TestBernstein5, DISABLED_CompareBatchScalar)

