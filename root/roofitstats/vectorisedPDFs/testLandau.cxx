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
#include "RooLandau.h"

class TestLandauEvil: public PDFTest
{
  protected:
    TestLandauEvil() :
      PDFTest("Landau_evil", 300000)
    {
      // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
        auto x = new RooRealVar("x", "x", -4000, 10000);
        auto mean = new RooRealVar("mean", "mean of landau", 1, -1000, 1000);
        auto sigma = new RooRealVar("sigma", "width of landau", 100, 0.1, 1000);

        // Build landau p.d.f in terms of x,mean and sigma
        _pdf = std::make_unique<RooLandau>("landau", "landau PDF", *x, *mean, *sigma);

      _variables.addOwned(*x);

      _variablesToPlot.add(*x);

      for (auto par : {mean, sigma}) {
        _parameters.addOwned(*par);
      }

     
    }
};

COMPARE_FIXED_VALUES_UNNORM(TestLandauEvil, CompareFixedValuesUnnorm)
COMPARE_FIXED_VALUES_NORM(TestLandauEvil, CompareFixedValuesNorm)
COMPARE_FIXED_VALUES_NORM_LOG(TestLandauEvil, CompareFixedNormLog)
FIT_TEST_SCALAR(TestLandauEvil, DISABLED_RunScalar)
FIT_TEST_BATCH(TestLandauEvil, DISABLED_RunBatch)
FIT_TEST_BATCH_VS_SCALAR(TestLandauEvil, DISABLED_CompareBatchScalar)


class TestLandau: public PDFTest
{
  protected:
    TestLandau() :
      PDFTest("Landau_convenient", 300000)
    {
      // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
        auto x = new RooRealVar("x", "x", -100, 400);
        auto mean = new RooRealVar("mean", "mean of landau", 1, -100, 100);
        auto sigma = new RooRealVar("sigma", "width of landau", 100, 0.1, 1000);

        // Build landau p.d.f in terms of x,mean and sigma
        _pdf = std::make_unique<RooLandau>("landau", "landau PDF", *x, *mean, *sigma);

      _variables.addOwned(*x);

//      _variablesToPlot.add(x);

      for (auto par : {mean, sigma}) {
        _parameters.addOwned(*par);
      }

     
    }
};

COMPARE_FIXED_VALUES_UNNORM(TestLandau, CompareFixedValuesUnnorm)
COMPARE_FIXED_VALUES_NORM(TestLandau, CompareFixedValuesNorm)
COMPARE_FIXED_VALUES_NORM_LOG(TestLandau, DISABLED_CompareFixedNormLog)
FIT_TEST_SCALAR(TestLandau, DISABLED_RunScalar)
FIT_TEST_BATCH(TestLandau, DISABLED_RunBatch)
FIT_TEST_BATCH_VS_SCALAR(TestLandau, DISABLED_CompareBatchScalar)

