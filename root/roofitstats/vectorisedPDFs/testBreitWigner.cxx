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

#include "RooBreitWigner.h"

class TestBreitWigner : public PDFTest
{
  protected:
    TestBreitWigner() :
      PDFTest("BreitWigner", 300000)
  {
      // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
        auto x = new RooRealVar("x", "x", -10, 10);
        auto mean = new RooRealVar("mean", "mean", 1, -7, 7);
        auto width = new RooRealVar("a2", "a2", 1, 0.5, 2.5);

        // Build gaussian p.d.f in terms of x,mean and sigma
        _pdf = std::make_unique<RooBreitWigner>("breitWigner", "breitWigner PDF", *x, *mean, *width);


      _variables.addOwned(*x);

      _variablesToPlot.add(*x);

      for (auto par : {mean, width}) {
        _parameters.addOwned(*par);
      }
  }
};

COMPARE_FIXED_VALUES_UNNORM(TestBreitWigner, CompareFixedValuesUnnorm)
COMPARE_FIXED_VALUES_NORM(TestBreitWigner, CompareFixedValuesNorm)
COMPARE_FIXED_VALUES_NORM_LOG(TestBreitWigner, DISABLED_CompareFixedNormLog)
FIT_TEST_SCALAR(TestBreitWigner, DISABLED_RunScalar)
FIT_TEST_BATCH(TestBreitWigner, DISABLED_RunBatch)
FIT_TEST_BATCH_VS_SCALAR(TestBreitWigner, DISABLED_CompareBatchScalar)

