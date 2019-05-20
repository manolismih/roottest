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


class TestPoisson : public PDFTest
{
  protected:
    TestPoisson() :
      PDFTest("Poisson", 100000)
  {
      auto x = new RooRealVar("x", "x", -10, 10);
      auto mean = new RooRealVar("mean", "Mean of Poisson", 2., 0., 10);
      _pdf = std::make_unique<RooPoisson>("Pois", "Poisson PDF", *x, *mean);

      _variables.addOwned(x);

      for (auto par : {mean}) {
        _parameters.addOwned(par);
      }

  }
};

RUN_BATCH(TestPoisson, RunBatch)
RUN_BATCH_VS_SCALAR(TestPoisson, CompareBatchScalar)
