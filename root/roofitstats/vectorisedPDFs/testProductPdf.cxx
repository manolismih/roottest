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

#include "RooProdPdf.h"
#include "RooGaussian.h"

class TestProdPdf : public PDFFitTest
{
  protected:
    TestProdPdf() :
      PDFFitTest("Gauss(x) * Gauss(y)")
  {
      auto x = new RooRealVar("x", "x", 0, -5, 5);
      auto m1 = new RooRealVar("m1", "m1", 0 , -5, 5);
      auto s1 = new RooRealVar("s1", "s1", 1, 0.7, 1.3);
      auto y = new RooRealVar("y", "y", 0, -5., 5.);
      auto m2 = new RooRealVar("m2", "m2", 1.5 , -5., 5.);
      auto s2 = new RooRealVar("s2", "s2", 2, 0.1, 5);

      //Make a 2D PDF
      auto g1 = new RooGaussian("gaus1", "gaus1", *x, *m1, *s1);
      auto g2 = new RooGaussian("gaus2", "gaus2", *y, *m2, *s2);
      _pdf = std::make_unique<RooProdPdf>("prod", "prod", RooArgSet(*g1, *g2));

      for (auto var : {x, y}) {
        _variables.addOwned(*var);
      }

      for (auto par : {m1, s1, m2, s2}) {
        _parameters.addOwned(*par);
      }

      for (auto obj : {g1, g2}) {
        _otherObjects.addOwned(*obj);
      }
  }
};

FIT_TEST_BATCH_VS_SCALAR(TestProdPdf, CompareBatchScalar)
