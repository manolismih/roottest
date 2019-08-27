//~ // Author: Stephan Hageboeck, CERN  26 Apr 2019

//~ /*****************************************************************************
 //~ * RooFit
 //~ * Authors:                                                                  *
 //~ *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 //~ *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 //~ *                                                                           *
 //~ * Copyright (c) 2000-2019, Regents of the University of California          *
 //~ *                          and Stanford University. All rights reserved.    *
 //~ *                                                                           *
 //~ * Redistribution and use in source and binary forms,                        *
 //~ * with or without modification, are permitted according to the terms        *
 //~ * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 //~ *****************************************************************************/

//~ #include "VectorisedPDFTests.h"

//~ #include "RooLegendre.h"

//~ class TestLegendre : public PDFTest
//~ {
  //~ protected:
    //~ TestLegendre() :
      //~ PDFTest("Legendre", 30000000)
  //~ {
    //~ // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
    //~ auto x = new RooRealVar("x", "x", 0, -1.2, 1.2);
    
    //~ // Build gaussian p.d.f in terms of x,mean and sigma
    //~ _pdf = std::make_unique<RooLegendre>("Legendre", "Legendre PDF", *x, 1, 2, 2, 3);
    
    //~ for (auto var : {x}) {
    //~ _variables.addOwned(*var);      
    //~ }

    //~ //for (auto par : {legendre, beta, mu}) {
      //~ //_parameters.addOwned(*par);
    //~ //}
    
    //~ _variablesToPlot.add(*x);
  //~ }
//~ };

//~ COMPARE_FIXED_VALUES_UNNORM(TestLegendre, CompareFixedValuesUnnorm)
//~ //COMPARE_FIXED_VALUES_NORM(TestLegendre, CompareFixedValuesNorm)
//~ //COMPARE_FIXED_VALUES_NORM_LOG(TestLegendre, CompareFixedNormLog)
//~ //FIT_TEST_SCALAR(TestLegendre, RunScalar)
//~ //FIT_TEST_BATCH(TestLegendre, RunBatch)
//~ //FIT_TEST_BATCH_VS_SCALAR(TestLegendre, CompareBatchScalar)

