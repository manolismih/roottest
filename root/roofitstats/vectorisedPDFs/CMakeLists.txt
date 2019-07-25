add_library(VectorisedPDFTests STATIC VectorisedPDFTests.cxx)
link_libraries(VectorisedPDFTests Gpad RooFitCore RooFit)

ROOT_ADD_GTEST(testCompatMode testCompatMode.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testGauss testGauss.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testPoisson testPoisson.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testAddPdf testAddPdf.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testExponential testExponential.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testNestedPDFs testNestedPDFs.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testProductPdf testProductPdf.cxx
              LIBRARIES Gpad RooFitCore RooFit)
ROOT_ADD_GTEST(testJohnson testJohnson.cxx
              LIBRARIES Gpad RooFitCore RooFit)

