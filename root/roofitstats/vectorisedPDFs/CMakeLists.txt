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

#  ROOTTEST_ADD_TEST(TestVectorisedPDFs
#                    MACRO ${CMAKE_CURRENT_SOURCE_DIR}/testVectorisedPDFs.cxx
#                    MACROARG "100000, true"
#                    )
                    #PASSREGEX "objfile.*libRooFitCore.so pretty-printers:")

