if(ROOT_roofit_FOUND)
  ROOTTEST_ADD_TEST(read-scientificnotation
                    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/read_scientific_notation.py
                    OUTCNVCMD grep -v -e "Wouter"
                    OUTREF read_scientific_notation.ref)

  ROOTTEST_ADD_TEST(RooDataSet_ASCII_in
                    MACRO ${CMAKE_CURRENT_SOURCE_DIR}/ASCII-in-out.C
                    MACROARG "false, \"${CMAKE_CURRENT_SOURCE_DIR}/ASCII-in-out_data.txt\""
                    OUTCNVCMD grep -v -e "Wouter" -e "\.txt" -e "reserved"
                    OUTREF ASCII-in-out.ref)
  ROOTTEST_ADD_TEST(RooDataSet_ASCII_out
                    MACRO ${CMAKE_CURRENT_SOURCE_DIR}/ASCII-in-out.C
                    MACROARG "true, \"./ASCII-in-out_result.txt\""
                    POSTCMD diff ${CMAKE_CURRENT_SOURCE_DIR}/ASCII-in-out_data.txt ./ASCII-in-out_result.txt)
if(CMAKE_BUILD_TYPE MATCHES "Debug|RelWithDebInfo")
  ROOTTEST_ADD_TEST(ListPrettyPrinters
                    COMMAND ls -l "${ROOTSYS}/lib/*-gdb.py"
                    PASSREGEX "libRooFitCore.so-gdb.py")
  ROOTTEST_ADD_TEST(TestLoadingPrettyPrinters
                    COMMAND gdb root.exe -ex "add-auto-load-safe-path ${ROOTSYS}/lib" -ex "run -b -q -e \"RooRealVar x\"" -ex "info pretty-printer" -ex "quit"
                    PASSREGEX "objfile.*libRooFitCore.so pretty-printers:")
endif(CMAKE_BUILD_TYPE MATCHES "Debug|RelWithDebInfo")

endif()
