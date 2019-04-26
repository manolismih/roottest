#include "PDFTests.h"
#include "RooFitResult.h"

#include <fenv.h>

#include "gtest/gtest.h"


#include <ctime>
#ifdef __INTEL_COMPILER
#include "ittnotify.h"
#else
void __itt_resume() {}
void __itt_pause() {}
#endif


class MyTimer {
   public:
   MyTimer(std::string&& name)
   : m_name(name), m_startTime(clock()), m_endTime(clock()) {
      m_startTime = clock();
   }
   
   ~MyTimer() {
      print(std::cout);
   }
   
   clock_t diffTime() const {
      return clock() - m_startTime;
   }
   
   void interval() {
      m_endTime = clock();
   }
   
   void print(std::ostream& str) {
      clock_t diff = m_endTime - m_startTime;
      str << "\n" << "Timer '" << m_name << "':\t" << double(diff)/CLOCKS_PER_SEC << "s" << std::endl;
   }
   
   private:
   std::string m_name;
   clock_t m_startTime;
   clock_t m_endTime;
};

std::ostream& operator<<(std::ostream& str, MyTimer& timer) {
   timer.interval();
   timer.print(str);
   return str;
}

std::vector<MyTimer> g_timers;


class PDFTest : public ::testing::TestWithParam<TestPDFs::TestConfig>
{
  protected:
  PDFTest() {
    auto& testData = GetParam();

    _pdf = testData._pdf.get();
    _data.reset(testData._pdf->generate(*testData._variables, 100000));

    setupParameters();

    //Shut up integration messages
    auto& msg = RooMsgService::instance();
    msg.getStream(0).minLevel = RooFit::WARNING;
    msg.getStream(1).minLevel = RooFit::WARNING;
  }

  ~PDFTest() {
    auto& testData = GetParam();

    //Reset paramters
    for (auto param : *testData._parameters) {
      auto lval = static_cast<RooAbsRealLValue*>(param);
      auto orig = static_cast<RooAbsRealLValue*>(testData._origParameters->find(param->GetName()));
      *lval = *orig;
    }

    for (auto obj : *testData._otherObjects) {
      auto lvalue = dynamic_cast<RooAbsRealLValue*>(obj);
      if (lvalue)
        lvalue->setConstant(false);
    }
  }

  void checkParams() {
    auto& testData = GetParam();

    for (auto param : *testData._parameters) {
        auto postFit = static_cast<RooRealVar*>(param);
        auto preFit  = static_cast<RooRealVar*>(testData._origParameters->find(param->GetName()));
        ASSERT_NE(preFit, nullptr) << "for parameter " << param->GetName();
        EXPECT_NEAR(postFit->getVal(), preFit->getVal(), fabs(postFit->getVal())*1.E-2)
            << "for parameter " << param->GetName();
        EXPECT_LE(fabs(postFit->getVal() - preFit->getVal()), postFit->getError())
            << "for parameter " << param->GetName();
      }
  }

  void setupParameters() {
    auto& testData = GetParam();
    //Kick parameters away from best-fit value
    for (auto param : *testData._parameters) {
      auto lval = static_cast<RooAbsRealLValue*>(param);
      auto orig = static_cast<RooAbsRealLValue*>(testData._origParameters->find(param->GetName()));
      *lval = orig->getVal() * 1.3;
    }

    for (auto obj : *testData._otherObjects) {
      auto lvalue = dynamic_cast<RooAbsRealLValue*>(obj);
      if (lvalue)
        lvalue->setConstant(true);
    }
  }

  RooAbsPdf* _pdf;
  std::unique_ptr<RooDataSet> _data;
};


TEST_P(PDFTest, DISABLED_FitSingle) {
  auto& testData = GetParam();

  MyTimer fitTimer("Fitting single mode " + testData._name);
  __itt_resume();
  auto result = _pdf->fitTo(*_data, RooFit::BatchMode(false), RooFit::PrintLevel(-1), RooFit::Save());
  __itt_pause();
  std::cout << fitTimer << std::endl;

  g_timers.push_back(std::move(fitTimer));

  SCOPED_TRACE("FitNoBatch for " + testData._name);
  EXPECT_EQ(result->status(), 0) << " for fit result " << *result;
  checkParams();
}

TEST_P(PDFTest, DISABLED_FitBatch) {
  auto& testData = GetParam();

  MyTimer fitTimer("Fitting batch mode " + testData._name);
  __itt_resume();
  auto result = _pdf->fitTo(*_data, RooFit::BatchMode(true), RooFit::PrintLevel(-1), RooFit::Save());
  __itt_pause();
  std::cout << fitTimer << std::endl;

  g_timers.push_back(std::move(fitTimer));


  SCOPED_TRACE("FitBatch for " + testData._name);
  EXPECT_EQ(result->status(), 0) << " for fit result " << *result;
  checkParams();
}

TEST_P(PDFTest, CompareBatchVsNormal) {
  auto& testData = GetParam();

  setupParameters();
  MyTimer singleTimer("Fitting single mode " + testData._name);
  auto resultSingle = _pdf->fitTo(*_data, RooFit::BatchMode(false), RooFit::PrintLevel(-1), RooFit::Save());
  std::cout << singleTimer << std::endl;

  setupParameters();
  MyTimer batchTimer("Fitting batch mode " + testData._name);
  auto resultBatch = _pdf->fitTo(*_data, RooFit::BatchMode(true), RooFit::PrintLevel(-1), RooFit::Save());
  std::cout << batchTimer << std::endl;



  g_timers.push_back(std::move(singleTimer));
  g_timers.push_back(std::move(batchTimer));

  EXPECT_TRUE(resultSingle->isIdentical(*resultBatch));
}


INSTANTIATE_TEST_CASE_P(FitPDFs,
    PDFTest,
    ::testing::ValuesIn(TestPDFs::makeTestPDFs()));
