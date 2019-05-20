#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooMsgService.h"
#include "RooFitResult.h"
#include "RooAbsRealLValue.h"
#include "RooGaussian.h"
#include "ROOT/RMakeUnique.hxx"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooRandom.h"

#include <fenv.h>
#include <numeric>

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


class PDFTest : public ::testing::Test
{
  protected:
  PDFTest(std::string&& name, std::size_t nEvt = 100000) :
    _name(name),
    _nEvents(nEvt)
  {
    //Shut up integration messages
    auto& msg = RooMsgService::instance();
    msg.getStream(0).minLevel = RooFit::WARNING;
    msg.getStream(1).minLevel = RooFit::WARNING;

    RooRandom::randomGenerator()->SetSeed(1337);
  }

  void SetUp() override {
    _origParameters.addClone(_parameters);
    _origYields.addClone(_yields);

    setUpData();

    makePlots(::testing::UnitTest::GetInstance()->current_test_info()->name()+std::string("_prefit"));

    resetParameters();
  }

  void TearDown() override {
    setValuesConstant(_otherObjects, false);
    makePlots(::testing::UnitTest::GetInstance()->current_test_info()->name()+std::string("_postfit"));
  }

  virtual void setUpData() {
    _data.reset(_pdf->generate(_variables, _nEvents));
  }

  void makePlots(std::string&& fitStage) const{
    for (auto elm : _variablesToPlot) {
      auto var = static_cast<RooRealVar*>(elm);
      auto canv = std::make_unique<TCanvas>();
      auto frame = std::unique_ptr<RooPlot>(var->frame());
      _data->plotOn(frame.get());
      _pdf->plotOn(frame.get(), RooFit::Precision(-1.));
      _pdf->paramOn(frame.get());
      frame->Draw();
      canv->Draw();
      std::string filename = _plotDirectory + _name + "_";
      filename += var->GetName();
      filename += "_" + fitStage + ".png";
      std::replace(filename.begin(), filename.end(), ' ', '_');
      canv->SaveAs(filename.c_str());
    }
  }

  ~PDFTest() {

  }

  void checkParams() {
    ASSERT_FALSE(_parameters.overlaps(_otherObjects)) << "Collections of parameters and other objects "
        << "cannot overlap. This will lead to wrong results, as parameters get kicked before the fit, "
        << "other objects are set constant. Hence, the fit cannot change them.";
    ASSERT_FALSE(_yields.overlaps(_otherObjects)) << "Collections of yields and other objects "
        << "cannot overlap. This will lead to wrong results, as parameters get kicked before the fit, "
        << "other objects are set constant. Hence, the fit cannot change them.";

    for (auto param : _parameters) {
        auto postFit = static_cast<RooRealVar*>(param);
        auto preFit  = static_cast<RooRealVar*>(_origParameters.find(param->GetName()));
        ASSERT_NE(preFit, nullptr) << "for parameter '" << param->GetName() << '\'';
        EXPECT_LE(fabs(postFit->getVal() - preFit->getVal()), 2.*postFit->getError())
                    << "[Within 2 std-dev: " << param->GetName()
                    << " (" << postFit->getVal() << " +- " << 2.*postFit->getError() << ")"
                    << " == " << preFit->getVal() << "]";

        EXPECT_LE(fabs(postFit->getVal() - preFit->getVal()), 1.5*postFit->getError())
                    << "[Within 1.5 std-dev: " << param->GetName()
                    << " (" << postFit->getVal() << " +- " << 1.5*postFit->getError() << ")"
                    << " == " << preFit->getVal() << "]";

        EXPECT_NEAR(postFit->getVal(), preFit->getVal(), fabs(postFit->getVal())*5.E-2)
            << "[Within 5% for parameter '" << param->GetName() << "']";

      }

    if (!_yields.empty()) {
      const double totalPre = std::accumulate(_origYields.begin(), _origYields.end(), 0.,
          [](double acc, const RooAbsArg* arg){
        return acc + static_cast<const RooAbsReal*>(arg)->getVal();
      });
      const double totalPost = std::accumulate(_yields.begin(), _yields.end(), 0.,
          [](double acc, const RooAbsArg* arg){
        return acc + static_cast<const RooAbsReal*>(arg)->getVal();
      });
      ASSERT_NE(totalPre, 0.);
      ASSERT_NE(totalPost, 0.);
      ASSERT_LE(fabs(totalPost - _nEvents) / _nEvents, 0.1) << "Total event yield not matching"
          << " number of generated events.";

      for (auto yield : _yields) {
        auto postFit = static_cast<RooRealVar*>(yield);
        auto preFit  = static_cast<RooRealVar*>(_origYields.find(yield->GetName()));
        ASSERT_NE(preFit, nullptr) << "for parameter '" << yield->GetName() << '\'';

        EXPECT_NEAR(postFit->getVal()/totalPost,
            preFit->getVal()/totalPre, 0.01) << "Yield " << yield->GetName()
              << " = " << postFit->getVal()
              << " does not match pre-fit ratios.";
      }
    }
  }

  void setValuesConstant(const RooAbsCollection& coll, bool constant) const {
    for (auto obj : coll) {
      auto lvalue = dynamic_cast<RooAbsRealLValue*>(obj);
      if (lvalue)
        lvalue->setConstant(constant);
    }
  }

  void resetParameters() const {
    //Kick parameters away from best-fit value
    for (auto param : _parameters) {
      auto lval = static_cast<RooAbsRealLValue*>(param);
      auto orig = static_cast<RooAbsRealLValue*>(_origParameters.find(param->GetName()));
      *lval = orig->getVal() * 1.3;
    }

    for (auto yield : _yields) {
      auto lval = static_cast<RooAbsRealLValue*>(yield);
      auto orig = static_cast<RooAbsRealLValue*>(_origYields.find(yield->GetName()));
      *lval = orig->getVal() * 1.3;
    }

    setValuesConstant(_otherObjects, true);
  }

  std::unique_ptr<RooAbsPdf> _pdf;
  std::unique_ptr<RooDataSet> _data;

  std::string _name;
  std::string _plotDirectory{"/tmp/"};
  RooArgSet _variables;
  RooArgSet _variablesToPlot;
  RooArgSet _parameters;
  RooArgSet _yields;
  RooArgSet _origYields;
  RooArgSet _origParameters;
  RooArgSet _otherObjects;
  const std::size_t _nEvents;
  double _toleranceParameter{5.E-5};
  double _toleranceCorrelation{1.E-4};
};



class PDFTestWeightedData : public PDFTest {
  protected:
    PDFTestWeightedData(const char* name, std::size_t events = 100000) :
      PDFTest(name, events) { }

    void setUpData() override {
      PDFTest::setUpData();
      RooRealVar var("gausWeight", "gausWeight", 0, 10);
      RooConstVar mean("meanWeight", "", 1.);
      RooConstVar sigma("sigmaWeight", "", 0.2);
      RooGaussian gausDistr("gausDistr", "gausDistr", var, mean, sigma);
      std::unique_ptr<RooDataSet> gaussData(gausDistr.generate(RooArgSet(var), _data->numEntries()));
      _data->merge(gaussData.get());

      auto wdata = new RooDataSet(_data->GetName(), _data->GetTitle(), *_data->get(),
          RooFit::Import(*_data), RooFit::WeightVar("gausWeight"));
      wdata->Print("V");
      wdata->get(0)->Print("V");
      std::cout << "Weight is " << wdata->weight() << std::endl;
      _data.reset(wdata);
    }
};



#define RUN_BATCH_VS_SCALAR(TEST_CLASS, TEST_NAME) \
TEST_F(TEST_CLASS, TEST_NAME) {\
  resetParameters();\
  MyTimer singleTimer("Fitting scalar mode " + _name);\
  auto resultSingle = _pdf->fitTo(*_data,\
      RooFit::BatchMode(false),\
      RooFit::SumW2Error(false),\
      RooFit::PrintLevel(-1), RooFit::Save());\
  std::cout << singleTimer << std::endl;\
  ASSERT_TRUE(resultSingle);\
  EXPECT_EQ(resultSingle->status(), 0) << "[Scalar fit did not converge.]";\
  \
  resetParameters();\
  MyTimer batchTimer("Fitting batch mode " + _name);\
  auto resultBatch = _pdf->fitTo(*_data,\
      RooFit::BatchMode(true), RooFit::Optimize(0),\
      RooFit::SumW2Error(false),\
      RooFit::PrintLevel(-1), RooFit::Save());\
  std::cout << batchTimer << std::endl;\
  ASSERT_TRUE(resultBatch);\
  EXPECT_EQ(resultBatch->status(), 0) << "[Batch fit did not converge.]";\
  \
  EXPECT_TRUE(resultSingle->isIdentical(*resultBatch, _toleranceParameter, _toleranceCorrelation));\
}


#define RUN_BATCH(TEST_CLASS, TEST_NAME) \
TEST_F(TEST_CLASS, TEST_NAME) {\
  resetParameters();\
  MyTimer batchTimer("Fitting batch mode " + _name);\
  auto resultBatch = _pdf->fitTo(*_data,\
      RooFit::BatchMode(true),\
      RooFit::SumW2Error(false),\
      RooFit::Optimize(0), RooFit::PrintLevel(-1), RooFit::Save());\
  std::cout << batchTimer << std::endl;\
  EXPECT_EQ(resultBatch->status(), 0) << "[Batch fit did not converge.]";\
  \
  checkParams();\
}


#define RUN_SCALAR(TEST_CLASS, TEST_NAME) \
TEST_F(TEST_CLASS, TEST_NAME) {\
  resetParameters();\
  MyTimer singleTimer("Fitting scalar mode " + _name);\
  auto resultSingle = _pdf->fitTo(*_data,\
      RooFit::BatchMode(false),\
      RooFit::SumW2Error(false),\
      RooFit::PrintLevel(-1), RooFit::Save());\
  EXPECT_EQ(resultSingle->status(), 0) << "[Scalar fit did not converge.]";\
  std::cout << singleTimer << std::endl;\
  \
  checkParams();\
}

//TEST_P(PDFTest, DISABLED_FitSingle) {
//  auto& testData = GetParam();
//
//  MyTimer fitTimer("Fitting single mode " + testData._name);
//  __itt_resume();
//  auto result = _pdf->fitTo(*_data, RooFit::BatchMode(false), RooFit::PrintLevel(-1), RooFit::Save());
//  __itt_pause();
//  std::cout << fitTimer << std::endl;
//
//  g_timers.push_back(std::move(fitTimer));
//
//  SCOPED_TRACE("FitNoBatch for " + testData._name);
//  EXPECT_EQ(result->status(), 0) << " for fit result " << *result;
//  checkParams();
//}
//
//TEST_P(PDFTest, DISABLED_FitBatch) {
//  auto& testData = GetParam();
//
//  MyTimer fitTimer("Fitting batch mode " + testData._name);
//  __itt_resume();
//  auto result = _pdf->fitTo(*_data, RooFit::BatchMode(true), RooFit::PrintLevel(-1), RooFit::Save());
//  __itt_pause();
//  std::cout << fitTimer << std::endl;
//
//  g_timers.push_back(std::move(fitTimer));
//
//
//  SCOPED_TRACE("FitBatch for " + testData._name);
//  EXPECT_EQ(result->status(), 0) << " for fit result " << *result;
//  checkParams();
//}
//

