#include <iostream>
#include <box/boxutils.hpp>
#include <oneobj/contboxconstr/dejong.hpp>
#include <oneobj/contboxconstr/rosenbrock.hpp>
#include <oneobj/contboxconstr/ackley1.hpp>
#include <funccnt.hpp>
#include <methods/lins/goldsec/goldsec.hpp>
#include <methods/lins/smartls/smartls.hpp>
#include "random_method.hpp"
#include "testfunc.hpp"

/*
 * 
 */
/*int main(int argc, char** argv) {
    const int n = 2;
    
    OPTITEST::TestProblemFactory fact(n, -4, 8);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    auto obj = std::make_shared<COMPI::FuncCnt<double>>(mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);
    
    LOCSEARCH::RandomMethod<double> desc(*mpp);   
    desc.getOptions().mDoTracing = true;
    desc.getOptions().mInc = 1.618;
    desc.getOptions().mDec = 0.618;

    double x[n];
    x[0] = 8.;
    x[1] = 9.;

    
    double v;
    bool rv = desc.search(x, v);

    std::cout << desc.about() << "\n";
    std::cout << "Found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
    return 0;
}*/
int main(int argc, char** argv) {
    const int n = 1000;
    
    OPTITEST::RosenbrockProblemFactory fact(n, -4, 8);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    auto obj = std::make_shared<COMPI::FuncCnt<double>>(mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);
    LOCSEARCH::RandomMethod<double> desc(*mpp);
#if 1    
    LOCSEARCH::GoldenSecLS<double>* locs = new LOCSEARCH::GoldenSecLS<double>(*mpp);
    locs->getOptions().mSInit = 0.1;
    locs->getOptions().mDelta = 0.02;
    locs->getOptions().mMaxBackSteps = 16;
    locs->getOptions().mDoTracing = true;
#endif
#if  0  
    LOCSEARCH::SmartLS<double> *locs = new LOCSEARCH::SmartLS<double>(*mpp);
    locs->getOptions().mSInit = 1e-4;
    locs->getOptions().mDoTracing = true;
    locs->getOptions().mMaxFailStepsBack = 0;
#endif    
    desc.getLineSearch().reset(locs);    
    desc.getOptions().mInc = 1.418;
    desc.getOptions().mDec = 0.368;
    //desc.getOptions().mDoTracing = false;
    desc.getOptions().numbOfPoints = 100;
    desc.getOptions().maxStepNumber = 2000;

    //desc.getOptions().mSearchType = LOCSEARCH::AdvancedCoordinateDescent<double>::SearchTypes::NO_DESCENT;
    //desc.getOptions().mVicinityAdaptation = LOCSEARCH::AdvancedCoordinateDescent<double>::VARIABLE_ADAPTATION;
    //desc.getOptions().mVicinityAdaptation = LOCSEARCH::AdvancedCoordinateDescent<double>::UNIFORM_ADAPTATION;
    
    double x[n];
    snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
    double v;
    bool rv = desc.search(x, v);
    std::cout << desc.about() << "\n";
    std::cout << "Found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
    SG_ASSERT(v <= 0.01);

    return 0;
}