/* 
 * File:   tesgfsdesc.cpp
 * Author: medved
 *
 * Created on February 22, 2016, 3:17 PM
 */

#include <iostream>
#include <oneobj/contboxconstr/dejong.hpp>
#include <funccnt.hpp>
#include <methods/lins/dichotls/dichotls.hpp>
#include <methods/lins/quadls/quadls.hpp>
#include "varcoorgrad.hpp"



/*
 * 
 */
int main(int argc, char** argv) {
    const int n = 10;
    OPTITEST::DejongProblemFactory fact(n, -4, 8);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);


    LOCSEARCH::VarCoorGrad<double> desc(*mpp);

    desc.getOptions().mHInit = .1;
    desc.getOptions().mGradStep = 1;
    desc.getOptions().mGradMaxSteps = 8;
    double x[n];

    for (int i = 0; i < n; i++)
        x[i] = i * 100 + 1;
    double v;
    bool rv = desc.search(x, v);
    std::cout << desc.about() << "\n";
    std::cout << "Found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
    SG_ASSERT(v <= 0.01);
    
    return 0;
}
