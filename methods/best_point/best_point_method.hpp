/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   best_point_method.hpp
 * Author: kate
 *
 * Created on 27 декабря 2018 г., 22:15
 */

#ifndef BEST_POINT_METHOD_HPP
#define BEST_POINT_METHOD_HPP

#include <sstream>
#include <stdlib.h>
#include <vector>
#include <functional>
#include <memory>
#include <solver.hpp>
#include <common/dummyls.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include <common/sgerrcheck.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>
#include <common/lineseach.hpp>
#include <math.h>
#include <string>
#include <chrono>
#include <random>

namespace LOCSEARCH {
    /**
     * Random method implementation
     */
    template <typename FT> class BestPointMethod : public COMPI::Solver<FT> {
    public:

        /**
         * Determines stopping conditions
         * @param fval current best value found
         * @param x current best point
         * @param stpn step number
         * @return true if the search should stop
         */
        using Stopper = std::function<bool(FT fval, const FT* x, int stpn) >;

        /**
         * Watches the current step
         * @param fval current best value found
         * @param current best point
         * @param stpn step number
         * @param gran - current granularity vector
         */
        using Watcher = std::function<void(FT fval, const FT* x, const std::vector<FT>& gran, int stpn) >;
        /**
         * Options for Gradient Box Descent method
         */
        struct Options {
            /**
             * Amount of points (better 3n) and max of unsuccessful steps
             */
            int numbOfPoints = 600;
            /**
             * Minimal value of step
             */
            FT minStep = 0.00001;
            /**
             * Initial value of granularity
             */
            FT* mHInit = nullptr;
            /**
             * Increase in the case of success
             */
            FT mInc = 1.618;
            /**
             * Decrease in the case of failure
             */
            FT mDec = 0.618;
            /**
             * Lower bound for granularity
             */
            FT mHLB = 1e-08;
            /**
             * Upper bound on granularity
             */
            FT mHUB = 1e+02;
            /**
             * Trace on/off
             */
            bool mDoTracing = true;
            /**
             * Max steps number
             */
            int maxStepNumber = 800;
        };

        /**
         * The constructor
         * @param prob - reference to the problem
         * @param stopper - reference to the stopper
         * @param ls - pointer to the line search
         */
        BestPointMethod(const COMPI::MPProblem<FT>& prob) :
        mProblem(prob) {
            unsigned int typ = COMPI::MPUtils::getProblemType(prob);
            SG_ASSERT(typ == COMPI::MPUtils::ProblemTypes::BOXCONSTR | COMPI::MPUtils::ProblemTypes::CONTINUOUS | COMPI::MPUtils::ProblemTypes::SINGLEOBJ);
        }

        /**
         * Retrieve the pointer to the line search
         * @return 
         */
        std::unique_ptr<LineSearch<FT>>&getLineSearch() {
            return mLS;
        }
        /**
         * Perform search
         * @param x start point and result
         * @param v  the resulting value
         * @return true if search converged and false otherwise
         */
        bool search(FT* x, FT& v) override {
            bool rv = false;
            auto obj = mProblem.mObjectives.at(0);
            int n = mProblem.mVarTypes.size();
            
            FT fcur = obj->func(x);
            int StepNumber = 0; 
            int Unsuccess = 0;
            bool br = false;
            
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator(seed);
            std::normal_distribution<FT> distribution(0.0,1.0);
            std::mt19937_64 rng;
            uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
            std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
            rng.seed(ss);
            std::uniform_real_distribution<double> unif(0, 1);
            
           
            FT* dirs;
            FT* main_dir;
            const snowgoose::Box<double>& box = *(mProblem.mBox);
            FT sft = 1.0;

            dirs = new FT[n * mOptions.numbOfPoints];
            //generator, based on normal distribution
            auto direction = [&] (int amount_of_points) {
                for (int j = 0; j < amount_of_points; j++)
                {
                    FT sum = 0.0;
                    for (int i = 0; i< n; i++)
                    {
                        dirs [n*j + i] = distribution(generator);
                        sum += (dirs [n*j + i]) * (dirs [n*j + i]);  
                    }
                    sum = sqrt(sum);
                    for (int i = 0; i< n; i++)
                    {
                        dirs [n*j + i] /= sum;  
                    }
                }
            };
            auto normalize = [&] () {
                        FT sum = 0.0;
                        for (int i = 0; i< n; i++)
                        {
                            sum += (main_dir [i]) * (main_dir [i]);  
                        }
                        sum = sqrt(sum);
                        for (int i = 0; i< n; i++)
                        {
                            main_dir [ i] /= sum;  
                        }
            };

            auto inc = [this] (FT h) {
                FT t = h;
                t = h * mOptions.mInc;
                t = SGMIN(t, mOptions.mHUB);
                return t;
            };

            auto dec = [this](FT h) {
                FT t = h;
                t = h * mOptions.mDec;
                t = SGMAX(t, mOptions.mHLB);
                return t;
            };

            auto step = [&] () {
                //std::cout << "\n*** Step " << StepNumber << " ***\n";
                bool isStepSuccessful = false;
                const FT h = sft;
                int Unsuccess = 0;
                const double e = 2.718281828;
                
                    FT delta_f[mOptions.numbOfPoints];
                    FT delta_x[n];
                    FT xtmp[n];
                    main_dir = new FT[n];
                    snowgoose::VecUtils::vecSet(n, 0., main_dir);

                    for (int i = 0; i < mOptions.numbOfPoints; i++) 
                    {
                        for (int j = 0; j < n; j++)
                        {
                            delta_x[j] = x[j] + dirs[i * n + j] * h;
                        }
                        delta_f[i] = obj->func(delta_x);
                        delta_f[i] = delta_f[i] - fcur;
                        for (int j = 0; j < n; j++)
                        {   
                            main_dir[j] += dirs[i * n + j] * delta_f[i];
                        }
                    }
                    snowgoose::VecUtils::vecMult(n, main_dir, (- 1.0 / h), main_dir);
                    normalize();
                    snowgoose::VecUtils::vecSaxpy(n, x, main_dir, h, xtmp);
                    FT fn = obj->func(xtmp);
                    if (fn < fcur)
                    {
                        isStepSuccessful = true;
                        snowgoose::VecUtils::vecCopy(n, xtmp, x);
                        fcur = fn;
                    }
                return isStepSuccessful;
            };
           
            while (!br) {

                direction(mOptions.numbOfPoints);
                
                bool success = step();
 
                StepNumber++;
                if (mOptions.mDoTracing) {
                    std::cout << (success ? "Success" : "Not success") << std::endl;
                    std::cout << "f =" << fcur << std::endl;
                    std::cout << "sft =" << sft << std::endl;
                }

                if (!success) {
                        if (sft > mOptions.minStep) 
                            {
                                sft = dec(sft);
                            } 
                            else
                            {
                              br = true;
                            }
                }    
                if (StepNumber >= mOptions.maxStepNumber) {
                    br = true;
                }
                
                for (auto s : mStoppers) {
                    if (s(fcur, x, StepNumber)) {
                        br = true;
                        break;
                    }
                }
            }
            v = fcur;
            delete [] dirs;
            return rv;
        }

        std::string about() const {
            std::ostringstream os;
            return os.str();
        }

        /**
         * Retrieve options
         * @return options
         */
        Options & getOptions() {
            return mOptions;
        }

        /**
         * Retrieve stoppers vector reference
         * @return stoppers vector reference
         */
        std::vector<Stopper>& getStoppers() {
            return mStoppers;
        }

        /**
         * Get watchers' vector
         * @return watchers vector
         
        std::vector<Watcher>& getWatchers() {
            return mWatchers;
        }*/

    private:

        const COMPI::MPProblem<FT>& mProblem;
        Options mOptions;
        std::vector<Stopper> mStoppers;
        //std::vector<Watcher> mWatchers;
        std::unique_ptr<LineSearch<FT>> mLS;
        
        void printArray(int n, FT * array) {
            std::cout << " dirs = ";
            std::cout << snowgoose::VecUtils::vecPrint(n, array) << std::endl;
        }
        
        void printVector(int n, std::vector<FT> vector) {
            std::cout << " dirs = ";
            for (int i = 0; i < n; i++) {
                std::cout << vector[i] << ", ";
            }
            std::cout << " ]" << std::endl;
        }
    };
}

#endif /* BEST_POINT_METHOD_HPP */

