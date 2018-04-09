#ifndef RANDOM_METHOD_HPP
#define RANDOM_METHOD_HPP

#include <sstream>
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
#include <random>


namespace LOCSEARCH {
    /**
     * Random method implementation
     */
    template <typename FT> class RandomMethod : public COMPI::Solver<FT> {
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
            FT minStep = 0.0001;
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
            FT mHLB = 1e-08;//-1e+02;
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
        RandomMethod(const COMPI::MPProblem<FT>& prob) :
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
            bool br = false;
            
            std::default_random_engine generator;
            std::normal_distribution<FT> distribution(0.0,1.0);
            
            const snowgoose::Box<double>& box = *(mProblem.mBox);
            FT sft = 1.0;
            FT sum = 0.0;

            FT* dirs = new FT[n * mOptions.numbOfPoints];
            auto direction = [&] () {
                for (int j = 0; j < mOptions.numbOfPoints; j++)
                {
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
                std::cout << "\n*** Step " << StepNumber << " ***\n";
                bool isStepSuccessful = false;
                int numb_of_best_vec = -1;
                FT best_f = fcur;
                FT x_best[n];
                bool trigger = false;
                const FT h = sft;

                for (int i = 0; i < mOptions.numbOfPoints; i++) {
                    
                    bool global_continue = false;
                    FT xtmp[n]; 
                    for (int j = 0; j < n; j++)
                    {
                        xtmp[j] = x[j] + dirs[i * n + j] * h;
                        if ((xtmp[j] != SGMAX(xtmp[j], box.mA[j])) || (xtmp[j] != SGMIN(xtmp[j], box.mB[j])))
                        {
                            global_continue = true;
                            break;
                        }  
                    }
                    if (global_continue) 
                    {
                        continue;
                    }
                    
                    FT fn = obj->func(xtmp);
                  
                    if (mOptions.mDoTracing) {
                        printArray("xtmp", n, xtmp);
                    }
                    
                    if (fn < fcur) {
                        FT x_continued[n];
                        for (int q = 0; q < n; q++)
                        {
                            x_continued[q] = x[q] + mOptions.mInc * (xtmp[q] - x[q]);
                            if((x_continued[q] != SGMAX(x_continued[q], box.mA[q])) || (x_continued[q] != SGMIN(x_continued[q], box.mB[q])))
                            {
                                global_continue = true;
                                break;
                            }
                        }
                        if (global_continue) 
                        {
                            continue;
                        }
                        FT f_continued = obj->func(x_continued);
                        if (f_continued < fcur) {
                            if (f_continued < best_f)
                                {
                                    
                                    best_f = f_continued;
                                    snowgoose::VecUtils::vecCopy(n, x_continued, x_best);
                                    numb_of_best_vec = i;
                                } 
                        }
                    } 
                }
                if (numb_of_best_vec != -1)
                {
                    isStepSuccessful = true;
                    snowgoose::VecUtils::vecCopy(n, x_best, x);
                    fcur = best_f;
                }
                return isStepSuccessful;
            };
           
            while (!br) {
                direction();
                bool success = step();
 
                StepNumber++;
                if (mOptions.mDoTracing) {
                    std::cout << (success ? "Success" : "Not success") << std::endl;
                    printArray("x", n, x);
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
                else 
                {
                   sft = inc(sft);
                }
                
                if (StepNumber >= mOptions.maxStepNumber) {
                    br = true;
                }
                /*for (auto w : mWatchers) {
                    w(fcur, x, sft, StepNumber);
                }*/
                
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
        
        void printArray(const char * name, int n, FT * array) {
            std::cout << name << " = ";
            std::cout << snowgoose::VecUtils::vecPrint(n, array) << std::endl;
        }
        
        void printVector(const char * name, int n, std::vector<FT> vector) {
            std::cout << name << " = ";
            std::cout << "[ ";
            for (int i = 0; i < n; i++) {
                std::cout << vector[i] << ", ";
            }
            std::cout << " ]" << std::endl;
        }
    };
}

#endif 