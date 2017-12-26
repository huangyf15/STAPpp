/*==============================================================================
                                    MPM3D++
           C++ code for Three Dimensional Material Point Method
  ==============================================================================

   Copyright (C) 2006 -

   Computational Dynamics Group
   Department of Engineering Mechanics
   Tsinghua University
   Beijing 100084, P. R. China

   Email: xzhang@tsinghua.edu.cn

  ==============================================================================
                      Implementation of class 'Clock'
  ==============================================================================*/

#include "Clock.h"

std::chrono::time_point<std::chrono::high_resolution_clock> getClock()
{
    return std::chrono::high_resolution_clock::now();
}

// Constructor
Clock::Clock()
{
    ct_ = 0;
    st0_ = st1_ = false;
}

// Start the clock
void Clock::Start()
{
    t0_ = getClock();
    st0_ = true;
}

// Stop the clock
void Clock::Stop()
{
    if (!st0_)
    {
        cerr << "\n*** Error *** In Clock :: Stop()";
        cerr << " : Method Start() must have been called before.\n";
    }

    if (!st1_)
    {
        t1_ = getClock();
        ct_ += ((std::chrono::duration<double>)(t1_ - t0_)).count();
        st1_ = true;
    }
}

// Resume the stopped clock
void Clock::Resume()
{
    if (!st0_)
    {
        cerr << "\n*** Error *** In Clock :: Resume()";
        cerr << " : Method Start() must have been called before.\n";
    }

    if (!st1_)
    {
        cerr << "\n*** Error *** In Clock::Resume()";
        cerr << " : Method Stop() must have been called before.\n";
    }
    else
    {
        t0_ = getClock();
        st1_ = false;
    }
}

// Clear the clock
void Clock::Clear()
{
    ct_ = 0;
    st0_ = st1_ = false;
}

// Return the elapsed time since the clock started
double Clock::ElapsedTime()
{
    double elapsed = 0.0f;

    if (!st0_)
    {
        cerr << "\n*** Error *** In Clock :: ElapsedTime()";
        cerr << " : Method Start() must have been called before.\n";
    }

    if (st1_) // Timer has been stopped.
        elapsed = ct_;
    else
    {
        t1_ = getClock();
        elapsed = ct_ + ((std::chrono::duration<double>)(t1_ - t0_)).count();
    }

    return elapsed;
}
