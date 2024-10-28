//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/vst/vsttypes.h"
// #include "pluginterfaces/base/futils.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <algorithm>
#include <numeric>

namespace yg331 {
//------------------------------------------------------------------------
using ParamValue    = Steinberg::Vst::ParamValue;
using SampleRate    = Steinberg::Vst::SampleRate;
using int32         = Steinberg::int32;
using uint32        = Steinberg::uint32;
using TBool         = Steinberg::TBool;

static SMTG_CONSTEXPR int overSample_1x  = 0;
static SMTG_CONSTEXPR int overSample_2x  = 1;
static SMTG_CONSTEXPR int overSample_4x  = 2;
static SMTG_CONSTEXPR int overSample_8x  = 3;
static SMTG_CONSTEXPR int overSample_num = 3;

static SMTG_CONSTEXPR int detectorPeak = 0;
static SMTG_CONSTEXPR int detectorRMS  = 1;
static SMTG_CONSTEXPR int detectorNum  = 1;

//------------------------------------------------------------------------
//  Class for converter
//------------------------------------------------------------------------
class DecibelConverter
{
public:
    static SMTG_CONSTEXPR double log2dB = 20.0 / M_LN10;
    static SMTG_CONSTEXPR double dB2log = M_LN10 / 20.0;
    static inline ParamValue ToGain (ParamValue dB)
    {
        // return std::pow(10.0, dB * 0.05);
        return std::exp(dB * dB2log);
    }
    
    static inline ParamValue ToDecibel (ParamValue gain)
    {
        return (gain > 0.0) ? (log2dB * std::log(gain)) : (-140.0);
    }
private:
    DecibelConverter() = delete;
};

class ParameterConverter
{
public:
    static constexpr int range = 0;
    static constexpr int log = 1;
    static constexpr int list = 2;
    
    ParameterConverter (ParamValue minValue, ParamValue maxValue, int type = range, int32 numSteps = -1)
    : minValue (minValue), maxValue (maxValue), type(type), numSteps (numSteps) {};
    
    ParamValue ToNormalized (ParamValue plain) const
    {
        switch (type) {
            case range: return ToNormalizedRange(plain); break;
            case log:   return ToNormalizedLog(plain);   break;
            case list:  return ToNormalizedList(plain);  break;
            default: break;
        }
        return ToNormalizedRange(plain);
    }
    
    ParamValue ToPlain (ParamValue normalized) const
    {
        switch (type) {
            case range: return ToPlainRange(normalized); break;
            case log:   return ToPlainLog(normalized);   break;
            case list:  return ToPlainList(normalized);  break;
            default: break;
        }
        return ToNormalizedRange(normalized);
    }
    
    ParamValue ToNormalizedRange (ParamValue plain) const
    {
        if (maxValue == minValue) return 0.0;
        return (plain - minValue) / (maxValue - minValue);
    }
    
    ParamValue ToPlainRange (ParamValue normalized) const
    {
        return normalized * (maxValue - minValue) + minValue;
    }
    
    /* Log Conversion
     
     x - x0    log(y) - log(y0)
     ------- = -----------------
     x1 - x0   log(y1) - log(y0)
     
     x0 = 0, x1 = 1 : normalized parameter
     
     norm  = ( log(plain) - log(minValue) ) / ( log(maxValue) - log(minValue) )
     = log(plain / minValue) / log(maxValue / minValue)
     
     log(plain / minValue) = norm * log(maxValue / minValue)
     plain / minValue = exp(norm * log(maxValue / minValue))
     plain = minValue * exp(norm * log(maxValue / minValue))
     
     */
    
    ParamValue ToNormalizedLog (ParamValue plain) const
    {
        if (minValue == 0.0) return 0.0;
        if (!((plain    / minValue) > 0.0)) return 0.0;
        if (!((maxValue / minValue) > 0.0)) return 0.0;
        return std::log(plain / minValue) / std::log(maxValue / minValue);
    }
    
    ParamValue ToPlainLog (ParamValue normalized) const
    {
        if (minValue == 0.0) return 0.0;
        if (!((maxValue / minValue) > 0.0)) return 0.0;
        return minValue * std::exp(normalized * std::log(maxValue / minValue));
    }

    ParamValue ToNormalizedList (int32 plain) const
    {
        // return (Steinberg::ToNormalized<type> (v, step));
        if (numSteps <= 0) return 0;
        return plain / ParamValue (numSteps);
    }
    
    int32 ToPlainList (ParamValue normalized) const
    {
        // return (Steinberg::FromNormalized<type> (v, step));
        if (numSteps <= 0) return 0;
        int32 cmp(normalized * (numSteps + 1));
        return numSteps < cmp ? numSteps : cmp;
    }
    
    void setRange (ParamValue _min, ParamValue _max)
    {
        minValue = _min;
        maxValue = _max;
    }
    
    void setSteps (int32 _steps)
    {
        numSteps = _steps;
    }
    
private:
    ParamValue minValue = 0.0;
    ParamValue maxValue = 1.0;
    int type = range;
    int32 numSteps = -1;
};


//------------------------------------------------------------------------
//  Min, Max, Default of Parameters
//------------------------------------------------------------------------
static SMTG_CONSTEXPR ParamValue minParamFreq = 20.0;
static SMTG_CONSTEXPR ParamValue maxParamFreq = 41000.0;
static SMTG_CONSTEXPR ParamValue dftBand1Freq = 80.0;
static SMTG_CONSTEXPR ParamValue dftBand2Freq = 200.0;
static SMTG_CONSTEXPR ParamValue dftBand3Freq = 2000.0;
static SMTG_CONSTEXPR ParamValue dftBand4Freq = 6000.0;
static SMTG_CONSTEXPR ParamValue dftBand5Freq = 16000.0;

static SMTG_CONSTEXPR ParamValue minParamGain = -12.0;
static SMTG_CONSTEXPR ParamValue maxParamGain = 12.0;
static SMTG_CONSTEXPR ParamValue dftParamGain = 0.0;

static SMTG_CONSTEXPR ParamValue minParamQlty = 0.1;
static SMTG_CONSTEXPR ParamValue maxParamQlty = 25.6;
static SMTG_CONSTEXPR ParamValue dftParamQlty = M_SQRT2;

static SMTG_CONSTEXPR int32      dftParamType = SVF::tBell;
static SMTG_CONSTEXPR int32      dftParamOrdr = SVF::o12dBoct;

static const ParameterConverter paramGain      (minParamGain, maxParamGain, ParameterConverter::range);
static const ParameterConverter paramFreq      (minParamFreq, maxParamFreq, ParameterConverter::log);
static const ParameterConverter paramQlty      (minParamQlty, maxParamQlty, ParameterConverter::log);
static const ParameterConverter paramType      (0,            0,            ParameterConverter::list, SVF::tNum);
static const ParameterConverter paramOrdr      (0,            0,            ParameterConverter::list, SVF::oNum);

static const double nrmBand1Freq = paramFreq.ToNormalized(dftBand1Freq);
static const double nrmBand2Freq = paramFreq.ToNormalized(dftBand2Freq);
static const double nrmBand3Freq = paramFreq.ToNormalized(dftBand3Freq);
static const double nrmBand4Freq = paramFreq.ToNormalized(dftBand4Freq);
static const double nrmBand5Freq = paramFreq.ToNormalized(dftBand5Freq);

static const double nrmParamGain = paramGain.ToNormalized(dftParamGain);
static const double nrmParamQlty = paramQlty.ToNormalized(dftParamQlty);
static const double nrmParamType = paramQlty.ToNormalizedList(dftParamType);
static const double nrmParamOrdr = paramOrdr.ToNormalizedList(dftParamOrdr);

static SMTG_CONSTEXPR bool dftBypass          = false;
static SMTG_CONSTEXPR bool dftSoftBypass      = false;

//------------------------------------------------------------------------
//  Kaiser-Bessel Filter designser for Oversampling
//------------------------------------------------------------------------
class Kaiser {
public:
    static constexpr int maxTap = 512;
    
    static inline double Ino(double x)
    {
        double d = 0, ds = 1, s = 1;
        do
        {
            d += 2;
            ds *= x * x / (d * d);
            s += ds;
        } while (ds > s * 1e-100);
        return s;
    };
    
    static void calcFilter(double Fs, double Fa, double Fb, int M, double Att, double dest [])
    {
        // Kaiser windowed FIR filter "DIGITAL SIGNAL PROCESSING, II" IEEE Press pp 123-126.
        
        int Np = (M - 1) / 2;
        double A[maxTap] = { 0, };
        double Alpha; // == pi * alpha (wikipedia)
        double Inoalpha;
        
        A[0] = 2 * (Fb - Fa) / Fs;
        
        for (int j = 1; j <= Np; j++)
            A[j] = (sin(2.0 * j * M_PI * Fb / Fs) - sin(2.0 * j * M_PI * Fa / Fs)) / (j * M_PI);
        
        if (Att < 21.0)
            Alpha = 0;
        else if (Att > 50.0)
            Alpha = 0.1102 * (Att - 8.7);
        else
            Alpha = 0.5842 * pow((Att - 21), 0.4) + 0.07886 * (Att - 21);
        
        Inoalpha = Ino(Alpha);
        
        for (int j = 0; j <= Np; j++)
        {
            dest[Np + j] = A[j] * Ino(Alpha * std::sqrt(1.0 - ((double)(j * j) / (double)(Np * Np)))) / Inoalpha;
        }
        dest[Np + Np] = A[Np] * Ino(0.0) / Inoalpha; // ARM with optimizer level O3 returns NaN == sqrt(1.0 - n/n), while x64 does not...
        for (int j = 0; j < Np; j++)
        {
            dest[j] = dest[M - 1 - j];
        }
    }
    
    static void calcFilter2(int length, double alpha, double dest [])
    {
        int N = length - 1;
        float Alpha = M_PI * alpha;
        float Inoalpha;
        
        Inoalpha = Ino(Alpha);
        
        for (int n = 0; n <= N; n++)
        {
            dest[n] = Ino(Alpha * sqrt(1.0 - (2.0 * (double)n / (double)N - 1.0) * (2.0 * (double)n / (double)N - 1.0))) / Inoalpha;
        }
        dest[0] = Ino(0.0) / Inoalpha; // ARM with optimizer level O3 returns NaN == sqrt(1.0 - n/n), while x64 does not...
        dest[N] = Ino(0.0) / Inoalpha;
        
        float pwr = 0.0;
        for (int i = 0; i < length; i++) {
            pwr += dest[i];
        }
        pwr = 1.0 / pwr;
        for (int i = 0; i < length; i++) {
            dest[i] *= pwr;
        }
    }
};

}
