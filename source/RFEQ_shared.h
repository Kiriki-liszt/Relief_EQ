//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/vst/vsttypes.h"
// #include "pluginterfaces/base/futils.h"

#define _USE_MATH_DEFINES
#include <cmath>     // Decibel, Kaiser-Bessel
#include <vector>

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

static SMTG_CONSTEXPR int bandIn    = 0;
static SMTG_CONSTEXPR int bandHz    = 1;
static SMTG_CONSTEXPR int bandQ     = 2;
static SMTG_CONSTEXPR int banddB    = 3;
static SMTG_CONSTEXPR int bandType  = 4;
static SMTG_CONSTEXPR int bandOrder = 5;
static SMTG_CONSTEXPR int bandNum   = 5;
static SMTG_CONSTEXPR int bandSize  = 6;
typedef ParamValue ParamBand_Array[bandSize];

class SVF {
public:
    static constexpr int tBell         = 0;
    static constexpr int tLowShelf     = 1;
    static constexpr int tHighShelf    = 2;
    static constexpr int tLowShelfHiQ  = 3;
    static constexpr int tHighShelfHiQ = 4;
    static constexpr int tLowPass      = 5;
    static constexpr int tHighPass     = 6;
    static constexpr int tNum          = 6;

    static constexpr int o6dBoct  = 0;
    static constexpr int o12dBoct = 1;
    static constexpr int o18dBoct = 2;
    static constexpr int o24dBoct = 3;
    static constexpr int oNum     = 3;

    static constexpr double _24dBoct_1 = 1.08239220029239402443; // sqrt(2) * sqrt(2 - sqrt(2))
    static constexpr double _24dBoct_2 = 2.61312592975275315155; // sqrt(2) * sqrt(2 + sqrt(2))

    typedef struct ds{
        int    In = 0;
        double Hz = 1000.0;
        double Q = 1.0;
        double dB = 0.0;
        int    Type = tBell;
        int    Order = o12dBoct;
        double Fs = 48000.0;

        double w = Hz * M_PI / Fs;;
        double g = std::tan(w);
        double k = 2.0 / Q;
        double gt0 = 1 / (1 + g * (g + k));
        double gk0 = (g + k) * gt0;

        double m0 = 1.0, m1 = 0.0, m2 = 1.0;
        double v0 = 0.0, v1 = 0.0, v2 = 0.0;
        double t0 = 0.0, t1 = 0.0, t2 = 0.0;
        double ic1eq = 0.0;
        double ic2eq = 0.0;
    } dataset;

    SVF()
    {
        initSVF();
    };
    
    SVF(const SVF& svf)
    {
        setSVF (static_cast<double>(svf.filter_1.In),
                svf.filter_1.Hz,
                svf.filter_1.Q,
                svf.filter_1.dB,
                static_cast<double>(svf.filter_1.Type),
                static_cast<double>(svf.filter_1.Order),
                svf.filter_1.Fs);
        initSVF ();
    };
    
    void setIn (double v) {filter_1.In = v == 0.0 ? 0 : 1; filter_2.In = v == 0.0 ? 0 : 1;}
    void setHz (double v) {filter_1.Hz = v; filter_2.Hz = v;}
    void setQ  (double v) {filter_1.Q = v; filter_2.Q = v;}
    void setdB (double v) {filter_1.dB = v; filter_2.dB = v;}
    void setType (double v) {filter_1.Type = static_cast<int>(v); filter_2.Type = static_cast<int>(v);}
    void setOrder (double v) {filter_1.Order = static_cast<int>(v); filter_2.Order = static_cast<int>(v);}

    inline void initSVF() {
        filter_1.ic1eq = 0.0; filter_1.ic2eq = 0.0;
        filter_2.ic1eq = 0.0; filter_2.ic2eq = 0.0;
    };

    inline void setSVF (double plainIn, double plainHz, double plainQ, double plaindB, double plainType, double plainOrder, double plainFs)
    {
        filter_1.In    = plainIn == 0 ? 0 : 1;
        filter_1.Hz    = plainHz;
        filter_1.Q     = plainQ;
        filter_1.dB    = plaindB;
        filter_1.Type  = static_cast<int>(plainType);
        filter_1.Order = static_cast<int>(plainOrder);
        filter_1.Fs    = plainFs;

        filter_2.In    = filter_1.In;
        filter_2.Fs    = filter_1.Fs;
        filter_2.Hz    = filter_1.Hz;
        filter_2.Q     = filter_1.Q;
        filter_2.dB    = filter_1.dB;
        filter_2.Type  = filter_1.Type;
        filter_2.Order = filter_1.Order;

        if (static_cast<int>(plainType) == tLowPass ||
            static_cast<int>(plainType) == tHighPass)
        {
            switch (static_cast<int>(plainOrder)) {
                case o6dBoct:                                                        makeSVF(&filter_1); break;
                case o12dBoct: filter_1.Q = M_SQRT2;                                 makeSVF(&filter_1); break;
                case o18dBoct: filter_1.Q = (2.0);        filter_2.Order = o6dBoct;  makeSVF(&filter_1); makeSVF(&filter_2); break;
                case o24dBoct: filter_1.Q = (_24dBoct_1); filter_2.Q = (_24dBoct_2); makeSVF(&filter_1); makeSVF(&filter_2); break;
                default:                                                             makeSVF(&filter_1); break;
            }
        }
        else {
            makeSVF(&filter_1);
        }
    }

    static SMTG_CONSTEXPR double dBpow2log = M_LN10 / 40.0; // loge(10) / 40
    
    inline void makeSVF(dataset* filter)
    {
        if (filter->Hz > filter->Fs / 2.0) filter->Hz = filter->Fs / 2.0;
        filter->w = filter->Hz * M_PI / filter->Fs;
        filter->g = std::tan(filter->w);
        filter->k = 2.0 / filter->Q;
        double A = std::pow(10.0, filter->dB / 40.0);

        // 0.05756462732 = M_LN10 / 40.0;
        // why not M_LN10 / 20.0? because it's A = std::pow(10.0, filter->dB / 40.0), not pow(10, dB / 20.0)
        double mm = std::exp(-dBpow2log * std::abs(filter->dB));
        double bk = 1 / (filter->Q * mm);
        double s = M_SQRT2 / std::log2(filter->Q * 0.5 + 1);

        double kdA = bk / A;
        double kmA = bk * A;
        double smA = s * A;
        //double gdA = filter->g / A;
        //double gmA = filter->g * A;
        double gdSA = filter->g / std::sqrt(A);
        double gmSA = filter->g * std::sqrt(A);
        double AmA = A * A;

        switch (filter->Type)
        {
            case tLowPass:      filter->m0 = 0;   filter->m1 = 0;   filter->m2 = 1;   break;
            case tHighPass:     filter->m0 = 1;   filter->m1 = 0;   filter->m2 = 0;   break;
            case tBell:         filter->m0 = 1;   filter->m1 = kmA; filter->m2 = 1;   filter->g = filter->g; filter->k = kdA;   break;
            case tLowShelf:     filter->m0 = 1;   filter->m1 = 0;   filter->m2 = AmA; filter->g = filter->g; filter->k = 1 - filter->g; break;
            case tHighShelf:    filter->m0 = AmA; filter->m1 = 0;   filter->m2 = 1;   filter->g = filter->g; filter->k = 1 - filter->g; break;
            case tLowShelfHiQ:  filter->m0 = 1;   filter->m1 = smA; filter->m2 = AmA; filter->g = gdSA;      filter->k = s;     break;
            case tHighShelfHiQ: filter->m0 = AmA; filter->m1 = smA; filter->m2 = 1;   filter->g = gmSA;      filter->k = s;     break;
            default: break;
        }

        filter->gt0 = 1.0 / (1.0 + filter->g * (filter->g + filter->k));
        filter->gk0 = (filter->g + filter->k) * filter->gt0;
        return;
    };

    inline double _6dBoct_tick(dataset* filter, double vin) {
        // disable v1 stage
        filter->t0 = vin - filter->ic2eq;
        filter->v0 = filter->t0 / (1.0 + filter->g);// gt0 * t0;
        filter->t2 = filter->g * filter->v0;
        filter->v2 = filter->ic2eq + filter->t2;
        filter->ic2eq += 2.0 * filter->t2;

        if (filter->In != 1) return vin;
        return filter->m0 * filter->v0 + filter->m2 * filter->v2;
    }
    inline double _12dBoct_tick(dataset* filter, double vin) {
        // tick serial(possibly quicker on cpus with low latencies)
        filter->t0 = vin - filter->ic2eq;
        filter->v0 = filter->gt0 * filter->t0 - filter->gk0 * filter->ic1eq; // high
        filter->t1 = filter->g * filter->v0;
        filter->v1 = filter->ic1eq + filter->t1; // band
        filter->t2 = filter->g * filter->v1;
        filter->v2 = filter->ic2eq + filter->t2; // low
        filter->ic1eq += 2.0 * filter->t1;
        filter->ic2eq += 2.0 * filter->t2;

        if (filter->In != 1) return vin;
        return filter->m0 * filter->v0 + filter->m1 * filter->v1 + filter->m2 * filter->v2;
    }

    inline double computeSVF (double vin)
    {
        if (filter_1.Type == tLowPass ||
            filter_1.Type == tHighPass)
        {
            switch (filter_1.Order) {
                case o6dBoct:  return  _6dBoct_tick(&filter_1, vin); break;
                case o12dBoct: return _12dBoct_tick(&filter_1, vin); break;
                case o18dBoct: return  _6dBoct_tick(&filter_2, _12dBoct_tick(&filter_1, vin)); break;
                case o24dBoct: return _12dBoct_tick(&filter_2, _12dBoct_tick(&filter_1, vin)); break;
                default:       return  _6dBoct_tick(&filter_1, vin); break;
            }
        }

        if (filter_1.Type == tLowShelf  ||
            filter_1.Type == tHighShelf)
        {
            return _6dBoct_tick(&filter_1, vin);
        }

        return _12dBoct_tick(&filter_1, vin);
    };

    inline double mag_response(dataset* filter, double freq) {
        if (!filter->In) return 1.0;

        double ONE_OVER_SAMPLE_RATE = 1.0 / filter->Fs;

        // exp(complex(0.0, -2.0 * pi) * frequency / sampleRate)
        double _zr = (0.0) * freq * ONE_OVER_SAMPLE_RATE;
        double _zi = (-2.0 * M_PI) * freq * ONE_OVER_SAMPLE_RATE;

        // z = zr + zi;
        double zr = std::exp(_zr) * std::cos(_zi);
        double zi = std::exp(_zr) * std::sin(_zi);

        double nr = 0, ni = 0;
        double dr = 0, di = 0;

        if ((filter->Type == tLowShelf || filter->Type == tHighShelf) ||
            (filter->Type == tLowPass  || filter->Type == tHighPass) && (filter->Order == o6dBoct))
        {
            // Numerator complex
            nr = zr * (-filter->m0 /* + m1 * (g - 1) */ + filter->m2 * filter->g) + (filter->m0 /* + m1 * (g + 1) */ + filter->m2 * filter->g);
            ni = zi * (-filter->m0 /* + m1 * (g - 1) */ + filter->m2 * filter->g);

            // Denominator complex
            dr = zr * (filter->g - 1) + (filter->g + 1);
            di = zi * (filter->g - 1);
        }
        else {
            // z * z
            double zsq_r = zr * zr - zi * zi;
            double zsq_i = zi * zr + zr * zi;
            double gsq = filter->g * filter->g;

            // Numerator complex
            double c_nzsq = (filter->m0 + filter->m1 * filter->g + filter->m2 * gsq);
            double c_nz = (filter->m0 * -2.0 + filter->m2 * 2.0 * gsq);
            double c_n = (filter->m0 + filter->m1 * -filter->g + filter->m2 * gsq);
            nr = zsq_r * c_nzsq + zr * c_nz + c_n;
            ni = zsq_i * c_nzsq + zi * c_nz;

            // Denominator complex
            double c_dzsq = (1.0 + filter->k * filter->g + gsq);
            double c_dz = (-2.0 + 2.0 * gsq);
            double c_d = (1.0 + filter->k * -filter->g + gsq);
            dr = zsq_r * c_dzsq + zr * c_dz + c_d;
            di = zsq_i * c_dzsq + zi * c_dz;
        }

        // Numerator / Denominator
        double norm = dr * dr + di * di;
        double ddr = (nr * dr + ni * di) / norm;
        double ddi = (ni * dr - nr * di) / norm;

        return std::sqrt(ddr * ddr + ddi * ddi);
    }

    inline double mag_response(double freq) {
        if (filter_1.Type == tLowPass ||
            filter_1.Type == tHighPass)
        {
            switch (filter_1.Order) {
            case o6dBoct:  return mag_response(&filter_1, freq); break;
            case o12dBoct: return mag_response(&filter_1, freq); break;
            case o18dBoct: return mag_response(&filter_1, freq) * mag_response(&filter_2, freq); break;
            case o24dBoct: return mag_response(&filter_1, freq) * mag_response(&filter_2, freq); break;
            default:       return mag_response(&filter_1, freq); break;
            }
        }
        return mag_response(&filter_1, freq);
    }

// private: // will it inline?
    dataset filter_1;
    dataset filter_2;
};


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
            A[j] = (std::sin(2.0 * j * M_PI * Fb / Fs) - std::sin(2.0 * j * M_PI * Fa / Fs)) / (j * M_PI);
        
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
            dest[n] = Ino(Alpha * std::sqrt(1.0 - (2.0 * (double)n / (double)N - 1.0) * (2.0 * (double)n / (double)N - 1.0))) / Inoalpha;
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

// zoom title-value struct
struct ZoomFactor
{
    const Steinberg::tchar* title;
    double factor;

    ZoomFactor(const Steinberg::tchar* title, double factor) : title(title), factor(factor) {}
};
typedef std::vector<ZoomFactor> ZoomFactorVector;
static const ZoomFactorVector zoomFactors {ZoomFactor(STR("50%"),   0.50),
                                           ZoomFactor(STR("75%"),   0.75),
                                           ZoomFactor(STR("100%"),  1.00),
                                           ZoomFactor(STR("125%"),  1.25),
                                           ZoomFactor(STR("150%"),  1.50),
                                           ZoomFactor(STR("175%"),  1.75),
                                           ZoomFactor(STR("200%"),  2.00)};
static SMTG_CONSTEXPR int32 zoomNum = 6;
static SMTG_CONSTEXPR int32 dftZoom = 2;

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

}
