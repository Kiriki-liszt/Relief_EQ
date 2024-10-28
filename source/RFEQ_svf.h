#pragma once
#include "pluginterfaces/base/futils.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>

static constexpr int bandIn    = 0;
static constexpr int bandHz    = 1;
static constexpr int bandQ     = 2;
static constexpr int banddB    = 3;
static constexpr int bandType  = 4;
static constexpr int bandOrder = 5;
static constexpr int bandNum   = bandOrder+1;
typedef double ParamBand_Array[bandNum];

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
        double g = tan(w);
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
        setSVF(static_cast<double>(svf.filter_1.In), svf.filter_1.Hz, svf.filter_1.Q, svf.filter_1.dB, static_cast<double>(svf.filter_1.Type), static_cast<double>(svf.filter_1.Order), svf.filter_1.Fs);
        initSVF();
    };
    
    void setIn (double v) {filter_1.In = v == 0.0 ? 0 : 1; filter_2.In = v == 0.0 ? 0 : 1;}
    void setHz (double v) {filter_1.Hz = v; filter_2.Hz = v;}
    void setQ  (double v) {filter_1.Q = v; filter_2.Q = v;}
    void setdB (double v) {filter_1.dB = v; filter_2.dB = v;}
    void setType (double v) {filter_1.Type = static_cast<int>(v); filter_2.Type = static_cast<int>(v);}
    void setOrder (double v) {filter_1.Order = static_cast<int>(v); filter_2.Order = static_cast<int>(v);}

    void initSVF() {
        filter_1.ic1eq = 0.0; filter_1.ic2eq = 0.0;
        filter_2.ic1eq = 0.0; filter_2.ic2eq = 0.0;
    };

    void setSVF (double plainIn, double plainHz, double plainQ, double plaindB, double plainType, double plainOrder, double plainFs)
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

    void makeSVF(dataset* filter)
    {
        if (filter->Hz > filter->Fs / 2.0) filter->Hz = filter->Fs / 2.0;
        filter->w = filter->Hz * M_PI / filter->Fs;
        filter->g = tan(filter->w);
        filter->k = 2.0 / filter->Q;
        double A = pow(10.0, filter->dB / 40.0);

        double mm = exp(-0.0575 * abs(filter->dB));
        double bk = 1 / (filter->Q * mm);
        double s = M_SQRT2 / log2(filter->Q * 0.5 + 1);

        double kdA = bk / A;
        double kmA = bk * A;
        double smA = s * A;
        //double gdA = filter->g / A;
        //double gmA = filter->g * A;
        double gdSA = filter->g / sqrt(A);
        double gmSA = filter->g * sqrt(A);
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

    double _6dBoct_tick(dataset* filter, double vin) {
        // disable v1 stage
        filter->t0 = vin - filter->ic2eq;
        filter->v0 = filter->t0 / (1.0 + filter->g);// gt0 * t0;
        filter->t2 = filter->g * filter->v0;
        filter->v2 = filter->ic2eq + filter->t2;
        filter->ic2eq += 2.0 * filter->t2;

        if (filter->In != 1) return vin;
        return filter->m0 * filter->v0 + filter->m2 * filter->v2;
    }
    double _12dBoct_tick(dataset* filter, double vin) {
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

    double computeSVF (double vin)
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

    double mag_response(dataset* filter, double freq) {
        if (!filter->In) return 1.0;

        double ONE_OVER_SAMPLE_RATE = 1.0 / filter->Fs;

        // exp(complex(0.0, -2.0 * pi) * frequency / sampleRate)
        double _zr = (0.0) * freq * ONE_OVER_SAMPLE_RATE;
        double _zi = (-2.0 * M_PI) * freq * ONE_OVER_SAMPLE_RATE;

        // z = zr + zi;
        double zr = exp(_zr) * cos(_zi);
        double zi = exp(_zr) * sin(_zi);

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

        return sqrt(ddr * ddr + ddi * ddi);
    }

    double mag_response(double freq) {
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

private:
    dataset filter_1;
    dataset filter_2;
};
