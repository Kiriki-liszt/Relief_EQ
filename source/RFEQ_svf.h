#pragma once
#include "pluginterfaces/base/futils.h"
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi             */
#endif
#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880   // sqrt(2)
#endif


typedef enum {
	overSample_1x,
	overSample_2x,
	overSample_4x,
	overSample_8x,
	overSample_num = 3
} overSample;

typedef enum {
	ParamArray_In,
	ParamArray_Hz,
	ParamArray_Q,
	ParamArray_dB,
	ParamArray_Type,
	ParamArray_Order,
	ParamArray_size
} ParamArray;

// In, Hz, Q, dB, Type, Order
typedef double ParamBand_Array[ParamArray_size];

/*

Filters will have six parameters;
1. In/out
2. Frequency
3. Q    - 25.6 to 0.1
4. Gain - +15dB to -15dB
5. Type
6. Order

For Pass filters, Q parameter will be replaced with Order parameter.



There will be seven types of filters;
1. Bell
2. Low Shelf
3. High Shelf
4. Low Shelf High Q
5. High Shelf High Q
6. Low Pass
7. High Pass




Each SVF filter will have these values stored;

* Fs - Sample rate working on

* Hz
* Q
* dB
* Type

* gt0 - Precalculated values
* gt1
* gt2 
* gk0
* gk1

* m0 - coefficients for multipliing with each v outs
* m1
* m2

* v0 - outout signal for each LP, BP, HP
* v1
* v2

* t0 - biquad coefficients
* t1
* t2

* ic1eq - states 
* ic2eq



* Each SVF filter will have these functions;

void initSVF() - init states
void copySVF(SVF* src)
void setSVF() - parameters to values
void makeSVF() - values to coefficients
double computeSVF(double input) 
double mag_response(double freq)

*/
class SVF {
public:
	enum filter_Type
	{
		kBell,
		kLowShelf,
		kHighShelf,
		kLowShelfHiQ,
		kHighShelfHiQ,
		kLowPass,
		kHighPass,

		kFltNum = 6
	};

	enum filter_Order
	{
		_6dBoct,
		_12dBoct,
		_18dBoct,
		_24dBoct,
		kOrderNum = 3
	};

	typedef struct ds{
		int    In = 0;
		double dB = 0.0;
		double Hz = 1000.0;
		double Q = 1.0;
		filter_Type Type = kBell;
		filter_Order Order = _12dBoct;
		double Fs = 48000.0;

		double w = Hz * M_PI / Fs;;
		double g = tan(w);
		double k = 2.0 / Q;
		double gt0 = 1 / (1 + g * (g + k));
		double gk0 = (g + k) * gt0;

		double m0 = 1.0, m1 = 1.0, m2 = 1.0;
		double v0 = 0.0, v1 = 0.0, v2 = 0.0;
		double t0 = 0.0, t1 = 0.0, t2 = 0.0;
		double ic1eq = 0.0;
		double ic2eq = 0.0;
	} dataset;

#define _24dBoct_1 1.08239220029239402443 // sqrt(2) * sqrt(2 - sqrt(2))
#define _24dBoct_2 2.61312592975275315155 // sqrt(2) * sqrt(2 + sqrt(2))

	SVF()
	{
		initSVF();
	};

	void initSVF() { 
		filter_1.ic1eq = 0.0; filter_1.ic2eq = 0.0; 
		filter_2.ic1eq = 0.0; filter_2.ic2eq = 0.0; 
	};

	void copySVF(SVF* src)
	{
		this->filter_1.In    = src->filter_1.In;
		this->filter_1.Fs    = src->filter_1.Fs;
		this->filter_1.Hz    = src->filter_1.Hz;
		this->filter_1.Q     = src->filter_1.Q;
		this->filter_1.dB    = src->filter_1.dB;
		this->filter_1.Type  = src->filter_1.Type;
		this->filter_1.Order = src->filter_1.Order;

		this->filter_2.In    = src->filter_2.In;
		this->filter_2.Fs    = src->filter_2.Fs;
		this->filter_2.Hz    = src->filter_2.Hz;
		this->filter_2.Q     = src->filter_2.Q;
		this->filter_2.dB    = src->filter_2.dB;
		this->filter_2.Type  = src->filter_2.Type;
		this->filter_2.Order = src->filter_2.Order;

		if (filter_1.Type == kLowPass ||
			filter_1.Type == kHighPass)
		{
			switch (filter_1.Order) {
			case _6dBoct:                                                        makeSVF(&filter_1); break;
			case _12dBoct: filter_1.Q = M_SQRT2;                                 makeSVF(&filter_1); break;
			case _18dBoct: filter_1.Q = (2.0);        filter_2.Order = _6dBoct;  makeSVF(&filter_1); makeSVF(&filter_2); break;
			case _24dBoct: filter_1.Q = (_24dBoct_1); filter_2.Q = (_24dBoct_2); makeSVF(&filter_1); makeSVF(&filter_2); break;
			default:                                                             makeSVF(&filter_1); break;
			}
		}
		else {
			makeSVF(&filter_1);
		}

		return;
	}

	void setSVF(double fParamIn, double fParamHz, double fParamQ, double fParamdB, double fParamtype, double fParamOrder, double fParamFs)
	{
		filter_1.In    = fParamIn ? 1 : 0;
		filter_1.Fs    = fParamFs;
		filter_1.Hz    = _norm_to_Hz(fParamHz);
		filter_1.Q     = _norm_to_Q(fParamQ);
		filter_1.dB    = _norm_to_dB(fParamdB);
		filter_1.Type  = _norm_to_Type(fParamtype);
		filter_1.Order = _norm_to_Order(fParamOrder);

		filter_2.In    = filter_1.In;
		filter_2.Fs    = filter_1.Fs;
		filter_2.Hz    = filter_1.Hz;
		filter_2.Q     = filter_1.Q;
		filter_2.dB    = filter_1.dB;
		filter_2.Type  = filter_1.Type;
		filter_2.Order = filter_1.Order;

		if (filter_1.Type == kLowPass ||
			filter_1.Type == kHighPass)
		{
			switch (filter_1.Order) {
			case _6dBoct:                                                        makeSVF(&filter_1); break;
			case _12dBoct: filter_1.Q = M_SQRT2;                                 makeSVF(&filter_1); break;
			case _18dBoct: filter_1.Q = (2.0);        filter_2.Order = _6dBoct;  makeSVF(&filter_1); makeSVF(&filter_2); break;
			case _24dBoct: filter_1.Q = (_24dBoct_1); filter_2.Q = (_24dBoct_2); makeSVF(&filter_1); makeSVF(&filter_2); break;
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
		case kLowPass:      filter->m0 = 0;   filter->m1 = 0;   filter->m2 = 1;   break;
		case kHighPass:     filter->m0 = 1;   filter->m1 = 0;   filter->m2 = 0;   break;
		case kBell:         filter->m0 = 1;   filter->m1 = kmA; filter->m2 = 1;   filter->g = filter->g; filter->k = kdA;   break;
		case kLowShelf:     filter->m0 = 1;   filter->m1 = 0;   filter->m2 = AmA; filter->g = filter->g; filter->k = 1 - filter->g; break;
		case kHighShelf:    filter->m0 = AmA; filter->m1 = 0;   filter->m2 = 1;   filter->g = filter->g; filter->k = 1 - filter->g; break;
		case kLowShelfHiQ:  filter->m0 = 1;   filter->m1 = smA; filter->m2 = AmA; filter->g = gdSA;      filter->k = s;     break;
		case kHighShelfHiQ: filter->m0 = AmA; filter->m1 = smA; filter->m2 = 1;   filter->g = gmSA;      filter->k = s;     break;
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
		if (filter_1.Type == kLowPass || 
			filter_1.Type == kHighPass) 
		{
			switch (filter_1.Order) {
			case _6dBoct:  return  _6dBoct_tick(&filter_1, vin); break;
			case _12dBoct: return _12dBoct_tick(&filter_1, vin); break;
			case _18dBoct: return  _6dBoct_tick(&filter_2, _12dBoct_tick(&filter_1, vin)); break;
			case _24dBoct: return _12dBoct_tick(&filter_2, _12dBoct_tick(&filter_1, vin)); break;
			default:       return  _6dBoct_tick(&filter_1, vin); break;
			}
		}

		if (filter_1.Type == kLowShelf  ||
			filter_1.Type == kHighShelf) 
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

		if ((filter->Type == kLowShelf || filter->Type == kHighShelf) ||
			(filter->Type == kLowPass || filter->Type == kHighPass) && (filter->Order == _6dBoct))
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
		if (filter_1.Type == kLowPass ||
			filter_1.Type == kHighPass)
		{
			switch (filter_1.Order) {
			case _6dBoct:  return mag_response(&filter_1, freq); break;
			case _12dBoct: return mag_response(&filter_1, freq); break;
			case _18dBoct: return mag_response(&filter_1, freq) * mag_response(&filter_2, freq); break;
			case _24dBoct: return mag_response(&filter_1, freq) * mag_response(&filter_2, freq); break;
			default:       return mag_response(&filter_1, freq); break;
			}
		}
		return mag_response(&filter_1, freq);
	}


	static double getFreqMax() { return 41000.0; };
	static double getFreqMin() { return 20.0; };
	static double getdBMax()   { return  12.0; };
	static double getdBMin()   { return -12.0; };
	static double getQMax()    { return 25.6; };
	static double getQMin()    { return 0.1; };

	static double Init_Band1_Hz() { return 80.0; };
	static double Init_Band2_Hz() { return 200.0; };
	static double Init_Band3_Hz() { return 2000.0; };
	static double Init_Band4_Hz() { return 6000.0; };
	static double Init_Band5_Hz() { return 16000.0; };

	static filter_Type _norm_to_Type(double value)
	{
		return static_cast<filter_Type>(Steinberg::FromNormalized<double>(value, kFltNum));
	};
	static double _Type_to_norm(filter_Type type)
	{
		return Steinberg::ToNormalized<double>(type, kFltNum);
	};

	static filter_Order _norm_to_Order(double value)
	{
		return static_cast<filter_Order>(Steinberg::FromNormalized<double>(value, kOrderNum));
	};
	static double _Order_to_norm(filter_Order order)
	{
		return Steinberg::ToNormalized<double>(order, kOrderNum);
	};


	static double _dB_to_norm(double _dB)
	{
		return (_dB + getdBMax()) / (getdBMax() - getdBMin());
	};
	static double _norm_to_dB(double paramValue)
	{
		return ((getdBMax() - getdBMin()) * paramValue) + getdBMin();
	};

	static double _Hz_to_norm(double _Hz)
	{
		return log(_Hz / getFreqMin()) / log(getFreqMax() / getFreqMin());
	};
	static double _norm_to_Hz(double paramValue)
	{
		double FREQ_LOG_MAX = log(getFreqMax() / getFreqMin());
		double tmp = getFreqMin() * exp(FREQ_LOG_MAX * paramValue);
		return std::max(std::min(tmp, getFreqMax()), getFreqMin());
	};

	static double _Q_to_norm(double _Q)
	{
		return log(_Q / getQMin()) / log(getQMax() / getQMin());
	};
	static double _norm_to_Q(double paramValue)
	{
		double Q_LOG_MAX = log(getQMax() / getQMin());
		double tmq = getQMin() * exp(Q_LOG_MAX * paramValue);
		return std::max(std::min(tmq, getQMax()), getQMin());
	};

private:
	dataset filter_1;
	dataset filter_2;
};
