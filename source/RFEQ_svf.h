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

	SVF() :
		dB(0.5), Hz(1.0), Q(1.0), Fs(48000.0),
		gt0(0.0), gk0(0.0),
		m0(0.0), m1(0.0), m2(0.0),
		v0(0.0), v1(0.0), v2(0.0),
		t0(0.0), t1(0.0), t2(0.0),
		ic1eq(0.0), ic2eq(0.0), Type(kBell)
	{
		//setSVF();
		initSVF();
	};

	SVF(double _Hz) :
		dB(0.5), Hz(_Hz), Q(1.0), Fs(192000.0),
		gt0(0.0), gk0(0.0),
		m0(0.0), m1(0.0), m2(0.0),
		v0(0.0), v1(0.0), v2(0.0),
		t0(0.0), t1(0.0), t2(0.0),
		ic1eq(0.0), ic2eq(0.0), Type(kBell)
	{
		//fParamHz = _Hz_to_norm(_Hz);
		//setSVF();
		initSVF();
	};

	void initSVF() { ic1eq = 0.0; ic2eq = 0.0; };

	void copySVF(SVF* src)
	{
		this->In = src->In;
		this->dB = src->dB;
		this->Hz = src->Hz;
		this->Q = src->Q;
		this->Type = src->Type;
		this->Fs = src->Fs;

		this->gt0 = src->gt0;
		this->gk0 = src->gk0;
		this->m0 = src->m0;
		this->m1 = src->m1;
		this->m2 = src->m2;
		return;
	}

	void setSVF(double fParamIn, double fParamFs, double fParamdB, double fParamHz, double fParamQ, double fParamtype, double fParamOrder)
	{
		In = fParamIn ? 1 : 0;

		Fs = fParamFs;

		Hz = _norm_to_Hz(fParamHz);

		Q = _norm_to_Q(fParamQ);

		dB = _norm_to_dB(fParamdB);

		Type = _norm_to_Type(fParamtype);

		Order = _norm_to_Order(fParamOrder);
	}

	void makeSVF()
	{
		
		w = Hz * M_PI / Fs;
		g = tan(w);
		k = 2.0 / Q;
		double A = pow(10.0, dB / 40.0);

		double mm = exp(-0.0575 * abs(dB));
		double bk = 1 / (Q * mm);
		double s = M_SQRT2 / log2(Q * 0.5 + 1);

		double kdA = bk / A;
		double kmA = bk * A;
		double smA = s * A;
		double gdA = g / A;
		double gmA = g * A;
		double gdSA = g / sqrt(A);
		double gmSA = g * sqrt(A);
		double AmA = A * A;

		switch (Type)
		{
		case kLowPass:      m0 = 0;   m1 = 0;   m2 = 1;   break;
		case kHighPass:     m0 = 1;   m1 = 0;   m2 = 0;   break;
		case kBell:         m0 = 1;   m1 = kmA; m2 = 1;   g = g;    k = kdA;   break;
		case kLowShelf:     m0 = 1;   m1 = 0;   m2 = AmA; g = gdA;  break;
		case kHighShelf:    m0 = AmA; m1 = 0;   m2 = 1;   g = gmA;  break;
		case kLowShelfHiQ:  m0 = 1;   m1 = smA; m2 = AmA; g = gdSA; k = s;     break;
		case kHighShelfHiQ: m0 = AmA; m1 = smA; m2 = 1;   g = gmSA; k = s;     break;
		default: break;
		}

		gt0 = 1 / (1 + g * (g + k));
		gk0 = (g + k) * gt0;
		return;
	};


	double computeSVF
	(double vin)
	{
		if (Type == kLowShelf  || 
			Type == kHighShelf || 
			(Type == kLowPass  || Order == _6dBoct) || 
			(Type == kHighPass || Order == _6dBoct)) {

			// disable v1 stage
			t0 = vin - ic2eq;
			v0 = gt0 * t0; 
			t2 = g * v0;
			v2 = ic2eq + t2;
			ic2eq += 2.0 * t2;

			return m0 * v0 + m2 * v2;
		}

		// tick serial(possibly quicker on cpus with low latencies)
		t0 = vin - ic2eq;
		v0 = gt0 * t0 - gk0 * ic1eq; // high
		t1 = g * v0;
		v1 = ic1eq + t1; // band
		t2 = g * v1;
		v2 = ic2eq + t2; // low
		ic1eq += 2.0 * t1;
		ic2eq += 2.0 * t2;

		return m0 * v0 + m1 * v1 + m2 * v2;	
	};

	double mag_response(double freq) {
		if (!In) return 1.0;

		double ONE_OVER_SAMPLE_RATE = 1.0 / Fs;

		// exp(complex(0.0, -2.0 * pi) * frequency / sampleRate)
		double _zr = (0.0) * freq * ONE_OVER_SAMPLE_RATE;
		double _zi = (-2.0 * M_PI) * freq * ONE_OVER_SAMPLE_RATE;

		// z = zr + zi;
		double zr = exp(_zr) * cos(_zi);
		double zi = exp(_zr) * sin(_zi);

		double nr = 0, ni = 0;
		double dr = 0, di = 0;

		if (false/*_6dB != 0*/) {
			// Numerator complex
			nr = zr * (-m0 /* + m1 * (g - 1) */ + m2 * g) + (m0 /* + m1 * (g + 1) */ + m2 * g);
			ni = zi * (-m0 /* + m1 * (g - 1) */ + m2 * g);

			// Denominator complex
			dr = zr * (g - 1) + (g + 1);
			di = zi * (g - 1);
		}
		else {
			// z * z
			double zsq_r = zr * zr - zi * zi;
			double zsq_i = zi * zr + zr * zi;
			double gsq = g * g;

			// Numerator complex
			double c_nzsq = (m0 + m1 * g + m2 * gsq);
			double c_nz = (m0 * -2 + m2 * 2.0 * gsq);
			double c_n = (m0 + m1 * -g + m2 * gsq);
			nr = zsq_r * c_nzsq + zr * c_nz + c_n;
			ni = zsq_i * c_nzsq + zi * c_nz;

			// Denominator complex
			double c_dzsq = (1 + k * g + gsq);
			double c_dz = (-2 + 2.0 * gsq);
			double c_d = (1 + k * -g + gsq);
			dr = zsq_r * c_dzsq + zr * c_dz + c_d;
			di = zsq_i * c_dzsq + zi * c_dz;
		}

		// Numerator / Denominator
		double norm = dr * dr + di * di;
		double ddr = (nr * dr + ni * di) / norm;
		double ddi = (ni * dr - nr * di) / norm;

		return sqrt(ddr * ddr + ddi * ddi);
	}


	static double getFreqMax() { return 41000.0; };
	static double getFreqMin() { return 20.0; };
	static double getdBMax()   { return  12.0; };
	static double getdBMin()   { return -12.0; };
	static double getQMax()    { return 25.6; };
	static double getQMin()    { return 0.5; };

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

	/*
	Steinberg::Vst::ParamValue fParamIn = 1.0;
	Steinberg::Vst::ParamValue fParamdB = _dB_to_norm(0.0);
	Steinberg::Vst::ParamValue fParamHz = _Hz_to_norm(1000.0);
	Steinberg::Vst::ParamValue fParamQ = _Q_to_norm(1.414);
	Steinberg::Vst::ParamValue fParam_6dB = 0.0;
	Steinberg::Vst::ParamValue fParamtype = _type_to_norm(kBell);
	Steinberg::Vst::Sample64   fParamFs = 192000.0;
	*/

	filter_Type Type;
	filter_Order Order;
	int In = 1;

private:
	double dB;
	double Hz;
	double Q;
	double Fs;

	double w, g, k;
	double gt0;
	double gk0;

	double m0, m1, m2;
	double v0, v1, v2;
	double t0, t1, t2;
	double ic1eq;
	double ic2eq;
};
/*
class PassFilter {
public:

	enum filter_degree
	{
		_6dBoct,
		_12dBoct,
		_18dBoct,
		_24dBoct,
		_36dBoct,
		kDegreeNum = 4
	};

	PassFilter(double _Hz, SVF::filter_type _type)
	{
		if (_type == SVF::kHP) fParamHz = _Hz_to_norm_HP(_Hz);
		else fParamHz = _Hz_to_norm_LP(_Hz);
		_6dB_1.fParamHz = fParamHz;
		_12dB_1.fParamHz = fParamHz;
		_12dB_2.fParamHz = fParamHz;
		_12dB_3.fParamHz = fParamHz;

		fParamtype = SVF::_type_to_norm(_type);
		_6dB_1.fParamtype = fParamtype;
		_12dB_1.fParamtype = fParamtype;
		_12dB_2.fParamtype = fParamtype;
		_12dB_3.fParamtype = fParamtype;

		_6dB_1.fParam_6dB = 1;
		_12dB_1.fParam_6dB = 0;
		_12dB_2.fParam_6dB = 0;
		_12dB_3.fParam_6dB = 0;

		fParamFs = 192000.0;
	};

	void initPassFilter() {
		_6dB_1.initSVF();
		_12dB_1.initSVF();
		_12dB_2.initSVF();
		_12dB_3.initSVF();
	};

	void copyPassFilter(PassFilter* src)
	{
		this->fParamdegree = src->fParamdegree;
		this->fParamHz = src->fParamHz;
		this->fParamIn = src->fParamIn;
		this->fParamtype = src->fParamtype;
		this->fParamFs = src->fParamFs;
		this->degree = src->degree;
		this->In = src->In;

		_6dB_1.copySVF(&(src->_6dB_1));
		_12dB_1.copySVF(&(src->_12dB_1));
		_12dB_2.copySVF(&(src->_12dB_2));
		_12dB_3.copySVF(&(src->_12dB_3));

		//setPassFilter();
	}

	void setPassFilter()
	{
		In = fParamIn ? 1 : 0;
		degree = _norm_to_degree(fParamdegree);

		_6dB_1.fParamIn = fParamIn;
		_12dB_1.fParamIn = fParamIn;
		_12dB_2.fParamIn = fParamIn;
		_12dB_3.fParamIn = fParamIn;

		double __Hz = 0;
		if (SVF::_norm_to_type(fParamtype) == SVF::kHP) __Hz = SVF::_Hz_to_norm(_norm_to_Hz_HP(fParamHz));
		else __Hz = SVF::_Hz_to_norm(_norm_to_Hz_LP(fParamHz));
		_6dB_1.fParamHz = __Hz;
		_12dB_1.fParamHz = __Hz;
		_12dB_2.fParamHz = __Hz;
		_12dB_3.fParamHz = __Hz;

		switch (degree)
		{
		case PassFilter::_6dBoct:
			break;
		case PassFilter::_12dBoct:
			_12dB_1.fParamQ = SVF::_Q_to_norm(M_SQRT2);
			break;
		case PassFilter::_18dBoct:
			_12dB_1.fParamQ = SVF::_Q_to_norm(2.0);
			break;
		case PassFilter::_24dBoct:
#define _24dBoct_1 1.08239220029239402443 // sqrt(2) * sqrt(2 - sqrt(2))
#define _24dBoct_2 2.61312592975275315155 // sqrt(2) * sqrt(2 + sqrt(2))
			_12dB_1.fParamQ = SVF::_Q_to_norm(_24dBoct_1);
			_12dB_2.fParamQ = SVF::_Q_to_norm(_24dBoct_2);
			break;
		case PassFilter::_36dBoct:
#define _36dBoct_1 1.03527618041008273586 // sqrt(6.0) - sqrt(2.0)
#define _36dBoct_2 M_SQRT2                // sqrt(2.0)
#define _36dBoct_3 3.86370330515627280477 // sqrt(6.0) + sqrt(2.0)
			_12dB_1.fParamQ = SVF::_Q_to_norm(_36dBoct_1);
			_12dB_2.fParamQ = SVF::_Q_to_norm(_36dBoct_2);
			_12dB_3.fParamQ = SVF::_Q_to_norm(_36dBoct_3);
			break;
		default:
			break;
		}

		_6dB_1.setSVF();
		_12dB_1.setSVF();
		_12dB_2.setSVF();
		_12dB_3.setSVF();
	}

	void makePassFilter()
	{
		_6dB_1.makeSVF();
		_12dB_1.makeSVF();
		_12dB_2.makeSVF();
		_12dB_3.makeSVF();
	}

	Steinberg::Vst::Sample64 computePassFilter
	(Steinberg::Vst::Sample64 input)
	{
		switch (degree)
		{
		case PassFilter::_6dBoct:
			return _6dB_1.computeSVF(input);
			break;
		case PassFilter::_12dBoct:
			return _12dB_1.computeSVF(input);
			break;
		case PassFilter::_18dBoct:
			return _12dB_1.computeSVF(_6dB_1.computeSVF(input));
			break;
		case PassFilter::_24dBoct:
			return _12dB_1.computeSVF(_12dB_2.computeSVF(input));
			break;
		case PassFilter::_36dBoct:
			return _12dB_1.computeSVF(_12dB_2.computeSVF(_12dB_3.computeSVF(input)));
			break;
		default:
			return input;
			break;
		}
	};

	double mag_response(double freq) {
		switch (degree)
		{
		case PassFilter::_6dBoct:
			return _6dB_1.mag_response(freq);
			break;
		case PassFilter::_12dBoct:
			return _12dB_1.mag_response(freq);
			break;
		case PassFilter::_18dBoct:
			return _12dB_1.mag_response(freq) * _6dB_1.mag_response(freq);
			break;
		case PassFilter::_24dBoct:
			return _12dB_1.mag_response(freq) * _12dB_2.mag_response(freq);
			break;
		case PassFilter::_36dBoct:
			return _12dB_1.mag_response(freq) * _12dB_2.mag_response(freq) * _12dB_3.mag_response(freq);
			break;
		default:
			return 1.0;
			break;
		}
	}

	static double Init_HP_Hz() { return 20.0; };
	static double Init_LP_Hz() { return 20000.0; };

	static double getFreqMax_HP() { return 400.0; };
	static double getFreqMin_HP() { return 20.0; };

	static double getFreqMax_LP() { return 41000.0; };
	static double getFreqMin_LP() { return 1000.0; };

	static Steinberg::Vst::ParamValue _Hz_to_norm_HP(double _Hz)
	{
		return log(_Hz / getFreqMin_HP()) / log(getFreqMax_HP() / getFreqMin_HP());
	};
	static double _norm_to_Hz_HP(Steinberg::Vst::ParamValue paramValue)
	{
		double FREQ_LOG_MAX = log(getFreqMax_HP() / getFreqMin_HP());
		double tmp = getFreqMin_HP() * exp(FREQ_LOG_MAX * paramValue);
		return std::max(std::min(tmp, getFreqMax_HP()), getFreqMin_HP());
	};

	static Steinberg::Vst::ParamValue _Hz_to_norm_LP(double _Hz)
	{
		return log(_Hz / getFreqMin_LP()) / log(getFreqMax_LP() / getFreqMin_LP());
	};
	static double _norm_to_Hz_LP(Steinberg::Vst::ParamValue paramValue)
	{
		double FREQ_LOG_MAX_LP = log(getFreqMax_LP() / getFreqMin_LP());
		double tmp_LP = getFreqMin_LP() * exp(FREQ_LOG_MAX_LP * paramValue);
		return tmp_LP;//std::max(std::min(tmp_LP, getFreqMax_LP()), getFreqMin_LP());
	};
	static filter_degree _norm_to_degree(Steinberg::Vst::ParamValue value)
	{
		return static_cast<filter_degree>(Steinberg::FromNormalized<Steinberg::Vst::ParamValue>(value, kDegreeNum));
	};
	static Steinberg::Vst::ParamValue _degree_to_norm(filter_degree type)
	{
		return Steinberg::ToNormalized<Steinberg::Vst::ParamValue>(type, kDegreeNum);
	};

	SVF _6dB_1;
	SVF _12dB_1, _12dB_2, _12dB_3;
	Steinberg::Vst::ParamValue fParamIn = 0.0;
	Steinberg::Vst::ParamValue fParamHz;
	Steinberg::Vst::ParamValue fParamtype;
	Steinberg::Vst::ParamValue fParamFs;
	Steinberg::Vst::ParamValue fParamdegree = _degree_to_norm(_12dBoct);
	filter_degree degree = _12dBoct;
	int In = 0;
};
*/