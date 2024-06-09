//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once
#include "RFEQ_svf.h"
#include "RFEQ_fft.h"
#include "RFEQ_dataexchange.h"
#include "public.sdk/source/vst/vstaudioeffect.h"

#include <math.h>
#include <queue>

namespace yg331 {

class Kaiser {
public:

#define maxTap 512

    static inline double Ino(double x)
    {
        double d = 0, ds = 1, s = 1;
        do
        {
            d += 2;
            ds *= x * x / (d * d);
            s += ds;
        } while (ds > s * 1e-6);
        return s;
    };

    static void calcFilter(double Fs, double Fa, double Fb, int M, double Att, double* dest)
    {
        // Kaiser windowed FIR filter "DIGITAL SIGNAL PROCESSING, II" IEEE Press pp 123-126.

        int Np = (M - 1) / 2;
        double A[maxTap] = { 0, };
        double Alpha;
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
        for (int j = 0; j < Np; j++)
        {
            dest[j] = dest[M - 1 - j];
        }

    };
};

// Buffers ------------------------------------------------------------------
typedef struct _Flt {
	double coef alignas(16)[maxTap] = { 0, };
	double buff alignas(16)[maxTap] = { 0, };
} Flt;

//------------------------------------------------------------------------
static constexpr Steinberg::Vst::DataExchangeBlock InvalidDataExchangeBlock = {
    nullptr, 0, Steinberg::Vst::InvalidDataExchangeBlockID
};

//------------------------------------------------------------------------
//  RFEQ_Processor
//------------------------------------------------------------------------
class RFEQ_Processor : public Steinberg::Vst::AudioEffect
{
public:
	RFEQ_Processor ();
	~RFEQ_Processor () SMTG_OVERRIDE;

    // Create function
	static Steinberg::FUnknown* createInstance (void* /*context*/) 
	{ 
		return (Steinberg::Vst::IAudioProcessor*)new RFEQ_Processor; 
	}

	//------------------------------------------------------------------------
	// AudioEffect overrides:
	//------------------------------------------------------------------------
	/** Called at first after constructor */
	Steinberg::tresult PLUGIN_API initialize (Steinberg::FUnknown* context) SMTG_OVERRIDE;
	
	/** Called at the end before destructor */
	Steinberg::tresult PLUGIN_API terminate () SMTG_OVERRIDE;

	/** Switch the Plug-in on/off */
	Steinberg::tresult PLUGIN_API setActive (Steinberg::TBool state) SMTG_OVERRIDE;

	/** Will be called before any process call */
	Steinberg::tresult PLUGIN_API setupProcessing (Steinberg::Vst::ProcessSetup& newSetup) SMTG_OVERRIDE;
	
	/** Asks if a given sample size is supported see SymbolicSampleSizes. */
	Steinberg::tresult PLUGIN_API canProcessSampleSize (Steinberg::int32 symbolicSampleSize) SMTG_OVERRIDE;

	/** Gets the current Latency in samples. */
	Steinberg::uint32 PLUGIN_API getLatencySamples() SMTG_OVERRIDE;

	/** Here we go...the process call */
	Steinberg::tresult PLUGIN_API process (Steinberg::Vst::ProcessData& data) SMTG_OVERRIDE;
		
	/** For persistence */
	Steinberg::tresult PLUGIN_API setState (Steinberg::IBStream* state) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API getState (Steinberg::IBStream* state) SMTG_OVERRIDE;

	//---from IAudioProcessor-------
	Steinberg::tresult PLUGIN_API setBusArrangements(
		Steinberg::Vst::SpeakerArrangement* inputs,  Steinberg::int32 numIns,
		Steinberg::Vst::SpeakerArrangement* outputs, Steinberg::int32 numOuts
	) SMTG_OVERRIDE;

	//------------------------------------------------------------------------
	// IConnectionPoint overrides:
	//------------------------------------------------------------------------
	/** Connects this instance with another connection point. */
	Steinberg::tresult PLUGIN_API connect(Steinberg::Vst::IConnectionPoint* other) SMTG_OVERRIDE;
	/** Disconnects a given connection point from this. */
	Steinberg::tresult PLUGIN_API disconnect(Steinberg::Vst::IConnectionPoint* other) SMTG_OVERRIDE;
	/** Called when a message has been sent from the connection point to this. */
	//Steinberg::tresult PLUGIN_API notify(Steinberg::Vst::IMessage* message) SMTG_OVERRIDE;

//------------------------------------------------------------------------
protected:

	template <typename SampleType>
	void processSVF
	(
		SampleType** inputs,
		SampleType** outputs,
		Steinberg::int32 numChannels,
		Steinberg::Vst::SampleRate getSampleRate,
		Steinberg::int32 sampleFrames
	);

	bool                       bBypass = false;
	Steinberg::Vst::ParamValue fLevel  = 0.5;
	Steinberg::Vst::ParamValue fOutput = 0.0;
	Steinberg::Vst::ParamValue fZoom   = 2.0 / 6.0;
	overSample                 fParamOS = overSample_2x;

	ParamBand_Array fParamBand1_Array = { 1.0, SVF::_Hz_to_norm(80.0),    SVF::_Q_to_norm(1.0), SVF::_dB_to_norm(0.0), SVF::_Type_to_norm(SVF::kBell), SVF::_Order_to_norm(SVF::_12dBoct) };
	ParamBand_Array fParamBand2_Array = { 1.0, SVF::_Hz_to_norm(200.0),   SVF::_Q_to_norm(1.0), SVF::_dB_to_norm(0.0), SVF::_Type_to_norm(SVF::kBell), SVF::_Order_to_norm(SVF::_12dBoct) };
	ParamBand_Array fParamBand3_Array = { 1.0, SVF::_Hz_to_norm(2000.0),  SVF::_Q_to_norm(1.0), SVF::_dB_to_norm(0.0), SVF::_Type_to_norm(SVF::kBell), SVF::_Order_to_norm(SVF::_12dBoct) };
	ParamBand_Array fParamBand4_Array = { 1.0, SVF::_Hz_to_norm(6000.0),  SVF::_Q_to_norm(1.0), SVF::_dB_to_norm(0.0), SVF::_Type_to_norm(SVF::kBell), SVF::_Order_to_norm(SVF::_12dBoct) };
	ParamBand_Array fParamBand5_Array = { 1.0, SVF::_Hz_to_norm(16000.0), SVF::_Q_to_norm(1.0), SVF::_dB_to_norm(0.0), SVF::_Type_to_norm(SVF::kBell), SVF::_Order_to_norm(SVF::_12dBoct) };

	SVF Band1[2];
	SVF Band2[2];
	SVF Band3[2];
	SVF Band4[2];
	SVF Band5[2];

	// Oversampling and Latency

	Flt OS_filter_x2[2];
	const int fir_size = 69;
	const int tap_hm = (fir_size - 1) / 2;

	void Fir_x2_dn(Steinberg::Vst::Sample64* in, Steinberg::Vst::Sample64* out, Steinberg::int32 channel);

	const Steinberg::int32 latency_Fir_x2 = 17;

	std::queue<double> latency_q[2];

	// DataExchange

	void acquireNewExchangeBlock();

	std::unique_ptr<Steinberg::Vst::DataExchangeHandler> dataExchange;
	Steinberg::Vst::DataExchangeBlock currentExchangeBlock{ InvalidDataExchangeBlock };
	uint16_t numChannels{ 0 };

	// FFT

	FFTProcessor FFT;
    alignas(16) std::vector<float> fft_in = {0.0, }, fft_out = { 0.0, };
};

//------------------------------------------------------------------------
} // namespace yg331
