//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once
#include "RFEQ_shared.h"
#include "RFEQ_fft.h"

#include "public.sdk/source/vst/vstaudioeffect.h"

#include <deque>     // std::deque
#include <numeric>   // std::transform_reduce
#include <algorithm> // std::fill, std::for_each

namespace yg331 {
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
    /** Called when a message has been sent from the connection point to this. */
    // Steinberg::tresult PLUGIN_API notify(Steinberg::Vst::IMessage* message) SMTG_OVERRIDE;

//------------------------------------------------------------------------
protected:
    static SMTG_CONSTEXPR int32 maxChannel = 2;
    static SMTG_CONSTEXPR int32 fir_size = 69;
    static SMTG_CONSTEXPR int32 tap_hm = (fir_size - 1) / 2;
    static SMTG_CONSTEXPR int32 latency_Fir_x2 = (fir_size - 1) / 4;
    static SMTG_CONSTEXPR size_t sizeofDouble = sizeof(double);
    static SMTG_CONSTEXPR size_t sizeofOsMove = sizeof(double) * (fir_size - 2);
    
    template <typename SampleType>
    void processSVF ( SampleType** inputs, SampleType** outputs, int32 numChannels, SampleRate getSampleRate, int32 sampleFrames );
    
    void call_after_SR_changed ();
    void call_after_parameter_changed ();

    bool       bBypass = false;
    ParamValue fLevel  = 0.5;
    ParamValue fOutput = 0.5; // UNUSED
    ParamValue fZoom   = 2.0 / 6.0; // UNUSED
    int32      fParamOS = overSample_2x;
    
    ParamBand_Array band1 = {1.0, dftBand1Freq, dftParamQlty, dftParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array band2 = {1.0, dftBand2Freq, dftParamQlty, dftParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array band3 = {1.0, dftBand3Freq, dftParamQlty, dftParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array band4 = {1.0, dftBand4Freq, dftParamQlty, dftParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array band5 = {1.0, dftBand5Freq, dftParamQlty, dftParamGain, nrmParamType, nrmParamOrdr};

    SVF band1_svf[maxChannel];
    SVF band2_svf[maxChannel];
    SVF band3_svf[maxChannel];
    SVF band4_svf[maxChannel];
    SVF band5_svf[maxChannel];

    // plugin enviroment
    SampleRate projectSR = 48000.0;
    SampleRate targetSR  = 96000.0;
    int32 currLatency = latency_Fir_x2;
    
    // Oversampling and Latency
    std::deque<ParamValue> latencyDelayLine[maxChannel];
    double OS_coef alignas(16)[Kaiser::maxTap] = {};
    double OS_buff alignas(16)[maxChannel][Kaiser::maxTap] = {};
};

//------------------------------------------------------------------------
} // namespace yg331
