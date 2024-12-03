//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#include "RFEQ_processor.h"
#include "RFEQ_cids.h"

#include "base/source/fstreamer.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"

#include "public.sdk/source/vst/vstaudioprocessoralgo.h"
#include "public.sdk/source/vst/vsthelpers.h"

using namespace Steinberg;

namespace yg331 {
//------------------------------------------------------------------------
// RFEQ_Processor
//------------------------------------------------------------------------
RFEQ_Processor::RFEQ_Processor ()
{
    //--- set the wanted controller for our processor
    setControllerClass (kRFEQ_ControllerUID);
}

//------------------------------------------------------------------------
RFEQ_Processor::~RFEQ_Processor ()
{}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::initialize (FUnknown* context)
{
    // Here the Plug-in will be instantiated

    //---always initialize the parent-------
    tresult result = AudioEffect::initialize (context);
    // if everything Ok, continue
    if (result != kResultOk)
    {
        return result;
    }

    //--- create Audio IO ------
    addAudioInput  (STR16 ("Stereo In"),  Steinberg::Vst::SpeakerArr::kStereo);
    addAudioOutput (STR16 ("Stereo Out"), Steinberg::Vst::SpeakerArr::kStereo);

    /* If you don't need an event bus, you can remove the next line */
    //addEventInput (STR16 ("Event In"), 1);

    std::fill_n(&OS_buff[0][0], maxChannel * Kaiser::maxTap, 0.0);
    std::fill(OS_coef, OS_coef + maxChannel, 0.0);

    return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::terminate ()
{
    // Here the Plug-in will be de-instantiated, last possibility to remove some memory!

    //---do not forget to call parent ------
    return AudioEffect::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::setBusArrangements(
    Vst::SpeakerArrangement* inputs, int32 numIns,
    Vst::SpeakerArrangement* outputs, int32 numOuts)
{
    if (numIns == 1 && numOuts == 1)
    {
        // the host wants Mono => Mono (or 1 channel -> 1 channel)
        if (Vst::SpeakerArr::getChannelCount(inputs[0]) == 1 &&
            Vst::SpeakerArr::getChannelCount(outputs[0]) == 1)
        {
            auto* bus = FCast<Vst::AudioBus>(audioInputs.at(0));
            if (bus)
            {
                // check if we are Mono => Mono, if not we need to recreate the busses
                if (bus->getArrangement() != inputs[0])
                {
                    getAudioInput(0)->setArrangement(inputs[0]);
                    getAudioInput(0)->setName(STR16("Mono In"));
                    getAudioOutput(0)->setArrangement(outputs[0]);
                    getAudioOutput(0)->setName(STR16("Mono Out"));
                }
                return kResultOk;
            }
        }
        // the host wants something else than Mono => Mono,
        // in this case we are always Stereo => Stereo
        else
        {
            auto* bus = FCast<Vst::AudioBus>(audioInputs.at(0));
            if (bus)
            {
                tresult result = kResultFalse;

                // the host wants 2->2 (could be LsRs -> LsRs)
                if (Vst::SpeakerArr::getChannelCount(inputs[0]) == 2 &&
                    Vst::SpeakerArr::getChannelCount(outputs[0]) == 2)
                {
                    getAudioInput(0)->setArrangement(inputs[0]);
                    getAudioInput(0)->setName(STR16("Stereo In"));
                    getAudioOutput(0)->setArrangement(outputs[0]);
                    getAudioOutput(0)->setName(STR16("Stereo Out"));
                    result = kResultTrue;
                }
                // the host want something different than 1->1 or 2->2 : in this case we want stereo
                else if (bus->getArrangement() != Vst::SpeakerArr::kStereo)
                {
                    getAudioInput(0)->setArrangement(Vst::SpeakerArr::kStereo);
                    getAudioInput(0)->setName(STR16("Stereo In"));
                    getAudioOutput(0)->setArrangement(Vst::SpeakerArr::kStereo);
                    getAudioOutput(0)->setName(STR16("Stereo Out"));
                    result = kResultFalse;
                }

                return result;
            }
        }
    }
    return kResultFalse;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::setActive (TBool state)
{
    //--- called when the Plug-in is enable/disable (On/Off) -----
    return AudioEffect::setActive (state);
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::process (Vst::ProcessData& data)
{
    Vst::IParameterChanges* paramChanges = data.inputParameterChanges;

    if (paramChanges)
    {
        int32 numParamsChanged = paramChanges->getParameterCount();

        for (int32 index = 0; index < numParamsChanged; index++)
        {
            Vst::IParamValueQueue* paramQueue = paramChanges->getParameterData(index);

            if (paramQueue)
            {
                Vst::ParamValue value;
                int32 sampleOffset;
                int32 numPoints = paramQueue->getPointCount();

                /*/*/
                if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue)
                {
                    switch (paramQueue->getParameterId()) {
                        case kParamBypass: bBypass = (value > 0.5f); break;
                        // case kParamZoom:   fZoom = value; break;
                        case kParamLevel:  fLevel = value; break;
                        // case kParamOutput: fOutput = value; break;
                            
                        case kParamBand1_In: band1[bandIn] = value; break;
                        case kParamBand2_In: band2[bandIn] = value; break;
                        case kParamBand3_In: band3[bandIn] = value; break;
                        case kParamBand4_In: band4[bandIn] = value; break;
                        case kParamBand5_In: band5[bandIn] = value; break;

                        case kParamBand1_Hz: band1[bandHz] = paramFreq.ToPlain(value); break;
                        case kParamBand2_Hz: band2[bandHz] = paramFreq.ToPlain(value); break;
                        case kParamBand3_Hz: band3[bandHz] = paramFreq.ToPlain(value); break;
                        case kParamBand4_Hz: band4[bandHz] = paramFreq.ToPlain(value); break;
                        case kParamBand5_Hz: band5[bandHz] = paramFreq.ToPlain(value); break;

                        case kParamBand1_Q: band1[bandQ] = paramQlty.ToPlain(value); break;
                        case kParamBand2_Q: band2[bandQ] = paramQlty.ToPlain(value); break;
                        case kParamBand3_Q: band3[bandQ] = paramQlty.ToPlain(value); break;
                        case kParamBand4_Q: band4[bandQ] = paramQlty.ToPlain(value); break;
                        case kParamBand5_Q: band5[bandQ] = paramQlty.ToPlain(value); break;

                        case kParamBand1_dB: band1[banddB] = paramGain.ToPlain(value); break;
                        case kParamBand2_dB: band2[banddB] = paramGain.ToPlain(value); break;
                        case kParamBand3_dB: band3[banddB] = paramGain.ToPlain(value); break;
                        case kParamBand4_dB: band4[banddB] = paramGain.ToPlain(value); break;
                        case kParamBand5_dB: band5[banddB] = paramGain.ToPlain(value); break;

                        case kParamBand1_Type: band1[bandType] = paramType.ToPlain(value); break;
                        case kParamBand2_Type: band2[bandType] = paramType.ToPlain(value); break;
                        case kParamBand3_Type: band3[bandType] = paramType.ToPlain(value); break;
                        case kParamBand4_Type: band4[bandType] = paramType.ToPlain(value); break;
                        case kParamBand5_Type: band5[bandType] = paramType.ToPlain(value); break;

                        case kParamBand1_Order: band1[bandOrder] = paramOrdr.ToPlain(value); break;
                        case kParamBand2_Order: band2[bandOrder] = paramOrdr.ToPlain(value); break;
                        case kParamBand3_Order: band3[bandOrder] = paramOrdr.ToPlain(value); break;
                        case kParamBand4_Order: band4[bandOrder] = paramOrdr.ToPlain(value); break;
                        case kParamBand5_Order: band5[bandOrder] = paramOrdr.ToPlain(value); break;
                    }
                }
            }
        }
        call_after_parameter_changed ();
    }

    if (data.numSamples <= 0)
        return kResultOk; // nothing to do
    
    if (data.numInputs == 0 || data.numOutputs == 0)
        return kResultOk; // nothing to do

    // (simplification) we suppose in this example that we have the same input channel count than
    // the output
    int32 numChannels = data.inputs[0].numChannels;

    //---get audio buffers----------------
    uint32 sampleFramesSize = getSampleFramesSizeInBytes(processSetup, data.numSamples);
    void** in  = getChannelBuffersPointer(processSetup, data.inputs[0]);
    void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);
    Vst::SampleRate getSampleRate = processSetup.sampleRate;

    //---check if silence---------------
    // check if all channel are silent then process silent
    if (data.inputs[0].silenceFlags == Vst::getChannelMask(data.inputs[0].numChannels))
    {
        // mark output silence too (it will help the host to propagate the silence)
        data.outputs[0].silenceFlags = data.inputs[0].silenceFlags;

        // the plug-in has to be sure that if it sets the flags silence that the output buffer are
        // clear
        for (int32 i = 0; i < numChannels; i++)
        {
            // do not need to be cleared if the buffers are the same (in this case input buffer are
            // already cleared by the host)
            if (in[i] != out[i])
            {
                memset(out[i], 0, sampleFramesSize);
            }
        }
    }
    else {

        data.outputs[0].silenceFlags = data.inputs[0].silenceFlags;
    
        if (data.symbolicSampleSize == Vst::kSample32) {
            processSVF<Vst::Sample32>((Vst::Sample32**)in, (Vst::Sample32**)out, numChannels, getSampleRate, data.numSamples);
        }
        else {
            processSVF<Vst::Sample64>((Vst::Sample64**)in, (Vst::Sample64**)out, numChannels, getSampleRate, data.numSamples);
        }
    }
    
    {
        if (IPtr<Vst::IMessage> message = owned (allocateMessage ()))
        {
            message->setMessageID ("GUI");
            message->getAttributes()->setFloat ("projectSR", projectSR);
            sendMessage (message);
        }
    }
    {
        if (IPtr<Vst::IMessage> message = owned (allocateMessage ()))
        {
            message->setMessageID ("GUI");
            message->getAttributes()->setFloat ("targetSR", targetSR);
            sendMessage (message);
        }
    }
    {
        auto output = data.outputs[0];
        if (IPtr<Vst::IMessage> message = owned (allocateMessage ()))
        {
            message->setMessageID ("GUI");
            std::vector<float> c;
            if (numChannels == 2) { // Stereo
                std::transform (output.channelBuffers32[0],
                                output.channelBuffers32[0] + data.numSamples,
                                output.channelBuffers32[1],
                                std::back_inserter(c),
                                []( const auto &x, const auto &y ) { return std::max( x, y ); });
            }
            else {
                std::copy(output.channelBuffers32[0], output.channelBuffers32[0] + data.numSamples, std::back_inserter(c));
            }
            c.resize(data.numSamples, 0.0);
            if (auto attributes = message->getAttributes())
                attributes->setBinary ("sample", c.data(), sampleFramesSize);
            sendMessage (message);
        }
    }

    return kResultOk;
}

//------------------------------------------------------------------------
uint32 PLUGIN_API RFEQ_Processor::getLatencySamples()
{
    // fprintf (stdout, "getLatencySamples\n");

    return currLatency;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::setupProcessing (Vst::ProcessSetup& newSetup)
{
    //--- called before any processing ----
    
    /*
     createInstance
     initialize
     Processor::getState
     setupProcessing
     setState
     getLatencySamples
     */
    
    // fprintf (stdout, "setupProcessing\n");
    
    projectSR = newSetup.sampleRate;
    
    //fft_in.resize(newSetup.maxSamplesPerBlock+1, 0.0);
    
    call_after_SR_changed (); // includes call_after_parameter_changed ()

    return AudioEffect::setupProcessing (newSetup);
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::canProcessSampleSize (int32 symbolicSampleSize)
{
    // by default kSample32 is supported
    if (symbolicSampleSize == Vst::kSample32)
        return kResultTrue;

    // disable the following comment if your processing support kSample64
    if (symbolicSampleSize == Vst::kSample64)
        return kResultTrue;

    return kResultFalse;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::setState (IBStream* state)
{
    // fprintf (stdout, "Processor::setState\n");
    
    // called when we load a preset, the model has to be reloaded
    IBStreamer streamer(state, kLittleEndian);
    
    // I started saving in Normalized, so have to keep saving in Normalized

    int32           savedBypass = 0;
    Vst::ParamValue savedZoom   = 0.0; // UNUSED, left for compatibility
    Vst::ParamValue savedLevel  = 0.0;
    Vst::ParamValue savedOutput = 0.0; // UNUSED, left for compatibility

    ParamBand_Array savedBand1_Array = {1.0, dftBand1Freq, dftParamQlty, dftParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array savedBand2_Array = {1.0, dftBand2Freq, dftParamQlty, dftParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array savedBand3_Array = {1.0, dftBand3Freq, dftParamQlty, dftParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array savedBand4_Array = {1.0, dftBand4Freq, dftParamQlty, dftParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array savedBand5_Array = {1.0, dftBand5Freq, dftParamQlty, dftParamGain, nrmParamType, nrmParamOrdr};

    if (streamer.readInt32 (savedBypass) == false) return kResultFalse;
    if (streamer.readDouble(savedZoom  ) == false) return kResultFalse;
    if (streamer.readDouble(savedLevel ) == false) return kResultFalse;
    if (streamer.readDouble(savedOutput) == false) return kResultFalse;

    if (streamer.readDoubleArray(savedBand1_Array, bandSize) == false) return kResultFalse;
    if (streamer.readDoubleArray(savedBand2_Array, bandSize) == false) return kResultFalse;
    if (streamer.readDoubleArray(savedBand3_Array, bandSize) == false) return kResultFalse;
    if (streamer.readDoubleArray(savedBand4_Array, bandSize) == false) return kResultFalse;
    if (streamer.readDoubleArray(savedBand5_Array, bandSize) == false) return kResultFalse;

    bBypass = savedBypass > 0;
    // fZoom   = savedZoom;
    fLevel  = savedLevel;
    // fOutput = savedOutput;
    
    auto copyNormToPlain = [this](double Plain[], double Norm[])
    {
        Plain[bandIn]    = Norm[bandIn];
        Plain[bandHz]    = paramFreq.ToPlain(Norm[bandHz]);
        Plain[bandQ]     = paramQlty.ToPlain(Norm[bandQ]);
        Plain[banddB]    = paramGain.ToPlain(Norm[banddB]);
        Plain[bandType]  = paramType.ToPlain(Norm[bandType]);
        Plain[bandOrder] = paramOrdr.ToPlain(Norm[bandOrder]);
    };
    
    copyNormToPlain(band1, savedBand1_Array);
    copyNormToPlain(band2, savedBand2_Array);
    copyNormToPlain(band3, savedBand3_Array);
    copyNormToPlain(band4, savedBand4_Array);
    copyNormToPlain(band5, savedBand5_Array);
    
    call_after_parameter_changed ();

    return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::getState (IBStream* state)
{
    // fprintf (stdout, "Processor::getState\n");
    
    // here we need to save the model
    IBStreamer streamer(state, kLittleEndian);
    
    // I started saving in Normalized, so have to keep saving in Normalized
    
    streamer.writeInt32(bBypass ? 1 : 0);
    streamer.writeDouble(fZoom);   // UNUSED, left for compatibility
    streamer.writeDouble(fLevel);
    streamer.writeDouble(fOutput); // UNUSED, left for compatibility
    
    auto writeDoubleArrayNorm = [&streamer](double Array[])
    {
        ParamBand_Array normArray = {
            Array[bandIn],
            paramFreq.ToNormalized(Array[bandHz]),
            paramQlty.ToNormalized(Array[bandQ]),
            paramGain.ToNormalized(Array[banddB]),
            paramType.ToNormalized(Array[bandType]),
            paramOrdr.ToNormalized(Array[bandOrder])
        };
        streamer.writeDoubleArray(normArray, bandSize);
    };
    
    writeDoubleArrayNorm(band1);
    writeDoubleArrayNorm(band2);
    writeDoubleArrayNorm(band3);
    writeDoubleArrayNorm(band4);
    writeDoubleArrayNorm(band5);
    
    return kResultOk;
}

template <typename SampleType>
void RFEQ_Processor::processSVF
 (
    SampleType** inputs,
    SampleType** outputs,
    Steinberg::int32 numChannels,
    Steinberg::Vst::SampleRate getSampleRate,
    Steinberg::int32 sampleFrames
  )
{
    Vst::Sample64 level = DecibelConverter::ToGain(paramGain.ToPlain(fLevel));
    double div_by_channels = 1.0 / numChannels;

    int32 oversampling = (fParamOS == overSample_2x) ? 2 : 1;

    for (int32 channel = 0; channel < numChannels; channel++)
    {
        int32 samples = 0;
        while (samples < sampleFrames)
        {
            Vst::Sample64 inputSample = inputs[channel][samples];
            Vst::Sample64 drySample = inputSample;
            inputSample *= level;

            double up_x[4] = { 0.0, };
            double up_y[4] = { 0.0, };

            up_x[0] = inputSample;

            // Process
            for (int i = 0; i < oversampling; i++)
            {
                Vst::Sample64 overSampled = up_x[i];

                double v1 = band1_svf[channel].computeSVF(overSampled);
                double v2 = band2_svf[channel].computeSVF(v1);
                double v3 = band3_svf[channel].computeSVF(v2);
                double v4 = band4_svf[channel].computeSVF(v3);
                double v5 = band5_svf[channel].computeSVF(v4);

                up_y[i] = v5; // * gain;
            }

            // Downsampling
            if (fParamOS == overSample_1x) inputSample = up_y[0];
            else if (fParamOS == overSample_2x)
            {
                memmove(OS_buff[channel] + 2, OS_buff[channel], sizeofOsMove);
                OS_buff[channel][1] = up_y[0];
                OS_buff[channel][0] = up_y[1];
                // transform_reduce works faster in double[], and slow in std::deque<double>
                // but if loop order channel->sample, cache miss happens, and std::deque<double> works faster
                // Well, it just depends case-by-case.
                inputSample = std::transform_reduce(OS_coef, OS_coef + fir_size, OS_buff[channel] + 1, 0.0);
            }

            // Latency compensate
            latencyDelayLine[channel].push_back(drySample);
            Vst::Sample64 delayed = *(latencyDelayLine[channel].end() - 1 - currLatency);
            latencyDelayLine[channel].pop_front();

            if (bBypass)
                inputSample = delayed;

            outputs[channel][samples] = (SampleType)inputSample;

            samples++;
        }
    }
    return;
}

void RFEQ_Processor::call_after_SR_changed ()
{
    if (projectSR <= 48000.0)
    {
        targetSR = 2 * projectSR;
        fParamOS = overSample_2x;
        currLatency = latency_Fir_x2;
    }
    else
    {
        targetSR = projectSR;
        fParamOS = overSample_1x;
        currLatency = 0;
    }
    
    Kaiser::calcFilter(96000.0, 0.0, 24000.0, fir_size, 110.0, OS_coef); // half band filter
    std::for_each(OS_coef, OS_coef + fir_size, [](double &n) { n *= 2.0; });
    std::fill_n(&OS_buff[0][0], maxChannel * Kaiser::maxTap, 0.0);
    
    for (auto& iter : latencyDelayLine) { iter.resize(latency_Fir_x2, 0.0); std::fill(iter.begin(), iter.end(), 0.0); }
    
    //FFT.reset();
    
    call_after_parameter_changed ();
};
void RFEQ_Processor::call_after_parameter_changed ()
{
    for(int ch = 0; ch < maxChannel; ch++)
    {
        band1_svf[ch].setSVF(band1[bandIn], band1[bandHz], band1[bandQ], band1[banddB], band1[bandType], band1[bandOrder], targetSR);
        band2_svf[ch].setSVF(band2[bandIn], band2[bandHz], band2[bandQ], band2[banddB], band2[bandType], band2[bandOrder], targetSR);
        band3_svf[ch].setSVF(band3[bandIn], band3[bandHz], band3[bandQ], band3[banddB], band3[bandType], band3[bandOrder], targetSR);
        band4_svf[ch].setSVF(band4[bandIn], band4[bandHz], band4[bandQ], band4[banddB], band4[bandType], band4[bandOrder], targetSR);
        band5_svf[ch].setSVF(band5[bandIn], band5[bandHz], band5[bandQ], band5[banddB], band5[bandType], band5[bandOrder], targetSR);
    }
};

//------------------------------------------------------------------------
} // namespace yg331
