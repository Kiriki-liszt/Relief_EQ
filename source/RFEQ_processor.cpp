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

    return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::terminate ()
{
    // Here the Plug-in will be de-instantiated, last possibility to remove some memory!
    fft_in.clear();
    fft_in.shrink_to_fit();
    // fft_out.clear();
    // fft_out.shrink_to_fit();
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
tresult PLUGIN_API RFEQ_Processor::connect(Vst::IConnectionPoint* other)
{
    auto result = Vst::AudioEffect::connect(other);
    if (result == kResultTrue)
    {
        auto configCallback = [this] (Vst::DataExchangeHandler::Config& config,
                                      const Vst::ProcessSetup& setup) {
                // Vst::SpeakerArrangement arr;
                // getBusArrangement(Vst::BusDirections::kInput, 0, arr);
                // numChannels = static_cast<uint16_t> (Vst::SpeakerArr::getChannelCount(arr));
                auto sampleSize = sizeof(float);

                config.blockSize = fftSize * sampleSize + sizeof(DataBlock);
                config.numBlocks = 1;
                // config.alignment = 32;
                // config.userContextID = 0;
                return true;
            };

        dataExchange = std::make_unique<Vst::DataExchangeHandler>(this, configCallback);
        dataExchange->onConnect(other, getHostContext());
    }
    return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::disconnect(Vst::IConnectionPoint* other)
{
    if (dataExchange)
    {
        dataExchange->onDisconnect(other);
        dataExchange.reset();
    }
    return AudioEffect::disconnect(other);
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::setActive (TBool state)
{
    //--- called when the Plug-in is enable/disable (On/Off) -----
    if (state) {
        if (dataExchange)
            dataExchange->onActivate(processSetup);
    }
    else {
        if (dataExchange)
            dataExchange->onDeactivate();
    }

    return AudioEffect::setActive (state);
}

//------------------------------------------------------------------------
void RFEQ_Processor::acquireNewExchangeBlock()
{
    currentExchangeBlock = dataExchange->getCurrentOrNewBlock();
    if (auto block = toDataBlock(currentExchangeBlock))
    {
        block->filterSampleRate = static_cast<uint32_t> (processSetup.sampleRate);
        block->FFTSampleRate    = static_cast<uint32_t> (processSetup.sampleRate);
        block->FFTDataAvail     = 0;
    }
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

    if (data.numInputs == 0 || data.numOutputs == 0)
    {
        // nothing to do
        return kResultOk;
    }

    // (simplification) we suppose in this example that we have the same input channel count than
    // the output
    // int32 numChannels = data.inputs[0].numChannels;
    numChannels = data.inputs[0].numChannels;

    //---get audio buffers----------------
    uint32 sampleFramesSize = getSampleFramesSizeInBytes(processSetup, data.numSamples);
    void** in  = getChannelBuffersPointer(processSetup, data.inputs[0]);
    void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);
    Vst::SampleRate getSampleRate = processSetup.sampleRate;

    std::fill(fft_in.begin(), fft_in.end(), 0.0);
    // std::fill(fft_out.begin(), fft_out.end(), 0.0);

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

    // int data_avail = FFT.getData(fft_out.data());
    
    if (int data_avail = FFT.getData(fft_out); data_avail)
    {
        //--- send data ----------------
        if (currentExchangeBlock.blockID == Vst::InvalidDataExchangeBlockID)
            acquireNewExchangeBlock();
        if (auto block = toDataBlock(currentExchangeBlock))
        {
            memcpy(&block->samples[0], fft_out, numBins * sizeof(float));
            // memcpy(&block->samples[0], fft_out.data(), numBins * sizeof(float));
            // std::copy(fft_out.begin(), fft_out.end(), &block->samples[0]);
            block->FFTSampleRate = projectSR;
            block->FFTDataAvail = data_avail;
            block->numSamples = data.numSamples;
            block->filterSampleRate = targetSR;
            dataExchange->sendCurrentBlock();
            acquireNewExchangeBlock();
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

    call_after_SR_changed (); // includes call_after_parameter_changed ()

    fft_in.resize(newSetup.maxSamplesPerBlock+1, 0.0);
    // fft_out.resize(numBins, 0.0);

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
    
    auto copyNormToPlain = [this](double Norm[], double Plain[])
    {
        Plain[bandIn] = Norm[bandIn];
        Plain[bandHz] = paramFreq.ToPlain(Norm[bandHz]);
        Plain[bandQ]  = paramQlty.ToPlain(Norm[bandQ]);
        Plain[banddB] = paramGain.ToPlain(Norm[banddB]);
        Plain[bandType] = paramType.ToPlain(Norm[bandType]);
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

    Vst::Sample64 _db = (24.0 * fLevel - 12.0);
    Vst::Sample64 level = exp(log(10.0) * _db / 20.0);
    double div_by_channels = 1.0 / numChannels;

    int32 oversampling = 1;
    if (fParamOS == overSample_2x) oversampling = 2;

    int32 latency = 0;
    if (fParamOS == overSample_2x) latency = latency_Fir_x2;

    for (int32 channel = 0; channel < numChannels; channel++)
    {
        SampleType* ptrIn  = (SampleType*)inputs[channel];
        SampleType* ptrOut = (SampleType*)outputs[channel];

        float* fft_in_begin = fft_in.data();
        
        int32 samples = sampleFrames;
        while (--samples >= 0)
        {
            Vst::Sample64 inputSample = *ptrIn; ptrIn++;
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
                // transform works faster in double[], and slow in std::deque<double>
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

            *fft_in_begin += div_by_channels * (inputSample);
            fft_in_begin++;

            *ptrOut = (SampleType)inputSample;
            ptrOut++;
        }
    }

    FFT.processBlock(fft_in.data(), sampleFrames, 0);

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
    for (auto& iter : latencyDelayLine) iter.resize(latency_Fir_x2, 0.0);
    
    FFT.reset();
    
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
