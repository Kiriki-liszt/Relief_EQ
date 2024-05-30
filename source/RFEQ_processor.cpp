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

	for (int ch = 0; ch < 2; ch++) {
		Kaiser::calcFilter(96000.0, 0.0, 24000.0, fir_size, 110.0, OS_filter_x2[ch].coef);
		for (int i = 0; i < fir_size; i++) OS_filter_x2[ch].coef[i] *= 2.0;
	}

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
tresult PLUGIN_API RFEQ_Processor::connect(Vst::IConnectionPoint* other)
{
	auto result = Vst::AudioEffect::connect(other);
	if (result == kResultTrue)
	{
		auto configCallback = [this](Vst::DataExchangeHandler::Config& config,
			const Vst::ProcessSetup& setup) {
				Vst::SpeakerArrangement arr;
				getBusArrangement(Vst::BusDirections::kInput, 0, arr);
				numChannels = static_cast<uint16_t> (Vst::SpeakerArr::getChannelCount(arr));
				//auto sampleSize = sizeof(float);

				config.blockSize = sizeof(DataBlock);
				config.numBlocks = 2;
				//config.alignment = 32;
				//config.userContextID = 0;
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
	if (state)
		dataExchange->onActivate(processSetup);
	else
		dataExchange->onDeactivate();

	return AudioEffect::setActive (state);
}

//------------------------------------------------------------------------
void RFEQ_Processor::acquireNewExchangeBlock()
{
	currentExchangeBlock = dataExchange->getCurrentOrNewBlock();
	if (auto block = toDataBlock(currentExchangeBlock))
	{
		ParamBand_Array Clean = { 0.0, };

		memcpy(block->Band1, Clean, ParamArray::ParamArray_size * sizeof(double));
		memcpy(block->Band2, Clean, ParamArray::ParamArray_size * sizeof(double));
		memcpy(block->Band3, Clean, ParamArray::ParamArray_size * sizeof(double));
		memcpy(block->Band4, Clean, ParamArray::ParamArray_size * sizeof(double));
		memcpy(block->Band5, Clean, ParamArray::ParamArray_size * sizeof(double));

		block->Fs = 0.0;
		block->byPass = 0;
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
				if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue) {
					switch (paramQueue->getParameterId()) {
					case kParamBypass: bBypass = (value > 0.5f); break;
					case kParamZoom:   fZoom = value; break;
					case kParamLevel:  fLevel = value; break;
					case kParamOutput: fOutput = value; break;
						
					case kParamBand1_In: fParamBand1_Array[ParamArray_In] = value; break;
					case kParamBand2_In: fParamBand2_Array[ParamArray_In] = value; break;
					case kParamBand3_In: fParamBand3_Array[ParamArray_In] = value; break;
					case kParamBand4_In: fParamBand4_Array[ParamArray_In] = value; break;
					case kParamBand5_In: fParamBand5_Array[ParamArray_In] = value; break;

					case kParamBand1_Hz: fParamBand1_Array[ParamArray_Hz] = value; break;
					case kParamBand2_Hz: fParamBand2_Array[ParamArray_Hz] = value; break;
					case kParamBand3_Hz: fParamBand3_Array[ParamArray_Hz] = value; break;
					case kParamBand4_Hz: fParamBand4_Array[ParamArray_Hz] = value; break;
					case kParamBand5_Hz: fParamBand5_Array[ParamArray_Hz] = value; break;

					case kParamBand1_Q: fParamBand1_Array[ParamArray_Q] = value; break;
					case kParamBand2_Q: fParamBand2_Array[ParamArray_Q] = value; break;
					case kParamBand3_Q: fParamBand3_Array[ParamArray_Q] = value; break;
					case kParamBand4_Q: fParamBand4_Array[ParamArray_Q] = value; break;
					case kParamBand5_Q: fParamBand5_Array[ParamArray_Q] = value; break;

					case kParamBand1_dB: fParamBand1_Array[ParamArray_dB] = value; break;
					case kParamBand2_dB: fParamBand2_Array[ParamArray_dB] = value; break;
					case kParamBand3_dB: fParamBand3_Array[ParamArray_dB] = value; break;
					case kParamBand4_dB: fParamBand4_Array[ParamArray_dB] = value; break;
					case kParamBand5_dB: fParamBand5_Array[ParamArray_dB] = value; break;

					case kParamBand1_Type: fParamBand1_Array[ParamArray_Type] = value; break;
					case kParamBand2_Type: fParamBand2_Array[ParamArray_Type] = value; break;
					case kParamBand3_Type: fParamBand3_Array[ParamArray_Type] = value; break;
					case kParamBand4_Type: fParamBand4_Array[ParamArray_Type] = value; break;
					case kParamBand5_Type: fParamBand5_Array[ParamArray_Type] = value; break;

					case kParamBand1_Order: fParamBand1_Array[ParamArray_Order] = value; break;
					case kParamBand2_Order: fParamBand2_Array[ParamArray_Order] = value; break;
					case kParamBand3_Order: fParamBand3_Array[ParamArray_Order] = value; break;
					case kParamBand4_Order: fParamBand4_Array[ParamArray_Order] = value; break;
					case kParamBand5_Order: fParamBand5_Array[ParamArray_Order] = value; break;
					}
				}
			}
		}
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
	return kResultOk;
}

//------------------------------------------------------------------------
uint32 PLUGIN_API RFEQ_Processor::getLatencySamples()
{
	if (fParamOS == overSample_2x) return latency_Fir_x2;
	return 0;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::setupProcessing (Vst::ProcessSetup& newSetup)
{
	//--- called before any processing ----

	Vst::ParamValue OS_target = 0.0;
	if (newSetup.sampleRate <= 48000.0) {
		fParamOS = overSample_2x;
		OS_target = 2 * newSetup.sampleRate;
	}
	else {
		fParamOS = overSample_1x;
		OS_target = newSetup.sampleRate;
	}

	return AudioEffect::setupProcessing (newSetup);
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::canProcessSampleSize (int32 symbolicSampleSize)
{
	// by default kSample32 is supported
	if (symbolicSampleSize == Vst::kSample32)
		return kResultTrue;

	// disable the following comment if your processing support kSample64
	/* if (symbolicSampleSize == Vst::kSample64)
		return kResultTrue; */

	return kResultFalse;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::setState (IBStream* state)
{
	// called when we load a preset, the model has to be reloaded
	IBStreamer streamer(state, kLittleEndian);

	int32           savedBypass = 0;
	Vst::ParamValue savedZoom   = 0.0;
	Vst::ParamValue savedLevel  = 0.0;
	Vst::ParamValue savedOutput = 0.0;

	ParamBand_Array savedBand1_Array = { 0.0, };
	ParamBand_Array savedBand2_Array = { 0.0, };
	ParamBand_Array savedBand3_Array = { 0.0, };
	ParamBand_Array savedBand4_Array = { 0.0, };
	ParamBand_Array savedBand5_Array = { 0.0, };

	if (streamer.readInt32 (savedBypass) == false) return kResultFalse;
	if (streamer.readDouble(savedZoom  ) == false) return kResultFalse;
	if (streamer.readDouble(savedLevel ) == false) return kResultFalse;
	if (streamer.readDouble(savedOutput) == false) return kResultFalse;

	if (streamer.readDoubleArray(savedBand1_Array, ParamArray_size) == false) return kResultFalse;
	if (streamer.readDoubleArray(savedBand2_Array, ParamArray_size) == false) return kResultFalse;
	if (streamer.readDoubleArray(savedBand3_Array, ParamArray_size) == false) return kResultFalse;
	if (streamer.readDoubleArray(savedBand4_Array, ParamArray_size) == false) return kResultFalse;
	if (streamer.readDoubleArray(savedBand5_Array, ParamArray_size) == false) return kResultFalse;

	bBypass = savedBypass > 0;
	fZoom   = savedZoom;
	fLevel  = savedLevel;
	fOutput = savedOutput;

#define Arrcpy(dst, src) \
	dst[ParamArray_In]    = src[ParamArray_In];   \
	dst[ParamArray_Hz]    = src[ParamArray_Hz];   \
	dst[ParamArray_Q]     = src[ParamArray_Q];    \
	dst[ParamArray_dB]    = src[ParamArray_dB];   \
	dst[ParamArray_Type]  = src[ParamArray_Type]; \
	dst[ParamArray_Order] = src[ParamArray_Order]; 

	Arrcpy(fParamBand1_Array, savedBand1_Array)
	Arrcpy(fParamBand2_Array, savedBand2_Array)
	Arrcpy(fParamBand3_Array, savedBand3_Array)
	Arrcpy(fParamBand4_Array, savedBand4_Array)
	Arrcpy(fParamBand5_Array, savedBand5_Array)

	if (Vst::Helpers::isProjectState(state) == kResultTrue)
	{
		// we are in project loading context...

		// Example of using the IStreamAttributes interface
		FUnknownPtr<Vst::IStreamAttributes> stream(state);
		if (stream)
		{
			if (Vst::IAttributeList* list = stream->getAttributes())
			{
				// get the full file path of this state
				Vst::TChar fullPath[1024];
				memset(fullPath, 0, 1024 * sizeof(Vst::TChar));
				if (list->getString(Vst::PresetAttributes::kFilePathStringType, fullPath,
					1024 * sizeof(Vst::TChar)) == kResultTrue)
				{
					// here we have the full path ...
				}
			}
		}
	}
	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::getState (IBStream* state)
{
	// here we need to save the model
	IBStreamer streamer(state, kLittleEndian);

	streamer.writeInt32(bBypass ? 1 : 0);
	streamer.writeDouble(fZoom);
	streamer.writeDouble(fLevel);
	streamer.writeDouble(fOutput);

	streamer.writeDoubleArray(fParamBand1_Array, ParamArray_size);
	streamer.writeDoubleArray(fParamBand2_Array, ParamArray_size);
	streamer.writeDoubleArray(fParamBand3_Array, ParamArray_size);
	streamer.writeDoubleArray(fParamBand4_Array, ParamArray_size);
	streamer.writeDoubleArray(fParamBand5_Array, ParamArray_size);

	return kResultOk;
}



template <typename SampleType>
void RFEQ_Processor::processSVF(
	SampleType** inputs,
	SampleType** outputs,
	Steinberg::int32 numChannels,
	Steinberg::Vst::SampleRate getSampleRate,
	Steinberg::int32 sampleFrames)
{

	Vst::Sample64 _db = (24.0 * fLevel - 12.0);
	Vst::Sample64 In_Atten = exp(log(10.0) * _db / 20.0);

	Vst::SampleRate currFs = getSampleRate;
	if (fParamOS == overSample_2x) currFs = getSampleRate * 2.0;

	int32 oversampling = 1;
	if (fParamOS == overSample_2x) oversampling = 2;

	int32 latency = 0;
	if (fParamOS == overSample_2x) latency = latency_Fir_x2;

	//--- send data ----------------
	if (currentExchangeBlock.blockID == Vst::InvalidDataExchangeBlockID)
		acquireNewExchangeBlock();
	if (auto block = toDataBlock(currentExchangeBlock))
	{
		memcpy(block->Band1, fParamBand1_Array, ParamArray::ParamArray_size * sizeof(double));
		memcpy(block->Band2, fParamBand2_Array, ParamArray::ParamArray_size * sizeof(double));
		memcpy(block->Band3, fParamBand3_Array, ParamArray::ParamArray_size * sizeof(double));
		memcpy(block->Band4, fParamBand4_Array, ParamArray::ParamArray_size * sizeof(double));
		memcpy(block->Band5, fParamBand5_Array, ParamArray::ParamArray_size * sizeof(double));
		block->Fs = currFs;
		block->byPass = bBypass;
		dataExchange->sendCurrentBlock();
		acquireNewExchangeBlock();
	}


#define SVF_set_func(Band, Array, channel) \
	Band[channel].setSVF( \
		Array[ParamArray_In],\
		Array[ParamArray_Hz],\
		Array[ParamArray_Q],\
		Array[ParamArray_dB],\
		Array[ParamArray_Type],\
		Array[ParamArray_Order],\
		currFs\
	);

	for (int32 channel = 0; channel < numChannels; channel++)
	{
		SVF_set_func(Band1, fParamBand1_Array, channel)
		SVF_set_func(Band2, fParamBand2_Array, channel)
		SVF_set_func(Band3, fParamBand3_Array, channel)
		SVF_set_func(Band4, fParamBand4_Array, channel)
		SVF_set_func(Band5, fParamBand5_Array, channel)

		int32 samples = sampleFrames;

		SampleType* ptrIn = (SampleType*)inputs[channel];
		SampleType* ptrOut = (SampleType*)outputs[channel];

		if (latency != latency_q[channel].size()) {
			int32 diff = latency - (int32)latency_q[channel].size();
			if (diff > 0) for (int i = 0; i <  diff; i++) latency_q[channel].push(0.0);
			else          for (int i = 0; i < -diff; i++) latency_q[channel].pop();
		}

		while (--samples >= 0)
		{
			Vst::Sample64 inputSample = *ptrIn; ptrIn++;
			Vst::Sample64 drySample = inputSample;
			inputSample *= In_Atten;

			double up_x[4] = { 0.0, };
			double up_y[4] = { 0.0, };

			up_x[0] = inputSample;

			// Process
			for (int i = 0; i < oversampling; i++)
			{
				Vst::Sample64 overSampled = up_x[i];

				double v1 = Band1[channel].computeSVF(overSampled);
				double v2 = Band2[channel].computeSVF(v1);
				double v3 = Band3[channel].computeSVF(v2);
				double v4 = Band4[channel].computeSVF(v3);
				double v5 = Band5[channel].computeSVF(v4);

				up_y[i] = v5; // * gain;
			}

			// Downsampling
			if (fParamOS == overSample_1x) inputSample = up_y[0];
			else if (fParamOS == overSample_2x) Fir_x2_dn(up_y, &inputSample, channel);

			// Latency compensate
			Vst::Sample64 delayed;
			if (fParamOS == overSample_1x) {
				delayed = drySample;
			}
			else {
				delayed = latency_q[channel].front();
				latency_q[channel].pop();
				latency_q[channel].push(drySample);
			}

			if (bBypass) {
				inputSample = delayed;
			}

			*ptrOut = (SampleType)inputSample;

			ptrOut++;
		}
	}
	return;
}

void RFEQ_Processor::Fir_x2_dn(
	Vst::Sample64* in,
	Vst::Sample64* out,
	Steinberg::int32 channel
)
{
	
	// OS - downsample
	memmove(OS_filter_x2[channel].buff + 1, OS_filter_x2[channel].buff, sizeof(double) * (fir_size - 1));
	OS_filter_x2[channel].buff[0] = in[0];
	
	double acc = 0.0;

	for (int i = 0; i < tap_hm; i++) {
		double a = OS_filter_x2[channel].buff[i];
		double b = OS_filter_x2[channel].buff[fir_size - 1 - i];
		acc += OS_filter_x2[channel].coef[i] * (a + b);
	}   
	acc += OS_filter_x2[channel].coef[tap_hm] * OS_filter_x2[channel].buff[tap_hm];

	memmove(OS_filter_x2[channel].buff + 1, OS_filter_x2[channel].buff, sizeof(double) * (fir_size - 1));
	OS_filter_x2[channel].buff[0] = in[1];

	*out = acc;
}

//------------------------------------------------------------------------
} // namespace yg331
