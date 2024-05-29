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
	
	//---do not forget to call parent ------
	return AudioEffect::terminate ();
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
	int32 numChannels = data.inputs[0].numChannels;

	//---get audio buffers----------------
	uint32 sampleFramesSize = getSampleFramesSizeInBytes(processSetup, data.numSamples);
	void** in = getChannelBuffersPointer(processSetup, data.inputs[0]);
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
		
#define SVF_set_func(Band, Array, channel) \
	Band[channel].setSVF( \
		Array[ParamArray_In],\
		Array[ParamArray_Hz],\
		Array[ParamArray_Q],\
		Array[ParamArray_dB],\
		Array[ParamArray_Type],\
		Array[ParamArray_Order],\
		getSampleRate\
	);

		for (int ch = 0; ch < numChannels; ch++) {
			SVF_set_func(Band1, fParamBand1_Array, ch)
			SVF_set_func(Band2, fParamBand2_Array, ch)
			SVF_set_func(Band3, fParamBand3_Array, ch)
			SVF_set_func(Band4, fParamBand4_Array, ch)
			SVF_set_func(Band5, fParamBand5_Array, ch)

			Vst::Sample32* ptrIn = (Vst::Sample32*)in[ch];
			Vst::Sample32* ptrOut = (Vst::Sample32*)out[ch];
			int32 samples = data.numSamples;
			while (--samples >= 0)
			{
				Vst::Sample64 inputSample = *ptrIn;
				double v1 = Band1[ch].computeSVF(inputSample);
				double v2 = Band2[ch].computeSVF(v1);
				double v3 = Band3[ch].computeSVF(v2);
				double v4 = Band4[ch].computeSVF(v3);
				double v5 = Band5[ch].computeSVF(v4);
				if (bBypass == 1) v5 = inputSample;
				*ptrOut = (Vst::Sample32)v5;

				ptrIn++;
				ptrOut++;
			}
		}
		
	}
	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::setupProcessing (Vst::ProcessSetup& newSetup)
{
	//--- called before any processing ----
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

	std::memcpy(fParamBand1_Array, savedBand1_Array, sizeof(ParamBand_Array));
	std::memcpy(fParamBand2_Array, savedBand2_Array, sizeof(ParamBand_Array));
	std::memcpy(fParamBand3_Array, savedBand3_Array, sizeof(ParamBand_Array));
	std::memcpy(fParamBand4_Array, savedBand4_Array, sizeof(ParamBand_Array));
	std::memcpy(fParamBand5_Array, savedBand5_Array, sizeof(ParamBand_Array));


	Band1[1].copySVF(&Band1[0]);
	Band2[1].copySVF(&Band2[0]);
	Band3[1].copySVF(&Band3[0]);
	Band4[1].copySVF(&Band4[0]);
	Band5[1].copySVF(&Band5[0]);

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

//------------------------------------------------------------------------
} // namespace yg331
