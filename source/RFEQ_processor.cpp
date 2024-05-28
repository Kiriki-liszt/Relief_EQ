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
		
		for (int ch = 0; ch < numChannels; ch++) {
			one[ch].setSVF(
				true, 
				48000.0, 
				SVF::_dB_to_norm(9.0), 
				SVF::_Hz_to_norm(500.0), 
				SVF::_Q_to_norm(25.0),
				SVF::_Type_to_norm(SVF::kBell), 
				SVF::_Order_to_norm(SVF::_12dBoct)
			);
			one[ch].makeSVF();
			Vst::Sample32* ptrIn = (Vst::Sample32*)in[ch];
			Vst::Sample32* ptrOut = (Vst::Sample32*)out[ch];
			int32 samples = data.numSamples;
			while (--samples >= 0)
			{
				Vst::Sample64 inputSample = *ptrIn;
				*ptrOut = (Vst::Sample32)one[ch].computeSVF(inputSample);

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
	IBStreamer streamer (state, kLittleEndian);
	
	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Processor::getState (IBStream* state)
{
	// here we need to save the model
	IBStreamer streamer (state, kLittleEndian);

	return kResultOk;
}

//------------------------------------------------------------------------
} // namespace yg331
