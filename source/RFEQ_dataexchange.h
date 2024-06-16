//------------------------------------------------------------------------
// Copyright(c) 2023 Steinberg Media Technologies.
//------------------------------------------------------------------------

#pragma once

// #include "RFEQ_svf.h"

#include "public.sdk/source/vst/utility/dataexchange.h"
#include "pluginterfaces/vst/vsttypes.h"
#include <cstdint>

//------------------------------------------------------------------------
namespace yg331 {
	using SampleRate = Steinberg::Vst::SampleRate;
	using ParamValue = Steinberg::Vst::ParamValue;
	using TBool = Steinberg::TBool;
	using uint32 = Steinberg::uint32;
	//------------------------------------------------------------------------
	struct DataBlock
	{

		ParamBand_Array Band1;
		ParamBand_Array Band2;
		ParamBand_Array Band3;
		ParamBand_Array Band4;
		ParamBand_Array Band5;
		SampleRate Fs; // SVF filter sample rate, oversampled
		TBool byPass;
		ParamValue level;

		SampleRate sampleRate; // FFT sample rate
		uint32 numSamples;
		float samples[0];
	};

	//------------------------------------------------------------------------
	inline DataBlock* toDataBlock(const Steinberg::Vst::DataExchangeBlock& block)
	{
		if (block.blockID != Steinberg::Vst::InvalidDataExchangeBlockID)
			return reinterpret_cast<DataBlock*> (block.data);
		return nullptr;
	}

	//------------------------------------------------------------------------
} // Steinberg::Tutorial