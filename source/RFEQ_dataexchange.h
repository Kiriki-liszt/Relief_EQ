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
using ParamValue    = Steinberg::Vst::ParamValue;
using SampleRate    = Steinberg::Vst::SampleRate;
using int32         = Steinberg::int32;
using uint32        = Steinberg::uint32;
using TBool         = Steinberg::TBool;
	//------------------------------------------------------------------------
	struct DataBlock
	{
		SampleRate filterSampleRate;
		SampleRate FFTSampleRate;
		uint32     FFTDataAvail;
		uint32     numSamples;
		float      samples[0];
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
