//------------------------------------------------------------------------
// Copyright(c) 2023 Steinberg Media Technologies.
//------------------------------------------------------------------------

#pragma once

// #include "RFEQ_svf.h"

#include "public.sdk/source/vst/utility/dataexchange.h"
#include <cstdint>

//------------------------------------------------------------------------
namespace yg331 {

	//------------------------------------------------------------------------
	struct DataBlock
	{
		/*
		uint32_t sampleRate;
		uint16_t sampleSize;
		uint16_t numChannels;
		uint32_t numSamples;
		float samples[0];
		*/
		ParamBand_Array Band1;
		ParamBand_Array Band2;
		ParamBand_Array Band3;
		ParamBand_Array Band4;
		ParamBand_Array Band5;
		double Fs;
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