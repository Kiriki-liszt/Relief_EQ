//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/base/funknown.h"
#include "pluginterfaces/vst/vsttypes.h"

namespace yg331 {
//------------------------------------------------------------------------
static const Steinberg::FUID kRFEQ_ProcessorUID (0xEB8EA8E4, 0xADF058F5, 0xA94A0C17, 0x9DFE6C8D);
static const Steinberg::FUID kRFEQ_ControllerUID (0x7F073B47, 0x0EFA53C7, 0xBFCE3213, 0xC8E59E95);

#define RFEQ_VST3Category "Fx|EQ"

enum {
	kParamBypass = 0,
	kParamZoom,

	kParamLevel,
	kParamOutput,

	kParamBand1_In,
	kParamBand2_In,
	kParamBand3_In,
	kParamBand4_In,
	kParamBand5_In,

	kParamBand1_Hz,
	kParamBand2_Hz,
	kParamBand3_Hz,
	kParamBand4_Hz,
	kParamBand5_Hz,

	kParamBand1_Q,
	kParamBand2_Q,
	kParamBand3_Q,
	kParamBand4_Q,
	kParamBand5_Q,

	kParamBand1_dB,
	kParamBand2_dB,
	kParamBand3_dB,
	kParamBand4_dB,
	kParamBand5_dB,

	kParamBand1_Type,
	kParamBand2_Type,
	kParamBand3_Type,
	kParamBand4_Type,
	kParamBand5_Type,

	kParamBand1_Order,
	kParamBand2_Order,
	kParamBand3_Order,
	kParamBand4_Order,
	kParamBand5_Order,
};
//------------------------------------------------------------------------
} // namespace yg331
