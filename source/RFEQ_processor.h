//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once
#include "RFEQ_svf.h"
#include "public.sdk/source/vst/vstaudioeffect.h"

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

	//--- ---------------------------------------------------------------------
	// AudioEffect overrides:
	//--- ---------------------------------------------------------------------
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

	/** Here we go...the process call */
	Steinberg::tresult PLUGIN_API process (Steinberg::Vst::ProcessData& data) SMTG_OVERRIDE;
		
	/** For persistence */
	Steinberg::tresult PLUGIN_API setState (Steinberg::IBStream* state) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API getState (Steinberg::IBStream* state) SMTG_OVERRIDE;

//------------------------------------------------------------------------
protected:

	SVF one[2];

	bool                       bBypass = false;
	Steinberg::Vst::ParamValue fLevel  = 0.5;
	Steinberg::Vst::ParamValue fOutput = 0.0;
	Steinberg::Vst::ParamValue fZoom  = 2.0 / 6.0;

	ParamBand_Array fParamBand1_Array = { 1.0, SVF::_Hz_to_norm(20.0), SVF::_Q_to_norm(1.1), SVF::_dB_to_norm(1.0), SVF::_Type_to_norm(SVF::kBell), SVF::_Order_to_norm(SVF::_6dBoct)};
	ParamBand_Array fParamBand2_Array = { 1.0, SVF::_Hz_to_norm(125.0), SVF::_Q_to_norm(2.2), SVF::_dB_to_norm(2.0), SVF::_Type_to_norm(SVF::kBell), SVF::_Order_to_norm(SVF::_12dBoct) };
	ParamBand_Array fParamBand3_Array = { 1.0, SVF::_Hz_to_norm(600.0), SVF::_Q_to_norm(3.3), SVF::_dB_to_norm(3.0), SVF::_Type_to_norm(SVF::kBell), SVF::_Order_to_norm(SVF::_12dBoct) };
	ParamBand_Array fParamBand4_Array = { 1.0, SVF::_Hz_to_norm(5000.0), SVF::_Q_to_norm(4.4), SVF::_dB_to_norm(4.0), SVF::_Type_to_norm(SVF::kBell), SVF::_Order_to_norm(SVF::_12dBoct) };
	ParamBand_Array fParamBand5_Array = { 1.0, SVF::_Hz_to_norm(15000.0), SVF::_Q_to_norm(5.5), SVF::_dB_to_norm(5.0), SVF::_Type_to_norm(SVF::kBell), SVF::_Order_to_norm(SVF::_12dBoct) };

	overSample      fParamOS = overSample_4x;

	SVF Band1[2] = { (SVF::Init_Band1_Hz()), (SVF::Init_Band1_Hz()) };
	SVF Band2[2] = { (SVF::Init_Band2_Hz()), (SVF::Init_Band2_Hz()) };
	SVF Band3[2] = { (SVF::Init_Band3_Hz()), (SVF::Init_Band3_Hz()) };
	SVF Band4[2] = { (SVF::Init_Band4_Hz()), (SVF::Init_Band4_Hz()) };
	SVF Band5[2] = { (SVF::Init_Band5_Hz()), (SVF::Init_Band5_Hz()) };

};

//------------------------------------------------------------------------
} // namespace yg331
