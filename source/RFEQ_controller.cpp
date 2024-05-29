//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#include "RFEQ_controller.h"
#include "RFEQ_cids.h"
#include "vstgui/plugin-bindings/vst3editor.h"
#include "pluginterfaces/base/ustring.h"
#include "base/source/fstreamer.h"

using namespace Steinberg;

namespace yg331 {


//------------------------------------------------------------------------
// LogRangeParameter Declaration
//------------------------------------------------------------------------
class LogRangeParameter : public Vst::RangeParameter
{
public:
	using RangeParameter::RangeParameter;
	Vst::ParamValue toPlain(Vst::ParamValue _valueNormalized) const SMTG_OVERRIDE;
	Vst::ParamValue toNormalized(Vst::ParamValue plainValue) const SMTG_OVERRIDE;
	void toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const SMTG_OVERRIDE;
};
//------------------------------------------------------------------------
// LogRangeParameter Implementation
//------------------------------------------------------------------------
Vst::ParamValue LogRangeParameter::toPlain(Vst::ParamValue _valueNormalized) const
{
	double FREQ_LOG_MAX = log(getMax() / getMin());
	double tmp = getMin() * exp(FREQ_LOG_MAX * _valueNormalized);
	double freq = (std::max)((std::min)(tmp, getMax()), getMin());
	return freq;
	//return _valueNormalized * (getMax() - getMin()) + getMin();
}

//------------------------------------------------------------------------
Vst::ParamValue LogRangeParameter::toNormalized(Vst::ParamValue plainValue) const
{
	SMTG_ASSERT(getMax() - getMin() != 0);
	double FREQ_LOG_MAX = log(getMax() / getMin());
	return log(plainValue / getMin()) / FREQ_LOG_MAX;
	//return (plainValue - getMin()) / (getMax() - getMin());
}

void LogRangeParameter::toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const
{
	{
		//Parameter::toString(toPlain(_valueNormalized), string);
		UString wrapper(string, str16BufferSize(Vst::String128));
		{
			if (!wrapper.printFloat(toPlain(_valueNormalized), precision))
				string[0] = 0;
			//wrapper.append(STR16(" "));
			//wrapper.append(getInfo().units);
		}
	}
}

//------------------------------------------------------------------------
// LinRangeParameter Declaration
//------------------------------------------------------------------------
class LinRangeParameter : public Vst::RangeParameter
{
public:
	using RangeParameter::RangeParameter;
	void toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const SMTG_OVERRIDE;
};
//------------------------------------------------------------------------
// LinRangeParameter Implementation
//------------------------------------------------------------------------
void LinRangeParameter::toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const
{
	{
		//Parameter::toString(toPlain(_valueNormalized), string);
		UString wrapper(string, str16BufferSize(Vst::String128));
		{
			if (!wrapper.printFloat(toPlain(_valueNormalized), precision))
				string[0] = 0;
			//wrapper.append(STR16(" "));
			//wrapper.append(getInfo().units);
		}
	}
}

//------------------------------------------------------------------------
// RFEQ_Controller Implementation
//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Controller::initialize (FUnknown* context)
{
	// Here the Plug-in will be instantiated

	//---do not forget to call parent ------
	tresult result = EditControllerEx1::initialize (context);
	if (result != kResultOk)
	{
		return result;
	}

	// Here you could register some parameters

	int32 stepCount;
	int32 flags;
	int32 tag;
	Vst::ParamValue defaultVal;
	Vst::ParamValue defaultPlain;
	Vst::ParamValue minPlain;
	Vst::ParamValue maxPlain;

	tag = kParamBypass;
	stepCount = 1;
	defaultVal = Init_Bypass ? 1 : 0;
	flags = Vst::ParameterInfo::kCanAutomate | Vst::ParameterInfo::kIsBypass;
	parameters.addParameter(STR16("Bypass"), nullptr, stepCount, defaultVal, flags, tag);

	if (zoomFactors.empty())
	{
		Vst::ParamValue zoom_coef = 1.0;
		zoomFactors.push_back(ZoomFactor(STR("50%"),  zoom_coef * 0.5));  // 0/6
		zoomFactors.push_back(ZoomFactor(STR("75%"),  zoom_coef * 0.75)); // 1/6
		zoomFactors.push_back(ZoomFactor(STR("100%"), zoom_coef * 1.0));  // 2/6
		zoomFactors.push_back(ZoomFactor(STR("125%"), zoom_coef * 1.25)); // 3/6
		zoomFactors.push_back(ZoomFactor(STR("150%"), zoom_coef * 1.5));  // 4/6
		zoomFactors.push_back(ZoomFactor(STR("175%"), zoom_coef * 1.75)); // 5/6
		zoomFactors.push_back(ZoomFactor(STR("200%"), zoom_coef * 2.0));  // 6/6
	}

	auto zoomParameter = new Vst::StringListParameter(STR("Zoom"), kParamZoom);
	for (auto it = zoomFactors.begin(), end = zoomFactors.end(); it != end; ++it)
	{
		zoomParameter->appendString(it->title);
	}
	zoomParameter->setNormalized(zoomParameter->toNormalized(Init_Zoom)); // toNorm(2) == 100%
	zoomParameter->addDependent(this);
	parameters.addParameter(zoomParameter);

	flags = Vst::ParameterInfo::kCanAutomate;

	minPlain = SVF::getdBMin();
	maxPlain = SVF::getdBMax();
	defaultPlain = 0.0;
	stepCount = 0;

	tag = kParamLevel;
	auto* ParamIn = new LinRangeParameter(STR16("Level"), tag, STR16("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
	ParamIn->setPrecision(2);
	parameters.addParameter(ParamIn);

	tag = kParamOutput;
	auto* ParamOut = new Vst::RangeParameter(STR16("Output"), tag, STR16("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
	ParamOut->setPrecision(2);
	parameters.addParameter(ParamOut);


	stepCount = 1;
	defaultVal = Init_Bypass ? 1 : 0;
	parameters.addParameter(STR16("Band1_In"), nullptr, stepCount, defaultVal, flags, kParamBand1_In);
	parameters.addParameter(STR16("Band2_In"), nullptr, stepCount, defaultVal, flags, kParamBand2_In);
	parameters.addParameter(STR16("Band3_In"), nullptr, stepCount, defaultVal, flags, kParamBand3_In);
	parameters.addParameter(STR16("Band4_In"), nullptr, stepCount, defaultVal, flags, kParamBand4_In);
	parameters.addParameter(STR16("Band5_In"), nullptr, stepCount, defaultVal, flags, kParamBand5_In);


	minPlain = SVF::getdBMin();
	maxPlain = SVF::getdBMax();
	defaultPlain = 0.0;
	stepCount = 0;

	auto* Band1_dB = new LinRangeParameter(STR("Band1_dB"), kParamBand1_dB, STR("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
	auto* Band2_dB = new LinRangeParameter(STR("Band2_dB"), kParamBand2_dB, STR("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
	auto* Band3_dB = new LinRangeParameter(STR("Band3_dB"), kParamBand3_dB, STR("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
	auto* Band4_dB = new LinRangeParameter(STR("Band4_dB"), kParamBand4_dB, STR("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
	auto* Band5_dB = new LinRangeParameter(STR("Band5_dB"), kParamBand5_dB, STR("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
	Band1_dB->setPrecision(2);
	Band2_dB->setPrecision(2);
	Band3_dB->setPrecision(2);
	Band4_dB->setPrecision(2);
	Band5_dB->setPrecision(2);
	parameters.addParameter(Band1_dB);
	parameters.addParameter(Band2_dB);
	parameters.addParameter(Band3_dB);
	parameters.addParameter(Band4_dB);
	parameters.addParameter(Band5_dB);

	minPlain = SVF::getFreqMin();
	maxPlain = SVF::getFreqMax();
	stepCount = 0;

	auto* Band1_Hz = new LogRangeParameter(STR("Band1_Hz"), kParamBand1_Hz, STR("Hz"), minPlain, maxPlain, SVF::Init_Band1_Hz(), stepCount, flags);
	auto* Band2_Hz = new LogRangeParameter(STR("Band2_Hz"), kParamBand2_Hz, STR("Hz"), minPlain, maxPlain, SVF::Init_Band2_Hz(), stepCount, flags);
	auto* Band3_Hz = new LogRangeParameter(STR("Band3_Hz"), kParamBand3_Hz, STR("Hz"), minPlain, maxPlain, SVF::Init_Band3_Hz(), stepCount, flags);
	auto* Band4_Hz = new LogRangeParameter(STR("Band4_Hz"), kParamBand4_Hz, STR("Hz"), minPlain, maxPlain, SVF::Init_Band4_Hz(), stepCount, flags);
	auto* Band5_Hz = new LogRangeParameter(STR("Band5_Hz"), kParamBand5_Hz, STR("Hz"), minPlain, maxPlain, SVF::Init_Band5_Hz(), stepCount, flags);
	Band1_Hz->setPrecision(0);
	Band2_Hz->setPrecision(0);
	Band3_Hz->setPrecision(0);
	Band4_Hz->setPrecision(0);
	Band5_Hz->setPrecision(0);
	parameters.addParameter(Band1_Hz);
	parameters.addParameter(Band2_Hz);
	parameters.addParameter(Band3_Hz);
	parameters.addParameter(Band4_Hz);
	parameters.addParameter(Band5_Hz);


	minPlain = SVF::getQMin();
	maxPlain = SVF::getQMax();
	defaultPlain = 1.414;
	stepCount = 0;

	auto* Band1_Q = new LogRangeParameter(STR16("Band1_Q"), kParamBand1_Q, STR16("Q"), minPlain, maxPlain, defaultPlain, stepCount, flags);
	auto* Band2_Q = new LogRangeParameter(STR16("Band2_Q"), kParamBand2_Q, STR16("Q"), minPlain, maxPlain, defaultPlain, stepCount, flags);
	auto* Band3_Q = new LogRangeParameter(STR16("Band3_Q"), kParamBand3_Q, STR16("Q"), minPlain, maxPlain, defaultPlain, stepCount, flags);
	auto* Band4_Q = new LogRangeParameter(STR16("Band4_Q"), kParamBand4_Q, STR16("Q"), minPlain, maxPlain, defaultPlain, stepCount, flags);
	auto* Band5_Q = new LogRangeParameter(STR16("Band5_Q"), kParamBand5_Q, STR16("Q"), minPlain, maxPlain, defaultPlain, stepCount, flags);
	Band1_Q->setPrecision(2);
	Band2_Q->setPrecision(2);
	Band3_Q->setPrecision(2);
	Band4_Q->setPrecision(2);
	Band5_Q->setPrecision(2);
	parameters.addParameter(Band1_Q);
	parameters.addParameter(Band2_Q);
	parameters.addParameter(Band3_Q);
	parameters.addParameter(Band4_Q);
	parameters.addParameter(Band5_Q);

	auto* Band1_Type = new Vst::StringListParameter(STR16("Band1_Type"), kParamBand1_Type, STR16(""), flags);
	auto* Band2_Type = new Vst::StringListParameter(STR16("Band2_Type"), kParamBand2_Type, STR16(""), flags);
	auto* Band3_Type = new Vst::StringListParameter(STR16("Band3_Type"), kParamBand3_Type, STR16(""), flags);
	auto* Band4_Type = new Vst::StringListParameter(STR16("Band4_Type"), kParamBand4_Type, STR16(""), flags);
	auto* Band5_Type = new Vst::StringListParameter(STR16("Band5_Type"), kParamBand5_Type, STR16(""), flags);

	for (int i = 0; i < SVF::kFltNum + 1; i++) {
		Band1_Type->appendString(Filter_Types[i]);
		Band2_Type->appendString(Filter_Types[i]);
		Band3_Type->appendString(Filter_Types[i]);
		Band4_Type->appendString(Filter_Types[i]);
		Band5_Type->appendString(Filter_Types[i]);
	}

	Band1_Type->setNormalized(SVF::_Type_to_norm(SVF::kBell));
	Band2_Type->setNormalized(SVF::_Type_to_norm(SVF::kBell));
	Band3_Type->setNormalized(SVF::_Type_to_norm(SVF::kBell));
	Band4_Type->setNormalized(SVF::_Type_to_norm(SVF::kBell));
	Band5_Type->setNormalized(SVF::_Type_to_norm(SVF::kBell));

	parameters.addParameter(Band1_Type);
	parameters.addParameter(Band2_Type);
	parameters.addParameter(Band3_Type);
	parameters.addParameter(Band4_Type);
	parameters.addParameter(Band5_Type);


	auto* Band1_Order = new Vst::StringListParameter(STR16("Band1_Order"), kParamBand1_Order, STR16(""), flags);
	auto* Band2_Order = new Vst::StringListParameter(STR16("Band2_Order"), kParamBand2_Order, STR16(""), flags);
	auto* Band3_Order = new Vst::StringListParameter(STR16("Band3_Order"), kParamBand3_Order, STR16(""), flags);
	auto* Band4_Order = new Vst::StringListParameter(STR16("Band4_Order"), kParamBand4_Order, STR16(""), flags);
	auto* Band5_Order = new Vst::StringListParameter(STR16("Band5_Order"), kParamBand5_Order, STR16(""), flags);
	for (int i = 0; i < SVF::kOrderNum + 1; i++) {
		Band1_Order->appendString(Filter_Order[i]);
		Band2_Order->appendString(Filter_Order[i]);
		Band3_Order->appendString(Filter_Order[i]);
		Band4_Order->appendString(Filter_Order[i]);
		Band5_Order->appendString(Filter_Order[i]);
	}
	Band1_Order->setNormalized(SVF::_Order_to_norm(SVF::_12dBoct));
	Band2_Order->setNormalized(SVF::_Order_to_norm(SVF::_12dBoct));
	Band3_Order->setNormalized(SVF::_Order_to_norm(SVF::_12dBoct));
	Band4_Order->setNormalized(SVF::_Order_to_norm(SVF::_12dBoct));
	Band5_Order->setNormalized(SVF::_Order_to_norm(SVF::_12dBoct));
	parameters.addParameter(Band1_Order);
	parameters.addParameter(Band2_Order);
	parameters.addParameter(Band3_Order);
	parameters.addParameter(Band4_Order);
	parameters.addParameter(Band5_Order);

	return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Controller::terminate ()
{
	// Here the Plug-in will be de-instantiated, last possibility to remove some memory!

	//---do not forget to call parent ------
	return EditControllerEx1::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Controller::setComponentState (IBStream* state)
{
	// Here you get the state of the component (Processor part)
	if (!state)
		return kResultFalse;

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

	setParamNormalized(kParamBypass, savedBypass ? 1 : 0);
	setParamNormalized(kParamZoom,   savedZoom);
	setParamNormalized(kParamLevel,  savedLevel);
	setParamNormalized(kParamOutput, savedOutput);

	setParamNormalized(kParamBand1_In, savedBand1_Array[ParamArray_In] ? 1 : 0);
	setParamNormalized(kParamBand2_In, savedBand2_Array[ParamArray_In] ? 1 : 0);
	setParamNormalized(kParamBand3_In, savedBand3_Array[ParamArray_In] ? 1 : 0);
	setParamNormalized(kParamBand4_In, savedBand4_Array[ParamArray_In] ? 1 : 0);
	setParamNormalized(kParamBand5_In, savedBand5_Array[ParamArray_In] ? 1 : 0);
	setParamNormalized(kParamBand1_Hz, savedBand1_Array[ParamArray_Hz]);
	setParamNormalized(kParamBand2_Hz, savedBand2_Array[ParamArray_Hz]);
	setParamNormalized(kParamBand3_Hz, savedBand3_Array[ParamArray_Hz]);
	setParamNormalized(kParamBand4_Hz, savedBand4_Array[ParamArray_Hz]);
	setParamNormalized(kParamBand5_Hz, savedBand5_Array[ParamArray_Hz]);
	setParamNormalized(kParamBand1_Q, savedBand1_Array[ParamArray_Q]);
	setParamNormalized(kParamBand2_Q, savedBand2_Array[ParamArray_Q]);
	setParamNormalized(kParamBand3_Q, savedBand3_Array[ParamArray_Q]);
	setParamNormalized(kParamBand4_Q, savedBand4_Array[ParamArray_Q]);
	setParamNormalized(kParamBand5_Q, savedBand5_Array[ParamArray_Q]);
	setParamNormalized(kParamBand1_dB, savedBand1_Array[ParamArray_dB]);
	setParamNormalized(kParamBand2_dB, savedBand2_Array[ParamArray_dB]);
	setParamNormalized(kParamBand3_dB, savedBand3_Array[ParamArray_dB]);
	setParamNormalized(kParamBand4_dB, savedBand4_Array[ParamArray_dB]);
	setParamNormalized(kParamBand5_dB, savedBand5_Array[ParamArray_dB]);
	setParamNormalized(kParamBand1_Type, savedBand1_Array[ParamArray_Type]);
	setParamNormalized(kParamBand2_Type, savedBand2_Array[ParamArray_Type]);
	setParamNormalized(kParamBand3_Type, savedBand3_Array[ParamArray_Type]);
	setParamNormalized(kParamBand4_Type, savedBand4_Array[ParamArray_Type]);
	setParamNormalized(kParamBand5_Type, savedBand5_Array[ParamArray_Type]);
	setParamNormalized(kParamBand1_Order, savedBand1_Array[ParamArray_Order]);
	setParamNormalized(kParamBand2_Order, savedBand2_Array[ParamArray_Order]);
	setParamNormalized(kParamBand3_Order, savedBand3_Array[ParamArray_Order]);
	setParamNormalized(kParamBand4_Order, savedBand4_Array[ParamArray_Order]);
	setParamNormalized(kParamBand5_Order, savedBand5_Array[ParamArray_Order]);

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Controller::setState (IBStream* state)
{
	/*
	// Here you get the state of the controller
	
	IBStreamer streamer(state, kLittleEndian);

	int8 byteOrder;
	if (streamer.readInt8(byteOrder) == false)
		return kResultFalse;
	if (streamer.readRaw(defaultMessageText, 128 * sizeof(TChar)) == false)
		return kResultFalse;

	// if the byteorder doesn't match, byte swap the text array ...
	if (byteOrder != BYTEORDER)
	{
		for (int32 i = 0; i < 128; i++)
		{
			SWAP_16(defaultMessageText[i])
		}
	}

	// update our editors
	for (auto& uiMessageController : uiMessageControllers)
		uiMessageController->setMessageText(defaultMessageText);
	*/
	return kResultTrue;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Controller::getState (IBStream* state)
{
	// Here you are asked to deliver the state of the controller (if needed)
	// Note: the real state of your plug-in is saved in the processor

	/*
	// here we can save UI settings for example

	// as we save a Unicode string, we must know the byteorder when setState is called

	IBStreamer streamer (state, kLittleEndian);

	int8 byteOrder = BYTEORDER;
	if (streamer.writeInt8 (byteOrder) == false)
		return kResultFalse;

	if (streamer.writeRaw (defaultMessageText, 128 * sizeof (TChar)) == false)
		return kResultFalse;
	*/

	return kResultTrue;
}

//------------------------------------------------------------------------
IPlugView* PLUGIN_API RFEQ_Controller::createView (FIDString name)
{
	// Here the Host wants to open your editor (if you have one)
	if (FIDStringsEqual (name, Vst::ViewType::kEditor))
	{
		// create your editor here and return a IPlugView ptr of it
		auto* view = new VSTGUI::VST3Editor (this, "view", "RFEQ_editor.uidesc");
		view->setZoomFactor(1.0);
		setKnobMode(Steinberg::Vst::KnobModes::kLinearMode);
		main_editor = view;
		return view;
	}
	return nullptr;
}

void PLUGIN_API RFEQ_Controller::update(FUnknown* changedUnknown, int32 message)
{
	EditControllerEx1::update(changedUnknown, message);

	// GUI Resizing
	// check 'zoomtest' code at
	// https://github.com/steinbergmedia/vstgui/tree/vstgui4_10/vstgui/tests/uidescription%20vst3/source

	Vst::Parameter* param = FCast<Vst::Parameter>(changedUnknown);
	if (!param)
		return;

	if (param->getInfo().id == kParamZoom)
	{
		size_t index = static_cast<size_t> (param->toPlain(param->getNormalized()));

		if (index >= zoomFactors.size())
			return;

		if (main_editor)
			main_editor->setZoomFactor(zoomFactors[index].factor);
		
		/*
		for (EditorVector::const_iterator it = editors.begin(), end = editors.end(); it != end; ++it)
		{
			VSTGUI::VST3Editor* editor = dynamic_cast<VSTGUI::VST3Editor*>(*it);
			if (editor)
				editor->setZoomFactor(zoomFactors[index].factor);
		}
		*/
	}
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Controller::setParamNormalized (Vst::ParamID tag, Vst::ParamValue value)
{
	// called by host to update your parameters
	tresult result = EditControllerEx1::setParamNormalized (tag, value);
	return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Controller::getParamStringByValue (Vst::ParamID tag, Vst::ParamValue valueNormalized, Vst::String128 string)
{
	// called by host to get a string for given normalized value of a specific parameter
	// (without having to set the value!)
	return EditControllerEx1::getParamStringByValue (tag, valueNormalized, string);
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Controller::getParamValueByString (Vst::ParamID tag, Vst::TChar* string, Vst::ParamValue& valueNormalized)
{
	// called by host to get a normalized value from a string representation of a specific parameter
	// (without having to set the value!)
	return EditControllerEx1::getParamValueByString (tag, string, valueNormalized);
}


//------------------------------------------------------------------------
// DataExchangeController Implementation
//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Controller::notify(Vst::IMessage* message)
{
	if (dataExchange.onMessage(message))
		return kResultTrue;

	if (!message)
		return kInvalidArgument;

	return EditControllerEx1::notify(message);
}


//------------------------------------------------------------------------
void PLUGIN_API RFEQ_Controller::queueOpened(Vst::DataExchangeUserContextID userContextID,
	uint32 blockSize,
	TBool& dispatchOnBackgroundThread)
{
	//FDebugPrint("Data Exchange Queue opened.\n");
}

//------------------------------------------------------------------------
void PLUGIN_API RFEQ_Controller::queueClosed(Vst::DataExchangeUserContextID userContextID)
{
	//FDebugPrint("Data Exchange Queue closed.\n");
}

//------------------------------------------------------------------------
void PLUGIN_API RFEQ_Controller::onDataExchangeBlocksReceived(
	Vst::DataExchangeUserContextID userContextID,
	uint32 numBlocks,
	Vst::DataExchangeBlock* blocks,
	TBool onBackgroundThread
)
{
	for (auto index = 0u; index < numBlocks; ++index)
	{
		auto dataBlock = toDataBlock(blocks[index]);
		/*
		InMeter[0] = dataBlock->inL;
		InMeter[1] = dataBlock->inR;
		OutMeter[0] = dataBlock->outL;
		OutMeter[1] = dataBlock->outR;
		GRMeter = dataBlock->gR;
		if (!vuMeterControllers.empty()) {
			for (auto iter = vuMeterControllers.begin(); iter != vuMeterControllers.end(); iter++) {
				(*iter)->setMeterValue(InMeter[0] * 0.5 + InMeter[1] * 0.5, kIn);
				(*iter)->setMeterValue(OutMeter[0] * 0.5 + OutMeter[1] * 0.5, kOut);
				(*iter)->setMeterValue(GRMeter, kGR);
				(*iter)->setVuMeterValue(
					InMeter[0], InMeter[1],
					OutMeter[0], OutMeter[1],
					GRMeter
				);
			}
		}
		*/
		/*
		FDebugPrint(
			"Received Data Block: SampleRate: %d, SampleSize: %d, NumChannels: %d, NumSamples: %d\n",
			dataBlock->sampleRate,
			static_cast<uint32_t> (dataBlock->sampleSize),
			static_cast<uint32_t> (dataBlock->numChannels),
			static_cast<uint32_t> (dataBlock->numSamples));
		
		FDebugPrint(\
			"Received Data Block: %f %f %f %f %f\n",\
			dataBlock->inL,\
			dataBlock->inR,\
			dataBlock->outL,\
			dataBlock->outR,\
			dataBlock->gR);
		*/
	}
}

//------------------------------------------------------------------------
} // namespace yg331
