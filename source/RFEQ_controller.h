//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "RFEQ_svf.h"
#include "RFEQ_dataexchange.h"

#include "public.sdk/source/vst/vsteditcontroller.h"
#include "vstgui/plugin-bindings/vst3editor.h"

namespace yg331 {

//------------------------------------------------------------------------
//  RFEQ_Controller
//------------------------------------------------------------------------
class RFEQ_Controller 
	: public Steinberg::Vst::EditControllerEx1
	, public VSTGUI::VST3EditorDelegate
	, public Steinberg::Vst::IDataExchangeReceiver
{
public:
//------------------------------------------------------------------------
	RFEQ_Controller () = default;
	~RFEQ_Controller () SMTG_OVERRIDE = default;

    // Create function
	static Steinberg::FUnknown* createInstance (void* /*context*/)
	{
		return (Steinberg::Vst::IEditController*)new RFEQ_Controller;
	}

	// IPluginBase
	Steinberg::tresult PLUGIN_API initialize (Steinberg::FUnknown* context) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API terminate () SMTG_OVERRIDE;

	// EditController
	Steinberg::tresult PLUGIN_API setComponentState (Steinberg::IBStream* state) SMTG_OVERRIDE;
	Steinberg::IPlugView* PLUGIN_API createView (Steinberg::FIDString name) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API setState (Steinberg::IBStream* state) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API getState (Steinberg::IBStream* state) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API setParamNormalized (Steinberg::Vst::ParamID tag,
                                                      Steinberg::Vst::ParamValue value) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API getParamStringByValue (Steinberg::Vst::ParamID tag,
                                                         Steinberg::Vst::ParamValue valueNormalized,
                                                         Steinberg::Vst::String128 string) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API getParamValueByString (Steinberg::Vst::ParamID tag,
                                                         Steinberg::Vst::TChar* string,
                                                         Steinberg::Vst::ParamValue& valueNormalized) SMTG_OVERRIDE;

	// EditController
	Steinberg::tresult PLUGIN_API notify(Steinberg::Vst::IMessage* message) override;
	void PLUGIN_API update(Steinberg::FUnknown* changedUnknown, Steinberg::int32 message) SMTG_OVERRIDE;
	//void editorAttached(Steinberg::Vst::EditorView* editor) SMTG_OVERRIDE;
	//void editorRemoved(Steinberg::Vst::EditorView* editor) SMTG_OVERRIDE;


	// IDataExchangeReceiver
	void PLUGIN_API queueOpened (Steinberg::Vst::DataExchangeUserContextID userContextID,
	                             Steinberg::uint32 blockSize,
	                             Steinberg::TBool& dispatchOnBackgroundThread) override;
	void PLUGIN_API queueClosed (Steinberg::Vst::DataExchangeUserContextID userContextID) override;
	void PLUGIN_API onDataExchangeBlocksReceived (Steinberg::Vst::DataExchangeUserContextID userContextID,
	                                              Steinberg::uint32 numBlocks, 
	                                              Steinberg::Vst::DataExchangeBlock* blocks,
	                                              Steinberg::TBool onBackgroundThread) override;

 	//---Interface---------
	DEFINE_INTERFACES
		// Here you can add more supported VST3 interfaces
		// DEF_INTERFACE (Vst::IXXX)
		DEF_INTERFACE(IDataExchangeReceiver)
	END_DEFINE_INTERFACES (EditController)
    DELEGATE_REFCOUNT (EditController)

//------------------------------------------------------------------------
protected:
	typedef std::vector<Steinberg::Vst::EditorView*> EditorVector;
	EditorVector editors;
	VSTGUI::VST3Editor* main_editor = nullptr;

	struct ZoomFactor {
		const Steinberg::tchar* title;
		double factor;

		ZoomFactor(const Steinberg::tchar* title, double factor) : title(title), factor(factor) {}
	};
	typedef std::vector<ZoomFactor> ZoomFactorVector;
	ZoomFactorVector zoomFactors;

	Steinberg::Vst::DataExchangeReceiverHandler dataExchange{ this };

	Steinberg::tchar* Filter_Types[SVF::kFltNum + 1] = {
		(Steinberg::tchar*)STR("Bell"),
		(Steinberg::tchar*)STR("Low Shelf"),
		(Steinberg::tchar*)STR("High Shelf"),
		(Steinberg::tchar*)STR("LowShelf +"),
		(Steinberg::tchar*)STR("HighShelf +"),
		(Steinberg::tchar*)STR("Low Pass"),
		(Steinberg::tchar*)STR("High Pass")
	};

	Steinberg::tchar* Filter_Order[SVF::kOrderNum + 1] = {
		(Steinberg::tchar*)STR("6"), 
		(Steinberg::tchar*)STR("12"),
		(Steinberg::tchar*)STR("18"), 
		(Steinberg::tchar*)STR("24")
	};
};

//------------------------------------------------------------------------
} // namespace yg331
