//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "RFEQ_svf.h"
#include "RFEQ_fft.h"
#include "RFEQ_dataexchange.h"

#include "public.sdk/source/vst/vsteditcontroller.h"
#include "vstgui/plugin-bindings/vst3editor.h"

#include "vstgui/lib/controls/cknob.h"


namespace VSTGUI {
	//------------------------------------------------------------------------
	//  TextEdit + Knob mouse control
	//------------------------------------------------------------------------
	class MyKnobText : public CTextEdit, protected CMouseWheelEditingSupport {
	public:
		
		MyKnobText(const CRect& size, IControlListener* listener, int32_t tag, UTF8StringPtr txt, CBitmap* background = nullptr, const int32_t style = 0 );
		MyKnobText(const MyKnobText& v);

		virtual void  setMinPlain(float val)     { minPlain = val; }
		virtual float getMinPlain() const        { return minPlain; }
		virtual void  setMaxPlain(float val)     { maxPlain = val; }
		virtual float getMaxPlain() const        { return maxPlain; }
		int32_t       getLogScale() const        { return logScale; }
		virtual void  setLogScale(int32_t style) { if (style != logScale) { logScale = style; setDirty(); } };

		//-----------------------------------------------------------------------------
		/// @name CKnobBase Methods
		//-----------------------------------------------------------------------------
		//@{
		virtual void  valueToPoint  (CPoint& point) const;
		virtual float valueFromPoint(CPoint& point) const;

		virtual void  setStartAngle(float val) { startAngle = val; compute(); };
		virtual float getStartAngle() const    { return startAngle; }

		virtual void  setRangeAngle(float val) { rangeAngle = val; compute(); };
		virtual float getRangeAngle() const    { return rangeAngle; }

		virtual void  setZoomFactor(float val) { zoomFactor = val; }
		virtual float getZoomFactor() const    { return zoomFactor; }

		virtual CCoord getInsetValue() const   { return inset; }
		virtual void  setInsetValue(CCoord val) { inset = val; }
		//@}

		//-----------------------------------------------------------------------------
		// overrides
		//-----------------------------------------------------------------------------
		void setText          (const UTF8String& txt)  override;
		void onMouseWheelEvent(MouseWheelEvent& event) override;
		void onKeyboardEvent  (KeyboardEvent& event)   override ;
		void setViewSize      (const CRect& rect, bool invalid = true) override ;
		bool sizeToFit        () override;
		void setMin(float val) override { CControl::setMin(val); if (getValue() < val) { setValue(val); } compute(); };
		void setMax(float val) override { CControl::setMax(val); if (getValue() > val) { setValue(val); } compute(); };

		CMouseEventResult onMouseDown  (CPoint& where, const CButtonState& buttons) override;
		CMouseEventResult onMouseUp    (CPoint& where, const CButtonState& buttons) override;
		CMouseEventResult onMouseMoved (CPoint& where, const CButtonState& buttons) override;
		CMouseEventResult onMouseCancel() override;

		CLASS_METHODS(MyKnobText, CTextEdit)

	protected:
		~MyKnobText() noexcept override {};

		void compute() {
			setDirty();
		};

		float startAngle, rangeAngle;
		float zoomFactor;
		float minPlain, maxPlain;
		int32_t logScale;
		CCoord inset;

	private:
		struct MouseEditingState;

		MouseEditingState& getMouseEditingState();
		void clearMouseEditingState();
	};

	//------------------------------------------------------------------------
	// EQ Curve Display
	//------------------------------------------------------------------------
	class EQCurveView : public CControl {
	public:
		EQCurveView(const CRect& size, IControlListener* listener, int32_t tag, CBitmap* background );
		EQCurveView(const EQCurveView& v);

		void setFracOct();
		void setFFTArray(float* array, double sampleRate);
		void setBandArray(double* array, double Fs, bool _byPass, int num);

		// get/set Attributes
		virtual void setBackColor(CColor color) { if (BackColor != color) { BackColor = color; setDirty(true); } }
		CColor getBackColor() const { return BackColor; }

		virtual void setBorderColor(CColor color) { if (BorderColor != color) { BorderColor = color; setDirty(true); } }
		CColor getBorderColor() const { return BorderColor; }

		virtual void setLineColor(CColor color) { if (LineColor != color) { LineColor = color; setDirty(true); } }
		CColor getLineColor() const { return LineColor; }

		virtual void setFFTLineColor(CColor color) { if (FFTLineColor != color) { FFTLineColor = color; setDirty(true); } }
		CColor getFFTLineColor() const { return FFTLineColor; }

		virtual void setFFTFillColor(CColor color) { if (FFTFillColor != color) { FFTFillColor = color; setDirty(true); } }
		CColor getFFTFillColor() const { return FFTFillColor; }

		// overrides
		void setDirty(bool state) override { CView::setDirty(state); };
		void draw(CDrawContext* pContext) override;
		void setViewSize(const CRect& newSize, bool invalid = true) override { CControl::setViewSize(newSize, invalid); };
		bool sizeToFit() override;

		/** called on idle when view wants idle */
		void onIdle() override {
			invalid();
		};

		CLASS_METHODS(EQCurveView, CControl)

		//------------------------------------------------------------------------

	protected:
		~EQCurveView() noexcept override
		{
		};

		CColor		BackColor;
		CColor		LineColor;
		CColor		BorderColor;
		CColor		FFTLineColor;
		CColor		FFTFillColor;

		SVF Band1;
		SVF Band2;
		SVF Band3;
		SVF Band4;
		SVF Band5;

		double filterSamplerate = 96000.0;
		bool byPass;

		static constexpr int fftOrder = yg331::_fftOrder;
		static constexpr int fftSize = 1 << fftOrder;      // 4096 samples
		static constexpr int numBins = fftSize / 2 + 1;    // 2049 bins

		float fft_linear[numBins] = { 0.0 };
		float fft_RMS[numBins] = { 0.0 };
		float fft_freq[numBins] = { 0.0 };

		// 11 octs, 12 per oct = 132 bands max.
#define MAX_BANDS 132
#define INTERVAL 12
		float bandsCenter[MAX_BANDS] = { 0.0, };
		float bandsLower[MAX_BANDS] = { 0.0, };
		float bandsUpper[MAX_BANDS] = { 0.0, };
		float bandsOutput[MAX_BANDS] = { 0.0, };
	};
}


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
	// VST3EditorDelegate
	/** verify a view after it was created */
	VSTGUI::CView* PLUGIN_API verifyView(VSTGUI::CView* view, 
	                                     const VSTGUI::UIAttributes& attributes,
	                                     const VSTGUI::IUIDescription* description, 
	                                     VSTGUI::VST3Editor* editor) SMTG_OVERRIDE;

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


	void viewIsBeingDeleted(VSTGUI::CView* view)
	{
		if (view == EQCurveView_saved)
			EQCurveView_saved = nullptr;
	}

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

	VSTGUI::EQCurveView* EQCurveView_saved = nullptr;

	Steinberg::tchar* Filter_Types[SVF::kFltNum + 1] = {
		(Steinberg::tchar*)STR("Bell"),
		(Steinberg::tchar*)STR("Low Shelf"),
		(Steinberg::tchar*)STR("High Shelf"),
		(Steinberg::tchar*)STR("L Shelf 12"),
		(Steinberg::tchar*)STR("H Shelf 12"),
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
