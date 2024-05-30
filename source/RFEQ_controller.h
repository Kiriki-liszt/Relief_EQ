//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "RFEQ_svf.h"
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
		
		MyKnobText::MyKnobText(const CRect& size, IControlListener* listener, int32_t tag, UTF8StringPtr txt, CBitmap* background = nullptr, const int32_t style = 0 );
		MyKnobText::MyKnobText(const MyKnobText& v);

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
		EQCurveView(
			const CRect& size, 
			IControlListener* listener, 
			int32_t tag, 
			CBitmap* background
		)
			: CControl(size, listener, tag, background)
		{
			BackColor = kWhiteCColor;
			LineColor = kBlackCColor;
			BorderColor = kBlackCColor;
			setWantsIdle(true);
		}
		EQCurveView(const EQCurveView& v)
			: CControl(v)
			, BackColor(v.BackColor)
			, LineColor(v.LineColor)
			, BorderColor(v.BorderColor)
		{
			setWantsIdle(true);
		}

#define SVF_set(Band, Array, SampleRate) \
			Band.setSVF( \
				Array[ParamArray_In],\
				Array[ParamArray_Hz],\
				Array[ParamArray_Q],\
				Array[ParamArray_dB],\
				Array[ParamArray_Type],\
				Array[ParamArray_Order],\
				SampleRate\
			);

		void setBandArray(double* array, double Fs, int num) {
			switch (num)
			{
			case 1: SVF_set(Band1, array, Fs) break;
			case 2: SVF_set(Band2, array, Fs) break;
			case 3: SVF_set(Band3, array, Fs) break;
			case 4: SVF_set(Band4, array, Fs) break;
			case 5: SVF_set(Band5, array, Fs) break;
			default:
				break;
			}
		}

		// get/set Attributes
		virtual void setBackColor(CColor color) { if (BackColor != color) { BackColor = color; setDirty(true); } }
		CColor getBackColor() const { return BackColor; }

		virtual void setBorderColor(CColor color) { if (BorderColor != color) { BorderColor = color; setDirty(true); } }
		CColor getBorderColor() const { return BorderColor; }

		virtual void setLineColor(CColor color) { if (LineColor != color) { LineColor = color; setDirty(true); } }
		CColor getLineColor() const { return LineColor; }

		// overrides
		void setDirty(bool state) override { CView::setDirty(state); };
		void draw(CDrawContext* pContext) override {

			pContext->setLineWidth(1);
			pContext->setFillColor(getBackColor());
			pContext->setFrameColor(getBorderColor());
			pContext->drawRect(getViewSize(), VSTGUI::kDrawFilledAndStroked);

			double MAX_FREQ = 22000.0;
			double MIN_FREQ = 10.0;
			double FREQ_LOG_MAX = log(MAX_FREQ / MIN_FREQ);
			double DB_EQ_RANGE = 15.0;

			auto border = getBorderColor();
			border.setNormAlpha(0.5);

			{
				VSTGUI::CRect r(getViewSize());
				pContext->setFrameColor(border);
				for (int x = 2; x < 10; x++) {
					VSTGUI::CCoord Hz_10 = r.getWidth() * log(10.0 * x / MIN_FREQ) / FREQ_LOG_MAX;
					const VSTGUI::CPoint _p1(r.left + Hz_10, r.bottom);
					const VSTGUI::CPoint _p2(r.left + Hz_10, r.top);
					pContext->drawLine(_p1, _p2);
				}
				for (int x = 2; x < 10; x++) {
					VSTGUI::CCoord Hz_100 = r.getWidth() * log(100.0 * x / MIN_FREQ) / FREQ_LOG_MAX;
					const VSTGUI::CPoint _p1(r.left + Hz_100, r.bottom);
					const VSTGUI::CPoint _p2(r.left + Hz_100, r.top);
					pContext->drawLine(_p1, _p2);
				}
				for (int x = 2; x < 10; x++) {
					VSTGUI::CCoord Hz_1000 = r.getWidth() * log(1000.0 * x / MIN_FREQ) / FREQ_LOG_MAX;
					const VSTGUI::CPoint _p1(r.left + Hz_1000, r.bottom);
					const VSTGUI::CPoint _p2(r.left + Hz_1000, r.top);
					pContext->drawLine(_p1, _p2);
				}

				for (int x = 2; x < 3; x++) {
					VSTGUI::CCoord Hz_10000 = r.getWidth() * log(10000.0 * x / MIN_FREQ) / FREQ_LOG_MAX;
					const VSTGUI::CPoint _p1(r.left + Hz_10000, r.bottom);
					const VSTGUI::CPoint _p2(r.left + Hz_10000, r.top);
					pContext->drawLine(_p1, _p2);
				}

				pContext->setFrameColor(border);
				{
					VSTGUI::CCoord Hz_100 = r.getWidth() * log(100.0 / MIN_FREQ) / FREQ_LOG_MAX;
					const VSTGUI::CPoint _p1(r.left + Hz_100, r.bottom);
					const VSTGUI::CPoint _p2(r.left + Hz_100, r.top);
					pContext->drawLine(_p1, _p2);
				}
				{
					VSTGUI::CCoord Hz_1000 = r.getWidth() * log(1000.0 / MIN_FREQ) / FREQ_LOG_MAX;
					const VSTGUI::CPoint _p1(r.left + Hz_1000, r.bottom);
					const VSTGUI::CPoint _p2(r.left + Hz_1000, r.top);
					pContext->drawLine(_p1, _p2);
				}
				{
					VSTGUI::CCoord Hz_10000 = r.getWidth() * log(10000.0 / MIN_FREQ) / FREQ_LOG_MAX;
					const VSTGUI::CPoint _p1(r.left + Hz_10000, r.bottom);
					const VSTGUI::CPoint _p2(r.left + Hz_10000, r.top);
					pContext->drawLine(_p1, _p2);
				}
			}

			{
				VSTGUI::CRect r(getViewSize());
				pContext->setFrameColor(border);
				{
					VSTGUI::CCoord dB_p15 = r.getHeight() * (1.0 - (((-15.0 / DB_EQ_RANGE) / 2) + 0.5));
					const VSTGUI::CPoint _p1(r.left, r.bottom - dB_p15);
					const VSTGUI::CPoint _p2(r.right, r.bottom - dB_p15);
					pContext->drawLine(_p1, _p2);
				}
				{
					VSTGUI::CCoord dB_p10 = r.getHeight() * (1.0 - (((-10.0 / DB_EQ_RANGE) / 2) + 0.5));
					const VSTGUI::CPoint _p1(r.left, r.bottom - dB_p10);
					const VSTGUI::CPoint _p2(r.right, r.bottom - dB_p10);
					pContext->drawLine(_p1, _p2);
				}
				{
					VSTGUI::CCoord dB_p5 = r.getHeight() * (1.0 - (((-5.0 / DB_EQ_RANGE) / 2) + 0.5));
					const VSTGUI::CPoint _p1(r.left, r.bottom - dB_p5);
					const VSTGUI::CPoint _p2(r.right, r.bottom - dB_p5);
					pContext->drawLine(_p1, _p2);
				}
				{
					VSTGUI::CCoord dB_0 = r.getHeight() * (1.0 - (((.0 / DB_EQ_RANGE) / 2) + 0.5));
					const VSTGUI::CPoint _p1(r.left, r.bottom - dB_0);
					const VSTGUI::CPoint _p2(r.right, r.bottom - dB_0);
					pContext->drawLine(_p1, _p2);
				}
				{
					VSTGUI::CCoord dB_m5 = r.getHeight() * (1.0 - (((5.0 / DB_EQ_RANGE) / 2) + 0.5));
					const VSTGUI::CPoint _p1(r.left, r.bottom - dB_m5);
					const VSTGUI::CPoint _p2(r.right, r.bottom - dB_m5);
					pContext->drawLine(_p1, _p2);
				}
				{
					VSTGUI::CCoord dB_m10 = r.getHeight() * (1.0 - (((10.0 / DB_EQ_RANGE) / 2) + 0.5));
					const VSTGUI::CPoint _p1(r.left, r.bottom - dB_m10);
					const VSTGUI::CPoint _p2(r.right, r.bottom - dB_m10);
					pContext->drawLine(_p1, _p2);
				}
				{
					VSTGUI::CCoord dB_m15 = r.getHeight() * (1.0 - (((15.0 / DB_EQ_RANGE) / 2) + 0.5));
					const VSTGUI::CPoint _p1(r.left, r.bottom - dB_m15);
					const VSTGUI::CPoint _p2(r.right, r.bottom - dB_m15);
					pContext->drawLine(_p1, _p2);
				}
			}


			VSTGUI::CCoord inset = 30;

			VSTGUI::CGraphicsPath* path = pContext->createGraphicsPath();
			if (path)
			{
				VSTGUI::CRect r(getViewSize());
				//r.inset(inset, 0);


				/*
				path->beginSubpath(VSTGUI::CPoint(r.left + r.getWidth() / 2, r.top));
				path->addLine(VSTGUI::CPoint(r.left, r.bottom));
				path->addLine(VSTGUI::CPoint(r.right, r.bottom));
				path->closeSubpath();
				*/

				VSTGUI::CCoord y_mid = r.bottom - (r.getHeight() / 2.0);
				path->beginSubpath(VSTGUI::CPoint(r.left - 1, y_mid));
				for (int x = -1; x <= r.getWidth() + 1; x++) {
					double tmp = MIN_FREQ * exp(FREQ_LOG_MAX * x / r.getWidth());
					double freq = (std::max)((std::min)(tmp, MAX_FREQ), MIN_FREQ);

					double dB_level = 0.0; // level;

					double dB_1 = 20 * log10(Band1.mag_response(freq));
					double dB_2 = 20 * log10(Band2.mag_response(freq));
					double dB_3 = 20 * log10(Band3.mag_response(freq));
					double dB_4 = 20 * log10(Band4.mag_response(freq));
					double dB_5 = 20 * log10(Band5.mag_response(freq));
					double dB = dB_level + dB_1 + dB_2 + dB_3 + dB_4 + dB_5;

					double m = 1.0 - (((dB / DB_EQ_RANGE) / 2) + 0.5);
					double scy = m * r.getHeight();

					// scy = 0.5 * r.getHeight();
					path->addLine(VSTGUI::CPoint(r.left + x, r.top + scy));
				}
				path->addLine(VSTGUI::CPoint(r.right + 1, r.bottom + 1));
				path->addLine(VSTGUI::CPoint(r.left - 1, r.bottom + 1));
				path->closeSubpath();

				pContext->setFrameColor(getLineColor());
				pContext->setDrawMode(VSTGUI::kAntiAliasing);
				pContext->setLineWidth(1.5);
				pContext->setLineStyle(VSTGUI::kLineSolid);
				pContext->drawGraphicsPath(path, VSTGUI::CDrawContext::kPathStroked);
				path->forget();
			}

			setDirty(false);
		};

		void setViewSize(const CRect& newSize, bool invalid = true) override
		{
			CControl::setViewSize(newSize, invalid);
		};

		bool sizeToFit() override {
			if (getDrawBackground())
			{
				CRect vs(getViewSize());
				vs.setWidth(getDrawBackground()->getWidth());
				vs.setHeight(getDrawBackground()->getHeight());
				setViewSize(vs);
				setMouseableArea(vs);
				return true;
			}
			return false;
		};

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

		SVF Band1;
		SVF Band2;
		SVF Band3;
		SVF Band4;
		SVF Band5;
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
