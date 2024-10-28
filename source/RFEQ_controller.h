//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "RFEQ_svf.h"
#include "RFEQ_fft.h"
#include "RFEQ_dataexchange.h"
#include "RFEQ_shared.h"

#include "public.sdk/source/vst/vsteditcontroller.h"
#include "vstgui/plugin-bindings/vst3editor.h"
#include "vstgui/uidescription/delegationcontroller.h"
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
    virtual void   setInsetValue(CCoord val) { inset = val; }
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

    void setFFTArray(float* array, int sampleBlockSize, double sampleRate);

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
    
    void setParamNorm(Steinberg::Vst::ParamID tag, Steinberg::Vst::ParamValue normValue);
    void setLevel(Steinberg::Vst::ParamValue normValue) {level = yg331::paramGain.ToPlain(normValue);}
    void setBypass(Steinberg::Vst::ParamValue normValue) {byPass = normValue == 0.0 ? false : true;}
    void setEQsampleRate(double SR) {EQ_SR = SR;}

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

    CColor  BackColor;
    CColor  LineColor;
    CColor  BorderColor;
    CColor  FFTLineColor;
    CColor  FFTFillColor;
    
    ParamBand_Array band1 = {1.0, yg331::dftBand1Freq, yg331::dftParamQlty, yg331::dftParamGain, yg331::nrmParamType, yg331::nrmParamOrdr};
    ParamBand_Array band2 = {1.0, yg331::dftBand2Freq, yg331::dftParamQlty, yg331::dftParamGain, yg331::nrmParamType, yg331::nrmParamOrdr};
    ParamBand_Array band3 = {1.0, yg331::dftBand3Freq, yg331::dftParamQlty, yg331::dftParamGain, yg331::nrmParamType, yg331::nrmParamOrdr};
    ParamBand_Array band4 = {1.0, yg331::dftBand4Freq, yg331::dftParamQlty, yg331::dftParamGain, yg331::nrmParamType, yg331::nrmParamOrdr};
    ParamBand_Array band5 = {1.0, yg331::dftBand5Freq, yg331::dftParamQlty, yg331::dftParamGain, yg331::nrmParamType, yg331::nrmParamOrdr};

    SVF     band1_svf;
    SVF     band2_svf;
    SVF     band3_svf;
    SVF     band4_svf;
    SVF     band5_svf;

    bool    byPass;
    double  level;
    double  EQ_SR = 96000.0;

    float fft_linear[yg331::numBins] = { 0.0 };
    float fft_RMS[yg331::numBins]    = { 0.0 };
    float fft_freq[yg331::numBins]   = { 0.0 };
};
}


namespace yg331 {

class EQCurveViewController;

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
    /*
    VSTGUI::CView* PLUGIN_API verifyView(VSTGUI::CView* view,
                                         const VSTGUI::UIAttributes& attributes,
                                         const VSTGUI::IUIDescription* description,
                                         VSTGUI::VST3Editor* editor) SMTG_OVERRIDE;
    */

    // EditController
    Steinberg::tresult PLUGIN_API notify(Steinberg::Vst::IMessage* message) SMTG_OVERRIDE;
    void PLUGIN_API update(Steinberg::FUnknown* changedUnknown, Steinberg::int32 message) SMTG_OVERRIDE;
    void editorAttached(Steinberg::Vst::EditorView* editor) SMTG_OVERRIDE; ///< called from EditorView if it was attached to a parent
    void editorRemoved (Steinberg::Vst::EditorView* editor) SMTG_OVERRIDE; ///< called from EditorView if it was removed from a parent

    //------------------------------------------------------------------------
    Steinberg::Vst::Parameter* getParameterObject(Steinberg::Vst::ParamID tag) SMTG_OVERRIDE
    {
        Steinberg::Vst::Parameter* param = EditControllerEx1::getParameterObject(tag);
        if (param == 0)
        {
            param = uiParameters.getParameter(tag);
        }
        return param;
    }
    bool isPrivateParameter(const Steinberg::Vst::ParamID paramID) SMTG_OVERRIDE
    {
        return uiParameters.getParameter(paramID) != 0 ? true : false;
    }

    // make sure that our UI only parameters doesn't call the following three EditController methods: beginEdit, endEdit, performEdit
    //------------------------------------------------------------------------
    Steinberg::tresult beginEdit(Steinberg::Vst::ParamID tag) SMTG_OVERRIDE
    {
        if (EditControllerEx1::getParameterObject(tag))
            return EditControllerEx1::beginEdit(tag);
        return Steinberg::kResultFalse;
    }

    //------------------------------------------------------------------------
    Steinberg::tresult performEdit(Steinberg::Vst::ParamID tag, Steinberg::Vst::ParamValue valueNormalized) SMTG_OVERRIDE
    {
        if (EditControllerEx1::getParameterObject(tag))
            return EditControllerEx1::performEdit(tag, valueNormalized);
        return Steinberg::kResultFalse;
    }

    //------------------------------------------------------------------------
    Steinberg::tresult endEdit(Steinberg::Vst::ParamID tag) SMTG_OVERRIDE
    {
        if (EditControllerEx1::getParameterObject(tag))
            return EditControllerEx1::endEdit(tag);
        return Steinberg::kResultFalse;
    }


    // IDataExchangeReceiver
    void PLUGIN_API queueOpened (Steinberg::Vst::DataExchangeUserContextID userContextID,
                                 Steinberg::uint32 blockSize,
                                 Steinberg::TBool& dispatchOnBackgroundThread) SMTG_OVERRIDE;
    void PLUGIN_API queueClosed (Steinberg::Vst::DataExchangeUserContextID userContextID) SMTG_OVERRIDE;
    void PLUGIN_API onDataExchangeBlocksReceived (Steinberg::Vst::DataExchangeUserContextID userContextID,
                                                  Steinberg::uint32 numBlocks,
                                                  Steinberg::Vst::DataExchangeBlock* blocks,
                                                  Steinberg::TBool onBackgroundThread) SMTG_OVERRIDE;

    //---from VST3EditorDelegate-----------
    VSTGUI::IController* createSubController (VSTGUI::UTF8StringPtr name,
                                              const VSTGUI::IUIDescription* description,
                                              VSTGUI::VST3Editor* editor) SMTG_OVERRIDE;

    //---Internal functions-------
    void addEQCurveViewController(EQCurveViewController* controller)
    {
        curveControllers.push_back(controller);
    };
    void removeEQCurveViewController(EQCurveViewController* controller)
    {
        auto it = std::find(curveControllers.begin(), curveControllers.end(), controller);
        if (it != curveControllers.end())
            curveControllers.erase(it);
    };

    //---Interface---------
    DEFINE_INTERFACES
        // Here you can add more supported VST3 interfaces
        // DEF_INTERFACE (Vst::IXXX)
        DEF_INTERFACE(IDataExchangeReceiver)
    END_DEFINE_INTERFACES (EditController)
    DELEGATE_REFCOUNT (EditController)

//------------------------------------------------------------------------
protected:
    // UI only parameter list
    Steinberg::Vst::ParameterContainer uiParameters;

    // editor list
    typedef std::vector<Steinberg::Vst::EditorView*> EditorVector;
    EditorVector editors;

    struct ZoomFactor {
        const Steinberg::tchar* title;
        double factor;

        ZoomFactor(const Steinberg::tchar* title, double factor) : title(title), factor(factor) {}
    };
    typedef std::vector<ZoomFactor> ZoomFactorVector;
    ZoomFactorVector zoomFactors;

    Steinberg::Vst::DataExchangeReceiverHandler dataExchange{ this };

    using UICurveControllerList = std::vector<EQCurveViewController*>;
    UICurveControllerList curveControllers;

    Steinberg::tchar* Filter_Types[SVF::tNum + 1] = {
        (Steinberg::tchar*)STR("Bell"),
        (Steinberg::tchar*)STR("Low Shelf"),
        (Steinberg::tchar*)STR("High Shelf"),
        (Steinberg::tchar*)STR("L Shelf 12"),
        (Steinberg::tchar*)STR("H Shelf 12"),
        (Steinberg::tchar*)STR("Low Pass"),
        (Steinberg::tchar*)STR("High Pass")
    };

    Steinberg::tchar* Filter_Order[SVF::oNum + 1] = {
        (Steinberg::tchar*)STR("6"),
        (Steinberg::tchar*)STR("12"),
        (Steinberg::tchar*)STR("18"),
        (Steinberg::tchar*)STR("24")
    };

    ParamValue fZoom   = 2.0 / 6.0;
};


//------------------------------------------------------------------------
// EQCurveViewController
//------------------------------------------------------------------------
class EQCurveViewController :
    public Steinberg::FObject,
    public VSTGUI::DelegationController
{
public:
    using Parameter          = Steinberg::Vst::Parameter;
    using ParameterContainer = Steinberg::Vst::ParameterContainer;
    
    EQCurveViewController ( IController* baseController, RFEQ_Controller* mainController ) :
        DelegationController(baseController),
        mainController(mainController),
        eqCurveView(nullptr)
    {
    }

    ~EQCurveViewController() override
    {
        for (int iter = 0; iter < pBand.getParameterCount(); iter++)
        {
            Parameter* param = pBand.getParameterByIndex(iter);
            if (param) param->removeDependent(this);
        }
        // pBand.removeAll();
        if (pLevel)  { pLevel ->removeDependent(this); pLevel  = nullptr; }
        if (pBypass) { pBypass->removeDependent(this); pBypass = nullptr; }
        
        mainController->removeEQCurveViewController(this);
    }

    enum Tags
    {
        kSendMessageTag = 1000
    };

    void setFFTArray(float* array, int sampleBlockSize, double sampleRate)
    {
        eqCurveView->setFFTArray(array, sampleBlockSize, sampleRate);
    }
    
    void setEQsampleRate(double SR) {eqCurveView->setEQsampleRate(SR);}
    
    void addBandParam(Parameter* pIn, Parameter* pHz, Parameter* pQ, Parameter* pdB, Parameter* pType, Parameter* pOrder)
    {
        if (Steinberg::FCast<Parameter> (pIn))    { pBand.addParameter(pIn);    pIn->addDependent(this); }
        if (Steinberg::FCast<Parameter> (pHz))    { pBand.addParameter(pHz);    pHz->addDependent(this); }
        if (Steinberg::FCast<Parameter> (pQ))     { pBand.addParameter(pQ);     pQ ->addDependent(this); }
        if (Steinberg::FCast<Parameter> (pdB))    { pBand.addParameter(pdB);    pdB->addDependent(this); }
        if (Steinberg::FCast<Parameter> (pType))  { pBand.addParameter(pType);  pType->addDependent(this); }
        if (Steinberg::FCast<Parameter> (pOrder)) { pBand.addParameter(pOrder); pOrder->addDependent(this); }
    }
    void addLevelParam(Parameter* p)
    {
        pLevel = p; p->addDependent(this);
    }
    void addBypassParam(Parameter* p)
    {
        pBypass = p; p->addDependent(this);
    }

private:
    using CControl       = VSTGUI::CControl;
    using CView          = VSTGUI::CView;
    using EQCurveView    = VSTGUI::EQCurveView;
    using UTF8String     = VSTGUI::UTF8String;
    using UIAttributes   = VSTGUI::UIAttributes;
    using IUIDescription = VSTGUI::IUIDescription;
    

    // FObject
    void PLUGIN_API update (Steinberg::FUnknown* changedUnknown, Steinberg::int32 message) SMTG_OVERRIDE
    {
        if (eqCurveView)
        {
            if (auto* p = Steinberg::FCast<Parameter>(changedUnknown); p)
            {
                if (message == kChanged)
                {
                    Steinberg::Vst::ParamID tag = p->getInfo().id;
                    if (Parameter* param = pBand.getParameter(tag); param)
                    {
                        eqCurveView->setParamNorm(tag, param->getNormalized());
                    }
                    if (p == pLevel  && pLevel)  eqCurveView->setLevel (p->getNormalized());
                    if (p == pBypass && pBypass) eqCurveView->setBypass(p->getNormalized());
                }
                else if (message == kWillDestroy)
                {
                    Steinberg::Vst::ParamID tag = p->getInfo().id;
                    if (Parameter* param = pBand.getParameter(tag); param ) { param->removeDependent(this); pBand.removeParameter(tag); }
                    if (p == pLevel  && pLevel)  { pLevel ->removeDependent(this); pLevel  = nullptr; }
                    if (p == pBypass && pBypass) { pBypass->removeDependent(this); pBypass = nullptr; }
                }
            }
        }
    }

    //--- from IControlListener ----------------------
    // void valueChanged(CControl* pControl) override { }
    // void controlBeginEdit(CControl* /*pControl*/) override {}
    // void controlEndEdit(CControl* pControl) override  { }

    //--- from IControlListener ----------------------
    //--- is called when a view is created -----
    CView* verifyView(
        CView* view,
        const UIAttributes& /*attributes*/,
        const IUIDescription* /*description*/) override
    {
        if (EQCurveView* control = dynamic_cast<EQCurveView*>(view); control)
        {
            eqCurveView = control;
            // eqCurveView->registerControlListener(this);
            // curveView->setValue(0.0); // do init thing
        }
        return view;
    }

    RFEQ_Controller* mainController;
    EQCurveView* eqCurveView;
    
    ParameterContainer pBand;
    Parameter* pLevel;
    Parameter* pBypass;
};

//------------------------------------------------------------------------
} // namespace yg331
