//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#include "RFEQ_controller.h"
#include "RFEQ_cids.h"
#include "vstgui/plugin-bindings/vst3editor.h"
#include "pluginterfaces/base/ustring.h"
#include "base/source/fstreamer.h"

//#include "vstgui/vstgui.h"
#include "vstgui/vstgui_uidescription.h"
#include "vstgui/uidescription/detail/uiviewcreatorattributes.h"

using namespace Steinberg;

namespace VSTGUI {
//------------------------------------------------------------------------
//  TextEdit with Knob mouse control
//------------------------------------------------------------------------
static constexpr CViewAttributeID kCKnobTextMouseStateAttribute = 'ktms';
static const std::string kAttrMinPlain = "min-plain";
static const std::string kAttrMaxPlain = "max-plain";
static const std::string kAttrLogScale = "Log-Scale";
//------------------------------------------------------------------------

static const float kCKnobTextRange = 200.f;

struct MyKnobText::MouseEditingState
{
    CPoint firstPoint;
    CPoint lastPoint;
    float startValue;
    float entryState;
    float range;
    float coef;
    CButtonState oldButton;
    bool modeLinear;
};
MyKnobText::MyKnobText(
    const CRect& size,
    IControlListener* listener,
    int32_t tag,
    UTF8StringPtr txt,
    CBitmap* background,
    const int32_t style
)
    : CTextEdit(size, listener, tag, txt, background)
{
    rangeAngle = 1.f;
    setStartAngle((float)(3.f * Constants::quarter_pi));
    setRangeAngle((float)(3.f * Constants::half_pi));
    zoomFactor = 1.5f;
    minPlain = 20.0;
    maxPlain = 22000.0;
    logScale = false;
}

MyKnobText::MyKnobText(const MyKnobText& v)
    : CTextEdit(v)
    , startAngle(v.startAngle)
    , rangeAngle(v.rangeAngle)
    , zoomFactor(v.zoomFactor)
    , inset(v.inset)
    , minPlain(v.minPlain)
    , maxPlain(v.maxPlain)
    , logScale(v.logScale)
{
}

void  MyKnobText::valueToPoint(CPoint& point) const {
    float alpha = (value - getMin()) / (getMax() - getMin());
    alpha = startAngle + alpha * rangeAngle;

    CPoint c(getViewSize().getWidth() / 2., getViewSize().getHeight() / 2.);
    double xradius = c.x - inset;
    double yradius = c.y - inset;

    point.x = (CCoord)(c.x + cosf(alpha) * xradius + 0.5f);
    point.y = (CCoord)(c.y + sinf(alpha) * yradius + 0.5f);
}

float MyKnobText::valueFromPoint(CPoint& point) const {
    float v;
    double d = rangeAngle * 0.5;
    double a = startAngle + d;

    CPoint c(getViewSize().getWidth() / 2., getViewSize().getHeight() / 2.);
    double xradius = c.x - inset;
    double yradius = c.y - inset;

    double dx = (point.x - c.x) / xradius;
    double dy = (point.y - c.y) / yradius;

    double alpha = atan2(dy, dx) - a;
    while (alpha >= Constants::pi)
        alpha -= Constants::double_pi;
    while (alpha < -Constants::pi)
        alpha += Constants::double_pi;

    if (d < 0.0)
        alpha = -alpha;

    if (alpha > d)
        v = getMax();
    else if (alpha < -d)
        v = getMin();
    else
    {
        v = float(0.5 + alpha / rangeAngle);
        v = getMin() + (v * getRange());
    }

    return v;
};

// overrides
void MyKnobText::setText(const UTF8String& txt) {
    float val = getValue();
    val = atof(txt.getString().c_str());
    if (getLogScale())
        val = log(val / getMinPlain()) / log(getMaxPlain() / getMinPlain());
    else
        val = (val - getMinPlain()) / (getMaxPlain() - getMinPlain());
    CTextLabel::setValue(val);
    CTextEdit::setText(txt);
};

void MyKnobText::onMouseWheelEvent(MouseWheelEvent& event) {
    onMouseWheelEditing(this);

    float v = getValueNormalized();
    if (buttonStateFromEventModifiers(event.modifiers) & kZoomModifier)
        v += 0.1f * static_cast<float> (event.deltaY) * getWheelInc();
    else
        v += static_cast<float> (event.deltaY) * getWheelInc();
    setValueNormalized(v);

    if (isDirty())
    {
        invalid();
        valueChanged();
    }
    event.consumed = true;
};

void MyKnobText::onKeyboardEvent(KeyboardEvent& event) {

    if (!platformControl || event.type != EventType::KeyDown)
        return;

    if (event.virt == VirtualKey::Escape)
    {
        bWasReturnPressed = false;
        platformControl->setText(text);
        getFrame()->setFocusView(nullptr);
        looseFocus();
        event.consumed = true;
    }
    else if (event.virt == VirtualKey::Return)
    {
        bWasReturnPressed = true;
        getFrame()->setFocusView(nullptr);
        looseFocus();
        event.consumed = true;
    }


    if (event.type != EventType::KeyDown)
        return;
    switch (event.virt)
    {
    case VirtualKey::Up:
    case VirtualKey::Right:
    case VirtualKey::Down:
    case VirtualKey::Left:
    {
        float distance = 1.f;
        if (event.virt == VirtualKey::Down || event.virt == VirtualKey::Left)
            distance = -distance;

        float v = getValueNormalized();
        if (buttonStateFromEventModifiers(event.modifiers) & kZoomModifier)
            v += 0.1f * distance * getWheelInc();
        else
            v += distance * getWheelInc();
        setValueNormalized(v);

        if (isDirty())
        {
            invalid();
            beginEdit();
            valueChanged();
            endEdit();
        }
        event.consumed = true;
    }
    case VirtualKey::Escape:
    {
        if (isEditing())
        {
            onMouseCancel();
            event.consumed = true;
        }
        break;
    }
    default: return;
    }
}

void MyKnobText::setViewSize(const CRect& rect, bool invalid)
{
    CControl::setViewSize(rect, invalid);
    compute();
}

bool MyKnobText::sizeToFit() {
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
}

CMouseEventResult MyKnobText::onMouseDown(CPoint& where, const CButtonState& buttons) {
    if (!buttons.isLeftButton())
        return kMouseEventNotHandled;

    if (getFrame()->getFocusView() != this)
    {
        if (isDoubleClickStyle() && (buttons & kDoubleClick))
        {
            takeFocus();
            return kMouseDownEventHandledButDontNeedMovedOrUpEvents;
        }

        // takeFocus();
    }

    CMouseWheelEditingSupport::invalidMouseWheelEditTimer(this);
    beginEdit();

    auto& mouseState = getMouseEditingState();
    mouseState.firstPoint = where;
    mouseState.lastPoint(-1, -1);
    mouseState.startValue = getOldValue();

    mouseState.modeLinear = false;
    mouseState.entryState = value;
    mouseState.range = kCKnobTextRange;
    mouseState.coef = (getMax() - getMin()) / mouseState.range;
    mouseState.oldButton = buttons;

    int32_t mode = kCircularMode;
    int32_t newMode = getFrame()->getKnobMode();
    if (kLinearMode == newMode)
    {
        if (!(buttons & kAlt))
            mode = newMode;
    }
    else if (buttons & kAlt)
    {
        mode = kLinearMode;
    }

    if (mode == kLinearMode)
    {
        if (buttons & kZoomModifier)
            mouseState.range *= zoomFactor;
        mouseState.lastPoint = where;
        mouseState.modeLinear = true;
        mouseState.coef = (getMax() - getMin()) / mouseState.range;
    }
    else
    {
        CPoint where2(where);
        where2.offset(-getViewSize().left, -getViewSize().top);
        mouseState.startValue = valueFromPoint(where2);
        mouseState.lastPoint = where;
    }

    return onMouseMoved(where, buttons);
};

CMouseEventResult MyKnobText::onMouseUp(CPoint& where, const CButtonState& buttons) {
    if (isEditing())
    {
        endEdit();
        clearMouseEditingState();
    }
    return kMouseEventHandled;
};

CMouseEventResult MyKnobText::onMouseMoved(CPoint& where, const CButtonState& buttons) {
    if (buttons.isLeftButton() && isEditing())
    {
        auto& mouseState = getMouseEditingState();

        float middle = (getMax() - getMin()) * 0.5f;

        if (where != mouseState.lastPoint)
        {
            mouseState.lastPoint = where;
            if (mouseState.modeLinear)
            {
                CCoord diff = (mouseState.firstPoint.y - where.y) + (where.x - mouseState.firstPoint.x);
                if (buttons != mouseState.oldButton)
                {
                    mouseState.range = kCKnobTextRange;
                    if (buttons & kZoomModifier)
                        mouseState.range *= zoomFactor;

                    float coef2 = (getMax() - getMin()) / mouseState.range;
                    mouseState.entryState += (float)(diff * (mouseState.coef - coef2));
                    mouseState.coef = coef2;
                    mouseState.oldButton = buttons;
                }
                //value = (float)(mouseState.entryState + diff * mouseState.coef);
                setValue((float)(mouseState.entryState + diff * mouseState.coef));
                bounceValue();
            }
            else
            {
                where.offset(-getViewSize().left, -getViewSize().top);
                //value = valueFromPoint(where);
                setValue(valueFromPoint(where));
                if (mouseState.startValue - value > middle)
                    setValue(getMax()); // value = getMax();
                else if (value - mouseState.startValue > middle)
                    setValue(getMin()); // value = getMin();
                else
                    mouseState.startValue = value;
            }
            if (value != getOldValue())
                valueChanged();
            if (isDirty())
                invalid();
        }
        return kMouseEventHandled;
    }
    return kMouseEventNotHandled;
};

CMouseEventResult MyKnobText::onMouseCancel() {
    if (isEditing())
    {
        auto& mouseState = getMouseEditingState();
        //value = mouseState.startValue;
        setValue(mouseState.startValue);
        if (isDirty())
        {
            valueChanged();
            invalid();
        }
        endEdit();
        clearMouseEditingState();
    }
    return kMouseEventHandled;
};

auto MyKnobText::getMouseEditingState() -> MouseEditingState& {
    MouseEditingState* state = nullptr;
    if (!getAttribute(kCKnobTextMouseStateAttribute, state))
    {
        state = new MouseEditingState;
        setAttribute(kCKnobTextMouseStateAttribute, state);
    }
    return *state;
};

void MyKnobText::clearMouseEditingState() {
    MouseEditingState* state = nullptr;
    if (!getAttribute(kCKnobTextMouseStateAttribute, state))
        return;
    delete state;
    removeAttribute(kCKnobTextMouseStateAttribute);
};


EQCurveView::EQCurveView(
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
    FFTLineColor = kBlackCColor;
    FFTFillColor = kBlackCColor;
    byPass = false;
    level = 0.0;
    idleRate = 60;
    setWantsIdle(true);
}
EQCurveView::EQCurveView(const EQCurveView& v)
    : CControl(v)
    , BackColor(v.BackColor)
    , LineColor(v.LineColor)
    , BorderColor(v.BorderColor)
    , FFTLineColor(v.FFTLineColor)
    , FFTFillColor(v.FFTFillColor)
    , level(v.level)
    , byPass(v.byPass)
{
    setWantsIdle(true);
}

void EQCurveView::setParamNorm(Steinberg::Vst::ParamID tag, Steinberg::Vst::ParamValue normValue)
{
    switch (tag) {
        case yg331::kParamBand1_In: band1[yg331::bandIn] = normValue; break;
        case yg331::kParamBand2_In: band2[yg331::bandIn] = normValue; break;
        case yg331::kParamBand3_In: band3[yg331::bandIn] = normValue; break;
        case yg331::kParamBand4_In: band4[yg331::bandIn] = normValue; break;
        case yg331::kParamBand5_In: band5[yg331::bandIn] = normValue; break;
            
        case yg331::kParamBand1_Hz: band1[yg331::bandHz] = yg331::paramFreq.ToPlain(normValue); break;
        case yg331::kParamBand2_Hz: band2[yg331::bandHz] = yg331::paramFreq.ToPlain(normValue); break;
        case yg331::kParamBand3_Hz: band3[yg331::bandHz] = yg331::paramFreq.ToPlain(normValue); break;
        case yg331::kParamBand4_Hz: band4[yg331::bandHz] = yg331::paramFreq.ToPlain(normValue); break;
        case yg331::kParamBand5_Hz: band5[yg331::bandHz] = yg331::paramFreq.ToPlain(normValue); break;
            
        case yg331::kParamBand1_Q:  band1[yg331::bandQ]  = yg331::paramQlty.ToPlain(normValue); break;
        case yg331::kParamBand2_Q:  band2[yg331::bandQ]  = yg331::paramQlty.ToPlain(normValue); break;
        case yg331::kParamBand3_Q:  band3[yg331::bandQ]  = yg331::paramQlty.ToPlain(normValue); break;
        case yg331::kParamBand4_Q:  band4[yg331::bandQ]  = yg331::paramQlty.ToPlain(normValue); break;
        case yg331::kParamBand5_Q:  band5[yg331::bandQ]  = yg331::paramQlty.ToPlain(normValue); break;
            
        case yg331::kParamBand1_dB: band1[yg331::banddB] = yg331::paramGain.ToPlain(normValue); break;
        case yg331::kParamBand2_dB: band2[yg331::banddB] = yg331::paramGain.ToPlain(normValue); break;
        case yg331::kParamBand3_dB: band3[yg331::banddB] = yg331::paramGain.ToPlain(normValue); break;
        case yg331::kParamBand4_dB: band4[yg331::banddB] = yg331::paramGain.ToPlain(normValue); break;
        case yg331::kParamBand5_dB: band5[yg331::banddB] = yg331::paramGain.ToPlain(normValue); break;
            
        case yg331::kParamBand1_Type: band1[yg331::bandType] = (yg331::paramType.ToPlain(normValue)); break;
        case yg331::kParamBand2_Type: band2[yg331::bandType] = (yg331::paramType.ToPlain(normValue)); break;
        case yg331::kParamBand3_Type: band3[yg331::bandType] = (yg331::paramType.ToPlain(normValue)); break;
        case yg331::kParamBand4_Type: band4[yg331::bandType] = (yg331::paramType.ToPlain(normValue)); break;
        case yg331::kParamBand5_Type: band5[yg331::bandType] = (yg331::paramType.ToPlain(normValue)); break;
            
        case yg331::kParamBand1_Order: band1[yg331::bandOrder] = (yg331::paramOrdr.ToPlain(normValue)); break;
        case yg331::kParamBand2_Order: band2[yg331::bandOrder] = (yg331::paramOrdr.ToPlain(normValue)); break;
        case yg331::kParamBand3_Order: band3[yg331::bandOrder] = (yg331::paramOrdr.ToPlain(normValue)); break;
        case yg331::kParamBand4_Order: band4[yg331::bandOrder] = (yg331::paramOrdr.ToPlain(normValue)); break;
        case yg331::kParamBand5_Order: band5[yg331::bandOrder] = (yg331::paramOrdr.ToPlain(normValue)); break;
        default: break;
    }
    band1_svf.setSVF(band1[yg331::bandIn], band1[yg331::bandHz], band1[yg331::bandQ], band1[yg331::banddB], band1[yg331::bandType], band1[yg331::bandOrder], EQ_SR);
    band2_svf.setSVF(band2[yg331::bandIn], band2[yg331::bandHz], band2[yg331::bandQ], band2[yg331::banddB], band2[yg331::bandType], band2[yg331::bandOrder], EQ_SR);
    band3_svf.setSVF(band3[yg331::bandIn], band3[yg331::bandHz], band3[yg331::bandQ], band3[yg331::banddB], band3[yg331::bandType], band3[yg331::bandOrder], EQ_SR);
    band4_svf.setSVF(band4[yg331::bandIn], band4[yg331::bandHz], band4[yg331::bandQ], band4[yg331::banddB], band4[yg331::bandType], band4[yg331::bandOrder], EQ_SR);
    band5_svf.setSVF(band5[yg331::bandIn], band5[yg331::bandHz], band5[yg331::bandQ], band5[yg331::banddB], band5[yg331::bandType], band5[yg331::bandOrder], EQ_SR);
}

#define cubic_hermite(A, B, C, D, t) \
    (/*a0*/(D - C - A + B) * t * t * t + /*a1*/(A - B - D + C + A - B) * t * t + /*a2*/(C - A) * t + /*a3*/(B) )

#define lanczos(x, a) \
    (a * sin(M_PI * (double)x) * sin(M_PI * (double)x / a) / (M_PI * M_PI * (double)x * (double)x))

#define lanczos_kernel(x, a) \
    ((x == 0) ? ( 1.0 ) : ((-a < x) && (x < a) ? (lanczos(x, a)) : (0.0)))

#define lanczos_interpolation(A, B, C, D, t) \
    (A * lanczos_kernel((t + 1.0), 2) + B * lanczos_kernel(t, 2) + C * lanczos_kernel((t - 1.0), 2) + D * lanczos_kernel((t - 2.0), 2))

#define lanczos_interpolation_6(Z, A, B, C, D, E, t) \
    (Z * lanczos_kernel((t + 2.0), 3) + (A * lanczos_kernel((t + 1.0), 3) + B * lanczos_kernel(t, 3) + C * lanczos_kernel((t - 1.0), 3) + D * lanczos_kernel((t - 2.0), 3))  + E * lanczos_kernel((t - 3.0), 3))

#define safe_bin(bin, x)   std::max(std::min((int)bin  + (int)x, (int)numBins   - 1) , 0)
#define safe_band(band, x) std::max(std::min((int)band + (int)x, (int)MAX_BANDS - 1) , 0)

void EQCurveView::setFFTArray(float* array, int sampleBlockSize, double sampleRate)
{
    // Unit frequency per bin, with sample rate
    double freqBin_width = sampleRate / yg331::fftSize;
    double coeff = exp(-1.0 / (0.3 * 0.001 * sampleRate));
    // double coeff = exp(-1.0 * (double)sampleBlockSize / (0.03 * 0.001 * sampleRate)); // gets faster as block size gets large
    // double coeff = exp(-1.0 / (0.03 * 0.001 * sampleRate * (double)sampleBlockSize)); // gets slower as block size gets large

    double icoef = 1.0 - coeff;

    for (int i = 0; i < yg331::numBins; ++i) {
        if (std::isnan(array[i])) array[i] = 0.00000001;
        fft_RMS[i] = (fft_RMS[i] * coeff) + (icoef * array[i]);
        fft_linear[i] = fft_RMS[i];
        fft_freq[i] = (i + 0.5) * freqBin_width;
    }
}


// overrides
void EQCurveView::draw(CDrawContext* pContext) {

    pContext->setLineWidth(1);
    pContext->setFillColor(getBackColor());
    pContext->setFrameColor(getBorderColor());
    pContext->drawRect(getViewSize(), VSTGUI::kDrawFilledAndStroked);

    // Given frequency, return screen x position
    auto freq_to_x = [this](double width, double freq) -> double {
        return width * log(freq / MIN_FREQ) / FREQ_LOG_MAX;
    };

    // Given screen x position, return frequency
    auto x_to_freq = [this](double width, double x) -> double {
        return std::max(std::min(MIN_FREQ * exp(FREQ_LOG_MAX * x / width), MAX_FREQ), MIN_FREQ);
    };

    // Given a magnitude, return y screen position as 0..1 with applied tilt
    auto mag_to_01 = [](double m, double freq) -> double {
        if (m == 0) m = 0.00001;
        if (freq == 0) freq = 1.0;
        return 1.0 - (( ((20.0 * log10(m)) + (4.5 * ((log(freq) / log(2.0)) - (log(1024.0) / log(2.0))))) - ceiling) / (noise_floor - ceiling));
    };

    // Given a magnitude (1.0 .... very small number), return y screen position
    auto mag_to_y = [](double height, double m) -> double {
        if (m == 0) m = 0.00001;
        return (((20.0 * log10(m)) - ceiling) / (noise_floor - ceiling)) * height;
    };

    // Given decibels, return screen y position
    auto db_to_y = [](double height, double dB) -> double {
        return (((dB - ceiling) / (noise_floor - ceiling)) * height);
    };

    // Given screen y position, return decibels
    auto y_to_db = [](double height, double y) -> double {
        return ceiling + ((y / height) * (noise_floor - ceiling));
    };

    auto dB_to_y_EQ = [](double height, double dB) -> double {
        return height * (1.0 - (((dB / DB_EQ_RANGE) / 2) + 0.5));
    };
    

    auto border = getBorderColor();
    border.setNormAlpha(0.5);
    
    const VSTGUI::CRect r(getViewSize());
    const double r_width = r.getWidth();
    const double r_height = r.getHeight();

    {
        pContext->setFrameColor(border);
        for (int x = 2; x < 10; x++) {
            VSTGUI::CCoord Hz_10 = freq_to_x(r_width, 10.0 * x);
            const VSTGUI::CPoint _p1(r.left + Hz_10, r.bottom);
            const VSTGUI::CPoint _p2(r.left + Hz_10, r.top);
            pContext->drawLine(_p1, _p2);
        }
        for (int x = 1; x < 10; x++) {
            VSTGUI::CCoord Hz_100 = freq_to_x(r_width, 100.0 * x);
            const VSTGUI::CPoint _p1(r.left + Hz_100, r.bottom);
            const VSTGUI::CPoint _p2(r.left + Hz_100, r.top);
            pContext->drawLine(_p1, _p2);
        }
        for (int x = 1; x < 10; x++) {
            VSTGUI::CCoord Hz_1000 = freq_to_x(r_width, 1000.0 * x);
            const VSTGUI::CPoint _p1(r.left + Hz_1000, r.bottom);
            const VSTGUI::CPoint _p2(r.left + Hz_1000, r.top);
            pContext->drawLine(_p1, _p2);
        }

        for (int x = 1; x < 3; x++) {
            VSTGUI::CCoord Hz_10000 = freq_to_x(r_width, 10000.0 * x);
            const VSTGUI::CPoint _p1(r.left + Hz_10000, r.bottom);
            const VSTGUI::CPoint _p2(r.left + Hz_10000, r.top);
            pContext->drawLine(_p1, _p2);
        }
    }

    {
        pContext->setFrameColor(border);
        for (int cnt = -(int)DB_EQ_RANGE; cnt < (int)DB_EQ_RANGE; cnt += 5)
        {
            VSTGUI::CCoord dB = dB_to_y_EQ(r_height, cnt);
            const VSTGUI::CPoint _p1(r.left,  r.bottom - dB);
            const VSTGUI::CPoint _p2(r.right, r.bottom - dB);
            pContext->drawLine(_p1, _p2);
        }
    }


    VSTGUI::CGraphicsPath* FFT_curve = pContext->createGraphicsPath();
    if (FFT_curve)
    {
        double y_start = mag_to_01(fft_linear[0], fft_freq[0]);
        y_start = (std::max)((std::min)(y_start, 1.0), 0.0);
        y_start *= r_height;
        FFT_curve->beginSubpath(VSTGUI::CPoint(r.left - 1, r.bottom - y_start));
        double x_last = 0.0;
        // RAW
        for (int bin = 0; bin < yg331::numBins; ++bin)
        {
            double x = freq_to_x(r_width, fft_freq[bin]);
            x = (std::max)((std::min)(x, r_width), 0.0);
            double y = mag_to_01(fft_linear[bin], fft_freq[bin]);
            y = (std::max)((std::min)(y, 1.0), 0.0);
            y *= r_height;
            if (x - x_last > 0.1)
            {
                x_last = x;
                FFT_curve->addLine(VSTGUI::CPoint(r.left + x, r.bottom - y));
            }
        }

        FFT_curve->addLine(VSTGUI::CPoint(r.right + 1, r.bottom + 1));
        FFT_curve->addLine(VSTGUI::CPoint(r.left - 1, r.bottom + 1));
        FFT_curve->closeSubpath();


        VSTGUI::CColor ff = getFFTFillColor();
        ff.setNormAlpha(0.8);
        pContext->setFrameColor(VSTGUI::kTransparentCColor);
        pContext->setFillColor(ff);
        pContext->setDrawMode(VSTGUI::kAntiAliasing);
        pContext->setLineWidth(0.0);
        pContext->setLineStyle(VSTGUI::kLineSolid);
        pContext->drawGraphicsPath(FFT_curve, VSTGUI::CDrawContext::kPathFilled);


        pContext->setFrameColor(getFFTLineColor());
        pContext->setDrawMode(VSTGUI::kAntiAliasing);
        pContext->setLineWidth(1.0);
        pContext->setLineStyle(VSTGUI::kLineSolid);
        pContext->drawGraphicsPath(FFT_curve, VSTGUI::CDrawContext::kPathStroked);

        FFT_curve->forget();
    }

    VSTGUI::CGraphicsPath* EQ_curve = pContext->createGraphicsPath();
    if (EQ_curve)
    {
        VSTGUI::CCoord y_mid = r.bottom - (r.getHeight() / 2.0);
        EQ_curve->beginSubpath(VSTGUI::CPoint(r.left - 1, y_mid));
        for (double x = -0.5; x <= r.getWidth() + 1; x+=0.5)
        {
            double tmp = MIN_FREQ * std::exp(FREQ_LOG_MAX * x / r.getWidth());
            double freq = (std::max)((std::min)(tmp, MAX_FREQ), MIN_FREQ);

            double dB_1 = 20 * log10(band1_svf.mag_response(freq));
            double dB_2 = 20 * log10(band2_svf.mag_response(freq));
            double dB_3 = 20 * log10(band3_svf.mag_response(freq));
            double dB_4 = 20 * log10(band4_svf.mag_response(freq));
            double dB_5 = 20 * log10(band5_svf.mag_response(freq));
            double dB = level + dB_1 + dB_2 + dB_3 + dB_4 + dB_5;

            double m = 1.0 - (((dB / DB_EQ_RANGE) / 2) + 0.5);
            double scy = m * r.getHeight();

            if (byPass) scy = 0.5 * r.getHeight();
            EQ_curve->addLine(VSTGUI::CPoint(r.left + x, r.top + scy));
        }
        EQ_curve->addLine(VSTGUI::CPoint(r.right + 1, r.bottom + 1));
        EQ_curve->addLine(VSTGUI::CPoint(r.left - 1, r.bottom + 1));
        EQ_curve->closeSubpath();

        pContext->setFrameColor(getLineColor());
        pContext->setDrawMode(VSTGUI::kAntiAliasing);
        pContext->setLineWidth(1.5);
        pContext->setLineStyle(VSTGUI::kLineSolid);
        pContext->drawGraphicsPath(EQ_curve, VSTGUI::CDrawContext::kPathStroked);
        EQ_curve->forget();
    }

    // box outline
    pContext->setLineWidth(1);
    pContext->setFrameColor(getBorderColor());
    pContext->drawRect(getViewSize(), VSTGUI::kDrawStroked);

    setDirty(false);
};

bool EQCurveView::sizeToFit() {
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


static const std::string kAttrBackColor = "back-color";
static const std::string kAttrBorderColor = "border-color";
static const std::string kAttrLineColor = "line-color";
static const std::string kAttrFFTLineColor = "FFT-line-color";
static const std::string kAttrFFTFillColor = "FFT-fill-color";

//------------------------------------------------------------------------
//  Factory for TextEdit
//------------------------------------------------------------------------
class MyKnobTextFactory : public ViewCreatorAdapter
{
public:
    //register this class with the view factory
    MyKnobTextFactory() { UIViewFactory::registerViewCreator(*this); }

    //return an unique name here
    IdStringPtr getViewName() const override { return "KnobText"; }

    //return the name here from where your custom view inherites.
    //	Your view automatically supports the attributes from it.
    IdStringPtr getBaseViewName() const override { return UIViewCreator::kCTextEdit; }

    //create your view here.
    //	Note you don't need to apply attributes here as
    //	the apply method will be called with this new view
    CView* create(const UIAttributes& attributes, const IUIDescription* description) const override
    {
        return new MyKnobText(CRect(0, 0, 100, 20), nullptr, -1, nullptr, nullptr);
    }
    bool apply(
        CView* view,
        const UIAttributes& attributes,
        const IUIDescription* description) const override
    {
        auto* KnobText = dynamic_cast<MyKnobText*> (view);

        if (!KnobText)
            return false;

        double d;
        if (attributes.getDoubleAttribute(UIViewCreator::kAttrAngleStart, d))
        {
            // convert from degree
            d = d / 180.f * static_cast<float> (Constants::pi);
            KnobText->setStartAngle(static_cast<float> (d));
        }
        if (attributes.getDoubleAttribute(UIViewCreator::kAttrAngleRange, d))
        {
            // convert from degree
            d = d / 180.f * static_cast<float> (Constants::pi);
            KnobText->setRangeAngle(static_cast<float> (d));
        }
        if (attributes.getDoubleAttribute(UIViewCreator::kAttrValueInset, d))
            KnobText->setInsetValue(d);
        if (attributes.getDoubleAttribute(UIViewCreator::kAttrZoomFactor, d))
            KnobText->setZoomFactor(static_cast<float> (d));
        if (attributes.getDoubleAttribute(kAttrMinPlain, d))
            KnobText->setMinPlain(static_cast<float> (d));
        if (attributes.getDoubleAttribute(kAttrMaxPlain, d))
            KnobText->setMaxPlain(static_cast<float> (d));

        bool b;
        if (attributes.getBooleanAttribute(kAttrLogScale, b))
            KnobText->setLogScale(b);

        return true;
    }

    bool getAttributeNames(StringList& attributeNames) const override
    {
        attributeNames.emplace_back(UIViewCreator::kAttrAngleStart);
        attributeNames.emplace_back(UIViewCreator::kAttrAngleRange);
        attributeNames.emplace_back(UIViewCreator::kAttrValueInset);
        attributeNames.emplace_back(UIViewCreator::kAttrZoomFactor);
        attributeNames.emplace_back(kAttrMinPlain);
        attributeNames.emplace_back(kAttrMaxPlain);
        attributeNames.emplace_back(kAttrLogScale);
        return true;
    }

    AttrType getAttributeType(const std::string& attributeName) const override
    {
        if (attributeName == UIViewCreator::kAttrAngleStart)
            return kFloatType;
        if (attributeName == UIViewCreator::kAttrAngleRange)
            return kFloatType;
        if (attributeName == UIViewCreator::kAttrValueInset)
            return kFloatType;
        if (attributeName == UIViewCreator::kAttrZoomFactor)
            return kFloatType;
        if (attributeName == kAttrMinPlain)
            return kFloatType;
        if (attributeName == kAttrMaxPlain)
            return kFloatType;
        if (attributeName == kAttrLogScale)
            return kBooleanType;
        return kUnknownType;
    }

    //------------------------------------------------------------------------
    bool getAttributeValue(
        CView* view,
        const string& attributeName,
        string& stringValue,
        const IUIDescription* desc) const override
    {
        auto* KnobText = dynamic_cast<MyKnobText*> (view);
        if (!KnobText)
            return false;

        if (attributeName == UIViewCreator::kAttrAngleStart)
        {
            stringValue =
                UIAttributes::doubleToString((KnobText->getStartAngle() / Constants::pi * 180.), 5);
            return true;
        }
        if (attributeName == UIViewCreator::kAttrAngleRange)
        {
            stringValue =
                UIAttributes::doubleToString((KnobText->getRangeAngle() / Constants::pi * 180.), 5);
            return true;
        }
        if (attributeName == UIViewCreator::kAttrValueInset)
        {
            stringValue = UIAttributes::doubleToString(KnobText->getInsetValue());
            return true;
        }
        if (attributeName == kAttrMinPlain)
        {
            stringValue = UIAttributes::doubleToString(KnobText->getMinPlain());
            return true;
        }
        if (attributeName == kAttrMaxPlain)
        {
            stringValue = UIAttributes::doubleToString(KnobText->getMaxPlain());
            return true;
        }
        if (attributeName == kAttrLogScale)
        {
            stringValue = KnobText->getLogScale() ? UIViewCreator::strTrue : UIViewCreator::strFalse;
            return true;
        }
        return false;
    }
};

class MyEQCurveViewFactory : public ViewCreatorAdapter
{
public:
    //register this class with the view factory
    MyEQCurveViewFactory() { UIViewFactory::registerViewCreator(*this); }

    //return an unique name here
    IdStringPtr getViewName() const override { return "EQ Curve View"; }

    //return the name here from where your custom view inherites.
    //	Your view automatically supports the attributes from it.
    IdStringPtr getBaseViewName() const override { return UIViewCreator::kCControl; }

    //create your view here.
    //	Note you don't need to apply attributes here as
    //	the apply method will be called with this new view
    CView* create(const UIAttributes& attributes, const IUIDescription* description) const override
    {
        return new EQCurveView(CRect(0, 0, 100, 20), nullptr, -1, nullptr);
    }
    bool apply(
        CView* view,
        const UIAttributes& attributes,
        const IUIDescription* description) const override
    {
        auto* v = dynamic_cast<EQCurveView*> (view);

        if (!v)
            return false;

        CColor color;
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrBackColor), color, description))
            v->setBackColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrBorderColor), color, description))
            v->setBorderColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrLineColor), color, description))
            v->setLineColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrFFTLineColor), color, description))
            v->setFFTLineColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrFFTFillColor), color, description))
            v->setFFTFillColor(color);

        return true;
    }

    bool getAttributeNames(StringList& attributeNames) const override
    {
        attributeNames.emplace_back(kAttrBackColor);
        attributeNames.emplace_back(kAttrBorderColor);
        attributeNames.emplace_back(kAttrLineColor);
        attributeNames.emplace_back(kAttrFFTLineColor);
        attributeNames.emplace_back(kAttrFFTFillColor);
        return true;
    }

    AttrType getAttributeType(const std::string& attributeName) const override
    {
        if (attributeName == kAttrBackColor)
            return kColorType;
        if (attributeName == kAttrBorderColor)
            return kColorType;
        if (attributeName == kAttrLineColor)
            return kColorType;
        if (attributeName == kAttrFFTLineColor)
            return kColorType;
        if (attributeName == kAttrFFTFillColor)
            return kColorType;
        return kUnknownType;
    }

    //------------------------------------------------------------------------
    bool getAttributeValue(
        CView* view,
        const string& attributeName,
        string& stringValue,
        const IUIDescription* desc) const override
    {
        auto* v = dynamic_cast<EQCurveView*> (view);

        if (!v)
            return false;

        if (attributeName == kAttrBackColor)
        {
            UIViewCreator::colorToString(v->getBackColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrBorderColor)
        {
            UIViewCreator::colorToString(v->getBorderColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrLineColor)
        {
            UIViewCreator::colorToString(v->getLineColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrFFTLineColor)
        {
            UIViewCreator::colorToString(v->getFFTLineColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrFFTFillColor)
        {
            UIViewCreator::colorToString(v->getFFTFillColor(), stringValue, desc);
            return true;
        }

        return false;
    }
};


//create a static instance so that it registers itself with the view factory
MyKnobTextFactory    __gMyMyKnobTextFactory;
MyEQCurveViewFactory __gMyEQCurveViewFactory;
} // namespace VSTGUI

namespace yg331 {


//------------------------------------------------------------------------
// LogRangeParameter Declaration
//------------------------------------------------------------------------
class LogRangeParameter_noUnit : public Vst::RangeParameter
{
public:
    using RangeParameter::RangeParameter;
    LogRangeParameter_noUnit (const Vst::TChar* title, Vst::ParamID tag, const Vst::TChar* units = nullptr,
                       Vst::ParamValue minPlain = 0., Vst::ParamValue maxPlain = 1.,
                       Vst::ParamValue defaultValuePlain = 0., int32 stepCount = 0,
                       int32 flags = Steinberg::Vst::ParameterInfo::kCanAutomate, Vst::UnitID unitID = Steinberg::Vst::kRootUnitId,
                       const Vst::TChar* shortTitle = nullptr)
    : Vst::RangeParameter(title, tag, units, minPlain, maxPlain, defaultValuePlain, stepCount, flags, unitID, shortTitle)
    {
        UString (info.title, str16BufferSize (Vst::String128)).assign (title);
        if (units)
            UString (info.units, str16BufferSize (Vst::String128)).assign (units);
        if (shortTitle)
            UString (info.shortTitle, str16BufferSize (Vst::String128)).assign (shortTitle);

        info.stepCount = stepCount;
        info.defaultNormalizedValue = valueNormalized = toNormalized (defaultValuePlain);
        info.flags = flags;
        info.id = tag;
        info.unitId = unitID;
    }
    
    /** Converts a normalized value to plain value (e.g. 0.5 to 10000.0Hz). */
    Vst::ParamValue toPlain(Vst::ParamValue _valueNormalized) const SMTG_OVERRIDE;
    
    /** Converts a plain value to a normalized value (e.g. 10000 to 0.5). */
    Vst::ParamValue toNormalized(Vst::ParamValue plainValue) const SMTG_OVERRIDE;
    
    /** Converts a normalized value to a string. */
    void toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const SMTG_OVERRIDE;
    
    OBJ_METHODS (LogRangeParameter_noUnit, RangeParameter)
};
//------------------------------------------------------------------------
// LogRangeParameter Implementation
//------------------------------------------------------------------------
Vst::ParamValue LogRangeParameter_noUnit::toPlain(Vst::ParamValue _valueNormalized) const
{
    double FREQ_LOG_MAX = std::log(getMax() / getMin());
    double tmp = getMin() * std::exp(FREQ_LOG_MAX * _valueNormalized);
    double freq = (std::max)((std::min)(tmp, getMax()), getMin());
    return freq;
    //return _valueNormalized * (getMax() - getMin()) + getMin();
}

//------------------------------------------------------------------------
Vst::ParamValue LogRangeParameter_noUnit::toNormalized(Vst::ParamValue plainValue) const
{
    SMTG_ASSERT(getMax() - getMin() != 0);
    double FREQ_LOG_MAX = std::log(getMax() / getMin());
    return std::log(plainValue / getMin()) / FREQ_LOG_MAX;
    //return (plainValue - getMin()) / (getMax() - getMin());
}

void LogRangeParameter_noUnit::toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const
{
    {
        //Parameter::toString(toPlain(_valueNormalized), string);
        UString wrapper(string, str16BufferSize(Vst::String128));
        {
            if (!wrapper.printFloat(toPlain(_valueNormalized), precision))
                string[0] = 0;
            // wrapper.append(STR16(" "));
            // wrapper.append(getInfo().units);
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
    LinRangeParameter (const Vst::TChar* title, Vst::ParamID tag, const Vst::TChar* units = nullptr,
                       Vst::ParamValue minPlain = 0., Vst::ParamValue maxPlain = 1.,
                       Vst::ParamValue defaultValuePlain = 0., int32 stepCount = 0,
                       int32 flags = Steinberg::Vst::ParameterInfo::kCanAutomate, Vst::UnitID unitID = Steinberg::Vst::kRootUnitId,
                       const Vst::TChar* shortTitle = nullptr)
    : Vst::RangeParameter(title, tag, units, minPlain, maxPlain, defaultValuePlain, stepCount, flags, unitID, shortTitle)
    {
        UString (info.title, str16BufferSize (Vst::String128)).assign (title);
        if (units)
            UString (info.units, str16BufferSize (Vst::String128)).assign (units);
        if (shortTitle)
            UString (info.shortTitle, str16BufferSize (Vst::String128)).assign (shortTitle);

        info.stepCount = stepCount;
        info.defaultNormalizedValue = valueNormalized = toNormalized (defaultValuePlain);
        info.flags = flags;
        info.id = tag;
        info.unitId = unitID;
    }
    
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
            wrapper.append(STR16(" "));
            wrapper.append(getInfo().units);
        }
    }
}

//------------------------------------------------------------------------
// LinRangeParameter Declaration
//------------------------------------------------------------------------
class LinRangeParameter_noUnit : public Vst::RangeParameter
{
public:
    using RangeParameter::RangeParameter;
    LinRangeParameter_noUnit (const Vst::TChar* title, Vst::ParamID tag, const Vst::TChar* units = nullptr,
                       Vst::ParamValue minPlain = 0., Vst::ParamValue maxPlain = 1.,
                       Vst::ParamValue defaultValuePlain = 0., int32 stepCount = 0,
                       int32 flags = Steinberg::Vst::ParameterInfo::kCanAutomate, Vst::UnitID unitID = Steinberg::Vst::kRootUnitId,
                       const Vst::TChar* shortTitle = nullptr)
    : Vst::RangeParameter(title, tag, units, minPlain, maxPlain, defaultValuePlain, stepCount, flags, unitID, shortTitle)
    {
        UString (info.title, str16BufferSize (Vst::String128)).assign (title);
        if (units)
            UString (info.units, str16BufferSize (Vst::String128)).assign (units);
        if (shortTitle)
            UString (info.shortTitle, str16BufferSize (Vst::String128)).assign (shortTitle);

        info.stepCount = stepCount;
        info.defaultNormalizedValue = valueNormalized = toNormalized (defaultValuePlain);
        info.flags = flags;
        info.id = tag;
        info.unitId = unitID;
    }
    void toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const SMTG_OVERRIDE;
};
//------------------------------------------------------------------------
// LinRangeParameter Implementation
//------------------------------------------------------------------------
void LinRangeParameter_noUnit::toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const
{
    {
        //Parameter::toString(toPlain(_valueNormalized), string);
        UString wrapper(string, str16BufferSize(Vst::String128));
        {
            if (!wrapper.printFloat(toPlain(_valueNormalized), precision))
                string[0] = 0;
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
    defaultVal = dftBypass ? 1 : 0;
    flags = Vst::ParameterInfo::kCanAutomate | Vst::ParameterInfo::kIsBypass;
    parameters.addParameter(STR16("Bypass"), nullptr, stepCount, defaultVal, flags, tag);

    flags = Vst::ParameterInfo::kCanAutomate;

    stepCount = 0;

    tag = kParamLevel;
    auto* ParamLevel = new LinRangeParameter(STR16("Level"), tag, STR16("dB"), minParamGain, maxParamGain, dftParamGain, stepCount, flags);
    ParamLevel->setPrecision(1);
    parameters.addParameter(ParamLevel);

    stepCount = 1;
    defaultVal = 1;
    parameters.addParameter(STR16("Band1_In"), nullptr, stepCount, defaultVal, flags, kParamBand1_In);
    parameters.addParameter(STR16("Band2_In"), nullptr, stepCount, defaultVal, flags, kParamBand2_In);
    parameters.addParameter(STR16("Band3_In"), nullptr, stepCount, defaultVal, flags, kParamBand3_In);
    parameters.addParameter(STR16("Band4_In"), nullptr, stepCount, defaultVal, flags, kParamBand4_In);
    parameters.addParameter(STR16("Band5_In"), nullptr, stepCount, defaultVal, flags, kParamBand5_In);

    stepCount = 0;

    auto* Band1_dB = new LinRangeParameter_noUnit(STR("Band1_dB"), kParamBand1_dB, STR("dB"), minParamGain, maxParamGain, dftParamGain, stepCount, flags);
    auto* Band2_dB = new LinRangeParameter_noUnit(STR("Band2_dB"), kParamBand2_dB, STR("dB"), minParamGain, maxParamGain, dftParamGain, stepCount, flags);
    auto* Band3_dB = new LinRangeParameter_noUnit(STR("Band3_dB"), kParamBand3_dB, STR("dB"), minParamGain, maxParamGain, dftParamGain, stepCount, flags);
    auto* Band4_dB = new LinRangeParameter_noUnit(STR("Band4_dB"), kParamBand4_dB, STR("dB"), minParamGain, maxParamGain, dftParamGain, stepCount, flags);
    auto* Band5_dB = new LinRangeParameter_noUnit(STR("Band5_dB"), kParamBand5_dB, STR("dB"), minParamGain, maxParamGain, dftParamGain, stepCount, flags);
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
    
    auto* Band1_Hz = new LogRangeParameter_noUnit(STR("Band1_Hz"), kParamBand1_Hz, STR("Hz"), minParamFreq, maxParamFreq, dftBand1Freq, stepCount, flags);
    auto* Band2_Hz = new LogRangeParameter_noUnit(STR("Band2_Hz"), kParamBand2_Hz, STR("Hz"), minParamFreq, maxParamFreq, dftBand2Freq, stepCount, flags);
    auto* Band3_Hz = new LogRangeParameter_noUnit(STR("Band3_Hz"), kParamBand3_Hz, STR("Hz"), minParamFreq, maxParamFreq, dftBand3Freq, stepCount, flags);
    auto* Band4_Hz = new LogRangeParameter_noUnit(STR("Band4_Hz"), kParamBand4_Hz, STR("Hz"), minParamFreq, maxParamFreq, dftBand4Freq, stepCount, flags);
    auto* Band5_Hz = new LogRangeParameter_noUnit(STR("Band5_Hz"), kParamBand5_Hz, STR("Hz"), minParamFreq, maxParamFreq, dftBand5Freq, stepCount, flags);
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

    auto* Band1_Q = new LogRangeParameter_noUnit(STR16("Band1_Q"), kParamBand1_Q, STR16("Q"), minParamQlty, maxParamQlty, dftParamQlty, stepCount, flags);
    auto* Band2_Q = new LogRangeParameter_noUnit(STR16("Band2_Q"), kParamBand2_Q, STR16("Q"), minParamQlty, maxParamQlty, dftParamQlty, stepCount, flags);
    auto* Band3_Q = new LogRangeParameter_noUnit(STR16("Band3_Q"), kParamBand3_Q, STR16("Q"), minParamQlty, maxParamQlty, dftParamQlty, stepCount, flags);
    auto* Band4_Q = new LogRangeParameter_noUnit(STR16("Band4_Q"), kParamBand4_Q, STR16("Q"), minParamQlty, maxParamQlty, dftParamQlty, stepCount, flags);
    auto* Band5_Q = new LogRangeParameter_noUnit(STR16("Band5_Q"), kParamBand5_Q, STR16("Q"), minParamQlty, maxParamQlty, dftParamQlty, stepCount, flags);
    Band1_Q->setPrecision(1);
    Band2_Q->setPrecision(1);
    Band3_Q->setPrecision(1);
    Band4_Q->setPrecision(1);
    Band5_Q->setPrecision(1);
    parameters.addParameter(Band1_Q);
    parameters.addParameter(Band2_Q);
    parameters.addParameter(Band3_Q);
    parameters.addParameter(Band4_Q);
    parameters.addParameter(Band5_Q);
    
    flags = Vst::ParameterInfo::kCanAutomate | Vst::ParameterInfo::kIsList;

    auto* Band1_Type = new Vst::StringListParameter(STR16("Band1_Type"), kParamBand1_Type, STR16(""), flags);
    auto* Band2_Type = new Vst::StringListParameter(STR16("Band2_Type"), kParamBand2_Type, STR16(""), flags);
    auto* Band3_Type = new Vst::StringListParameter(STR16("Band3_Type"), kParamBand3_Type, STR16(""), flags);
    auto* Band4_Type = new Vst::StringListParameter(STR16("Band4_Type"), kParamBand4_Type, STR16(""), flags);
    auto* Band5_Type = new Vst::StringListParameter(STR16("Band5_Type"), kParamBand5_Type, STR16(""), flags);

    for (int i = 0; i < SVF::tNum + 1; i++) {
        Band1_Type->appendString(Filter_Types[i]);
        Band2_Type->appendString(Filter_Types[i]);
        Band3_Type->appendString(Filter_Types[i]);
        Band4_Type->appendString(Filter_Types[i]);
        Band5_Type->appendString(Filter_Types[i]);
    }

    Band1_Type->setNormalized(nrmParamType);
    Band2_Type->setNormalized(nrmParamType);
    Band3_Type->setNormalized(nrmParamType);
    Band4_Type->setNormalized(nrmParamType);
    Band5_Type->setNormalized(nrmParamType);

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
    for (int i = 0; i < SVF::oNum + 1; i++) {
        Band1_Order->appendString(Filter_Order[i]);
        Band2_Order->appendString(Filter_Order[i]);
        Band3_Order->appendString(Filter_Order[i]);
        Band4_Order->appendString(Filter_Order[i]);
        Band5_Order->appendString(Filter_Order[i]);
    }
    Band1_Order->setNormalized(nrmParamOrdr);
    Band2_Order->setNormalized(nrmParamOrdr);
    Band3_Order->setNormalized(nrmParamOrdr);
    Band4_Order->setNormalized(nrmParamOrdr);
    Band5_Order->setNormalized(nrmParamOrdr);
    parameters.addParameter(Band1_Order);
    parameters.addParameter(Band2_Order);
    parameters.addParameter(Band3_Order);
    parameters.addParameter(Band4_Order);
    parameters.addParameter(Band5_Order);
    
    // GUI only parameter
    Vst::StringListParameter* zoomParameter = new Vst::StringListParameter(STR("Zoom"), kParamZoom);
    for (ZoomFactorVector::const_iterator it = zoomFactors.begin(), end = zoomFactors.end(); it != end; ++it)
    {
        zoomParameter->appendString(it->title);
    }
    zoomParameter->setNormalized(zoomParameter->toNormalized(dftZoom));
    zoomParameter->getInfo().defaultNormalizedValue = zoomParameter->toNormalized(dftZoom);
    zoomParameter->addDependent(this);
    uiParameters.addParameter(zoomParameter);
    
    return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Controller::terminate ()
{
    // Here the Plug-in will be de-instantiated, last possibility to remove some memory!

    getParameterObject(kParamZoom)->removeDependent(this);
    
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
    
    ParamBand_Array savedBand1_Array = {1.0, nrmBand1Freq, nrmParamQlty, nrmParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array savedBand2_Array = {1.0, nrmBand2Freq, nrmParamQlty, nrmParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array savedBand3_Array = {1.0, nrmBand3Freq, nrmParamQlty, nrmParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array savedBand4_Array = {1.0, nrmBand4Freq, nrmParamQlty, nrmParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array savedBand5_Array = {1.0, nrmBand5Freq, nrmParamQlty, nrmParamGain, nrmParamType, nrmParamOrdr};

    if (streamer.readInt32 (savedBypass) == false) return kResultFalse;
    if (streamer.readDouble(savedZoom  ) == false) return kResultFalse;
    if (streamer.readDouble(savedLevel ) == false) return kResultFalse;
    if (streamer.readDouble(savedOutput) == false) return kResultFalse;

    if (streamer.readDoubleArray(savedBand1_Array, bandSize) == false) return kResultFalse;
    if (streamer.readDoubleArray(savedBand2_Array, bandSize) == false) return kResultFalse;
    if (streamer.readDoubleArray(savedBand3_Array, bandSize) == false) return kResultFalse;
    if (streamer.readDoubleArray(savedBand4_Array, bandSize) == false) return kResultFalse;
    if (streamer.readDoubleArray(savedBand5_Array, bandSize) == false) return kResultFalse;

    auto setParamNormArray = [this](double Array[], int num)
    {
        setParamNormalized(kParamBand1_In    + num, Array[bandIn] ? 1 : 0);
        setParamNormalized(kParamBand1_Hz    + num, Array[bandHz]);
        setParamNormalized(kParamBand1_Q     + num, Array[bandQ]);
        setParamNormalized(kParamBand1_dB    + num, Array[banddB]);
        setParamNormalized(kParamBand1_Type  + num, Array[bandType]);
        setParamNormalized(kParamBand1_Order + num, Array[bandOrder]);
    };
    
    setParamNormalized(kParamBypass, savedBypass ? 1 : 0);
    // setParamNormalized(kParamZoom, savedZoom);
    setParamNormalized(kParamLevel, savedLevel);
    // setParamNormalized(kParamOutput, savedOutput);
    
    setParamNormArray(savedBand1_Array, 0);
    setParamNormArray(savedBand2_Array, 1);
    setParamNormArray(savedBand3_Array, 2);
    setParamNormArray(savedBand4_Array, 3);
    setParamNormArray(savedBand5_Array, 4);
    
    auto copyArray = [this](double dst[], double src[]) {
        dst[bandIn]    = src[bandIn];
        dst[bandHz]    = src[bandHz];
        dst[bandQ]     = src[bandQ];
        dst[banddB]    = src[banddB];
        dst[bandType]  = src[bandType];
        dst[bandOrder] = src[bandOrder];
    };

    pBypass = savedBypass > 0;
    // fZoom   = savedZoom;
    fLevel  = savedLevel;
    // fOutput = savedOutput;

    copyArray(fParamBand1_Array, savedBand1_Array);
    copyArray(fParamBand2_Array, savedBand2_Array);
    copyArray(fParamBand3_Array, savedBand3_Array);
    copyArray(fParamBand4_Array, savedBand4_Array);
    copyArray(fParamBand5_Array, savedBand5_Array);

    return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Controller::setState (IBStream* state)
{
    // Here you get the state of the controller

    if (!state)
        return kResultFalse;

    IBStreamer streamer(state, kLittleEndian);
    
    auto copyArray = [this](double dst[], double src[]) {
        dst[bandIn]    = src[bandIn];
        dst[bandHz]    = src[bandHz];
        dst[bandQ]     = src[bandQ];
        dst[banddB]    = src[banddB];
        dst[bandType]  = src[bandType];
        dst[bandOrder] = src[bandOrder];
    };

    Vst::ParamValue savedZoom   = dftZoom/zoomNum;
    Vst::ParamValue savedLevel  = 0.0;
    Vst::ParamValue savedOutput = 0.0;
    
    ParamBand_Array savedBand1_Array = {1.0, nrmBand1Freq, nrmParamQlty, nrmParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array savedBand2_Array = {1.0, nrmBand2Freq, nrmParamQlty, nrmParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array savedBand3_Array = {1.0, nrmBand3Freq, nrmParamQlty, nrmParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array savedBand4_Array = {1.0, nrmBand4Freq, nrmParamQlty, nrmParamGain, nrmParamType, nrmParamOrdr};
    ParamBand_Array savedBand5_Array = {1.0, nrmBand5Freq, nrmParamQlty, nrmParamGain, nrmParamType, nrmParamOrdr};

    if (streamer.readDouble(savedZoom  ) == false) savedZoom = dftZoom/zoomNum;
    if (streamer.readDouble(savedLevel ) == false) savedLevel = fLevel;
    if (streamer.readDouble(savedOutput) == false) savedOutput = 0.0;

    if (streamer.readDoubleArray(savedBand1_Array, bandSize) == false) copyArray(savedBand1_Array, fParamBand1_Array);
    if (streamer.readDoubleArray(savedBand2_Array, bandSize) == false) copyArray(savedBand2_Array, fParamBand2_Array);
    if (streamer.readDoubleArray(savedBand3_Array, bandSize) == false) copyArray(savedBand3_Array, fParamBand3_Array);
    if (streamer.readDoubleArray(savedBand4_Array, bandSize) == false) copyArray(savedBand4_Array, fParamBand4_Array);
    if (streamer.readDoubleArray(savedBand5_Array, bandSize) == false) copyArray(savedBand5_Array, fParamBand5_Array);

    auto setParamNormArray = [this](double Array[], int num)
    {
        setParamNormalized(kParamBand1_In    + num, Array[bandIn] ? 1 : 0);
        setParamNormalized(kParamBand1_Hz    + num, Array[bandHz]);
        setParamNormalized(kParamBand1_Q     + num, Array[bandQ]);
        setParamNormalized(kParamBand1_dB    + num, Array[banddB]);
        setParamNormalized(kParamBand1_Type  + num, Array[bandType]);
        setParamNormalized(kParamBand1_Order + num, Array[bandOrder]);
    };

    setParamNormalized(kParamZoom,   savedZoom);
    setParamNormalized(kParamLevel,  savedLevel);
    // setParamNormalized(kParamOutput, savedOutput);

    setParamNormArray(savedBand1_Array, 0);
    setParamNormArray(savedBand2_Array, 1);
    setParamNormArray(savedBand3_Array, 2);
    setParamNormArray(savedBand4_Array, 3);
    setParamNormArray(savedBand5_Array, 4);

    fZoom   = savedZoom;
    fLevel  = savedLevel;
    // fOutput = savedOutput;

    copyArray(fParamBand1_Array, savedBand1_Array);
    copyArray(fParamBand2_Array, savedBand2_Array);
    copyArray(fParamBand3_Array, savedBand3_Array);
    copyArray(fParamBand4_Array, savedBand4_Array);
    copyArray(fParamBand5_Array, savedBand5_Array);
    
    return kResultTrue;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RFEQ_Controller::getState (IBStream* state)
{
    // Here you are asked to deliver the state of the controller (if needed)
    // Note: the real state of your plug-in is saved in the processor

    if (!state)
        return kResultFalse;

    IBStreamer streamer(state, kLittleEndian);

    fZoom   = getParamNormalized(kParamZoom);
    fLevel  = getParamNormalized(kParamLevel);
    fOutput = getParamNormalized(kParamOutput);
    
    auto getParamNormArray = [this](double Array[], int num)
    {
        Array[bandIn]    = getParamNormalized(kParamBand1_In    + num);
        Array[bandHz]    = getParamNormalized(kParamBand1_Hz    + num);
        Array[bandQ]     = getParamNormalized(kParamBand1_Q     + num);
        Array[banddB]    = getParamNormalized(kParamBand1_dB    + num);
        Array[bandType]  = getParamNormalized(kParamBand1_Type  + num);
        Array[bandOrder] = getParamNormalized(kParamBand1_Order + num);
    };
    
    getParamNormArray(fParamBand1_Array, 0);
    getParamNormArray(fParamBand2_Array, 1);
    getParamNormArray(fParamBand3_Array, 2);
    getParamNormArray(fParamBand4_Array, 3);
    getParamNormArray(fParamBand5_Array, 4);

    streamer.writeDouble(fZoom);
    streamer.writeDouble(fLevel);
    streamer.writeDouble(fOutput);

    streamer.writeDoubleArray(fParamBand1_Array, bandSize);
    streamer.writeDoubleArray(fParamBand2_Array, bandSize);
    streamer.writeDoubleArray(fParamBand3_Array, bandSize);
    streamer.writeDoubleArray(fParamBand4_Array, bandSize);
    streamer.writeDoubleArray(fParamBand5_Array, bandSize);

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
        std::vector<double> _zoomFactors;
        _zoomFactors.push_back(0.50);
        _zoomFactors.push_back(0.75);
        _zoomFactors.push_back(1.00);
        _zoomFactors.push_back(1.25);
        _zoomFactors.push_back(1.50);
        _zoomFactors.push_back(1.75);
        _zoomFactors.push_back(2.00);
        view->setAllowedZoomFactors(_zoomFactors);
        view->setZoomFactor(1.0);
        view->setIdleRate(1000.0/60.0);
        setKnobMode(Steinberg::Vst::KnobModes::kLinearMode);
        return view;
    }
    return nullptr;
}

VSTGUI::IController* RFEQ_Controller::createSubController (VSTGUI::UTF8StringPtr name,
                                          const VSTGUI::IUIDescription* description,
                                          VSTGUI::VST3Editor* editor)
{
    if (VSTGUI::UTF8StringView(name) == "EQCurveViewController")
    {
        EQCurveViewController* controller = new EQCurveViewController(editor, this);
        controller->addBandParam(getParameterObject(kParamBand1_In),
                                 getParameterObject(kParamBand1_Hz),
                                 getParameterObject(kParamBand1_Q),
                                 getParameterObject(kParamBand1_dB),
                                 getParameterObject(kParamBand1_Type),
                                 getParameterObject(kParamBand1_Order));
        controller->addBandParam(getParameterObject(kParamBand2_In),
                                 getParameterObject(kParamBand2_Hz),
                                 getParameterObject(kParamBand2_Q),
                                 getParameterObject(kParamBand2_dB),
                                 getParameterObject(kParamBand2_Type),
                                 getParameterObject(kParamBand2_Order));
        controller->addBandParam(getParameterObject(kParamBand3_In),
                                 getParameterObject(kParamBand3_Hz),
                                 getParameterObject(kParamBand3_Q),
                                 getParameterObject(kParamBand3_dB),
                                 getParameterObject(kParamBand3_Type),
                                 getParameterObject(kParamBand3_Order));
        controller->addBandParam(getParameterObject(kParamBand4_In),
                                 getParameterObject(kParamBand4_Hz),
                                 getParameterObject(kParamBand4_Q),
                                 getParameterObject(kParamBand4_dB),
                                 getParameterObject(kParamBand4_Type),
                                 getParameterObject(kParamBand4_Order));
        controller->addBandParam(getParameterObject(kParamBand5_In),
                                 getParameterObject(kParamBand5_Hz),
                                 getParameterObject(kParamBand5_Q),
                                 getParameterObject(kParamBand5_dB),
                                 getParameterObject(kParamBand5_Type),
                                 getParameterObject(kParamBand5_Order));
        controller->addLevelParam(getParameterObject(kParamLevel));
        controller->addBypassParam(getParameterObject(kParamBypass));
        addEQCurveViewController(controller);
        return controller;
    }
    return nullptr;
};

//------------------------------------------------------------------------
void RFEQ_Controller::editorAttached(Vst::EditorView* editor)
{
    editors.push_back(editor);
}

//------------------------------------------------------------------------
void RFEQ_Controller::editorRemoved(Vst::EditorView* editor)
{
    editors.erase(std::find(editors.begin(), editors.end(), editor));
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

        for (EditorVector::const_iterator it = editors.begin(), end = editors.end(); it != end; ++it)
        {
            VSTGUI::VST3Editor* editor = dynamic_cast<VSTGUI::VST3Editor*>(*it);
            if (editor)
                editor->setZoomFactor(zoomFactors[index].factor);
        }
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
    if (!message)
        return kInvalidArgument;
    
    if (strcmp (message->getMessageID (), "GUI") == 0)
    {
        if (!curveControllers.empty())
        {
            for (auto iter = curveControllers.begin(); iter != curveControllers.end(); iter++)
            {
                if (auto attributes = message->getAttributes ())
                {
                    const void* data;
                    uint32 sizeInBytes;
                    ParamValue getValue = 0.0;
                    
                    if (attributes->getFloat  ("projectSR", getValue) == kResultTrue) projectSR = getValue;
                    if (attributes->getFloat  ("targetSR",  getValue) == kResultTrue) targetSR  = getValue;
                    if (attributes->getBinary ("sample", data, sizeInBytes) == kResultTrue)
                    {
                        auto numSamples = sizeInBytes / sizeof(float);
                        FFT.processBlock((float*)data, numSamples, 0);
                        std::fill(fft_out, fft_out + numBins, 0.0);
                        if (FFT.getData(fft_out))
                            (*iter)->setFFTArray(fft_out, numSamples, projectSR);
                        (*iter)->setEQsampleRate(targetSR);
                    }
                }
            }
        }
        return kResultOk;
    }

    return EditControllerEx1::notify(message);
}

//------------------------------------------------------------------------
} // namespace yg331
