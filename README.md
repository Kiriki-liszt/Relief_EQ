# Relief-EQ

Relief EQ is a EQ for everyday tasks, anywhere from mixing a track to master bus.  

Runs in double precision 64-bit internal processing. Also double precision input / output if supported.  
At 44.1 kHz and 48 kHz sampling rates, it upsamples to 88.2 kHz or 96 kHz respectively to do all EQ processing with 17 sample latency.  

Windows and Mac, VST3 and AU.  

[![GitHub Release](https://img.shields.io/github/v/release/kiriki-liszt/Relief_EQ?style=flat-square&label=Get%20latest%20Release)](https://github.com/Kiriki-liszt/Relief_EQ/releases/latest)
[![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/kiriki-liszt/Relief_EQ/total?style=flat-square&label=total%20downloads&color=blue)](https://tooomm.github.io/github-release-stats/?username=Kiriki-liszt&repository=Relief_EQ)  

[![Static Badge](https://img.shields.io/badge/coffee%20maybe%3F%20%3D%5D%20-gray?style=for-the-badge&logo=buy-me-a-coffee)](https://buymeacoffee.com/kirikiaris)  

<img src="https://github.com/Kiriki-liszt/Relief_EQ/blob/main/screenshot.png?raw=true"  width="600"/>  

## Compatibility  

VST3, AUv2  

## System Requirements

Audio Units  

* Mac OS X 10.13 or later (Intel or Apple Silicon Native)

VST3  

* Mac OS X 10.13 or later (Intel or Apple Silicon Native)
* Windows 10 or later

## How to use  

1. Windows

Unzip Win.zip from latest release and copy to "C:\Program Files\Common Files\VST3".  

2. MacOS(Intel tested, Apple Silicon not tested).  

Unzip MacOS.zip from latest release and copy vst3 to "/Library/Audio/Plug-Ins/VST3" and component to "/Library/Audio/Plug-Ins/Components".  

> If it doesn't go well, configure security options in console as  
>
> ``` console  
> sudo xattr -r -d com.apple.quarantine /Library/Audio/Plug-Ins/VST3/Relief_EQ.vst3  
> sudo xattr -r -d com.apple.quarantine /Library/Audio/Plug-Ins/Components/Relief_EQ.component
>
> sudo codesign --force --sign - /Library/Audio/Plug-Ins/VST3/Relief_EQ.vst3  
> sudo codesign --force --sign - /Library/Audio/Plug-Ins/Components/Relief_EQ.component
> ```  
>
> tested by @jonasborneland [here](https://github.com/Kiriki-liszt/JS_Inflator_to_VST2_VST3/issues/12#issuecomment-1616671177)

## Licensing  

Relief EQ is using GPL v3 license.  

### VST  

> Q: I would like to share the source code of my VST 3 plug-in/host on GitHub or other such platform.  
>
> - You can choose the GPLv3 license and feel free to share your plug-ins/host's source code including or referencing the VST 3 SDK's sources on GitHub.  
> - **You are allowed to provide a binary form of your plug-ins/host too, provided that you provide its source code as GPLv3 too.**  
> - Note that you have to follow the Steinberg VST usage guidelines.  
>  
> <https://steinbergmedia.github.io/vst3_dev_portal/pages/FAQ/Licensing.html>  

![VST Logo](https://github.com/Kiriki-liszt/Sky_Blue_EQ4/assets/107096260/142e3c12-cd5f-415d-9b72-8b4f04419633)

VSTSDK 3.7.9 used  
VSTGUI 4.12 used  

### PFFFT  

PFFFT is work of Julien Pommier, with Copyright (c) 2013  Julien Pommier ( <pommier@modartt.com> )  
It is under BSD-Like, and license is in the header at RFEQ_fft.h.  

### C++ Wrapper for pffft  

C++ Wrapper for pffft if work of Justin Johnson, obtained at [https://www.kvraudio.com/forum/viewtopic.php?p=8726913#p8726913](https://www.kvraudio.com/forum/viewtopic.php?p=8726913#p8726913)  
As other works of him are in MIT lisence, I assumed same.  
The lisence is in the header at RFEQ_fft.h.  

### FFT handler with FIFO managment  

It is work of Matthijs Hollemans, obtained at [https://github.com/hollance/fft-juce](https://github.com/hollance/fft-juce)  
It is under MIT lisence, and you can find the license in the header at RFEQ_fft.h.  

## Project Build  

Use CMake to build itself or make IDE project file.  
Check .github/workflows/Windows Build.yml.  
Remember to git clone VSTSDK, too.  
Supports Windows, Mac, Linux(same as VSTSDK).  

## Version logs

v1.0.0.b : intial try.  

v1.0.0 : Official release.  Changed shelf arch.  

v1.1.0 : FFT analyzer added.  

v1.1.1 : Crash after closing the plugin window fixed.  

v1.1.2 :

- Processor - data exchange check at setActive(Ableton crash).  
- Contoller - save state of Contoller for correct recalling(AUv2, VST with FL studio).

v1.1.3 : FFT graph more consistent over block size change.  

v1.1.4 : Bug in Apple Silicon Native build fixed.  

v1.1.5 : FFT analyzer tilt 4.5dB/oct added.  

v1.1.6 : Slight refactor, more accurate coef for constant gain-Q dependency, GUI default reset value corrected and cosmetic change to match Relief Compressor. Do not use anymore, as it turns out to be faulty.  

v1.1.6.1 : Hotfix for save-reload process.  

## How it started

When I started music, I was mixing live band and some PA stuffs.  
While doing that, I didn't really cared how EQ sounded.  
If it's not broken, and it can boost and cut frequencies, that's all I need!  

After I moved into DAW and digital mixing, there was so many many digital EQs to choose.  
At first, I used SSL Channel strips and it's EQs. It's similar to what I was doing.  
As I did more mixes, I felt a need for a surgical EQ to fix and patch a few spots.  

Since I major in EE, I thought there was no reason for digital EQs to sound different.  
It took a long path to do all researches of EQs and compare them.  

Doing some research, there was many things that makes EQ sound different.  

- Archiecture of filter implementaion
- Bit depth
- Cramping and De-cramping methods
- Oversampling
- Serial and Parallel topology
- Gain - Q dependency
- Q range and UI/UX
- Shelf definition(frequency point, filter order, symmetrical or asymmetrical)

Finally, I made a conclusion in these lines;  

1. If filter runs in 64bit double precision, and not Direct Form, the 'Quality' will be about same.  
2. Oversampling does change the 'Sound' of an EQ.  
3. Everything else matters in 'Feeling' of an EQ.

### 1. Quality  

There are many papers and researches made about Archiecture of filters, in aspects of some points;  

- stability
- rounding error causing noise
- transients

Long story short, a 64bit processing SVF architecture is best all-around fit.  
It's stable, robust in math, and handles transients well.  

### 2. How it sounds  

It starts from cramping.  
Nyquist Theorem lets every DSP possible, however, it cramps infinite frequency into limited range.  
It causes frequency and phase response of an EQ to be mismatched from analog, or theoretically perfect response.  
In result, high bells will be narrower and looese symmetry, and also high shlef will loose it's shape.  
So cramped EQ has it's own 'sound' to it, especially in high frequency and it's different then 'quality' of EQ.  

If this cramping is a issue for one, we can de-cramp it's frequency curve by adding shelf or adjusting Q and so on.  
It will make frequency response of a EQ more matched to a analog EQ.  
Compared to cramping EQ, decramped EQ will have more high frequency with same gain and Q, it 'sounds' more bright.  

However, it causes the mismatch between phase and frequency response.  
Since phase diffrence makes magnitude difference, I belive phase response matching frequency response affects how 'natural' a EQ sounds.  
For example, a bell curve is made from a zero crossing phase response at it's peak frequency, not the other way.  
In this manner, de-cramped EQ might feel more unnatural than cramping EQ.  

So, cramping EQ will sound less airy and less analog than decramped one, but matches frequency and phase response.  
decramped EQ will sound more like an analog EQ frequency-wise, but relation between frequency and phase response is lost.  

Now, the oversampling solves both issues.  
It decrampes both frequency and phase response, in cost of CPU and latency.  
Moreover, with clean digital EQs without nonlinearity, oversampling half band filters does not have to be strong.  
10 - 20 sample latency for x2 oversampling, 20 - 30 sample latency for x4 oversampling is enough for EQ.  
Oversampled EQ will sound analog due to it's matched curve in high frequency, and naturalness of EQ is also preserved as frequency-phase relation is not touched!  

### 3. How it behaves  

The 'Feeling' of an EQ was more important than I thought!

Professional mix engineers, they don't turn on Plugin Doctor while mixing.  
Also, they don't rely on EQ curve and frequency display either.  
They look at values what plugins show and listens how it sounds while turning knobs, and thats all.  

Then, shouldn't we focus on 'what values are showing' and 'what happens when we turn knobs?'  

Let's look at Gain - Q dependency first.  
It means turning gain knob affects Q value internally.  
EQ with minimal Gain - Q dependency will remain relativly high Q in low gains, feeling clinical.  
EQ with moderate Gain - Q dependency will have broad Q in low gains, feeling musical.  
They will sound exactly same if we match gain and Q in Plugin Doctor!
However, any engineer using these two EQs will say that they sound different, because the difference in Gain - Q dependency made them use EQ diffently!!  

Another thing that we should look into is the Shelf Definition.  
Many analog EQs and digital EQs with reputaion of nice highs have something in common.  
They use specific definition of shelf filter: 6dB/oct first order filter for shelves, with frequency centered at -3dB point.  
This is commonly missed point in digital EQs.  
They make shelves in 12dB/oct and frequency centered at mid point, because this is how RBJ style shelves are defined.  
So, what happens here?
6dB/oct filters are more shallow and much wider than 12dB/oct.  
Also, at same displayed frequency, -3dB point shelf makes much lower shelf then -6dB or mid point shelf.  
Since the magnitude changes gradually, starts early and ends late, it feels 'transparant' and 'smooth'.  

These choices in filter design and UI will guide a specific way of using an EQ, making the perception of 'they sound different'

## FFT  

The FFT in this plugin uses Kaiser Bessel windowing, RMS in frequency magnitudes, and no bands-per-oct smoothing.  
Also, it runs in 4096 sample size.  

## What I've learned  

### AUv2 and FL studio  

They call Controller state after Processer, overwriting states to default.  
Always save Contoller state too, but not bypass state.  

### Data Exchange  

Ableton calls setActive before connect/disconnect happens.  
This causes dataExchange pointer to be null while active process, so check for that in setActive.  

Also, Ableton passes Silence tag quite strongly, so call exchange data after silence is asserted(if you use data exchange for parameter to draw UI or STH).  

## TODO  

- [ ] Dark mode  

## Further lookings  

[https://gearspace.com/board/showpost.php?p=15864586&postcount=730](https://gearspace.com/board/showpost.php?p=15864586&postcount=730)  
[Vicanek, Martin. Matched Second Order Digital Filters. (2016).](https://www.vicanek.de/articles/BiquadFits.pdf)  
[John Flynn & Joshua D. Reiss (2018). Improving the frequency response magnitude and phase of analogue-matched digital filters](https://www.eecs.qmul.ac.uk/~josh/documents/2018/19412.pdf)  
[D. W. Gunness, O. S. Chauhan, “Optimizing the Magnitude Response of Matched z-Transform Filters (“MZTi”) for Loudspeaker Equalization”](https://www.khabdha.org/wp-content/uploads/2008/03/optimizing-the-magnitude-response-of-mzt-filters-mzti-2007.pdf)  

## Ref  

<https://dafx14.fau.de/papers/dafx14_aaron_wishnick_time_varying_filters_for_.pdf>  
<https://dafx2020.mdw.ac.at/proceedings/papers/DAFx2020_paper_52.pdf>  
<https://cytomic.com/files/dsp/SVF-vs-DF1.pdf>  
<https://www.researchgate.net/publication/282326563>  
<https://www.dsprelated.com/freebooks/filters/Implementation_Structures_Recursive_Digital.html>  
<https://forum.juce.com/t/dsp-module-discussion-iir-filter-and-statevariablefilter/23891>  
