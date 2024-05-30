# Relief-EQ
Relief EQ is a EQ for everyday tasks, anywhere from mixing a track to master bus.  

# How it started  

When I started music, I was mixing live band and some PA stuffs.  
While doing that, I didn't really cared how EQ sounded.  
If it's not broken, and it can boost and cut frequencies, that's all I need!  

After I moved into DAW and digital mixing, there was so many many digital EQs to choose.  
At first, I used SSL Channel strips and it's EQs. It's similar to what I was doing.  
As I did more mixes, I felt a need for a surgical EQ to fix and patch a few spots.  

Since I major in EE, I thought there was no reason for digital EQs to sound different.  
It took a long path to do all researches of EQs and compare them.  

Doing some research, there was many things that makes EQ sound different.  
* Archiecture of filter implementaion
* Bit depth
* Cramping and De-cramping methods
* Oversampling
* Serial and Parallel topology
* Gain - Q dependency
* Q range and UI/UX
* Shelf definition(frequency point, filter order, symmetrical or asymmetrical)

Finally, I made a conclusion in these lines; 
1. If filter runs in 64bit double precision, and not Direct Form, the 'Quality' will be about same.
2. Oversampling does change the 'Sound' of an EQ. 
3. Everything else matters in 'Feeling' of an EQ.

## 1. Quality  

There are many papers and researches made about Archiecture of filters, in aspects of some points;  

* stability
* rounding error causing noise
* transients

Long story short, a 64bit processing SVF architecture is best all-around fit.  
It's stable, robust in math, and handles transients well.  

## 2. How it sounds  

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

## 3. How it behaves    

The 'Feeling' of an EQ was more important than I thought!

Professional mix engineers, they don't turn on Plugin Doctor while mixing.  
Also, they don't rely on EQ curve and frequency display either.  
They look at values what plugins show and listens how it sounds while turning knobs, and thats all.  

Then, shouldn't we focus on 'what values are showing' and 'what happens when we turn knobs?'  

Let's look at Gain - Q dependency first.  
I means turning gain knob affects Q value internally.  
EQ with minimal Gain - Q dependency will remain relativly high Q in low gains, feeling clinical.  
EQ with moderate Gain - Q dependency will have borad Q in low gains, feeling musical.  
They will sound exactly same if we match gain and Q in Plugin Doctor!
However, any engineer using these two EQs will say that they sound different, because the difference in Gain - Q dependency made them use EQ diffently!!  

Another thing that we should look into is the Shelf Definition.  
Many analog EQs and digital EQs with reputaion of nice highs have something in common.  
They use specific definition of shelf filter: 6dB/oct first order filter for shelves, with frequency centered at -3dB point.  
This is commonly missed point in digital EQs.  
They make shelves in 12dB/oct and frequency centered at mid point, because this is how RBJ style shelves are defined.  
So, what happens here?
6dB/oct filters are more shallow and much wider than 12dB/oct.  
Also, at same displyed frequency, -3dB point shelf makes much lower shelf then -6dB or mid point shelf.  
Since the magnitude changes gradually, starts early and ends late, it feels 'transparant' and 'smooth'.  

These choices in filter design and UI will guide a specific way of using an EQ, making the perception of 'they sound different'

# Tech  
It uses SVF State-Space archiecture for implementation of biquad filters.  
It's better than Direct-Form implementaions in stability and error.  

It runs in 64bit double precision processing, and 64bit I/O if host supports it.  
It's to ensure there is no reason to blame sound quality of an EQ, and keep our focus in EQ decisions.  

It oversamples to 96/88kHz at 48/44kHz sample rates, and does not at 96/88kHz and higher.  
It's to keep both frequency and phase response of filters more natural and intuitive to our ears.  

# Ref  

<https://dafx14.fau.de/papers/dafx14_aaron_wishnick_time_varying_filters_for_.pdf>  
<https://dafx2020.mdw.ac.at/proceedings/papers/DAFx2020_paper_52.pdf>  
<https://cytomic.com/files/dsp/SVF-vs-DF1.pdf>  
<https://www.researchgate.net/publication/282326563>  
<https://www.dsprelated.com/freebooks/filters/Implementation_Structures_Recursive_Digital.html>  

<https://forum.juce.com/t/dsp-module-discussion-iir-filter-and-statevariablefilter/23891>  
