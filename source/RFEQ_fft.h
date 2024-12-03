#pragma once

//------------------------------------------------------------------------
namespace yg331 {

    static constexpr int fftOrder = 12;
    static constexpr int fftSize = 1 << fftOrder;      // 4096 samples
    static constexpr int numBins = fftSize / 2 + 1;    // 2049 bins

} // yg331

#include "complex"
#include "array"

// Definitions of useful mathematical constants
//
// Define _USE_MATH_DEFINES before including <math.h> to expose these macro
// definitions for common math constants.  These are placed under an #ifdef
// since these commonly-defined names are not part of the C or C++ standards
#ifndef M_E 
#define M_E        2.71828182845904523536   // e   
#endif  //M_E
#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340736   // log2(e)
#endif  //M_LOG2E
#ifndef M_LOG10E
#define M_LOG10E   0.434294481903251827651  // log10(e)
#endif  //M_LOG10E
#ifndef M_LN2
#define M_LN2      0.693147180559945309417  // ln(2)
#endif  //M_LN2
#ifndef M_LN10
#define M_LN10     2.30258509299404568402   // ln(10)
#endif  //M_LN10
#ifndef M_PI
#define M_PI       3.14159265358979323846   // pi
#endif  //M_PI
#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923   // pi/2
#endif  //M_PI_2
#ifndef M_PI_4
#define M_PI_4     0.785398163397448309616  // pi/4
#endif  //M_PI_4
#ifndef M_1_PI
#define M_1_PI     0.318309886183790671538  // 1/pi
#endif  //M_1_PI
#ifndef M_2_PI
#define M_2_PI     0.636619772367581343076  // 2/pi
#endif  //M_2_PI
#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390   // 2/sqrt(pi)
#endif  //M_2_SQRTPI
#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880   // sqrt(2)
#endif  //M_SQRT2
#ifndef M_SQRT1_2
#define M_SQRT1_2  0.707106781186547524401  // 1/sqrt(2)
#endif  //M_SQRT1_2

// #define PFFFT_SIMD_DISABLE

namespace yg331 {

    /* Copyright (c) 2013  Julien Pommier ( pommier@modartt.com )

       Based on original fortran 77 code from FFTPACKv4 from NETLIB,
       authored by Dr Paul Swarztrauber of NCAR, in 1985.

       As confirmed by the NCAR fftpack software curators, the following
       FFTPACKv5 license applies to FFTPACKv4 sources. My changes are
       released under the same terms.

       FFTPACK license:

       http://www.cisl.ucar.edu/css/software/fftpack5/ftpk.html

       Copyright (c) 2004 the University Corporation for Atmospheric
       Research ("UCAR"). All rights reserved. Developed by NCAR's
       Computational and Information Systems Laboratory, UCAR,
       www.cisl.ucar.edu.

       Redistribution and use of the Software in source and binary forms,
       with or without modification, is permitted provided that the
       following conditions are met:

       - Neither the names of NCAR's Computational and Information Systems
       Laboratory, the University Corporation for Atmospheric Research,
       nor the names of its sponsors or contributors may be used to
       endorse or promote products derived from this Software without
       specific prior written permission.

       - Redistributions of source code must retain the above copyright
       notices, this list of conditions, and the disclaimer below.

       - Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions, and the disclaimer below in the
       documentation and/or other materials provided with the
       distribution.

       THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
       EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO THE WARRANTIES OF
       MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
       NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT
       HOLDERS BE LIABLE FOR ANY CLAIM, INDIRECT, INCIDENTAL, SPECIAL,
       EXEMPLARY, OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY, WHETHER IN AN
       ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
       CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE
       SOFTWARE.
    */

    /*
       PFFFT : a Pretty Fast FFT.

       This is basically an adaptation of the single precision fftpack
       (v4) as found on netlib taking advantage of SIMD instruction found
       on cpus such as intel x86 (SSE1), powerpc (Altivec), and arm (NEON).

       For architectures where no SIMD instruction is available, the code
       falls back to a scalar version.

       Restrictions:

       - 1D transforms only, with 32-bit single precision.

       - supports only transforms for inputs of length N of the form
       N=(2^a)*(3^b)*(5^c), a >= 5, b >=0, c >= 0 (32, 48, 64, 96, 128,
       144, 160, etc are all acceptable lengths). Performance is best for
       128<=N<=8192.

       - all (float*) pointers in the functions below are expected to
       have an "simd-compatible" alignment, that is 16 bytes on x86 and
       powerpc CPUs.

       You can allocate such buffers with the functions
       pffft_aligned_malloc / pffft_aligned_free (or with stuff like
       posix_memalign..)

    */

#ifndef PFFFT_H
#define PFFFT_H

#include <stddef.h> // for size_t

#ifdef __cplusplus
    extern "C" {
#endif

        /* opaque struct holding internal stuff (precomputed twiddle factors)
           this struct can be shared by many threads as it contains only
           read-only data.
        */
        typedef struct PFFFT_Setup PFFFT_Setup;

        /* direction of the transform */
        typedef enum { PFFFT_FORWARD, PFFFT_BACKWARD } pffft_direction_t;

        /* type of transform */
        typedef enum { PFFFT_REAL, PFFFT_COMPLEX } pffft_transform_t;

        /*
          prepare for performing transforms of size N -- the returned
          PFFFT_Setup structure is read-only so it can safely be shared by
          multiple concurrent threads.
        */
        PFFFT_Setup* pffft_new_setup(int N, pffft_transform_t transform);
        void pffft_destroy_setup(PFFFT_Setup*);

        /*
           Perform a Fourier transform , The z-domain data is stored in the
           most efficient order for transforming it back, or using it for
           convolution. If you need to have its content sorted in the
           "usual" way, that is as an array of interleaved complex numbers,
           either use pffft_transform_ordered , or call pffft_zreorder after
           the forward fft, and before the backward fft.

           Transforms are not scaled: PFFFT_BACKWARD(PFFFT_FORWARD(x)) = N*x.
           Typically you will want to scale the backward transform by 1/N.

           The 'work' pointer should point to an area of N (2*N for complex
           fft) floats, properly aligned. If 'work' is NULL, then stack will
           be used instead (this is probably the best strategy for small
           FFTs, say for N < 16384).

           input and output may alias.
        */
        void pffft_transform(PFFFT_Setup* setup, const float* input, float* output, float* work, pffft_direction_t direction);

        /*
           Similar to pffft_transform, but makes sure that the output is
           ordered as expected (interleaved complex numbers).  This is
           similar to calling pffft_transform and then pffft_zreorder.

           input and output may alias.
        */
        void pffft_transform_ordered(PFFFT_Setup* setup, const float* input, float* output, float* work, pffft_direction_t direction);

        /*
           call pffft_zreorder(.., PFFFT_FORWARD) after pffft_transform(...,
           PFFFT_FORWARD) if you want to have the frequency components in
           the correct "canonical" order, as interleaved complex numbers.

           (for real transforms, both 0-frequency and half frequency
           components, which are real, are assembled in the first entry as
           F(0)+i*F(n/2+1). Note that the original fftpack did place
           F(n/2+1) at the end of the arrays).

           input and output should not alias.
        */
        void pffft_zreorder(PFFFT_Setup* setup, const float* input, float* output, pffft_direction_t direction);

        /*
           Perform a multiplication of the frequency components of dft_a and
           dft_b and accumulate them into dft_ab. The arrays should have
           been obtained with pffft_transform(.., PFFFT_FORWARD) and should
           *not* have been reordered with pffft_zreorder (otherwise just
           perform the operation yourself as the dft coefs are stored as
           interleaved complex numbers).

           the operation performed is: dft_ab += (dft_a * fdt_b)*scaling

           The dft_a, dft_b and dft_ab pointers may alias.
        */
        void pffft_zconvolve_accumulate(PFFFT_Setup* setup, const float* dft_a, const float* dft_b, float* dft_ab, float scaling);

        /*
          the float buffers must have the correct alignment (16-byte boundary
          on intel and powerpc). This function may be used to obtain such
          correctly aligned buffers.
        */
        void* pffft_aligned_malloc(size_t nb_bytes);
        void  pffft_aligned_free(void*);

        /* return 4 or 1 wether support SSE/Altivec instructions was enable when building pffft.c */
        int pffft_simd_size(void);

#ifndef PFFFT_SIMD_DISABLE
        void validate_pffft_simd(); // a small function inside pffft.c that will detect compiler bugs with respect to simd instruction 
#endif

#ifdef __cplusplus
    }
#endif

#endif // PFFFT_H



    // https://www.kvraudio.com/forum/viewtopic.php?p=8726913#p8726913
    /*
    ==============================================================================

    PFFT.h
    Created: 29 Aug 2022 1:08:12pm
    Author:  Justin Johnson

    ==============================================================================

    MIT License

    Copyright (c) 2021 Justin Johnson

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    */

    /**
     * C++ Wrapper for pffft, a reasonably fast FFT library.
     *  The class here reflects closely the Juce FFT class and is a drop
     *  in replacement.
     *  See: https://bitbucket.org/jpommier/pffft/src/master/
     */
    class PFFFT
    {
    public:
        PFFFT(int order)
        {
            size_ = 1 << order;
            scale_ = 1.f / size_;
            pSetup_ = pffft_new_setup(size_, PFFFT_REAL);
        }

        ~PFFFT()
        {
            pffft_destroy_setup(pSetup_);
        }

        void performRealOnlyForwardTransform(float* pBuffer, bool onlyCalculateNonNegativeFrequencies = false)
        {
            pffft_transform_ordered(pSetup_, pBuffer, pBuffer, NULL, PFFFT_FORWARD);
        }

        void performRealOnlyInverseTransform(float* pBuffer)
        {
            pffft_transform_ordered(pSetup_, pBuffer, pBuffer, NULL, PFFFT_BACKWARD);

            for (int i = 0; i < size_; ++i)
            {
                pBuffer[i] *= scale_;
            }
        }

        void performFrequencyOnlyForwardTransform(float* inputOutputData, bool ignoreNegativeFreqs = false) const noexcept
        {
            if (size_ == 1) return;

            pffft_transform_ordered(pSetup_, inputOutputData, inputOutputData, NULL, PFFFT_FORWARD);

            auto* out = reinterpret_cast<std::complex<float> *>(inputOutputData);

            const auto limit = ignoreNegativeFreqs ? (size_ / 2) + 1 : size_;

            for (int i = 0; i < limit; ++i)
            {
                inputOutputData[i] = std::abs(out[i]);
            }

            std::fill(inputOutputData + limit, inputOutputData + size_ * 2, 0.f);
        }

        [[nodiscard]] int getSize() const noexcept { return size_; }

    private:
        int size_;
        float scale_;

        PFFFT_Setup* pSetup_;
    };





    /*
    
    MIT License

    Copyright (c) 2023 Matthijs Hollemans

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    
    
    
    
    */
    /**
    STFT analysis and resynthesis of audio data.

    Each channel should have its own FFTProcessor.
    */
    class FFTProcessor
    {
    public:
        FFTProcessor();

        int getLatencyInSamples() const { return fftSize; }

        void reset();
        float processSample(float sample, bool bypassed);
        void processBlock(float* data, int numSamples, bool bypassed);
        void processBlock(double* data, int numSamples, bool bypassed);
        bool getData(float* out) {
            if (!data_avail)
                return false;

            auto* cdata = reinterpret_cast<std::complex<float>*>(&fftData);
            for (int i = 0; i < numBins; ++i) {
                float magnitude = std::abs(cdata[i]);
                out[i] = magnitude;
            }
            data_avail = 0;
            return true;
        };
        static void hannWindow(float* window, int length)
        {
            float delta = (M_PI * 2) / float(length);
            float phase = 0.0f;
            for (int i = 0; i < length; ++i) {
                window[i] = 0.5f * (1.0f - std::cos(phase));
                phase += delta;
                // phase = i * delta
            }
            float pwr = 0.0;
            for (int i = 0; i < length; ++i) {
                pwr += window[i];
            }
            pwr = 1.0 / pwr;
            for (int i = 0; i < length; ++i) {
                window[i] *= pwr;
            }
        }
        static void bkhsWindow(float* window, int length)
        {
            double dwindowpos = (M_PI * 2) / length;
            double pwr = 0.0;
            for (int i = 0; i < 0.5 * length + 1; i++) {
                double windowpos = i * dwindowpos;
                window[i] = (0.35875 - 0.48829 * cos(windowpos) + 0.14128 * cos(2.0 * windowpos) - 0.01168 * cos(3.0 * windowpos));
                pwr += window[i];
            }
            pwr = 0.5 / (pwr * 2 - window[(int)(0.5 * length)]);
            for (int i = 0; i < 0.5 * length + 1; i++) {
                window[i] *= pwr;
            }
            for (int i = 0; i < 0.5 * length; i++) {
                window[length - i - 1] = window[i];
            }
        }
        static inline float Ino(float x)
        {
            double d = 0, ds = 1, s = 1;
            do
            {
                d += 2;
                ds *= x * x / (d * d);
                s += ds;
            } while (ds > s * 1e-6);
            return s;
        };

        static void ksblWindow(float* window, int length)
        {
            int Np = (length - 1) / 2;
            float Alpha;
            float Inoalpha;

            Alpha = 3.0 * M_PI;

            Inoalpha = Ino(Alpha);

            for (int j = 0; j <= Np; j++)
            {
                window[Np + j] = Ino(Alpha * std::sqrt(1.0 - ((float)(j * j) / (float)(Np * Np)))) / Inoalpha;
            }
            for (int j = 0; j < Np; j++)
            {
                window[j] = window[length - 1 - j];
            }

            float pwr = 0.0;
            for (int i = 0; i < length; i++) {
                pwr += window[i];
            }
            pwr = 1.0 / pwr;
            for (int i = 0; i < length; i++) {
                window[i] *= pwr * 1.86; // Normalization & Amplitube correction
            }
        };

    private:
        void processFrame(bool bypassed);
        void processSpectrum(float* data, int numBins);

        // The FFT has 2^order points and fftSize/2 + 1 bins.
        static constexpr int fftOrder = yg331::fftOrder;
        static constexpr int fftSize = 1 << fftOrder;      // 4096 samples
        static constexpr int numBins = fftSize / 2 + 1;    // 2049 bins
        static constexpr int overlap = 4;                  // 75% overlap
        static constexpr int hopSize = fftSize / overlap;  // //256 samples

        // Gain correction for using Hann window with 75% overlap.
        static constexpr float windowCorrection = 2.0f / 3.0f;

        PFFFT fft;
        std::array<float, fftSize> window = { 0.0, };

        // Counts up until the next hop.
        int count = 0;

        int data_avail = 0;

        // Write position in input FIFO and read position in output FIFO.
        int pos = 0;

        // Circular buffers for incoming and outgoing audio data.
        /* SSE and co like 16-bytes aligned pointers */
        alignas(16) std::array<float, fftSize>  inputFifo;
        alignas(16) std::array<float, fftSize> outputFifo;

        // The FFT working space. Contains interleaved complex numbers.
        alignas(16) std::array<float, fftSize * 2>  fftData;
    };
}
