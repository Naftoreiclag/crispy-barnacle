/*
    Copyright 2016 James Fong

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MFCC_HPP
#define MFCC_HPP

#include <stdint.h>

#include "WaveformIO.hpp"

double melScale(double hertz);

double invMelScale(double mels);

struct MFCCParams {
    MFCCParams();
    
    double frameLengthMilliseconds;
    double frameStepMilliseconds;
    double filterMinFreq;
    double filterMaxFreq;
    int32_t numFilterbanks;
    int32_t numMfccs;
};

struct MFCC {
    MFCC(double** mfccs, int32_t windowCount, int32_t mfccCount);
    ~MFCC();
    
    double** data;
    
    int64_t numWindows;
    int32_t numMfccs;
    
    MFCCParams params;
};

double similarityIndex(MFCC* a, int64_t aw, MFCC* b, int64_t bw);

#ifdef NDEBUG
#define debugMostSimilarSamples(...)
#else
void debugMostSimilarSamples(MFCC* templat, MFCC* palette, Waveform waveform);
#endif

MFCC* generateMFCC(
    Waveform inputAudio, 
    const MFCCParams params,
    bool debugOutput = false,
    std::string debugPrefix = "debug");

#endif // MFCC_HPP
