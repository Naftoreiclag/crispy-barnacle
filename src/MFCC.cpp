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

#include "MFCC.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

#include "fftw3.h"
#ifndef NDEBUG
#include "stb_image_write.h"
#endif // !NDEBUG

#include "ComplexNumber.hpp"
#include "DebugOutput.hpp"

double melScale(double hertz) {
    return 1127 * std::log(hertz / 700 + 1);
}

double invMelScale(double mels) {
    return 700 * (std::exp(mels / 1127) - 1);
}

MFCCParams::MFCCParams()
: frameLengthMilliseconds(25)
, frameStepMilliseconds(10)
, filterMinFreq(300)
, filterMaxFreq(8000)
, numFilterbanks(26)
, numMfccs(12) {
}

MFCC::MFCC(double** mfccs, int32_t windowCount, int32_t mfccCount)
: data(mfccs)
, numWindows(windowCount)
, numMfccs(mfccCount) {
}

MFCC::~MFCC() {
    for(int64_t i = 0; i < numWindows; ++ i) {
        delete[] data[i];
    }
    delete[] data;
}

double similarityIndex(MFCC* a, int64_t aw, MFCC* b, int64_t bw) {
    assert(a->numMfccs == b->numMfccs && "Different MFCC counts!");
    
    double sum = 0;
    for(int32_t i = 1; i < a->numMfccs; ++ i) {
        double distSq = (b->data[bw][i] - a->data[aw][i]);
        distSq *= distSq;
        
        double weighted = distSq * std::pow(0.8, i - 1);
        sum += weighted;
    }
    
    return sqrt(sum);
}

void debugMostSimilarSamples(MFCC* templat, MFCC* palette, Waveform paletteAudio) {
    assert(palette->numMfccs > 0);
    int32_t* closestMatches = new int32_t[templat->numWindows];
    for(int64_t windowIndex = 0; windowIndex < templat->numWindows; ++ windowIndex) {
        
        double bestSim = similarityIndex(palette, 0, templat, windowIndex);
        int32_t bestPaletteIndex = 0;
        
        for(int32_t i = 1; i < palette->numWindows; ++ i) {
            double thisSim = similarityIndex(palette, i, templat, windowIndex);
            if(thisSim < bestSim) {
                bestSim = thisSim;
                bestPaletteIndex = i;
            }
        }
        
        closestMatches[windowIndex] = bestPaletteIndex;
    }
    
    {
        assert(palette->params.frameStepMilliseconds < palette->params.frameLengthMilliseconds);
        
        
        Waveform mixture;
        
        int32_t windowLength = (palette->params.frameLengthMilliseconds * paletteAudio.mSampleRate) / 1000;
        int32_t windowStep = (palette->params.frameStepMilliseconds * paletteAudio.mSampleRate) / 1000;
        
        mixture.mNumSamples = templat->numWindows * windowStep + windowLength;
        mixture.mFloatSamples = new double[mixture.mNumSamples];
        for(int64_t i = 0; i < mixture.mNumSamples; ++ i) {
            mixture.mFloatSamples[i] = 0;
        }
        mixture.mSampleRate = paletteAudio.mSampleRate;
        
        for(int64_t windowIndex = 0; windowIndex < templat->numWindows; ++ windowIndex) {
            for(int64_t windowSample = 0; windowSample < windowLength; ++ windowSample) {
                // Apply hanning window
                
                // Hooray for compiler optimizations
                double tau = 6.28318530717958647692528677;
                double numerator = tau * windowSample;
                double denominator = windowLength - 1;
                double hanning = 0.5 * (1.0 - std::cos(numerator / denominator));
                mixture.mFloatSamples[windowIndex * windowStep + windowSample] += 
                    hanning * paletteAudio.mFloatSamples[closestMatches[windowIndex] * windowStep + windowSample];
            }
        }
        
        writeWaveform("closest_matches.ogg", mixture);
    }
    
    {
        int width = templat->numWindows;
        int height = templat->numMfccs;
        
        if(width > 2048) width = 2048;
        char imageData[width * height * 3];
        
        for(int pixelY = 0; pixelY < height; ++ pixelY) {
            for(int pixelX = 0; pixelX < width; ++ pixelX) {
                
                int frame = pixelX;
                int spectrum = (height - pixelY) - 1;
                
                {
                    double power = normalized(palette->data[closestMatches[frame]][spectrum], -4, 18);
                    
                    RGB heat = colorrampSevenHeat(power);
                    
                    imageData[(pixelY * width + pixelX) * 3    ] = heat.RU8();
                    imageData[(pixelY * width + pixelX) * 3 + 1] = heat.GU8();
                    imageData[(pixelY * width + pixelX) * 3 + 2] = heat.BU8();
                }
            }
        }
        stbi_write_png("closest_matches.png", width, height, 3, imageData, width * 3);
    }
}

MFCC* generateMFCC(
    Waveform inputAudio, 
    const MFCCParams params,
    bool debugOutput,
    std::string debugPrefix) {
    
    // Yay compiler optimizations
    const double& frameLengthMilliseconds = params.frameLengthMilliseconds;
    const double& frameStepMilliseconds = params.frameStepMilliseconds;
    const double& filterMinFreq = params.filterMinFreq;
    const double& filterMaxFreq = params.filterMaxFreq;
    const int32_t& numFilterbanks = params.numFilterbanks;
    const int32_t& numMfccs = params.numMfccs;
    
    // (Rounded toward zero)
    int32_t windowLength = (frameLengthMilliseconds * inputAudio.mSampleRate) / 1000;
    int32_t spectrumLength = windowLength / 2;
    int32_t windowStep = (frameStepMilliseconds * inputAudio.mSampleRate) / 1000;
    
    std::cout << "\tWindow function: Hanning" << std::endl;
    std::cout << "\tWindow length: " << windowLength << " samples / " << frameLengthMilliseconds << "ms" << std::endl;
    if(windowLength < 0) {
        std::cerr << "Fatal error! Length cannot be negative!" << std::endl;
        return NULL;
    }
    
    std::cout << "\tWindow step: " << windowStep << " samples / " << frameStepMilliseconds << "ms" << std::endl;
    if(windowStep < 1) {
        std::cerr << "Fatal error! Step must be greater than 0!" << std::endl;
        return NULL;
    }
    
    int64_t numWindows = 0;
    // TODO: use constant time calculation
    for(int64_t windowIndex = 0; (windowIndex * windowStep) < inputAudio.mNumSamples; ++ windowIndex) {
        numWindows ++;
    }
    std::cout << "\tWindow count: " << numWindows << std::endl;
    
    std::cout << "\tFilterbank energy count: " << numFilterbanks << std::endl;
    std::cout << "\tMel frequency count:" << numMfccs << std::endl;
    if(numMfccs > numFilterbanks) {
        std::cerr << "Fatal error! Mel frequency count is greater than the number of filterbank energies!" << std::endl;
        return NULL;
    }
    
    std::cout << "\tLower filterbank range: " << filterMinFreq << "hz" << std::endl;
    std::cout << "\tUpper filterbank range: " << filterMaxFreq << "hz" << std::endl;
    if(filterMinFreq > filterMaxFreq) {
        std::cerr << "Fatal error! Invalid range!" << std::endl;
        return NULL;
    }
    
    // The frequency represented by fftwCompleteOutput[i][n] in hertz is equal to:
    // freq = n * inputAudio.mSampleRate / windowLength
    // Inverse is approximated by:
    // bin = floor((freq * (windowLength + 1)) / inputAudio.mSampleRate)

    // FFTW transform
    ComplexNumber** fftwCompleteOutput = new ComplexNumber*[numWindows];
    for(int64_t i = 0; i < numWindows; ++ i) {
        fftwCompleteOutput[i] = new ComplexNumber[spectrumLength];
    }
    {
        fftw_complex* fftwInput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * windowLength);
        fftw_complex* fftwOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * windowLength);
        
        std::cout << "\tOptimizing FFTW... ";
        fftw_plan fftwPlan = fftw_plan_dft_1d(windowLength, fftwInput, fftwOutput, FFTW_FORWARD, FFTW_MEASURE);
        std::cout << "\tdone" << std::endl;
        
        // Set imaginary components to be zero
        for(int64_t windowSample = 0; windowSample < windowLength; ++ windowSample) {
            fftwInput[windowSample][1] = 0;
        }
        
        std::cout << "\tPerforming DFT... ";
        for(int64_t windowIndex = 0; windowIndex < numWindows; ++ windowIndex) {
            
            int64_t beginningSample = windowIndex * windowStep;
            
            for(int64_t windowSample = 0; windowSample < windowLength; ++ windowSample) {
                int64_t currentSampleIndex = beginningSample + windowSample;
                
                double sample;
                if(currentSampleIndex >= inputAudio.mNumSamples) {
                    sample = 0;
                } else {
                    sample = inputAudio.mFloatSamples[currentSampleIndex];
                    
                    // Apply hanning window
                    
                    // Hooray for compiler optimizations
                    double tau = 6.28318530717958647692528677;
                    double numerator = tau * windowSample;
                    double denominator = windowLength - 1;
                    double hanning = 0.5 * (1.0 - std::cos(numerator / denominator));
                    
                    sample *= hanning;
                }
                
                fftwInput[windowSample][0] = sample;
            }
            
            fftw_execute(fftwPlan);
            
            for(int64_t windowSample = 0; windowSample < spectrumLength; ++ windowSample) {
                
                fftwCompleteOutput[windowIndex][windowSample].real = fftwOutput[windowSample][0];
                fftwCompleteOutput[windowIndex][windowSample].imag = fftwOutput[windowSample][1];
            }
        }
        if(debugOutput) writeFFTWOutputDebug(debugPrefix + "_fftw.png", fftwCompleteOutput, numWindows, spectrumLength);
        std::cout << "\tdone" << std::endl;
        
        fftw_destroy_plan(fftwPlan);
        fftw_free(fftwOutput);
        fftw_free(fftwInput);
    }
    
    // Power estmates
    double** powerEstimates = new double*[numWindows];
    for(int64_t i = 0; i < numWindows; ++ i) {
        powerEstimates[i] = new double[spectrumLength];
    }
    {
        std::cout << "\tComputing power estimates... ";
        for(int64_t windowIndex = 0; windowIndex < numWindows; ++ windowIndex) {
            for(int64_t windowSample = 0; windowSample < spectrumLength; ++ windowSample) {
                
                double real = fftwCompleteOutput[windowIndex][windowSample].real;
                double imaginary = fftwCompleteOutput[windowIndex][windowSample].imag;
                double absValSq = real * real + imaginary * imaginary;
                double denom = windowLength; // NOT spectrumLength!
                
                
                powerEstimates[windowIndex][windowSample] = absValSq / denom;
            }
        }
        if(debugOutput) writeGenericHeatOutput(debugPrefix + "_power.png", powerEstimates, numWindows, spectrumLength);
        std::cout << "\tdone" << std::endl;
    }
    
    for(int64_t i = 0; i < numWindows; ++ i) {
        delete[] fftwCompleteOutput[i];
    }
    delete[] fftwCompleteOutput;
    
    // Mel filterbank
    double* filterbank = new double[numFilterbanks + 2];
    {
        double filterMaxFreqMels = melScale(filterMaxFreq);
        double filterMinFreqMels = melScale(filterMinFreq);
        
        std::cout << "\tMel filterbank (hz): ";
        double melsStep = (filterMaxFreqMels - filterMinFreqMels) / (numFilterbanks + 1);
        for(int32_t i = 0; i < numFilterbanks + 2; ++ i) {
            double mels = melsStep * i + filterMinFreqMels;
            double hertz = invMelScale(mels);
            
            filterbank[i] = hertz;
            
            std::cout << hertz;
            if(i != numFilterbanks + 1) {
                std::cout << ", ";
            }
        }
        std::cout << std::endl;
    }
    
    //
    double** loggedFilterbankEnergies = new double*[numWindows];
    for(int64_t i = 0; i < numWindows; ++ i) {
        loggedFilterbankEnergies[i] = new double[numFilterbanks];
    }
    {
        #ifndef NDEBUG
        double** filterbankEnergies = new double*[numWindows];
        for(int64_t i = 0; i < numWindows; ++ i) {
            filterbankEnergies[i] = new double[numFilterbanks];
        }
        #endif // !NDEBUG
        std::cout << "\tComputing filterbank energies... ";
        for(int64_t windowIndex = 0; windowIndex < numWindows; ++ windowIndex) {
            for(int32_t filterbankIndex = 0; filterbankIndex < numFilterbanks; ++ filterbankIndex) {
                double filterBegin = filterbank[filterbankIndex];
                double filterMiddle = filterbank[filterbankIndex + 1];
                double filterEnd = filterbank[filterbankIndex + 2];
                
                
                double totalEnergy = 0;
                for(int64_t windowSample = 0; windowSample < spectrumLength; ++ windowSample) {
                    
                    double sampleFrequency = ((double) (inputAudio.mSampleRate * windowSample)) / ((double) windowLength);
                    
                    if(sampleFrequency < filterBegin) continue;
                    if(sampleFrequency > filterEnd) break;
                    
                    double filterY = filterMiddle - sampleFrequency;
                    if(filterY < 0) {
                        filterY /= filterMiddle - filterEnd;
                    } else {
                        filterY /= filterMiddle - filterBegin;
                    }
                    
                    totalEnergy += powerEstimates[windowIndex][windowSample] * filterY;
                }
                #ifndef NDEBUG
                filterbankEnergies[windowIndex][filterbankIndex] = totalEnergy;
                #endif // !NDEBUG
                loggedFilterbankEnergies[windowIndex][filterbankIndex] = std::log(totalEnergy); // Natural log, please
            }
        }
        if(debugOutput) writeGenericHeatOutput(debugPrefix + "_energies_log.png", loggedFilterbankEnergies, numWindows, numFilterbanks, -10, 1);
        #ifndef NDEBUG
        if(debugOutput) writeGenericHeatOutput(debugPrefix + "_energies.png", filterbankEnergies, numWindows, numFilterbanks);
        for(int64_t i = 0; i < numWindows; ++ i) {
            delete[] filterbankEnergies[i];
        }
        delete[] filterbankEnergies;
        #endif // !NDEBUG
        
        std::cout << "\tdone" << std::endl;
        
    }
    for(int64_t i = 0; i < numWindows; ++ i) {
        delete[] powerEstimates[i];
    }
    delete[] powerEstimates;
    delete[] filterbank;
    
    // MFCC (Discrete cosine transform and tossing out high-frequency data)
    double** mfccs = new double*[numWindows];
    for(int64_t i = 0; i < numWindows; ++ i) {
        mfccs[i] = new double[numMfccs];
    }
    {
        
        std::cout << "\tComputing mel frequency cepstral coefficients... ";
        for(int64_t windowIndex = 0; windowIndex < numWindows; ++ windowIndex) {
            for(int32_t mfccIndex = 0; mfccIndex < numMfccs; ++ mfccIndex) {
                
                double total = 0;
                for(int32_t filterbankIndex = 0; filterbankIndex < numFilterbanks; ++ filterbankIndex) {
                    total += loggedFilterbankEnergies[windowIndex][filterbankIndex] * std::cos((3.141592653589793 / numFilterbanks) * (0.5 + filterbankIndex) * mfccIndex);
                }
                
                // May or may not be needed here
                // Makes the transformation orthoganal
                if(mfccIndex == 0) total *= std::sqrt(0.5);
                total *= std::sqrt(2.0 / numFilterbanks);
                
                mfccs[windowIndex][mfccIndex] = total;
            }
            
        }
        if(debugOutput) writeGenericHeatOutput(debugPrefix + "_mfcc.png", mfccs, numWindows, numMfccs, -4, 18);
        std::cout << "\tdone" << std::endl;
    }
    for(int64_t i = 0; i < numWindows; ++ i) {
        delete[] loggedFilterbankEnergies[i];
    }
    delete[] loggedFilterbankEnergies;
    
    return new MFCC(mfccs, numWindows, numMfccs);
}
