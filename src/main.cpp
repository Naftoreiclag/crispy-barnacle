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

#include <iostream>
#include <string>
#include <cmath>
#include <stdint.h>

#include "fftw3.h"

#include "WaveformIO.hpp"
#include "DebugOutput.hpp"
#include "ComplexNumber.hpp"

double melScale(double hertz) {
    return 1127 * std::log(hertz / 700 + 1);
}

double invMelScale(double mels) {
    return 700 * (std::exp(mels / 1127) - 1);
}

struct MFCCParams {
    double frameLengthMilliseconds = 25;
    double frameStepMilliseconds = 10;
    double filterMinFreq = 300;
    double filterMaxFreq = 8000;
    int32_t numFilterbanks = 26;
    int32_t numMfccs = 12;
};

struct MFCC {
    MFCC(double** mfccs, int32_t windowCount, int32_t mfccCount)
    : data(mfccs)
    , numWindows(windowCount)
    , numMfccs(mfccCount) {
    }
    ~MFCC() {
        for(int64_t i = 0; i < numWindows; ++ i) {
            delete[] data[i];
        }
        delete[] data;
    }
    
    double** data;
    
    int64_t numWindows;
    int32_t numMfccs;
    
    MFCCParams params;
};

MFCC* generateMFCC(
    Waveform inputAudio, 
    const MFCCParams params,
    bool debugOutput = false) {
        
    std::cout << "Generating MFCC data for audio... " << std::endl;
    
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
        if(debugOutput) writeFFTWOutputDebug("debug_fftw.png", fftwCompleteOutput, numWindows, spectrumLength);
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
        if(debugOutput) writeGenericHeatOutput("debug_power.png", powerEstimates, numWindows, spectrumLength);
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
        if(debugOutput) writeGenericHeatOutput("debug_energies_log.png", loggedFilterbankEnergies, numWindows, numFilterbanks, -10, 1);
        #ifndef NDEBUG
        if(debugOutput) writeGenericHeatOutput("debug_energies.png", filterbankEnergies, numWindows, numFilterbanks);
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
        if(debugOutput) writeGenericHeatOutput("debug_mfcc.png", mfccs, numWindows, numMfccs, -4, 18);
        std::cout << "\tdone" << std::endl;
    }
    for(int64_t i = 0; i < numWindows; ++ i) {
        delete[] loggedFilterbankEnergies[i];
    }
    delete[] loggedFilterbankEnergies;
    
    return new MFCC(mfccs, numWindows, numMfccs);
}

int run(std::string inputAudioFilename) {
    
    std::cout << "Filename is " << inputAudioFilename << std::endl;
    
    std::cout << "Loading mimic waveform into memory... " << std::endl;
    Waveform inputAudio;
    int32_t errorCode = loadWaveform(inputAudioFilename, inputAudio);
    
    if(errorCode) {
        std::cerr << "Fatal error! Failed to load waveform!" << std::endl;
        return -1;
    }
    std::cout << "Done loading mimic waveform" << std::endl;
    
    MFCCParams params;
    
    params.frameLengthMilliseconds = 25;
    params.frameStepMilliseconds = 10;
    params.filterMinFreq = 300;
    params.filterMaxFreq = 8000;
    params.numFilterbanks = 26;
    params.numMfccs = 12;
    
    std::cout << "Generating mimic MFCC... " << std::endl;
    MFCC* mimicMFCCs = generateMFCC(inputAudio, params, true);
    if(!mimicMFCCs) {
        std::cerr << "Fatal error! Failed to generate MFCCs for mimic data!" << std::endl;
        return -1;
    }
    std::cout << "Done generating mimic MFCC" << std::endl;
    
    return 0;
}

int main(int argc, char** argv) {
    
    if(argc == 1) {
        std::cout << "Usage: " << argv[0] << std::endl;
        return 0;
    } else if(argc > 1) {
        return run(argv[1]);
    }
    
    return 0;
}
