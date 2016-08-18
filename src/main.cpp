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

int run(std::string inputAudioFilename) {
    
    std::cout << "Filename is " << inputAudioFilename << std::endl;
    
    Waveform inputAudio;
    int32_t errorCode = loadWaveform(inputAudioFilename, inputAudio);
    
    if(errorCode) {
        std::cerr << "Fatal error! Exiting application..." << std::endl;
        return -1;
    }
    
    int32_t frameLengthMilliseconds = 25;
    int32_t frameStepMilliseconds = 10;
    
    // Rounded toward zero
    int32_t frameLengthSamples = (frameLengthMilliseconds * inputAudio.mSampleRate) / 1000;
    int32_t frameStepSamples = (frameStepMilliseconds * inputAudio.mSampleRate) / 1000;
    
    std::cout << "Frame length: " << frameLengthSamples << " samples / " << frameLengthMilliseconds << "ms" << std::endl;
    std::cout << "Frame step: " << frameStepSamples << " samples / " << frameStepMilliseconds << "ms" << std::endl;
    int64_t numFrames = 0;
    for(int64_t frameIndex = 0; (frameIndex * frameStepSamples) < inputAudio.mNumSamples; ++ frameIndex) {
        numFrames ++;
    }
    std::cout << "Frame count: " << numFrames << std::endl;
    
    std::cout << "Allocating memory for power spectral estimates... ";
    
    double** powerEstimates = new double*[numFrames];
    for(int64_t i = 0; i < numFrames; ++ i) {
        powerEstimates[i] = new double[frameLengthSamples];
    }
    
    std::cout << "done" << std::endl;
    
    {
        fftw_complex* fftwInput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * frameLengthSamples);
        fftw_complex* fftwOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * frameLengthSamples);
        
        std::cout << "Optimizing FFTW... ";
        fftw_plan fftwPlan = fftw_plan_dft_1d(frameLengthSamples, fftwInput, fftwOutput, FFTW_FORWARD, FFTW_MEASURE);
        std::cout << "done" << std::endl;
        
        // Set imaginary components to be zero
        for(int64_t inFrameSample = 0; inFrameSample < frameLengthSamples; ++ inFrameSample) {
            fftwInput[inFrameSample][1] = 0;
        }
        
        std::cout << "Window function: Hanning" << std::endl;
        std::cout << "Performing DFT... ";
        for(int64_t frameIndex = 0; frameIndex < numFrames; ++ frameIndex) {
            int64_t frameStartSampleIndex = frameIndex * frameStepSamples;
            for(int64_t inFrameSample = 0; inFrameSample < frameLengthSamples; ++ inFrameSample) {
                int64_t currentSampleIndex = frameStartSampleIndex + inFrameSample;
                
                double sample;
                if(currentSampleIndex >= inputAudio.mNumSamples) {
                    sample = 0;
                } else {
                    sample = inputAudio.mFloatSamples[currentSampleIndex];
                    
                    // Apply hanning window
                    
                    // Horray for compiler optimizations
                    double tau = 6.28318530717958647692528677;
                    double numerator = tau * inFrameSample;
                    double denominator = frameLengthSamples - 1;
                    double hanning = 0.5 * (1.0 - std::cos(numerator / denominator));
                    
                    sample *= hanning;
                }
                
                fftwInput[inFrameSample][0] = sample;
            }
            
            fftw_execute(fftwPlan);
            
            for(int64_t inFrameSample = 0; inFrameSample < frameLengthSamples; ++ inFrameSample) {
                double real = fftwOutput[inFrameSample][0];
                double imaginary = fftwOutput[inFrameSample][1];
                double absValSq = real * real + imaginary * imaginary;
                double denom = frameLengthSamples;
                
                powerEstimates[frameIndex][inFrameSample] = absValSq / denom;
            }
        }
        std::cout << "done" << std::endl;
        
        
        fftw_destroy_plan(fftwPlan);
        fftw_free(fftwOutput);
        fftw_free(fftwInput);
    }
    
    for(int64_t i = 0; i < numFrames; ++ i) {
        delete[] powerEstimates[i];
    }
    delete[] powerEstimates;
    
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
