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

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "WaveformIO.hpp"

struct ComplexNumber {
    double real;
    double imag;
};

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
    
    ComplexNumber** fftwCompleteOutput = new ComplexNumber*[numFrames];
    for(int64_t i = 0; i < numFrames; ++ i) {
        fftwCompleteOutput[i] = new ComplexNumber[frameLengthSamples];
    }
    
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
                    
                    // Hooray for compiler optimizations
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
                
                fftwCompleteOutput[frameIndex][inFrameSample].real = fftwOutput[inFrameSample][0];
                fftwCompleteOutput[frameIndex][inFrameSample].imag = fftwOutput[inFrameSample][1];
                
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
    
    {
        std::cout << "Computing power estimates... ";
        for(int64_t frameIndex = 0; frameIndex < numFrames; ++ frameIndex) {
            for(int64_t inFrameSample = 0; inFrameSample < frameLengthSamples; ++ inFrameSample) {
                
                double real = fftwCompleteOutput[frameIndex][inFrameSample].real;
                double imaginary = fftwCompleteOutput[frameIndex][inFrameSample].imag;
                double absValSq = real * real + imaginary * imaginary;
                double denom = frameLengthSamples;
                
                powerEstimates[frameIndex][inFrameSample] = absValSq / denom;
            }
        }
        std::cout << "done" << std::endl;
    }
    
    // Debug output power estimates as image
    {
        std::cout << "Writing debug image... ";
        int width = numFrames;
        int height = frameLengthSamples;
        char debugPowerEstimates[width * height];
        
        for(int pixelY = height - 1; pixelY >= 0; -- pixelY) {
            for(int pixelX = 0; pixelX < width; ++ pixelX) {
                
                int frame = pixelX;
                int spectrum = height - pixelY;
                
                double power = powerEstimates[frame][spectrum];
                
                if(power > 1.0) {
                    power = 1.0;
                } else if(power < 0.0) {
                    power = 0.0;
                }
                
                /*
                double ay = spectrum;
                double ax = frame;
                
                power = (ay * width + ax) / ((double) (width * height));
                */
                
                debugPowerEstimates[(pixelY * width + pixelX)    ] = power * 255;
                //debugPowerEstimates[(pixelY * width + pixelX) * 2 + 1] = power * 255;
            }
        }
        
        stbi_write_png("debug.png", width, height, 1, debugPowerEstimates, width);
        
        std::cout << "done" << std::endl;
    }
    
    for(int64_t i = 0; i < numFrames; ++ i) {
        delete[] fftwCompleteOutput[i];
    }
    delete[] fftwCompleteOutput;
    
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
