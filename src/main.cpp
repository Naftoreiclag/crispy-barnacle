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
#include "DebugOutput.hpp"


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
    int32_t windowLength = (frameLengthMilliseconds * inputAudio.mSampleRate) / 1000;
    int32_t spectrumLength = windowLength / 2;
    int32_t windowStep = (frameStepMilliseconds * inputAudio.mSampleRate) / 1000;
    
    std::cout << "Frame length: " << windowLength << " samples / " << frameLengthMilliseconds << "ms" << std::endl;
    std::cout << "Frame step: " << windowStep << " samples / " << frameStepMilliseconds << "ms" << std::endl;
    int64_t numWindows = 0;
    for(int64_t windowIndex = 0; (windowIndex * windowStep) < inputAudio.mNumSamples; ++ windowIndex) {
        numWindows ++;
    }
    std::cout << "Frame count: " << numWindows << std::endl;
    
    std::cout << "Allocating additional memory... ";
    double** powerEstimates = new double*[numWindows];
    for(int64_t i = 0; i < numWindows; ++ i) {
        powerEstimates[i] = new double[spectrumLength];
    }
    ComplexNumber** fftwCompleteOutput = new ComplexNumber*[numWindows];
    for(int64_t i = 0; i < numWindows; ++ i) {
        fftwCompleteOutput[i] = new ComplexNumber[windowLength];
    }
    std::cout << "done" << std::endl;
    
    {
        fftw_complex* fftwInput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * windowLength);
        fftw_complex* fftwOutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * windowLength);
        
        std::cout << "Optimizing FFTW... ";
        fftw_plan fftwPlan = fftw_plan_dft_1d(windowLength, fftwInput, fftwOutput, FFTW_FORWARD, FFTW_MEASURE);
        std::cout << "done" << std::endl;
        
        // Set imaginary components to be zero
        for(int64_t windowSample = 0; windowSample < windowLength; ++ windowSample) {
            fftwInput[windowSample][1] = 0;
        }
        
        std::cout << "Window function: Hanning" << std::endl;
        std::cout << "Performing DFT... ";
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
            
            for(int64_t windowSample = 0; windowSample < windowLength; ++ windowSample) {
                
                fftwCompleteOutput[windowIndex][windowSample].real = fftwOutput[windowSample][0];
                fftwCompleteOutput[windowIndex][windowSample].imag = fftwOutput[windowSample][1];
            }
        }
        std::cout << "done" << std::endl;
        
        
        fftw_destroy_plan(fftwPlan);
        fftw_free(fftwOutput);
        fftw_free(fftwInput);
    }
    
    {
        std::cout << "Computing power estimates... ";
        for(int64_t windowIndex = 0; windowIndex < numWindows; ++ windowIndex) {
            for(int64_t windowSample = 0; windowSample < spectrumLength; ++ windowSample) {
                
                double real = fftwCompleteOutput[windowIndex][windowSample].real;
                double imaginary = fftwCompleteOutput[windowIndex][windowSample].imag;
                double absValSq = real * real + imaginary * imaginary;
                double denom = windowLength;
                
                
                powerEstimates[windowIndex][windowSample] = absValSq / denom;
            }
        }
        std::cout << "done" << std::endl;
    }
    
    // Debug output power estimates as image
    {
        
        std::cout << "Writing debug images... ";
        {
            int width = numWindows;
            int height = spectrumLength;
            char imageData[width * height * 3];
            
            for(int pixelY = height - 1; pixelY >= 0; -- pixelY) {
                for(int pixelX = 0; pixelX < width; ++ pixelX) {
                    
                    int frame = pixelX;
                    int spectrum = height - pixelY;
                    
                    {
                        double power = powerEstimates[frame][spectrum];
                        
                        if(power > 1.0) {
                            power = 1.0;
                        } else if(power < 0.0) {
                            power = 0.0;
                        }
                        
                        RGB heat = colorrampSevenHeat(power);
                        
                        imageData[(pixelY * width + pixelX) * 3    ] = heat.RU8();
                        imageData[(pixelY * width + pixelX) * 3 + 1] = heat.GU8();
                        imageData[(pixelY * width + pixelX) * 3 + 2] = heat.BU8();
                    }
                }
            }
            stbi_write_png("debug_power.png", width, height, 3, imageData, width * 3);
        }
        std::cout << "Writing debug images... ";
        {
            int width = numWindows;
            int height = spectrumLength;
            char imageData[width * height * 3];
            
            for(int pixelY = height - 1; pixelY >= 0; -- pixelY) {
                for(int pixelX = 0; pixelX < width; ++ pixelX) {
                    
                    int frame = pixelX;
                    int spectrum = height - pixelY;
                    
                    {
                        double power = powerEstimates[frame][spectrum];
                        power = std::sqrt(power);
                        
                        if(power > 1.0) {
                            power = 1.0;
                        } else if(power < 0.0) {
                            power = 0.0;
                        }
                        
                        RGB heat = colorrampSevenHeat(power);
                        
                        imageData[(pixelY * width + pixelX) * 3    ] = heat.RU8();
                        imageData[(pixelY * width + pixelX) * 3 + 1] = heat.GU8();
                        imageData[(pixelY * width + pixelX) * 3 + 2] = heat.BU8();
                    }
                }
            }
            stbi_write_png("debug_power_sqrt.png", width, height, 3, imageData, width * 3);
        }
        {
            int width = numWindows;
            int height = spectrumLength;
            char imageData[width * height * 3];
            
            for(int pixelY = height - 1; pixelY >= 0; -- pixelY) {
                for(int pixelX = 0; pixelX < width; ++ pixelX) {
                    
                    int frame = pixelX;
                    int spectrum = height - pixelY;
                    
                    double real = fftwCompleteOutput[frame][spectrum].real;
                    double imag = fftwCompleteOutput[frame][spectrum].imag;
                    
                    if(real < 0) real = -real;
                    if(imag < 0) imag = -imag;
                    
                    if(real > 1.0) {
                        real = 1.0;
                    } else if(real < 0.0) {
                        real = 0.0;
                    }
                    
                    if(imag > 1.0) {
                        imag = 1.0;
                    } else if(imag < 0.0) {
                        imag = 0.0;
                    }
                    
                    imageData[(pixelY * width + pixelX) * 3    ] = real * 255;
                    imageData[(pixelY * width + pixelX) * 3 + 1] = imag * 255;
                    imageData[(pixelY * width + pixelX) * 3 + 2] = 0;
                }
            }
            
            stbi_write_png("debug_fftw_output.png", width, height, 3, imageData, width * 3);
        }
        
        std::cout << "done" << std::endl;
    }
    
    for(int64_t i = 0; i < numWindows; ++ i) {
        delete[] fftwCompleteOutput[i];
    }
    delete[] fftwCompleteOutput;
    
    for(int64_t i = 0; i < numWindows; ++ i) {
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
