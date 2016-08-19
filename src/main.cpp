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
#include <stdint.h>

#include "WaveformIO.hpp"
#include "MFCC.hpp"

int run(std::string inputAudioFilename) {
    
    std::cout << "Filename is " << inputAudioFilename << std::endl;
    std::cout << std::endl;
    
    std::cout << "Loading mimic waveform into memory... " << std::endl;
    Waveform inputAudio;
    int32_t errorCode = loadWaveform(inputAudioFilename, inputAudio);
    
    if(errorCode) {
        std::cerr << "Fatal error! Failed to load waveform!" << std::endl;
        return -1;
    }
    std::cout << "Done loading mimic waveform" << std::endl;
    std::cout << std::endl;
    
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
    
    freeWaveform(inputAudio);
    
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
