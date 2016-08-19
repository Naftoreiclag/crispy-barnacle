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

int run(std::string paletteAudioFilename, std::string templateAudioFilename) {
    
    // PARAMETERS
    // should be kept globally constant
    
    MFCCParams params;
    
    params.frameLengthMilliseconds = 25;
    params.frameStepMilliseconds = 10;
    params.filterMinFreq = 300;
    params.filterMaxFreq = 8000;
    params.numFilterbanks = 26;
    params.numMfccs = 12;
    
    // ECHO
    
    std::cout << "Palette filename is " << paletteAudioFilename << std::endl;
    std::cout << "Template filename is " << templateAudioFilename << std::endl;
    std::cout << std::endl;
    
    // READ WAVEFORM
    int32_t errorCode;
    
    std::cout << "Loading palette waveform into memory... " << std::endl;
    Waveform paletteAudio;
    errorCode = loadWaveform(paletteAudioFilename, paletteAudio);
    if(errorCode) {
        std::cerr << "Fatal error! Failed to load waveform!" << std::endl;
        return -1;
    }
    std::cout << "Done loading palette waveform" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Loading template waveform into memory... " << std::endl;
    Waveform templateAudio;
    errorCode = loadWaveform(templateAudioFilename, templateAudio);
    if(errorCode) {
        std::cerr << "Fatal error! Failed to load waveform!" << std::endl;
        return -1;
    }
    std::cout << "Done loading template waveform" << std::endl;
    std::cout << std::endl;
    
    
    // MFCC GENERATION
    
    std::cout << "Generating palette MFCC... " << std::endl;
    MFCC* paletteMFCCs = generateMFCC(paletteAudio, params, true, "palette");
    if(!paletteMFCCs) {
        std::cerr << "Fatal error! Failed to generate MFCCs for palette data!" << std::endl;
        return -1;
    }
    std::cout << "Done generating palette MFCC" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Generating template MFCC... " << std::endl;
    MFCC* templateMFCCs = generateMFCC(templateAudio, params, true, "template");
    if(!templateMFCCs) {
        std::cerr << "Fatal error! Failed to generate MFCCs for template data!" << std::endl;
        return -1;
    }
    std::cout << "Done generating template MFCC" << std::endl;
    std::cout << std::endl;
    
    // CLEANUP
    
    std::cout << "Cleaning up... ";
    freeWaveform(paletteAudio);
    freeWaveform(templateAudio);
    std::cout << "\tdone" << std::endl;
    
    return 0;
}

int main(int argc, char** argv) {
    
    if(argc == 1) {
        std::cout << "Usage: " << argv[0] << std::endl;
        return 0;
    } else if(argc > 2) {
        return run(argv[1], argv[2]);
    }
    
    return 0;
}
