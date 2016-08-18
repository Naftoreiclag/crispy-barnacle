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

#ifndef WAVEFORMIO_HPP
#define WAVEFORMIO_HPP

#include <string>
#include <stdint.h>

struct Waveform {
    enum InputFileFormat {
        FLAC,
        OGG_VORBIS
    };
    
    InputFileFormat mInputFileFormat;
    
    // original samples in 64-bit words (max value is still determined by mSampleSize)
    int64_t* mOriginalSamples;
    
    // Converted into floats
    double* mFloatSamples;
    
    int32_t mSampleSize; // Size of individual samples in bytes
    int64_t mNumSamples; // Total number of samples
    int32_t mSampleRate; // Samples per second
    int64_t mRunningTotal; // Used only by FLAC decoding for write head persistence
};

// Loads a waveform into memory through a filename
int32_t loadWaveform(std::string filename, Waveform& returnVal);
// Delete a waveform from memory
void freeWaveform(Waveform& wavefrom);

// Internal
int32_t loadWaveformAsOggVorbis(std::string filename, Waveform& returnVal);
int32_t loadWaveformAsFLAC(std::string filename, Waveform& returnVal);

// Open a file just enough to tell what type it is
Waveform::InputFileFormat determineWaveformType(std::string filename);


#endif // WAVEFORMIO_HPP
