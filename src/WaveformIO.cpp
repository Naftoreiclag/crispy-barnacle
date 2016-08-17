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

#include "WaveformIO.hpp"

#include <iostream>

#include "vorbis/vorbisfile.h"
#include "vorbis/vorbisenc.h"
#include "vorbis/codec.h"
#include "FLAC/stream_decoder.h"


// Open a file just enough to tell what type it is
Waveform::InputFileFormat determineWaveformType(std::string inputFilename) {
    // Attempt ogg-vorbis
    {
        OggVorbis_File inputVorbisFile;
        int errorCode;
        errorCode = ov_fopen(inputFilename.c_str(), &inputVorbisFile);
        if(!errorCode) {
            std::cout << "File appears to be ogg-vorbis" << std::endl;
            return Waveform::OGG_VORBIS;
        }
    }
    
    // Assume FLAC
    std::cout << "File appears to be FLAC" << std::endl;
    return Waveform::FLAC;
}

// Loads a waveform into memory through a filename
uint32_t loadWaveform(std::string filename, Waveform& returnVal) {
    std::cout << "Loading waveform from " << filename << std::endl;
    std::cout << "Detecting file format..." << std::endl;
    Waveform::InputFileFormat format = determineWaveformType(filename);
    
    switch(format) {
        case Waveform::OGG_VORBIS: return loadWaveformAsOggVorbis(filename, returnVal);
        case Waveform::FLAC: return loadWaveformAsFLAC(filename, returnVal);
        default: return -1;
    }
}

uint32_t loadWaveformAsOggVorbis(std::string filename, Waveform& returnVal) {
    std::cout << "Reading file as ogg-vorbis..." << std::endl;
    
    returnVal.mInputFileFormat = Waveform::OGG_VORBIS;
    
    // Read
    OggVorbis_File inputVorbisFile;
    int errorCode;
    errorCode = ov_fopen(filename.c_str(), &inputVorbisFile);
    if(errorCode) {
        if(errorCode == OV_ENOTVORBIS) {
            // It is supposed to be ogg vorbis!
            std::cerr << "Fatal error, this particular case should not be possible! (ogg-vorbis)" << std::endl;
            ov_clear(&inputVorbisFile);
            return -1;
        } else {
            // Something else (fatal) went wrong
            std::cerr << "Fatal error while initializing ogg-vorbis stream!" << std::endl;
            ov_clear(&inputVorbisFile);
            return -1;
        }
    }
    
    const vorbis_info* inputVorbisInfo = ov_info(&inputVorbisFile, -1);
    returnVal.mSampleRate = inputVorbisInfo->rate;
    
    // Reading information
    int logicalBitstreamPtr;
    int bigEndian = 0; // 0 = little endian, 1 = big endian
    returnVal.mSampleSize = 2;
    int signage = 1; // 0 = unsigned, 1 = signed
    char inputBuffer[4096 * returnVal.mSampleSize];
    
    returnVal.mNumSamples = ov_pcm_total(&inputVorbisFile, -1);
    
    returnVal.mOriginalSamples = new uint8_t[returnVal.mNumSamples * returnVal.mSampleSize];
    
    returnVal.mRunningTotal = 0;

    std::cout << "Sample rate: " << returnVal.mSampleRate << std::endl;
    std::cout << "Sample size: " << returnVal.mSampleSize << " bytes" << std::endl;
    std::cout << "Sample count: " << returnVal.mNumSamples << std::endl;
    
    while(true) {
        int ovRet = ov_read(&inputVorbisFile, inputBuffer, sizeof(inputBuffer), bigEndian, returnVal.mSampleSize, signage, &logicalBitstreamPtr);
        
        // End of file
        if(ovRet == 0) {
            break;
        }
        // Error
        else if(ovRet < 0) {
            std::cerr << "Fatal error while reading ogg-vorbis stream!" << std::endl;
            ov_clear(&inputVorbisFile);
            return -1;
        }
        // Positive return values are the number of bytes read
        else {
            int numSamplesRead = ovRet / returnVal.mSampleSize;
            for(int i = 0; i < numSamplesRead; ++ i) {
                for(int32_t le = 0; le < returnVal.mSampleSize; ++ le) {
                    returnVal.mOriginalSamples[(returnVal.mRunningTotal + i) * returnVal.mSampleSize + le] = inputBuffer[i * returnVal.mSampleSize + le] & 0xff;
                }
            }
            returnVal.mRunningTotal += numSamplesRead;
            if(returnVal.mRunningTotal > returnVal.mNumSamples) {
                std::cerr << "Sample count mismatch! Header declares too few samples!" << std::endl;
                ov_clear(&inputVorbisFile);
                return -1;
            }
        }
    }
    
    if(returnVal.mRunningTotal < returnVal.mNumSamples) {
        std::cerr << "Sample count mismatch! Header declares too many samples!" << std::endl;
        ov_clear(&inputVorbisFile);
        return -1;
    }
    
    std::cout << "Successfully read ogg-vorbis file" << std::endl;
    
    ov_clear(&inputVorbisFile);
    
    return 0;
}

FLAC__StreamDecoderWriteStatus flacDecoderWriteCallback(
    const FLAC__StreamDecoder* inputFlacDecoder,
    const FLAC__Frame* inputFrame,
    const FLAC__int32* const inputBuffer[],
    void* userData) {
    
    Waveform& returnVal = *(reinterpret_cast<Waveform*>(userData));
    int errorCode;
    
    // Called only once
    if(inputFrame->header.number.sample_number == 0) {
        returnVal.mOriginalSamples = new uint8_t[returnVal.mNumSamples * returnVal.mSampleSize];
        returnVal.mRunningTotal = 0;
    }
    
    // (Past tense "read")
    unsigned int numSamplesRead = inputFrame->header.blocksize;
    
    for(unsigned int i = 0; i < numSamplesRead; ++ i) {
        for(int32_t le = 0; le < returnVal.mSampleSize; ++ le) {
            returnVal.mOriginalSamples[(returnVal.mRunningTotal + i) * returnVal.mSampleSize + le] = (inputBuffer[0][i] >> (le * 8)) & 0xff;
        }
        returnVal.mRunningTotal += numSamplesRead;
        if(returnVal.mRunningTotal > returnVal.mNumSamples) {
            std::cerr << "Sample count mismatch! Header declares too few samples!" << std::endl;
            std::cerr << "Discovered "<< returnVal.mRunningTotal << " samples!" << std::endl;
            return FLAC__STREAM_DECODER_WRITE_STATUS_ABORT;
        }
    }
    
    return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;
}

void flacDecoderMetadataCallback(
    const FLAC__StreamDecoder* inputFlacDecoder,
    const FLAC__StreamMetadata* inputMetadata,
    void* userData) {
    
    Waveform& returnVal = *(reinterpret_cast<Waveform*>(userData));
    
    // comments
    if(inputMetadata->type == FLAC__METADATA_TYPE_VORBIS_COMMENT) {
    }
    
    // other
    if(inputMetadata->type == FLAC__METADATA_TYPE_STREAMINFO) {
        returnVal.mSampleRate = inputMetadata->data.stream_info.sample_rate;
        returnVal.mSampleSize = inputMetadata->data.stream_info.bits_per_sample;
        returnVal.mNumSamples = inputMetadata->data.stream_info.total_samples;
        
        std::cout << "Sample rate: " << returnVal.mSampleRate << std::endl;
        std::cout << "Sample size: " << returnVal.mSampleSize << " bytes" << std::endl;
        std::cout << "Sample count: " << returnVal.mNumSamples << std::endl;
        // inputMetadata->data.stream_info.channels;
    }
    
}

void flacDecoderErrorCallback(
    const FLAC__StreamDecoder* inputFlacDecoder,
    const FLAC__StreamDecoderErrorStatus errorStatus,
    void* userData) {
    
    std::cout << "FLAC error callback: " << FLAC__StreamDecoderErrorStatusString[errorStatus];
}

uint32_t loadWaveformAsFLAC(std::string filename, Waveform& returnVal) {
    std::cout << "Reading file as FLAC..." << std::endl;
    returnVal.mInputFileFormat = Waveform::FLAC;
    
    FLAC__StreamDecoder* inputFlacDecoder = FLAC__stream_decoder_new();
    
    if(!inputFlacDecoder) {
        std::cerr << "Failed to allocate FLAC decoder!" << std::endl;
        return -1;
    }
    
    FLAC__stream_decoder_set_md5_checking(inputFlacDecoder, true); // Why not?
    
    FLAC__StreamDecoderInitStatus inputDecoderStatus = 
        FLAC__stream_decoder_init_file(inputFlacDecoder, filename.c_str(),
        flacDecoderWriteCallback, flacDecoderMetadataCallback, flacDecoderErrorCallback, &returnVal);
    
    if(inputDecoderStatus != FLAC__STREAM_DECODER_INIT_STATUS_OK) {
        std::cout << "Failure during FLAC decoder initialization: "
            << FLAC__StreamDecoderInitStatusString[inputDecoderStatus] << std::endl;
        FLAC__stream_decoder_delete(inputFlacDecoder);
        return -1;
    }
    
    FLAC__stream_decoder_process_until_end_of_stream(inputFlacDecoder);
    
    std::cout << "FLAC decoder status: " << FLAC__StreamDecoderStateString[FLAC__stream_decoder_get_state(inputFlacDecoder)] << std::endl;
    
    // Reading cleanup
    FLAC__stream_decoder_delete(inputFlacDecoder);
    return 0;
}

// Delete a waveform from memory
void freeWaveform(Waveform& wavefrom) {
    delete[] wavefrom.mOriginalSamples;
}