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

#include <fstream>
#include <iostream>
#include <cmath>

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
            std::cout << "\tFile appears to be ogg-vorbis" << std::endl;
            return Waveform::OGG_VORBIS;
        }
    }
    
    // Assume FLAC
    std::cout << "\tFile appears to be FLAC" << std::endl;
    return Waveform::FLAC;
}

// Loads a waveform into memory through a filename
int32_t loadWaveform(std::string filename, Waveform& returnVal) {
    std::cout << "\tLoading waveform from " << filename << std::endl;
    std::cout << "\tDetecting file format..." << std::endl;
    Waveform::InputFileFormat format = determineWaveformType(filename);
    
    int32_t errorVal;
    
    switch(format) {
        case Waveform::OGG_VORBIS: errorVal = loadWaveformAsOggVorbis(filename, returnVal); break;
        case Waveform::FLAC: errorVal = loadWaveformAsFLAC(filename, returnVal); break;
        default: return -1;
    }
    if(errorVal) {
        std::cerr << "Error while decoding!" << std::endl;
        
        delete[] returnVal.mOriginalSamples;
    } else {
        std::cout << "\tConverting into floats...";
        returnVal.mFloatSamples = new double[returnVal.mNumSamples];
        int64_t maxVal = 1 << (returnVal.mSampleSize * 8 - 1);
        double fmaxVal = maxVal;
        for(int64_t i = 0; i < returnVal.mNumSamples; ++ i) {
            returnVal.mFloatSamples[i] = ((double) returnVal.mOriginalSamples[i]) / fmaxVal;
        }
        std::cout << "\tdone" << std::endl;
        
        std::cout << "\tDone decoding ===" << std::endl;
        std::cout << "\tSample rate: " << returnVal.mSampleRate << std::endl;
        std::cout << "\tSample size: " << (returnVal.mSampleSize * 8) << " bits" << std::endl;
        std::cout << "\tSample count: " << returnVal.mNumSamples << std::endl;
        double lenSeconds = std::floor(1000.0 * returnVal.mNumSamples / returnVal.mSampleRate) / 1000.0;
        std::cout << "\tLength: " << lenSeconds << " seconds" << std::endl;
        
        
        int64_t allocByte = (sizeof(int64_t) + sizeof(double)) * returnVal.mNumSamples;
        double allocMB = std::floor(100.0 * allocByte / (1 << 20)) / 100.0;
        std::cout << "\tMemory allocated: " << allocByte << " bytes, (~ " << allocMB << "MB)" << std::endl;
    }
    
    return errorVal;
}

int32_t loadWaveformAsOggVorbis(std::string filename, Waveform& returnVal) {
    std::cout << "\tReading file as ogg-vorbis..." << std::endl;
    
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
    
    if(returnVal.mSampleSize > 2) {
        std::cerr << "Fatal error! Unsupported ogg-vorbis sample size (" << returnVal.mSampleSize << " bytes)!" << std::endl;
        ov_clear(&inputVorbisFile);
        return -1;
    }
    
    returnVal.mNumSamples = ov_pcm_total(&inputVorbisFile, -1);
    
    returnVal.mOriginalSamples = new int64_t[returnVal.mNumSamples];
    
    returnVal.mRunningTotal = 0;
    
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
                
                // Ogg vorbis sample sizes are either 1 or 2 bytes
                
                if(returnVal.mSampleSize == 2) {
                    int16_t sample = ((inputBuffer[i * 2] & 0x00ff) << 8) | (inputBuffer[i * 2 + 1] & 0x00ff);
                    returnVal.mOriginalSamples[returnVal.mRunningTotal + i] = sample;
                } else {
                    int8_t sample = inputBuffer[i] & 0xff;
                    returnVal.mOriginalSamples[returnVal.mRunningTotal + i] = sample;
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
    
    std::cout << "\tSuccessfully read ogg-vorbis file" << std::endl;
    
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
        returnVal.mOriginalSamples = new int64_t[returnVal.mNumSamples];
        returnVal.mRunningTotal = 0;
    }
    
    // (Past tense "read")
    unsigned int numSamplesRead = inputFrame->header.blocksize;
    
    for(unsigned int i = 0; i < numSamplesRead; ++ i) {
        returnVal.mOriginalSamples[returnVal.mRunningTotal] = inputBuffer[0][i];
        ++ returnVal.mRunningTotal;
        if(returnVal.mRunningTotal == returnVal.mNumSamples) {
            break;
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
        returnVal.mSampleSize = inputMetadata->data.stream_info.bits_per_sample >> 3;
        returnVal.mNumSamples = inputMetadata->data.stream_info.total_samples;
        // inputMetadata->data.stream_info.channels;
    }
    
}

void flacDecoderErrorCallback(
    const FLAC__StreamDecoder* inputFlacDecoder,
    const FLAC__StreamDecoderErrorStatus errorStatus,
    void* userData) {
    
    std::cout << "\tFLAC error callback: " << FLAC__StreamDecoderErrorStatusString[errorStatus];
}

int32_t loadWaveformAsFLAC(std::string filename, Waveform& returnVal) {
    std::cout << "\tReading file as FLAC..." << std::endl;
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
        std::cout << "\tFailure during FLAC decoder initialization: "
            << FLAC__StreamDecoderInitStatusString[inputDecoderStatus] << std::endl;
        FLAC__stream_decoder_delete(inputFlacDecoder);
        return -1;
    }
    
    FLAC__stream_decoder_process_until_end_of_stream(inputFlacDecoder);
    
    std::cout << "\tFLAC decoder final status: " << FLAC__StreamDecoderStateString[FLAC__stream_decoder_get_state(inputFlacDecoder)] << std::endl;
    
    // Reading cleanup
    FLAC__stream_decoder_delete(inputFlacDecoder);
    return 0;
}

int32_t writeWaveform(std::string filename, Waveform waveform) {
    int errorCode;
    
    // Write intialization
    vorbis_info outputVorbisInfo;
    vorbis_info_init(&outputVorbisInfo);
    int outputNumChannels = 1;
    errorCode = vorbis_encode_init_vbr(&outputVorbisInfo, outputNumChannels, waveform.mSampleRate, 1.0);
    if(errorCode) {
        std::cout << "\tFailed to initialize encoding for VBR!" << std::endl;
        return -1;
    }

    vorbis_dsp_state outputVorbisDspState;
    errorCode = vorbis_analysis_init(&outputVorbisDspState, &outputVorbisInfo);
    if(errorCode) {
        std::cout << "\tFailed to initialize vorbis dsp state!" << std::endl;
        return -1;
    }

    //
    vorbis_block outputVorbisBlock;
    vorbis_block_init(&outputVorbisDspState, &outputVorbisBlock);

    int cerealNumber = 1337; // Why?

    ogg_stream_state outputOggStream;
    ogg_stream_init(&outputOggStream, cerealNumber);

    ogg_page outputOggPage;

    std::ofstream outputData(filename.c_str(), std::ios::out | std::ios::binary);

    {
        ogg_packet streamIdentity;
        ogg_packet serializedMetadata;
        ogg_packet codebooks;
        
        vorbis_comment outputVorbisMetadata;
        vorbis_comment_init(&outputVorbisMetadata);
        vorbis_analysis_headerout(&outputVorbisDspState, &outputVorbisMetadata, &streamIdentity, &serializedMetadata, &codebooks);
        
        ogg_stream_packetin(&outputOggStream, &streamIdentity);
        ogg_stream_packetin(&outputOggStream, &serializedMetadata);
        ogg_stream_packetin(&outputOggStream, &codebooks);
        
        while(true) {
            int ret = ogg_stream_flush(&outputOggStream, &outputOggPage);
            if(ret == 0) {
                break;
            }
            outputData.write(reinterpret_cast<char*>(outputOggPage.header), outputOggPage.header_len);
            outputData.write(reinterpret_cast<char*>(outputOggPage.body), outputOggPage.body_len);
        }
        vorbis_comment_clear(&outputVorbisMetadata);
    }
    
    int64_t remainingSamples = waveform.mNumSamples;
    int64_t writtenSamples = 0;

    while(remainingSamples) {
        // Get a buffer which we must write to
        int numSamplesWrite = remainingSamples > 4096 ? 4096 : remainingSamples;
        
        float** outputBuffers = vorbis_analysis_buffer(&outputVorbisDspState, numSamplesWrite);
        
        for(int i = 0; i < numSamplesWrite; ++ i) {
            outputBuffers[0][i] = waveform.mFloatSamples[writtenSamples + i];
        }
        writtenSamples += numSamplesWrite;
        remainingSamples -= numSamplesWrite;
        
        // Tell vorbis that we are done writing to the buffer gained previously
        // (Necessary because we technically could have submitted less data than was requested,
        // and also because we need to notify vorbis that we are done writing for now)
        vorbis_analysis_wrote(&outputVorbisDspState, numSamplesWrite);
        
        //
        while(vorbis_analysis_blockout(&outputVorbisDspState, &outputVorbisBlock) == 1) {
            vorbis_analysis(&outputVorbisBlock, NULL);
            vorbis_bitrate_addblock(&outputVorbisBlock);
            
            ogg_packet outputOggPacket;
            while(vorbis_bitrate_flushpacket(&outputVorbisDspState, &outputOggPacket)) {
                ogg_stream_packetin(&outputOggStream, &outputOggPacket);
                while(true) {
                    int ret = ogg_stream_pageout(&outputOggStream, &outputOggPage);
                    if(ret == 0) {
                        break;
                    }
                    outputData.write(reinterpret_cast<char*>(outputOggPage.header), outputOggPage.header_len);
                    outputData.write(reinterpret_cast<char*>(outputOggPage.body), outputOggPage.body_len);
                }
            }
        }
        
        
        
    }
    vorbis_analysis_wrote(&outputVorbisDspState, 0);

    // Writing cleanup
    ogg_stream_clear(&outputOggStream);
    vorbis_block_clear(&outputVorbisBlock);
    vorbis_dsp_clear(&outputVorbisDspState);
    vorbis_info_clear(&outputVorbisInfo);
    return 0;
}
// Delete a waveform from memory
void freeWaveform(Waveform& wavefrom) {
    delete[] wavefrom.mOriginalSamples;
    delete[] wavefrom.mFloatSamples;
}
