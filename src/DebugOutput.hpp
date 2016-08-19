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

#ifndef DEBUGOUTPUT_HPP
#define DEBUGOUTPUT_HPP

#include <stdint.h>

#include "ComplexNumber.hpp"

struct RGB {
    RGB(float r, float g, float b);
    
    float r;
    float g;
    float b;
    
    unsigned char RU8();
    unsigned char GU8();
    unsigned char BU8();
};

RGB colorrampSevenHeat(float normalizedIntensity);

RGB interp(RGB a, RGB b, float amnt);

double normalized(double value, double min, double max);
double normalizedClamp(double value, double min, double max);

void writeFFTWOutputDebug(char* filename, ComplexNumber** fftwCompleteOutput, int32_t numWindows, int32_t spectrumLength);
void writeGenericHeatOutput(char* filename, double** data, int width, int height, double minRange = 0, double maxRange = 1);
#endif // DEBUGOUTPUT_HPP
