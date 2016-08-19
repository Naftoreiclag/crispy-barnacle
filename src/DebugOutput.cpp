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

#include "DebugOutput.hpp"

double normalized(double value, double min, double max) {
    return (value - min) / (max - min);
}
double normalizedClamp(double value, double min, double max) {
    double asdf = normalized(value, min, max);
    if(asdf < 0) asdf = 0;
    if(asdf > 1) asdf = 1;
    return asdf;
}

RGB::RGB(float ar, float ag, float ab)
: r(ar)
, g(ag)
, b(ab) {
}

unsigned char RGB::RU8() {
    if(r > 1.f) return 255;
    if(r < 0.f) return 0;
    return (unsigned char) (r * 255.f);
}
unsigned char RGB::GU8() {
    if(g > 1.f) return 255;
    if(g < 0.f) return 0;
    return (unsigned char) (g * 255.f);
}
unsigned char RGB::BU8() {
    if(b > 1.f) return 255;
    if(b < 0.f) return 0;
    return (unsigned char) (b * 255.f);
}

RGB colorrampSevenHeat(float normalizedIntensity) {
    if(normalizedIntensity < 0.f) {
        return RGB(0.3, 0.0, 0.3);
    }
    if(normalizedIntensity > 1.f) {
        return RGB(0.9, 0.0, 0.9);
    }
    if(normalizedIntensity < 1.f / 255.f) {
        return RGB(0.6, 0.6, 0.6);
    }
    
    for(int i = 1; i < 6; ++ i) {
        float upperBound = 0.16666666667f * i;
        if(normalizedIntensity < upperBound) {
            switch(i) {
                case 1: return interp(RGB(0, 0, 0), RGB(0, 0, 1), (upperBound - normalizedIntensity) * 6.f);
                case 2: return interp(RGB(0, 0, 1), RGB(0, 1, 1), (upperBound - normalizedIntensity) * 6.f);
                case 3: return interp(RGB(0, 1, 1), RGB(0, 1, 0), (upperBound - normalizedIntensity) * 6.f);
                case 4: return interp(RGB(0, 1, 0), RGB(1, 1, 0), (upperBound - normalizedIntensity) * 6.f);
                case 5: return interp(RGB(1, 1, 0), RGB(1, 0, 0), (upperBound - normalizedIntensity) * 6.f);
            }
        }
    }
    return interp(RGB(1, 0, 0), RGB(1, 1, 1), (1.0 - normalizedIntensity) * 6.f);
}

RGB interp(RGB a, RGB b, float aWeight) {
    float bWeight = 1 - aWeight;
    return RGB(a.r * aWeight + b.r * bWeight, a.g * aWeight + b.g * bWeight, a.b * aWeight + b.b * bWeight);
}
