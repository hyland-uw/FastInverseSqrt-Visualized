#include "sampling-harness.h"

//// Implementations of various inverse square roots.
//// all should accept:
////     an input float x, assumed to be positive and normal
////     an input int NR, assumed to be > 0
//// all should return:
////     a float which is an approximation of 1/sqrt(x)

//Kahan's "Magic" square root
// See http://www.arithmazium.org/library/lib/wk_square_root_aug80.pdf (1980)
// may have a predecessor in the 1970s in UNIX Version 5
// https://minnie.tuhs.org/cgi-bin/utree.pl?file=V5/usr/source/s3/sqrt.s
float MagicISR(float x, int NR) {
    int magic = 0x1FBF8300;
    union { float f; uint32_t u; } y = {x};
    y.u = magic - (y.u >> 1);
    while (NR > 0) {
        y.f = (y.f + (x / y.f))/2;
        NR--;
    }
    // Kahan's magic sqrt is for a square root, not an inverse sqrt
    // We want to compare the bit shifting properly so we reciprocate at the end
    return  1.0f / y.f;
}


// FISR demonstrating the logaritmic nature of the floating-point numbers
// The constant chosen is not "magic" but that which replaces the lost
// exponents from the shift of the number as an integer.
// 1997
float BlinnISR(float x, int NR) {
    int magic = 0x5F400000;
    union { float f; uint32_t u; } y = {x};
    // 0x5F400000 = 1598029824
    y.u = magic - (y.u >> 1);
    while (NR > 0) {
        // See https://doi.org/10.1109/38.595279
        y.f = y.f * (1.47f - 0.47f * x * y.f * y.f);
        NR--;
    }
    return y.f;
}

// See https://en.wikipedia.org/wiki/Fast_inverse_square_root
// Published ~1999
float QuakeISR(float x, int NR) {
    int magic = 0x5f3759df;
    union { float f; uint32_t u; } y = {x};
    y.u = magic - (y.u >> 1);

    while (NR > 0) {
        y.f = y.f * (1.5f - (0.5f * x * y.f * y.f));
        NR--;
    }
    return  y.f;
}

// "Square root without division" From Kahan 1999
// See http://www.arithmazium.org/library/lib/wk_square_root_without_division_feb99.pdf
float withoutDivISR(float x, int NR) {
    // 0x5f39d015 is the hex value for
    // the float which has an eponent of 190
    // and a fraction of ~.451
    int magic = 0x5f39d015;
    union { float f; uint32_t u; } y = {x};
    y.u = magic - (y.u >> 1);
    while (NR > 0) {
        y.f = y.f * (1.5f - (0.5f * x * y.f * y.f));
        NR--;
    }
    return y.f;
}


// See https://web.archive.org/web/20220826232306/http://rrrola.wz.cz/inv_sqrt.html
// Optimal Newton-Raphson / magic constant combination
// found viabrute force search of the 32-bit integer space
// for magic constant and coefficients in Newton-Raphson step
// post 2010
float optimalFISR (float x, int NR) {
    union { float f; uint32_t u; } y = {x};
    int magic = 0x5F1FFFF9ul;
    y.u = magic - (y.u >> 1);
    while (NR > 0) {
        y.f = 0.703952253f * y.f * (2.38924456f - x * y.f * y.f);
        NR--;
    }
    return y.f;
}

// Constant choice from Moroz et al (2018)
// https://doi.org/10.48550/arXiv.1603.04483 is original 2016 paper
// NR values and different constant chosen from 2018 paper
// https://doi.org/10.48550/arXiv.1802.06302
// See page 15 for InvSqrt2 constants
float MorozISR(float x, int NR) {
    float halfx = x * 0.5f;
    int magic = 0x5F376908;
    union { float f; uint32_t u; } y = {x};
    y.u = magic - (y.u >> 1);
    while (NR > 0) {
        y.f = y.f * (1.5008789f - halfx * y.f * y.f);
        NR--;
    }
    return y.f;
}

/*
LaLonde and Dawson's 1990 method was previously checked into this repository,
it is now at https://gist.github.com/Protonk/dfbcab17986777ff997f24dcdd8e3bbc
*/

// This construction allows easy dispatch for each of the FISR methods
// we are testing.

ISREntry isr_table[] = {
    {"MagicISR", MagicISR},
    {"BlinnISR", BlinnISR},
    {"QuakeISR", QuakeISR},
    {"withoutDivISR", withoutDivISR},
    {"optimalFISR", optimalFISR},
    {"MorozISR", MorozISR},
    {NULL, NULL} // Sentinel to mark the end of the array
};

float FISR(const char *name, float x, int NR) {
    for (int i = 0; isr_table[i].name != NULL; i++) {
        if (strcmp(name, isr_table[i].name) == 0) {
            return isr_table[i].func(x, NR);
        }
    }
    // Function not found
    fprintf(stderr, "Error: ISR function '%s' not found\n", name);
    return NAN;
}

int main() {
    // Draw parameters from the sampling harness
    int num_draws = NUM_FLOATS;
    float min_input = FLOAT_START;  // This will be adjusted to FLT_MIN if it's smaller
    float max_input = FLOAT_END;

    if (min_input >= max_input) {
        fprintf(stderr, "Error: min_input must be less than max_input\n");
        return 1;
    }
    srand(time(NULL));
    // Prints in tidy format, e.g.:
    //
    // ISR_function, input, reference, NR_0, NR_1
    // BlinnISR, 1.25, 0.8944272, 0.9, 0.895
    printf("ISR_function, input, reference, NR_0, NR_1\n");

    for (int draw = 0; draw < num_draws; draw++) {
        float input = reciprocalRange(min_input, max_input);
        float reference = 1.0f / sqrtf(input);

        for (int i = 0; isr_table[i].name != NULL; i++) {
            float result_nr0 = FISR(isr_table[i].name, input, 0);
            float result_nr1 = FISR(isr_table[i].name, input, 1);

            printf("%s, %.6e, %.7e, %.6e, %.6e\n",
                   isr_table[i].name, input, reference, result_nr0, result_nr1);
        }
    }

    return 0;
}