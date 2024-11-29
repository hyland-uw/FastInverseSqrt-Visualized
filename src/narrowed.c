#include "util-harness.h"

// 0x5f37642f +/- 4000
#define INT_LO 1597461647
#define INT_HI 1597469647

// This is a temporary solution to generate a plot comparing very close
// FISR approximation magic constants. It works by copying out the
// deconstructed flow with some changes for a much narrower search.
// Ideally this would be a variation within the interface provided by deconstructed.c


/*
This harness allows us to deconstruct the working and output of an idealized FISR.

The main purpose of this is to visualize the result of selecting both the
input and the magic constant at random, so the space of possible choices is seen.

Unlike in the individual methods, the FISR is allowed to converge and the
iterations required for convergence is stored.

Unions are used to avoid type punning, so this code should work at all
levels of compiler optimization.
*/

narrowHarness narrow_decon_rsqrt(float x, uint32_t magic, float tol) {
    narrowHarness result;

    // Compute a reference inverse square root
    // Used to compute error
    result.reference = 1.0f / sqrtf(x);

    // float tol should be somewhere around or above 0.000125f


    // Track if we reach a state which won't plot well
    result.invalid_float_reached = false;

    // The input is given two simultaneous representations:
    // a 32 bit float y.f and a 32 bit unsigned integer y.u,
    // both from the same bitfield associated with x
    // For an integer, the position of the bit does not have
    // special meaning save that the more significant bits
    // represent larger numbers. By contrast, bit position
    // in floats is significant,with bit 0 being the sign,
    // 1-8 being exponents, and the remaining 23 the fraction.
    union { float f; uint32_t u; } y = {x};

    // The unisgned integer representation y.u is right shifted
    // dividing it by two. It is then subtracted from a
    // magic constant. The constant's main purpose is to restore
    // the exponents lost in that bit shift.
    // A choice of "0x5F400000" is exactly will
    // only restore the exponents in high 8 bits where
    // the float exponents are stored.
    y.u = magic - (y.u >> 1);

    // Now that we have manipulated the bitfield as an integer
    // and restored the bits in the exponent, we extract the
    // floating point representation. Treating the value
    // has the effect of removing us from the logarithmic domain
    result.initial_approx = y.f;

    // All of the FISRs we see in this library use at least
    // one iteration of Newton-Raphson.
    // Because different approximations choose differing
    // values for halfthree and halfone, we can select them
    int iters = 0;
    while (iters < 2) {
        // Hardcode 1.5 and 0.5 for this version
        y.f = y.f * (1.5f - 0.5f * x * y.f * y.f);
        iters++;
        // terminate NR iteration when we are close
        // rather than after 1 or 2 to better show
        // the possibility space
        if (fabs(y.f - result.reference) < tol) {
            break;
        }
    }
    // Record output after the while loop, then check
    // validity
    result.output = y.f;
    result.NR_iters = iters;

    if (!isnormal(result.initial_approx)) {
        result.invalid_float_reached = true;
    }
    // A poor choice of restoring constant can make the
    // resulting float invalid. isnormal() is chosen to
    // exclude subnormal numbers, which won't work with
    // the trick
    // c.f. https://stackoverflow.com/q/75772363/1188479
    //
    // We may also reach an invalid float through
    // overflow with (very) poor choices of
    // three or half
    if (!isnormal(result.output) || !isnormal(result.initial_approx)) {
        result.invalid_float_reached = true;
    }

    return result;
}

typedef struct {
    float input;
    float reference;
    float initial_approx;
    float output;
    unsigned int NR_iters;
    uint32_t partial_fraction;
    uint32_t magic;
} SampleResult;

void total_decon_rsqrt(uint32_t int_min, uint32_t int_max, float min, float max, float tol) {
    uint32_t int_width = int_max - int_min;
    SampleResult* results = malloc(sizeof(SampleResult) * int_width);
    uint32_t valid_results = 0;
    float* inputrange = malloc(sizeof(float) * int_width);

    for (uint32_t i = 0; i < int_width; i++) {
        inputrange[i] = logStratifiedSampler(min, max);
    }

    #pragma omp parallel
    {
        // Allocate local results with sufficient size
        SampleResult* local_results = malloc(sizeof(SampleResult) * int_width);
        uint32_t local_valid_results = 0;

        #pragma omp for
        for (uint32_t i = int_min; i < int_max; i++) {
            float x = inputrange[i - int_min];  // Adjust index to match inputrange
            uint32_t magic = i;

            narrowHarness result = narrow_decon_rsqrt(x, magic, tol);

            if (!result.invalid_float_reached) {
                SampleResult sample = {
                    .input = x,
                    .reference = result.reference,
                    .initial_approx = result.initial_approx,
                    .output = result.output,
                    .NR_iters = result.NR_iters,
                    .partial_fraction = extract_top10_fraction(x),
                    .magic = magic
                };
                local_results[local_valid_results++] = sample;
            }
        }

        // Merge local results into global results
        #pragma omp critical
        {
            for (uint32_t i = 0; i < local_valid_results && valid_results < int_width; i++) {
                results[valid_results++] = local_results[i];
            }
        }

        free(local_results);
    }

    // Print results
    printf("input,reference,initial,final,iters,partial_fraction,magic\n");
    for (uint32_t i = 0; i < valid_results; i++) {
        printf("%f,%f,%f,%f,%u,%d,0x%08X\n",
               results[i].input,
               results[i].reference,
               results[i].initial_approx,
               results[i].output,
               results[i].NR_iters,
               results[i].partial_fraction,
               results[i].magic);
    }

    free(inputrange);
    free(results);
}

// Main function calls the sampler which
// outputs via printf for a csv
int main() {
    srand(time(NULL));
    total_decon_rsqrt(INT_LO, INT_HI, FLOAT_START, FLOAT_END, 0.00024f);
    return 0;
}
