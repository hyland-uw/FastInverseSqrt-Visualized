#include "util-harness.h"


// For deconstruction
// experiments show that if it does not converge after
// 94 it probably will not converge (tested up to 2000)
#define MAX_NR 95
#define INT_DRAWS 262144
#define INT_DRAW_MIN 1592500000
#define INT_DRAW_MAX 1602000000

/*
This harness allows us to deconstruct the working and output of an idealized FISR.

The main purpose of this is to visualize the result of selecting both the
input and the magic constant at random, so the space of possible choices is seen.

Unlike in the individual methods, the FISR is allowed to converge and the
iterations required for convergence is stored.

Unions are used to avoid type punning, so this code should work at all
levels of compiler optimization.
*/

GeneralizedHarness decon_rsqrt(float x, int NRmax, uint32_t magic, float tol) {
    return generalized_rsqrt(x, NRmax, magic, tol, true);
}

typedef struct {
    float input;
    float reference;
    float initial_approx;
    float after_one;
    float output;
    unsigned int NR_iters;
    uint32_t magic;
} SampleResult;

void sample_decon_rsqrt(int draws, int NRmax, float min, float max, float tol) {
    SampleResult* results = malloc(sizeof(SampleResult) * draws);
    int valid_results = 0;
    float inputrange[draws];
    for (int i = 0; i < draws; i++) {
         inputrange[i] = logStratifiedSampler(min, max);
    }

    #pragma omp parallel
    {
        // Thread-local storage for results
        int threads = omp_get_num_threads();
        SampleResult* local_results = malloc(sizeof(SampleResult) * draws / threads);
        int local_valid_results = 0;

        #pragma omp for
        for (int i = 0; i < draws; i++) {
            float x = inputrange[i];
            // Wider integer sample range than the others. we are looking for
            // poor fits
            uint32_t magic = sample_integer_range(INT_DRAW_MIN, INT_DRAW_MAX);

            GeneralizedHarness result = decon_rsqrt(x, NRmax, magic, tol);

            if (!result.invalid_float_reached) {
                SampleResult sample = {
                    .input = x,
                    .reference = result.reference,
                    .initial_approx = result.initial_approx,
                    .after_one = result.after_one,
                    .output = result.output,
                    .NR_iters = result.NR_iters,
                    .magic = magic
                };
                local_results[local_valid_results++] = sample;
            }
        }

        // Merge local results into global results
        #pragma omp critical
        {
            for (int i = 0; i < local_valid_results && valid_results < draws; i++) {
                results[valid_results++] = local_results[i];
            }
        }

        free(local_results);
    }

    // Print results
    printf("input,reference,initial,after_one,final,iters,magic\n");
    for (int i = 0; i < valid_results; i++) {
        printf("%f,%f,%f,%f,%f,%u,0x%08X\n",
               results[i].input,
               results[i].reference,
               results[i].initial_approx,
               results[i].after_one,
               results[i].output,
               results[i].NR_iters,
               results[i].magic);
    }

    free(results);
}

// Main function calls the sampler which
// outputs via printf for a csv
int main() {
    srand(time(NULL));
    sample_decon_rsqrt(INT_DRAWS, MAX_NR, FLOAT_START, FLOAT_END, FLOAT_TOL);
    return 0;
}
