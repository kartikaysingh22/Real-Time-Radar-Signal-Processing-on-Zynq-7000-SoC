// ============================================================
// Triggering Based Real Time Radar Signal Processing on Edge
// IIIT Delhi
// ============================================================
// Status:
// parameters check                     - done
// range axis computation               - done
// golay generation                     - done
// zadoff chu generation                - done
// fifo timing parameter calculation    - done
// velocity axis generation             - done
// cross correlation functionality      - done
// fifo loop                            - done
// fft and ifft functionality           - done
// timing instrumentation               - done
// ============================================================

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "xil_printf.h"
#include "platform.h"
#include "xparameters.h"
#include "xtime_l.h"


/***************************************************************************
 * SYSTEM PARAMETERS
 ***************************************************************************/

#define C                 3e8f        // Speed of light (m/s)
#define FC                30e9f       // Carrier Frequency (Hz) — 30 GHz
#define FS                300e6f      // Sampling Frequency (Hz) — 300 MHz

#define NSAMP             512         // Base sequence length (samples per pulse)
#define FIFO_SIZE_SAMPLES 32768       // FIFO memory size in samples
#define NUM_BATCHES       16          // Number of FIFO refills / batches
#define NUM_MC            50          // Number of Monte Carlo runs per SNR level

#define SNR_START         -45         // Start of SNR sweep (dB)
#define SNR_END           -30         // End of SNR sweep (dB)
#define NUM_SNR_LEVELS    (SNR_END - SNR_START + 1)  // Total SNR levels (16)

#define MAX_PULSES        512         // Maximum number of slow-time pulses (array bound)

#define PI                3.14159265358979323846f  // Pi constant


/***************************************************************************
 * HELPER MACRO
 * Converts XTime tick difference to microseconds.
 *
 * COUNTS_PER_SECOND is provided by the Xilinx BSP (xtime_l.h).
 * Dividing by (COUNTS_PER_SECOND / 1e6) converts raw tick counts
 * into floating-point microseconds.
 ***************************************************************************/
#define TIME_US(start, end) \
    ((float)((end) - (start)) / (COUNTS_PER_SECOND / 1000000))


/***************************************************************************
 * GLOBAL MEMORY
 * All large arrays declared globally to avoid stack overflow on Zedboard
 * (ARM Cortex-A9 default stack is small).
 ***************************************************************************/

// Range axis: r_axis[i] = i * (c / (2 * fs)) — maps bin index to metres
float r_axis[NSAMP];

// Radar data cubes: [range bin][pulse index]
// cube_g — Golay complementary pair accumulation
// cube_z — Zadoff-Chu accumulation
float complex cube_g[NSAMP][MAX_PULSES]; // Golay radar cube
float complex cube_z[NSAMP][MAX_PULSES]; // ZC radar cube

// Range-Doppler maps (output of slow-time FFT across cube columns)
float complex rd_map_g[NSAMP][MAX_PULSES];  // Golay Range-Doppler Map
float complex rd_map_z[NSAMP][MAX_PULSES];  // ZC Range-Doppler Map

// Golay complementary sequences (length NSAMP each)
float complex ga[NSAMP];  // Golay sequence A
float complex gb[NSAMP];  // Golay sequence B

// Zadoff-Chu base sequence (length NSAMP)
float complex zc[NSAMP];

// Zero-padded Golay Tx signals (sequence + zero guard interval)
float complex tx_ga[2 * NSAMP];   // Zero-padded Golay A for transmission
float complex tx_gb[2 * NSAMP];   // Zero-padded Golay B for transmission
int L_golay;                       // Actual Golay pulse length = 2 * NSAMP

// ZC Tx signal with cyclic prefix (first copy = CP, second = data)
float complex tx_zc[2 * NSAMP];   // ZC signal for transmission
int L_zc;                          // Actual ZC pulse length = 2 * NSAMP

// Doppler / velocity axes for Golay
float doppler_bins_g[MAX_PULSES];  // Symmetric Doppler bin indices (-N/2 … N/2-1)
float doppler_axis_g[MAX_PULSES];  // Doppler frequency axis (Hz)
float v_axis_g[MAX_PULSES];        // Velocity axis (m/s)

// Doppler / velocity axes for ZC
float doppler_bins_z[MAX_PULSES];  // Symmetric Doppler bin indices
float doppler_axis_z[MAX_PULSES];  // Doppler frequency axis (Hz)
float v_axis_z[MAX_PULSES];        // Velocity axis (m/s)

float lambda = C / FC;   // Wavelength of Tx Wave (metres)

// RMSE storage arrays (one value per SNR level)
float rmse_r_g[NUM_SNR_LEVELS];   // Golay range RMSE (metres)
float rmse_v_g[NUM_SNR_LEVELS];   // Golay velocity RMSE (m/s)
float rmse_r_z[NUM_SNR_LEVELS];   // ZC range RMSE (metres)
float rmse_v_z[NUM_SNR_LEVELS];   // ZC velocity RMSE (m/s)

// Shared processing buffers (global to avoid stack overflow)
float complex rx[2 * NSAMP];      // Received signal buffer (max pulse length)
float complex corr[2 * NSAMP];    // Full cross-correlation output buffer
float complex rx_win[NSAMP];      // ZC windowed receive (after cyclic-prefix removal)
float complex circ_out[NSAMP];    // ZC circular correlation output


// --- ADDED FOR FFT OPTIMIZATION ---
// Maximum FFT size needs to handle the padded length in correlate_fft (4096)
#define MAX_FFT_SIZE 4096
float complex twiddle_LUT[MAX_FFT_SIZE / 2];
int lut_initialized = 0;

void init_fft_luts() {
    for (int i = 0; i < MAX_FFT_SIZE / 2; i++) {
        float angle = -2.0f * PI * i / MAX_FFT_SIZE;
        twiddle_LUT[i] = cosf(angle) + I * sinf(angle);
    }
    lut_initialized = 1;
}
// ----------------------------------


/***************************************************************************
 * TIMING STRUCTURE
 *
 * This version uses a consolidated two-block timing model per waveform:
 *
 * TX  (one-time build, measured once before MC loop):
 * GOLAY TX = generate_golay() + build_golay_tx()
 * ZC    TX = generate_zc()    + build_zc_tx()
 *
 * CHANNEL (per MC run, summed over ALL pulses):
 * GOLAY CH = delay_signal()           [propagation delay]
 * + phase/shift computation  [Doppler prep]
 * + rx[i] *= shift loop      [Doppler apply]
 * + awgn_measured()          [AWGN injection]
 * -- measured for Golay A and Golay B combined --
 *
 * ZC CH    = same four steps as above for ZC pulses
 *
 * RX (per MC run, covers everything after channel injection):
 * GOLAY RX = correlate_fft() + cube store       [matched filter]
 * + compute_range_doppler()            [2-D FFT / RD map]
 * + find_peak()                        [peak detection]
 * + r_axis[] + v_axis_g[] lookup       [range + vel estimation]
 * -- measured for Golay A and Golay B combined + post-processing --
 *
 * ZC RX    = rx_win[] = rx[] loop               [remove cyclic prefix]
 * + correlate_circular_zc() + cube store[circular correlation]
 * + compute_range_doppler()            [2-D FFT / RD map]
 * + find_peak()                        [peak detection]
 * + r_axis[] + v_axis_z[] lookup       [range + vel estimation]
 *
 * All per-SNR values = total time across NUM_MC runs / NUM_MC (per-MC average)
 ***************************************************************************/

// TX build times — measured once as one-time setup cost
float golay_tx_build_time = 0.0f;  // generate_golay() + build_golay_tx()
float zc_tx_build_time    = 0.0f;  // generate_zc()    + build_zc_tx()

// Per-SNR averages (avg per MC run at each SNR level), all in microseconds
float snr_g_ch   [NUM_SNR_LEVELS];   // Golay channel avg per MC run
float snr_g_rx   [NUM_SNR_LEVELS];   // Golay RX avg per MC run
float snr_g_total[NUM_SNR_LEVELS];   // Golay total (CH + RX) avg per MC run

float snr_z_ch   [NUM_SNR_LEVELS];   // ZC channel avg per MC run
float snr_z_rx   [NUM_SNR_LEVELS];   // ZC RX avg per MC run
float snr_z_total[NUM_SNR_LEVELS];   // ZC total (CH + RX) avg per MC run


/***************************************************************************
 * RANGE AXIS COMPUTATION
 *
 * Fills r_axis[] so that r_axis[i] = i * (c / (2 * fs)).
 * Maps each fast-time sample index to the corresponding range in metres.
 * Resolution = c / (2 * fs) ≈ 0.5 m for fs = 300 MHz.
 ***************************************************************************/
void compute_range_axis()
{
    for(int i = 0; i < NSAMP; i++)
    {
        r_axis[i] = i * (C / (2.0f * FS));
    }
}


/***************************************************************************
 * GOLAY SEQUENCE GENERATION
 *
 * Constructs a Golay complementary pair (ga, gb) of length NSAMP
 * using the iterative doubling construction:
 *
 * ga_new = [ga_old | gb_old]
 * gb_new = [ga_old | -gb_old]
 *
 * Starting seeds:  ga = [1, 1],  gb = [1, -1]
 * After log2(NSAMP) iterations the sequences reach length NSAMP.
 *
 * Key property: |FFT(ga)|^2 + |FFT(gb)|^2 = 2*N (flat spectral energy).
 ***************************************************************************/
void generate_golay()
{
    // Seed the two-element Golay pair
    ga[0] = 1.0f + 0.0f*I;
    ga[1] = 1.0f + 0.0f*I;
    gb[0] = 1.0f + 0.0f*I;
    gb[1] = -1.0f + 0.0f*I;
    int length = 2;

    // Iteratively double until length reaches NSAMP
    while(length < NSAMP)
    {
        float complex ga_old[NSAMP];
        float complex gb_old[NSAMP];

        // Copy old values before overwriting
        for(int i = 0; i < length; i++)
        {
            ga_old[i] = ga[i];
            gb_old[i] = gb[i];
        }

        // Build new sequences using the doubling rule
        for(int i = 0; i < length; i++)
        {
            ga[i]          = ga_old[i];   // first half  of ga_new = ga_old
            ga[i + length] = gb_old[i];   // second half of ga_new = gb_old
            gb[i]          = ga_old[i];   // first half  of gb_new = ga_old
            gb[i + length] = -gb_old[i];  // second half of gb_new = -gb_old
        }

        length *= 2;
    }
}


/***************************************************************************
 * BUILD GOLAY TX SIGNAL
 *
 * Builds transmission buffers tx_ga and tx_gb by appending NSAMP zeros
 * after each Golay sequence (guard interval):
 *
 * tx_ga = [ga (NSAMP)] [zeros (NSAMP)]   total L_golay = 2*NSAMP samples
 * tx_gb = [gb (NSAMP)] [zeros (NSAMP)]
 *
 * Equivalent to MATLAB: tx_ga = [ga; zeros(NSAMP,1)]
 ***************************************************************************/
void build_golay_tx()
{
    // First half = original sequence
    for(int i = 0; i < NSAMP; i++)
    {
        tx_ga[i] = ga[i];
        tx_gb[i] = gb[i];
    }

    // Second half = zeros (guard interval)
    for(int i = NSAMP; i < 2 * NSAMP; i++)
    {
        tx_ga[i] = 0.0f + 0.0f * I;
        tx_gb[i] = 0.0f + 0.0f * I;
    }

    // Store total pulse length
    L_golay = 2 * NSAMP;
}


/***************************************************************************
 * ZADOFF-CHU GENERATION
 *
 * Generates a Zadoff-Chu (ZC) sequence of length NSAMP with root index u:
 *
 * zc[n] = exp(-j * pi * u * n^2 / NSAMP),  n = 0, 1, …, NSAMP-1
 *
 * ZC sequences have constant amplitude and ideal periodic autocorrelation.
 * Root index u must be coprime with NSAMP (typically u = 1).
 ***************************************************************************/
void generate_zc(int u)
{
    for(int n = 0; n < NSAMP; n++)
    {
        float phase = -PI * u * n * n / NSAMP;
        zc[n] = cosf(phase) + I * sinf(phase);
    }
}


/***************************************************************************
 * BUILD ZC TX SIGNAL
 *
 * Builds the ZC transmission buffer by prepending a full cyclic prefix:
 *
 * tx_zc = [zc (NSAMP)] [zc (NSAMP)]   total L_zc = 2*NSAMP samples
 *
 * The receiver discards the first NSAMP samples (cyclic prefix) and
 * processes only the second NSAMP samples as a circular convolution,
 * which is efficiently computed via FFT-based circular correlation.
 *
 * Equivalent to MATLAB: tx_zc = [zc; zc]   (second copy = cyclic prefix)
 ***************************************************************************/
void build_zc_tx()
{
    // First half = ZC sequence
    for(int i = 0; i < NSAMP; i++)
        tx_zc[i] = zc[i];

    // Second half (cyclic prefix — copy of ZC sequence)
    for(int i = 0; i < NSAMP; i++)
        tx_zc[i + NSAMP] = zc[i];

    L_zc = 2 * NSAMP;
}


/***************************************************************************
 * VELOCITY AXIS COMPUTATION
 *
 * Computes three arrays for each waveform:
 * 1. doppler_bins   — symmetric bin indices: (-N/2 to N/2-1)
 * equivalent to MATLAB: (-N/2 : N/2-1)
 * 2. doppler_axis   — Doppler frequency:  fd = bin * (PRF / N)  [Hz]
 * 3. v_axis         — radial velocity:    v  = fd * (lambda / 2) [m/s]
 *
 * Uses MAX_PULSES sized arrays globally;
 * fills up to TOTAL_PULSES_GOLAY or TOTAL_PULSES_ZC entries only.
 ***************************************************************************/
void compute_velocity_axes(int TOTAL_PULSES_GOLAY,
                           int TOTAL_PULSES_ZC,
                           float PRF_GOLAY,
                           float PRF_ZC)
{
    float lam = C / FC;  // Wavelength = c / fc

    /**************** GOLAY DOPPLER AXIS ****************/
    for(int i = 0; i < TOTAL_PULSES_GOLAY; i++)
    {
        // Symmetric Doppler bins (-N/2 to N/2-1)
        doppler_bins_g[i] = i - (TOTAL_PULSES_GOLAY / 2);

        // Doppler frequency: bin * PRF / N
        doppler_axis_g[i] = doppler_bins_g[i] * (PRF_GOLAY / TOTAL_PULSES_GOLAY);

        // Radial velocity: v = (lambda/2) * fd
        v_axis_g[i] = doppler_axis_g[i] * (lam / 2.0f);
    }

    /**************** ZC DOPPLER AXIS ****************/
    for(int i = 0; i < TOTAL_PULSES_ZC; i++)
    {
        doppler_bins_z[i] = i - (TOTAL_PULSES_ZC / 2);
        doppler_axis_z[i] = doppler_bins_z[i] * (PRF_ZC / TOTAL_PULSES_ZC);
        v_axis_z[i]       = doppler_axis_z[i] * (lam / 2.0f);
    }
}


/***************************************************************************
 * DELAY FUNCTION  (equivalent to MATLAB delayseq)
 *
 * Simulates propagation delay by shifting tx[] by delay_samples positions
 * and writing the result to rx[].
 *
 * rx[i] = tx[i - delay_samples]  for i >= delay_samples
 * rx[i] = 0                      for i <  delay_samples
 ***************************************************************************/
void delay_signal(const float complex *tx,
                  float complex *rx,
                  int length,
                  int delay_samples)
{
    int i;

    // Shift valid samples toward the end
    for(i = length - 1; i >= delay_samples; i--)
        rx[i] = tx[i - delay_samples];

    // Zero-fill the beginning (propagation delay region)
    for(i = 0; i < delay_samples; i++)
        rx[i] = 0.0f + 0.0f * I;
}


/***************************************************************************
 * AWGN GENERATION
 *
 * Generates a standard normal (Gaussian) random number using the
 * Box-Muller transform:
 *
 * Z = sqrt(-2 * ln(u1)) * cos(2*pi*u2)
 *
 * where u1, u2 are uniform(0,1] random numbers drawn from rand().
 ***************************************************************************/
float randn()
{
    float u1 = ((float)rand() + 1.0f) / ((float)RAND_MAX + 1.0f);
    float u2 = ((float)rand() + 1.0f) / ((float)RAND_MAX + 1.0f);

    return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * PI * u2);
}


/***************************************************************************
 * AWGN FUNCTION  (equivalent to MATLAB awgn(x, snr, 'measured'))
 *
 * Adds complex AWGN to signal[] so that the output SNR equals snr_db.
 *
 * Steps:
 * 1. Measure actual signal power: P_s = mean(|x|^2)
 * 2. Convert SNR from dB:         SNR_lin = 10^(snr_db/10)
 * 3. Compute noise power:         P_n = P_s / SNR_lin
 * 4. Per-component std deviation: sigma = sqrt(P_n / 2)
 * 5. Add independent Gaussian noise to real and imaginary parts.
 ***************************************************************************/
void awgn_measured(float complex *signal, int length, float snr_db)
{
    float signal_power = 0.0f;

    // Measure signal power
    for(int i = 0; i < length; i++)
    {
        signal_power += crealf(signal[i]) * crealf(signal[i])
                      + cimagf(signal[i]) * cimagf(signal[i]);
    }
    signal_power /= length;

    // Convert SNR from dB to linear
    float snr_linear  = powf(10.0f, snr_db / 10.0f);

    // Compute noise power and per-component standard deviation
    float noise_power = signal_power / snr_linear;
    float sigma       = sqrtf(noise_power / 2.0f); // each I/Q component gets half

    // Add independent complex Gaussian noise to each sample
    for(int i = 0; i < length; i++)
    {
        signal[i] += (sigma * randn()) + I * (sigma * randn());
    }
}


/***************************************************************************
 * FFT FUNCTION
 *
 * In-place Cooley-Tukey radix-2 DIT (Decimation-In-Time) FFT.
 * N must be a power of 2.
 *
 * Steps:
 * 1. Bit-reversal permutation of input array.
 * 2. Butterfly stages: for each stage length (2, 4, 8, … N),
 * apply DFT butterfly utilizing precomputed twiddle factors.
 ***************************************************************************/
void fft_complex(float complex *x, int N)
{
    // Ensure LUT is populated before calculating
    if (!lut_initialized) init_fft_luts();

    // Bit-reversal permutation
    int j = 0;
    for(int i = 1; i < N; i++)
    {
        int bit = N >> 1;
        while(j & bit) { j ^= bit; bit >>= 1; }
        j |= bit;
        if(i < j)
        {
            float complex temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
    }

    // FFT butterfly stages optimized with Twiddle LUT
    for(int len = 2; len <= N; len <<= 1)
    {
        int step = MAX_FFT_SIZE / len;
        for(int i = 0; i < N; i += len)
        {
            for(int j = 0; j < len/2; j++)
            {
                float complex u = x[i+j];
                float complex v = x[i+j+len/2] * twiddle_LUT[j * step];
                x[i+j]         = u + v;   // Butterfly sum
                x[i+j+len/2]   = u - v;   // Butterfly difference
            }
        }
    }
}


/***************************************************************************
 * IFFT FUNCTION
 *
 * In-place IFFT using the conjugate symmetry trick:
 * IFFT(x) = conj(FFT(conj(x))) / N
 *
 * Reuses the forward FFT implementation — no separate IFFT kernel needed.
 ***************************************************************************/
void ifft_complex(float complex *x, int N)
{
    // Conjugate input
    for(int i = 0; i < N; i++)
        x[i] = conjf(x[i]);

    // Forward FFT on conjugated input
    fft_complex(x, N);

    // Conjugate again and scale by 1/N
    for(int i = 0; i < N; i++)
        x[i] = conjf(x[i]) / N;
}


/***************************************************************************
 * FFT-BASED CROSS CORRELATION  (equivalent to MATLAB xcorr(a, b))
 *
 * Computes the full cross-correlation of a[] with b[] using:
 * out = IFFT( FFT(a) * conj(FFT(b)) )
 *
 * Output length = len_a + len_b - 1.
 * Both inputs are zero-padded to the next power-of-2 size N for efficiency.
 * Result is reordered to match MATLAB xcorr lag ordering (lag 0 at centre).
 ***************************************************************************/
void correlate_fft(const float complex *a, int len_a,
                   const float complex *b, int len_b,
                   float complex *out)
{
    int out_len = len_a + len_b - 1;

    // Find next power-of-2 >= out_len
    int N = 1;
    while(N < out_len) N <<= 1;

    static float complex A[4096];  // Zero-padded FFT buffer for signal a
    static float complex B[4096];  // Zero-padded FFT buffer for signal b

    // Zero padding
    for(int i = 0; i < N; i++)
    {
        A[i] = (i < len_a) ? a[i] : 0.0f + 0.0f*I;
        B[i] = (i < len_b) ? b[i] : 0.0f + 0.0f*I;
    }

    // Forward FFT of both signals
    fft_complex(A, N);
    fft_complex(B, N);

    // Multiply with conjugate (matched filter in frequency domain)
    for(int i = 0; i < N; i++)
        A[i] *= conjf(B[i]);

    // Inverse FFT to obtain correlation in time domain
    ifft_complex(A, N);

    // Reorder output to match MATLAB xcorr lag ordering
    int shift = len_b - 1;
    for(int i = 0; i < out_len; i++)
        out[i] = A[(i - shift + N) % N];
}


/***************************************************************************
 * FFTSHIFT  (equivalent to MATLAB fftshift for 1D)
 *
 * Circular shift of a 1-D complex array by N/2 positions (swap halves).
 * Centres the zero-Doppler bin in the middle of the Doppler axis
 * after taking the slow-time FFT.
 ***************************************************************************/
void fftshift_1d(float complex *x, int N)
{
    int half = N / 2;
    for(int i = 0; i < half; i++)
    {
        float complex temp = x[i];
        x[i]       = x[i + half];
        x[i + half] = temp;
    }
}


/***************************************************************************
 * RANGE-DOPPLER MAP COMPUTATION
 *
 * Computes the Range-Doppler map from a radar data cube by applying a
 * slow-time (Doppler) FFT across each range bin:
 *
 * For each range bin r:
 * 1. Extract slow-time vector: temp[] = cube[r][0 … num_pulses-1]
 * 2. Apply FFT along the pulse (slow-time) dimension.
 * 3. Scale by 1/num_pulses (match MATLAB convention).
 * 4. fftshift to centre zero-Doppler.
 * 5. Store result in rd[r][].
 *
 * Equivalent to MATLAB: fftshift(fft(cube,[],2),2) / TOTAL_PULSES
 ***************************************************************************/
void compute_range_doppler(float complex cube[][MAX_PULSES],
                             float complex rd[][MAX_PULSES],
                             int num_range_bins,
                             int num_pulses)
{
    float complex temp[MAX_PULSES];  // Temporary slow-time buffer for one range bin

    for(int r = 0; r < num_range_bins; r++)
    {
        // Copy slow-time vector for this range bin
        for(int p = 0; p < num_pulses; p++)
            temp[p] = cube[r][p];

        // Doppler FFT along the pulse dimension
        fft_complex(temp, num_pulses);

        // Scaling (match MATLAB — divides by N)
        for(int p = 0; p < num_pulses; p++)
            temp[p] /= num_pulses;

        // fftshift: centre zero-Doppler bin in the middle
        fftshift_1d(temp, num_pulses);

        // Store result back into the Range-Doppler map
        for(int p = 0; p < num_pulses; p++)
            rd[r][p] = temp[p];
    }
}


/***************************************************************************
 * PEAK DETECTION IN RANGE-DOPPLER MAP
 *
 * Locates the peak magnitude cell in the Range-Doppler map rd[][].
 * Returns the (range_idx, doppler_idx) of the maximum |rd[r][d]|.
 *
 * Used to estimate target range (via r_axis[range_idx]) and
 * radial velocity (via v_axis[doppler_idx]).
 ***************************************************************************/
void find_peak(float complex rd[][MAX_PULSES],
               int num_range,
               int num_pulses,
               int *range_idx,
               int *doppler_idx)
{
    float max_val = 0.0f;

    for(int r = 0; r < num_range; r++)
    {
        for(int d = 0; d < num_pulses; d++)
        {
            float mag = cabsf(rd[r][d]);
            if(mag > max_val)
            {
                max_val      = mag;
                *range_idx   = r;
                *doppler_idx = d;
            }
        }
    }
}


/***************************************************************************
 * ZC CIRCULAR CORRELATION
 *
 * Computes the circular cross-correlation of rx_in[] with the ZC
 * reference tx_in[] using the FFT:
 *
 * out = IFFT( FFT(rx_in) * conj(FFT(tx_in)) )
 *
 * Equivalent to MATLAB: ifft(fft(rx_win) .* conj(fft(zc))) / nsamp
 *
 * The cyclic structure of ZC sequences ensures an ideal impulse response
 * (Kronecker delta) at the correct delay bin.
 * Note: ifft_complex already divides by N, so no extra scaling is needed.
 ***************************************************************************/
void correlate_circular_zc(const float complex *rx_in,
                            const float complex *tx_in,
                            float complex *out, int N)
{
    // Local FFT buffers to avoid corrupting global arrays
    float complex RX[NSAMP], TX[NSAMP];

    for(int i = 0; i < N; i++) { RX[i] = rx_in[i]; TX[i] = tx_in[i]; }

    // Forward FFT of received and reference signals
    fft_complex(RX, N);
    fft_complex(TX, N);

    // Matched filter: multiply received spectrum by conjugate of reference
    for(int i = 0; i < N; i++)
        RX[i] *= conjf(TX[i]);

    // Inverse FFT to obtain circular correlation output
    ifft_complex(RX, N);   // ifft_complex already divides by N

    // Copy result to output buffer
    for(int i = 0; i < N; i++)
        out[i] = RX[i];
}


/***************************************************************************
 * MAIN
 ***************************************************************************/
int main()
{
    init_platform();
    init_fft_luts(); // Precalculate FFT twiddle factors once

    // -------------------------------------------------------
    // TX BUILD TIMING
    //
    // The TX waveforms are fixed for the entire simulation.
    // They are built once here and timed as a one-time setup cost,
    // separate from the per-MC channel and RX processing times.
    //
    //   GOLAY TX build = generate_golay() + build_golay_tx()
    //                    (range axis also computed here as setup)
    //   ZC    TX build = generate_zc()    + build_zc_tx()
    // -------------------------------------------------------
    XTime t0, t1;

    // Golay TX build time (includes range axis setup)
    XTime_GetTime(&t0);
    compute_range_axis();   // Fill r_axis[] with range bin values (metres)
    generate_golay();       // Build Golay complementary pair ga[], gb[]
    build_golay_tx();       // Append zero-guard to build tx_ga[], tx_gb[]
    XTime_GetTime(&t1);
    golay_tx_build_time = TIME_US(t0, t1);

    // ZC TX build time
    XTime_GetTime(&t0);
    generate_zc(1);         // Build ZC sequence with root index u = 1
    build_zc_tx();          // Append cyclic prefix to build tx_zc[]
    XTime_GetTime(&t1);
    zc_tx_build_time = TIME_US(t0, t1);

    printf("Hello World\n\r");

    // -------------------------------------------------------
    // FIFO / TIMING PARAMETERS
    //
    // How many pulses fit per FIFO batch:
    //   PULSES_PER_BATCH = floor(FIFO_SIZE_SAMPLES / SAMPLES_PER_PULSE)
    //
    // Total coherent processing interval (CPI) pulses:
    //   TOTAL_PULSES = PULSES_PER_BATCH * NUM_BATCHES
    //
    // PRI (Pulse Repetition Interval) and PRF (Pulse Repetition Freq):
    //   PRI = pulse_length / fs  [seconds]
    //   PRF = 1 / PRI            [Hz]
    // -------------------------------------------------------

    // Samples per pulse (same as pulse length for each waveform)
    int SAMPLES_PER_PULSE_GOLAY = L_golay;
    int SAMPLES_PER_PULSE_ZC    = L_zc;

    // Pulses per FIFO batch — floor division automatic in integer division
    int PULSES_PER_BATCH_GOLAY  = FIFO_SIZE_SAMPLES / SAMPLES_PER_PULSE_GOLAY;
    int PULSES_PER_BATCH_ZC     = FIFO_SIZE_SAMPLES / SAMPLES_PER_PULSE_ZC;

    // Total pulses across all batches (full CPI)
    int TOTAL_PULSES_GOLAY      = PULSES_PER_BATCH_GOLAY * NUM_BATCHES;
    int TOTAL_PULSES_ZC         = PULSES_PER_BATCH_ZC    * NUM_BATCHES;

    // PRI and PRF for each waveform
    float PRI_GOLAY = (float)L_golay / FS;   // Golay PRI in seconds
    float PRF_GOLAY = 1.0f / PRI_GOLAY;      // Golay PRF in Hz
    float PRI_ZC    = (float)L_zc    / FS;   // ZC PRI in seconds
    float PRF_ZC    = 1.0f / PRI_ZC;         // ZC PRF in Hz

    // Compute Doppler and velocity axes for both waveforms
    compute_velocity_axes(TOTAL_PULSES_GOLAY, TOTAL_PULSES_ZC,
                          PRF_GOLAY, PRF_ZC);

    // -------------------------------------------------------
    // INITIALISE RMSE ARRAYS
    // -------------------------------------------------------
    for(int i = 0; i < NUM_SNR_LEVELS; i++)
    {
        rmse_r_g[i] = 0.0f; rmse_v_g[i] = 0.0f;
        rmse_r_z[i] = 0.0f; rmse_v_z[i] = 0.0f;
    }

    // -------------------------------------------------------
    // INITIALISE PER-SNR TIMING ARRAYS
    // -------------------------------------------------------
    for(int i = 0; i < NUM_SNR_LEVELS; i++)
    {
        snr_g_ch[i]    = 0.0f;
        snr_g_rx[i]    = 0.0f;
        snr_g_total[i] = 0.0f;
        snr_z_ch[i]    = 0.0f;
        snr_z_rx[i]    = 0.0f;
        snr_z_total[i] = 0.0f;
    }

    printf("\n\n ================= MONTE CARLO SIMULATION =================\n\n");
    printf("-------------------------------------------------------------------\n");
    printf("| SNR | Golay R RMSE | Golay V RMSE | ZC R RMSE | ZC V RMSE |\n");
    printf("-------------------------------------------------------------------\n");

    // Wall-clock timer wrapping the entire SNR / MC double loop
    XTime PS_Start, PS_End;
    XTime_SetTime(0);
    XTime_GetTime(&PS_Start);

    // =========================================================
    // MONTE CARLO SIMULATION
    // Outer loop: SNR levels from SNR_START to SNR_END
    // Inner loop: NUM_MC independent Monte Carlo trials per SNR
    // =========================================================
    for(int s = 0; s < NUM_SNR_LEVELS; s++)
    {
        // Current SNR value in dB
        int snr_db = SNR_START + s;

        // Error accumulators for RMSE at this SNR level
        float err_rg = 0.0f, err_vg = 0.0f;  // Golay range/velocity squared errors
        float err_rz = 0.0f, err_vz = 0.0f;  // ZC range/velocity squared errors

        /*
         * Per-SNR timing accumulators.
         * Sum of all MC run times at this SNR level.
         * Divided by NUM_MC at end of SNR loop → per-MC average.
         */
        float acc_g_ch = 0.0f;   // Golay channel total across all MC runs at this SNR
        float acc_g_rx = 0.0f;   // Golay RX total across all MC runs at this SNR
        float acc_z_ch = 0.0f;   // ZC channel total across all MC runs at this SNR
        float acc_z_rx = 0.0f;   // ZC RX total across all MC runs at this SNR

        // ---- MC LOOP ----
        for(int mc = 0; mc < NUM_MC; mc++)
        {
            // Generate random target range (50 to 200 meters)
            float true_r = 50.0f   + 150.0f * ((float)rand() / RAND_MAX);

            // Generate random target velocity (-100 to +100 m/s)
            float true_v = -100.0f + 200.0f * ((float)rand() / RAND_MAX);

            /*
             * Per-MC-run timing accumulators.
             * Accumulate over all pulses within this single MC run.
             * Added into acc_xxx at end of MC run.
             *
             * mc_g_ch: total Golay channel time for this MC run
             * (summed over all Golay A + Golay B pulses)
             * mc_g_rx: total Golay RX time for this MC run
             * (correlation for all pulses + post-processing once)
             * mc_z_ch: total ZC channel time for this MC run
             * mc_z_rx: total ZC RX time for this MC run
             * (CP removal + correlation for all pulses + post-processing once)
             */
            float mc_g_ch = 0.0f;
            float mc_g_rx = 0.0f;
            float mc_z_ch = 0.0f;
            float mc_z_rx = 0.0f;

            // Clear Golay radar data cube for this MC run
            for(int i = 0; i < NSAMP; i++)
                for(int j = 0; j < TOTAL_PULSES_GOLAY; j++)
                    cube_g[i][j] = 0.0f + I*0.0f;

            // Clear ZC radar data cube for this MC run
            for(int i = 0; i < NSAMP; i++)
                for(int j = 0; j < TOTAL_PULSES_ZC; j++)
                    cube_z[i][j] = 0.0f + I*0.0f;

            // Slow-time counters (track absolute time for Doppler phase)
            float trigger_time_g = 0.0f;  // Golay slow-time counter (seconds)
            float trigger_time_z = 0.0f;  // ZC slow-time counter (seconds)

            // Pulse indices into the radar data cubes
            int pulse_idx_g = 0;
            int pulse_idx_z = 0;

            // Compute Doppler frequency for this target (constant per MC run)
            float fd = (2.0f * true_v * FC) / C;

            // ==================================================
            // FIFO / BATCH LOOP
            //
            // Each batch processes (PULSES_PER_BATCH_GOLAY / 2) Golay
            // complementary pairs and the same number of ZC pulses,
            // filling one FIFO's worth of data at a time.
            //
            // Channel and RX pulse-level times accumulate into the
            // mc_xxx accumulators across all batches.
            // ==================================================
            for(int b = 0; b < NUM_BATCHES; b++)
            {
                for(int p = 0; p < (PULSES_PER_BATCH_GOLAY / 2); p++)
                {
                    /******************* GOLAY A *******************/

                    float t_now       = trigger_time_g;

                    // Instantaneous range (moving target: r = r0 + v*t)
                    float r_now       = true_r + true_v * t_now;

                    // Round-trip delay in samples: tau = 2*r/c, delay = tau * fs
                    float tau         = (2.0f * r_now) / C;
                    int delay_samples = (int)roundf(tau * FS);

                    // CHANNEL GOLAY A:
                    // All four steps timed as one unbroken block:
                    //   delay_signal + doppler prep + doppler apply + awgn
                    XTime_GetTime(&t0);
                    delay_signal(tx_ga, rx, L_golay, delay_samples);
                    float phase = 2.0f * PI * fd * t_now;
                    float complex shift = cosf(phase) + I * sinf(phase);
                    for(int i = 0; i < L_golay; i++)
                        rx[i] *= shift;
                    awgn_measured(rx, L_golay, snr_db);
                    XTime_GetTime(&t1);
                    mc_g_ch += TIME_US(t0, t1);

                    // RX GOLAY A:
                    // correlate_fft + store NSAMP valid range bins into cube
                    // (valid lags start at index L_golay - 1 in full xcorr output)
                    XTime_GetTime(&t0);
                    correlate_fft(rx, L_golay, tx_ga, L_golay, corr);
                    for(int i = 0; i < NSAMP; i++)
                        cube_g[i][pulse_idx_g] = corr[L_golay - 1 + i];
                    XTime_GetTime(&t1);
                    mc_g_rx += TIME_US(t0, t1);

                    pulse_idx_g++;
                    trigger_time_g += PRI_GOLAY;  // Advance slow-time clock


                    /******************* GOLAY B *******************/

                    t_now         = trigger_time_g;
                    r_now         = true_r + true_v * t_now;
                    tau           = (2.0f * r_now) / C;
                    delay_samples = (int)roundf(tau * FS);

                    // CHANNEL GOLAY B:
                    // delay_signal + doppler prep + doppler apply + awgn
                    XTime_GetTime(&t0);
                    delay_signal(tx_gb, rx, L_golay, delay_samples);
                    phase = 2.0f * PI * fd * t_now;
                    shift = cosf(phase) + I * sinf(phase);
                    for(int i = 0; i < L_golay; i++)
                        rx[i] *= shift;
                    awgn_measured(rx, L_golay, snr_db);
                    XTime_GetTime(&t1);
                    mc_g_ch += TIME_US(t0, t1);

                    // RX GOLAY B:
                    // correlate_fft + cube store
                    XTime_GetTime(&t0);
                    correlate_fft(rx, L_golay, tx_gb, L_golay, corr);
                    for(int i = 0; i < NSAMP; i++)
                        cube_g[i][pulse_idx_g] = corr[L_golay - 1 + i];
                    XTime_GetTime(&t1);
                    mc_g_rx += TIME_US(t0, t1);

                    pulse_idx_g++;
                    trigger_time_g += PRI_GOLAY;


                    /******************* ZC *******************/

                    t_now         = trigger_time_z;
                    r_now         = true_r + true_v * t_now;
                    tau           = (2.0f * r_now) / C;
                    delay_samples = (int)roundf(tau * FS);

                    // CHANNEL ZC:
                    // delay_signal + doppler prep + doppler apply + awgn
                    XTime_GetTime(&t0);
                    delay_signal(tx_zc, rx, L_zc, delay_samples);
                    phase = 2.0f * PI * fd * t_now;
                    shift = cosf(phase) + I * sinf(phase);
                    for(int i = 0; i < L_zc; i++)
                        rx[i] *= shift;
                    awgn_measured(rx, L_zc, snr_db);
                    XTime_GetTime(&t1);
                    mc_z_ch += TIME_US(t0, t1);

                    // RX ZC:
                    // remove cyclic prefix (keep second NSAMP samples only)
                    // + circular correlation with ZC reference + cube store
                    XTime_GetTime(&t0);
                    for(int i = 0; i < NSAMP; i++)
                        rx_win[i] = rx[NSAMP + i];
                    correlate_circular_zc(rx_win, zc, circ_out, NSAMP);
                    for(int i = 0; i < NSAMP; i++)
                        cube_z[i][pulse_idx_z] = circ_out[i];
                    XTime_GetTime(&t1);
                    mc_z_rx += TIME_US(t0, t1);

                    pulse_idx_z++;
                    trigger_time_z += PRI_ZC;  // Advance slow-time clock

                } // end inner pulse loop (p)
            } // end batch loop (b)


            // ==================================================
            // POST-PROCESSING  (once per MC run)
            //
            // compute_range_doppler + find_peak + axis lookups
            // are all part of the RX timing block (timed together
            // in one unbroken block per waveform).
            // ==================================================

            // RX GOLAY POST-PROCESSING:
            // 2-D RD map (slow-time FFT) + peak detection + range and velocity estimation
            XTime_GetTime(&t0);
            compute_range_doppler(cube_g, rd_map_g, NSAMP, TOTAL_PULSES_GOLAY);
            int r_g_idx, v_g_idx;
            find_peak(rd_map_g, NSAMP, TOTAL_PULSES_GOLAY, &r_g_idx, &v_g_idx);
            float est_r_g = r_axis[r_g_idx];      // Map range bin index to metres
            float est_v_g = v_axis_g[v_g_idx];    // Map Doppler bin index to m/s
            XTime_GetTime(&t1);
            mc_g_rx += TIME_US(t0, t1);

            // RX ZC POST-PROCESSING:
            // 2-D RD map (slow-time FFT) + peak detection + range and velocity estimation
            XTime_GetTime(&t0);
            compute_range_doppler(cube_z, rd_map_z, NSAMP, TOTAL_PULSES_ZC);
            int r_z_idx, v_z_idx;
            find_peak(rd_map_z, NSAMP, TOTAL_PULSES_ZC, &r_z_idx, &v_z_idx);
            float est_r_z = r_axis[r_z_idx];      // Map range bin index to metres
            float est_v_z = v_axis_z[v_z_idx];    // Map Doppler bin index to m/s
            XTime_GetTime(&t1);
            mc_z_rx += TIME_US(t0, t1);

            // Flush this MC run's totals into the per-SNR accumulators
            acc_g_ch += mc_g_ch;
            acc_g_rx += mc_g_rx;
            acc_z_ch += mc_z_ch;
            acc_z_rx += mc_z_rx;

            /*
             * Accumulate squared estimation errors for RMSE:
             * error = (estimate - truth)^2
             * Range error in metres^2; velocity error in (m/s)^2.
             */
            err_rg += powf(est_r_g - true_r, 2);
            err_vg += powf(est_v_g - true_v, 2);
            err_rz += powf(est_r_z - true_r, 2);
            err_vz += powf(est_v_z - true_v, 2);

        } // end MC loop

        // After all Monte Carlo runs for this SNR level, compute RMSE
        // RMSE = sqrt(mean squared error over all MC runs)
        rmse_r_g[s] = sqrtf(err_rg / NUM_MC);  // Golay range RMSE (metres)
        rmse_v_g[s] = sqrtf(err_vg / NUM_MC);  // Golay velocity RMSE (m/s)
        rmse_r_z[s] = sqrtf(err_rz / NUM_MC);  // ZC range RMSE (metres)
        rmse_v_z[s] = sqrtf(err_vz / NUM_MC);  // ZC velocity RMSE (m/s)

        // Print results for this SNR level
        printf("| %3d | %12.4f | %12.4f | %10.4f | %10.4f |\n",
               snr_db, rmse_r_g[s], rmse_v_g[s], rmse_r_z[s], rmse_v_z[s]);

        // Divide by NUM_MC to get per-MC averages for this SNR level
        snr_g_ch[s]    = acc_g_ch / NUM_MC;
        snr_g_rx[s]    = acc_g_rx / NUM_MC;
        snr_g_total[s] = snr_g_ch[s] + snr_g_rx[s];  // CH + RX = total per MC

        snr_z_ch[s]    = acc_z_ch / NUM_MC;
        snr_z_rx[s]    = acc_z_rx / NUM_MC;
        snr_z_total[s] = snr_z_ch[s] + snr_z_rx[s];

    } // end SNR loop

    printf("-------------------------------------------------------------------\n");

    XTime_GetTime(&PS_End);
    float time_total = TIME_US(PS_Start, PS_End);  // Total wall-clock time in us

    // =========================================================
    // COMPUTE OVERALL AVERAGES ACROSS ALL SNR LEVELS
    // Sum per-SNR averages, then divide by NUM_SNR_LEVELS
    // to get the grand average per MC run across all SNR levels.
    // =========================================================
    float overall_g_ch    = 0.0f;
    float overall_g_rx    = 0.0f;
    float overall_g_total = 0.0f;
    float overall_z_ch    = 0.0f;
    float overall_z_rx    = 0.0f;
    float overall_z_total = 0.0f;

    // Sum across all SNR levels
    for(int s = 0; s < NUM_SNR_LEVELS; s++)
    {
        overall_g_ch    += snr_g_ch[s];
        overall_g_rx    += snr_g_rx[s];
        overall_g_total += snr_g_total[s];
        overall_z_ch    += snr_z_ch[s];
        overall_z_rx    += snr_z_rx[s];
        overall_z_total += snr_z_total[s];
    }

    // Divide by NUM_SNR_LEVELS for overall average across all SNR levels
    overall_g_ch    /= NUM_SNR_LEVELS;
    overall_g_rx    /= NUM_SNR_LEVELS;
    overall_g_total /= NUM_SNR_LEVELS;
    overall_z_ch    /= NUM_SNR_LEVELS;
    overall_z_rx    /= NUM_SNR_LEVELS;
    overall_z_total /= NUM_SNR_LEVELS;

    // =========================================================
    // PRINT: TX BUILD TIMES  (one-time setup cost)
    // =========================================================
    printf("\n=================================================================\n");
    printf("  TX BUILD TIMES  (one-time setup)\n");
    printf("=================================================================\n");
    printf("| Golay TX (generate_golay + build_golay_tx) : %10.4f us |\n",
           golay_tx_build_time);
    printf("| ZC    TX (generate_zc   + build_zc_tx)     : %10.4f us |\n",
           zc_tx_build_time);
    printf("=================================================================\n");

    // =========================================================
    // PRINT: PER-SNR TIMING TABLE
    // CH    = delay + doppler + awgn  (all pulses in one MC run)
    // RX    = correlation + RD map + peak + range/vel estimation
    // =========================================================
    printf("\n=================================================================\n");
    printf("  GOLAY: PER-SNR TIMING  (avg per MC run, us)\n");
    printf("  CH    = delay + doppler + awgn  (all pulses in MC run)\n");
    printf("  RX    = correlation + RD map + peak + range/vel estimation\n");
    printf("=================================================================\n");
    printf("| %4s | %12s | %12s | %12s |\n",
           "SNR", "CH (us)", "RX (us)", "TOTAL (us)");
    printf("|------|--------------|--------------|-------------- |\n");
    for(int s = 0; s < NUM_SNR_LEVELS; s++)
    {
        printf("| %4d | %12.4f | %12.4f | %12.4f |\n",
               SNR_START + s,
               snr_g_ch[s], snr_g_rx[s], snr_g_total[s]);
    }

    printf("\n=================================================================\n");
    printf("  ZADOFF-CHU: PER-SNR TIMING  (avg per MC run, us)\n");
    printf("  CH    = delay + doppler + awgn  (all pulses in MC run)\n");
    printf("  RX    = cp removal + correlation + RD map + peak + range/vel\n");
    printf("=================================================================\n");
    printf("| %4s | %12s | %12s | %12s |\n",
           "SNR", "CH (us)", "RX (us)", "TOTAL (us)");
    printf("|------|--------------|--------------|---------------|\n");
    for(int s = 0; s < NUM_SNR_LEVELS; s++)
    {
        printf("| %4d | %12.4f | %12.4f | %12.4f |\n",
               SNR_START + s,
               snr_z_ch[s], snr_z_rx[s], snr_z_total[s]);
    }

    // =========================================================
    // PRINT: OVERALL TIMING REPORT
    // =========================================================
    printf("\n=================================================================\n");
    printf("  OVERALL TIMING REPORT  (avg per MC run, across all SNR, us)\n");
    printf("=================================================================\n");

    printf("\n+---------------------------------------------------------------+\n");
    printf("|                     GOLAY TIMING                               |\n");
    printf("+---------------------------------------------------------------+\n");
    printf("| TX  generate_golay + build_golay_tx                          |\n");
    printf("|                          : %10.4f us  (one-time)       |\n",
           golay_tx_build_time);
    printf("| CH  delay + doppler + awgn                                   |\n");
    printf("|     (all pulses per MC)  : %10.4f us  (avg per MC)     |\n",
           overall_g_ch);
    printf("| RX  corr + RD map + peak + range/vel                         |\n");
    printf("|     (all pulses + post)  : %10.4f us  (avg per MC)     |\n",
           overall_g_rx);
    printf("+---------------------------------------------------------------+\n");
    printf("| GOLAY TOTAL (CH + RX)    : %10.4f us  (avg per MC)     |\n",
           overall_g_total);
    printf("+---------------------------------------------------------------+\n");

    printf("\n+---------------------------------------------------------------+\n");
    printf("|                   ZADOFF-CHU TIMING                          |\n");
    printf("+---------------------------------------------------------------+\n");
    printf("| TX  generate_zc + build_zc_tx                                |\n");
    printf("|                          : %10.4f us  (one-time)       |\n",
           zc_tx_build_time);
    printf("| CH  delay + doppler + awgn                                   |\n");
    printf("|     (all pulses per MC)  : %10.4f us  (avg per MC)     |\n",
           overall_z_ch);
    printf("| RX  cp removal + corr + RD map + peak + range/vel            |\n");
    printf("|     (all pulses + post)  : %10.4f us  (avg per MC)     |\n",
           overall_z_rx);
    printf("+---------------------------------------------------------------+\n");
    printf("| ZC TOTAL (CH + RX)       : %10.4f us  (avg per MC)     |\n",
           overall_z_total);
    printf("+---------------------------------------------------------------+\n");

    printf("\n+---------------------------------------------------------------+\n");
    printf("|                     OVERALL SUMMARY                          |\n");
    printf("+---------------------------------------------------------------+\n");
    printf("| Golay TX Build           : %10.4f us  (one-time)       |\n",
           golay_tx_build_time);
    printf("| ZC    TX Build           : %10.4f us  (one-time)       |\n",
           zc_tx_build_time);
    printf("| Golay  Total (avg/MC)    : %10.4f us                   |\n",
           overall_g_total);
    printf("| ZC     Total (avg/MC)    : %10.4f us                   |\n",
           overall_z_total);
    printf("| Combined Golay+ZC(avg/MC): %10.4f us                   |\n",
           overall_g_total + overall_z_total);
    printf("| Total Run Time           : %.2f us                     |\n",
           time_total);
    printf("+---------------------------------------------------------------+\n");

    cleanup_platform();
    return 0;
}