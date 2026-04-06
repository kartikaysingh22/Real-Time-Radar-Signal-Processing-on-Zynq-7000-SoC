//parametrs check
//range axis comp check
// golay generation check
// zdoff chu check
// fifo timing parameter claculation check
// velocity axis generation done
// cross correlation functionality achieved
// fifo loop written
// for zc we need ifft and fft functions, update the zc part of fifo/batch loop
// now next step: ifft and fft functionality




#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
//#include "xil_printf.h"
//#include "platform.h"
//#include "xaxidma.h"
//#include "xparameters.h"
//#include "xtime_l.h"


/***************************************************************************
 * SYSTEM PARAMETERS
 ***************************************************************************/

#define C 3e8f        // Speed of light
#define FC 30e9f        // Carrier Frequency
#define FS 300e6f       // Samplig Frequency

#define NSAMP 512       //base sequence length
#define FIFO_SIZE_SAMPLES 32768     // Fifo memory Size in samples
#define NUM_BATCHES 16          //no of FIFO refills/batches 
#define NUM_MC 50           //no of Monte Carlo runs

#define SNR_START -45       //SNR Range
#define SNR_END   -30
#define NUM_SNR_LEVELS (SNR_END - SNR_START + 1)


#define MAX_PULSES 512


#define PI 3.14159265358979323846f

//GLOBAL MEMORY 

// Range axis
float r_axis[NSAMP];

// Radar data cube (Range × Pulses)
float complex cube_g[NSAMP][MAX_PULSES]; // Initialize Golay radar cube
float complex cube_z[NSAMP][MAX_PULSES]; // Initialize ZC radar cube


// Range-Doppler Map
float complex rd_map_g[NSAMP][MAX_PULSES];
float complex rd_map_z[NSAMP][MAX_PULSES];

// Golay sequences
float complex ga[NSAMP];
float complex gb[NSAMP];

// Zadoff-Chu sequence
float complex zc[NSAMP];

float complex tx_ga[2 * NSAMP];   //Zero-padded Golay A for transmission
float complex tx_gb[2 * NSAMP];   //Zero-padded Golay B for transmission
int L_golay; 

float complex tx_zc[2 * NSAMP]; // Zadoff chu signal for transmission
int L_zc;

float doppler_bins_g[MAX_PULSES];  // Doppler bins Golay
float doppler_axis_g[MAX_PULSES];   // Doppler frequency axis Golay
float v_axis_g[MAX_PULSES];         // Velocity axis Golay

float doppler_bins_z[MAX_PULSES];       // Doppler bins ZC
float doppler_axis_z[MAX_PULSES];       // Doppler frequency axis ZC
float v_axis_z[MAX_PULSES];             // Velocity axis ZC

float lambda = C / FC;   // Wavelength of Tx Wave


// RMSE storage arrays (one value per SNR level)
float rmse_r_g[NUM_SNR_LEVELS];
float rmse_v_g[NUM_SNR_LEVELS];
float rmse_r_z[NUM_SNR_LEVELS];
float rmse_v_z[NUM_SNR_LEVELS];

// Buffers for processing (to avoid dynamic allocation and stack overflow in C)
float complex rx[2 * NSAMP];
float complex corr[2 * NSAMP];
float complex rx_win[NSAMP];
float complex circ_out[NSAMP];



// RANGE AXIS COMPUTATION
void compute_range_axis()
{
    for(int i = 0; i < NSAMP; i++)
    {
        r_axis[i] = i * (C / (2.0f * FS));
    }
}


// GOLAY SEQUENCE GENERATION
void generate_golay()
{
    ga[0] = 1.0f + 0.0f*I;
    ga[1] = 1.0f + 0.0f*I;
    gb[0] = 1.0f + 0.0f*I;
    gb[1] = -1.0f + 0.0f*I;
    int length = 2; 

    while(length < NSAMP)
    {
        float complex ga_old[NSAMP];
        float complex gb_old[NSAMP];

        // Copy old values
        for(int i = 0; i < length; i++)
        {
            ga_old[i] = ga[i];
            gb_old[i] = gb[i];
        }

        // Build new sequences
        for(int i = 0; i < length; i++)
        {
            ga[i] = ga_old[i];
            ga[i + length] = gb_old[i];

            gb[i] = ga_old[i];
            gb[i + length] = -gb_old[i];
        }

        length *= 2;
    }
}


// building the Golay Signal for Tx
void build_golay_tx()
{
    // First half = original sequence
    for(int i = 0; i < NSAMP; i++)
    {
        tx_ga[i] = ga[i];
        tx_gb[i] = gb[i];

    }

    // Second half = zeros
    for(int i = NSAMP; i < 2 * NSAMP; i++)
    {
        tx_ga[i] = 0.0f + 0.0f * I;
        tx_gb[i] = 0.0f + 0.0f * I;

    }

    // Length of pulse
    L_golay = 2 * NSAMP;
}



// ZADOFF-CHU GENERATION
void generate_zc(int u)
{
    for(int n = 0; n < NSAMP; n++)
    {
        float phase = -PI * u * n * n / NSAMP;
        zc[n] = cosf(phase) + I * sinf(phase);
    }
}

// building the Zadoof Chu Signal for Tx
void build_zc_tx()
{
    // First half
    for(int i = 0; i < NSAMP; i++)
    {
        tx_zc[i] = zc[i];
    }

    // Second half (cyclic prefix)
    for(int i = 0; i < NSAMP; i++)
    {
        tx_zc[i + NSAMP] = zc[i];
    }

    L_zc = 2 * NSAMP;
}



/***************************************************************************
 * VELOCITY AXIS COMPUTATION
 *
 * This function computes:
 * 1) Doppler bin indices  (-N/2 to N/2-1)
 * 2) Doppler frequency axis
 * 3) Velocity axis
 *
 * We use MAX_PULSES sized arrays globally,
 * but only fill up to TOTAL_PULSES_GOLAY or TOTAL_PULSES_ZC.
 ***************************************************************************/
void compute_velocity_axes(int TOTAL_PULSES_GOLAY,
                           int TOTAL_PULSES_ZC,
                           float PRF_GOLAY,
                           float PRF_ZC)
{
    // Wavelength = c / fc
    float lambda = C / FC;

    /**************** GOLAY DOPPLER AXIS ****************/

    for(int i = 0; i < TOTAL_PULSES_GOLAY; i++)
    {
        // Create symmetric Doppler bins:
        // Equivalent to MATLAB:
        // (-N/2 : N/2-1)
        doppler_bins_g[i] = i - (TOTAL_PULSES_GOLAY / 2);

        // Convert bin index to Doppler frequency
        // // Delta_f = PRF / N
        doppler_axis_g[i] =
            doppler_bins_g[i] * (PRF_GOLAY / TOTAL_PULSES_GOLAY);

        // Convert Doppler frequency to velocity
        // v = (lambda/2) * fd
        v_axis_g[i] =
            doppler_axis_g[i] * (lambda / 2.0f);
    }

    /**************** ZC DOPPLER AXIS ****************/

    for(int i = 0; i < TOTAL_PULSES_ZC; i++)
    {
        doppler_bins_z[i] = i - (TOTAL_PULSES_ZC / 2);

        doppler_axis_z[i] =
            doppler_bins_z[i] * (PRF_ZC / TOTAL_PULSES_ZC);

        v_axis_z[i] =
            doppler_axis_z[i] * (lambda / 2.0f);
    }
}



/* ------------------------------------------------ */
/* Delay Function (like MATLAB delayseq)           */
/* rx = delayseq(tx, delay_samples)                */
/* ------------------------------------------------ */
void delay_signal(const float complex *tx,
                  float complex *rx,
                  int length,
                  int delay_samples)
{
    int i;

    /* Shift samples */
    for (i = length - 1; i >= delay_samples; i--)
    {
        rx[i] = tx[i - delay_samples];
    }

    /* Zero-fill beginning */
    for (i = 0; i < delay_samples; i++)
    {
        rx[i] = 0.0f + 0.0f * I;
    }
}



//////// AWGN Generation function
/* Generate standard normal random number using Box-Muller */
float randn()
{
    float u1 = ((float)rand() + 1.0f) / ((float)RAND_MAX + 1.0f);
    float u2 = ((float)rand() + 1.0f) / ((float)RAND_MAX + 1.0f);

    return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * PI * u2);
}


/* ------------------------------------------------ */
/* AWGN Function (like MATLAB awgn(x, snr, 'measured')) */
/* ------------------------------------------------ */
void awgn_measured(float complex *signal,
                   int length,
                   float snr_db)
{
    float signal_power = 0.0f;

    /* Measure signal power */
    for (int i = 0; i < length; i++)
    {
        float mag = crealf(signal[i]) * crealf(signal[i]) +
                    cimagf(signal[i]) * cimagf(signal[i]);
        signal_power += mag;
    }

    signal_power /= length;

    /* Convert SNR from dB to linear */
    float snr_linear = powf(10.0f, snr_db / 10.0f);

    /* Compute noise power */
    float noise_power = signal_power / snr_linear;

    /* Each complex component gets half */
    float sigma = sqrtf(noise_power / 2.0f);

    /* Add noise */
    for (int i = 0; i < length; i++)
    {
        float noise_real = sigma * randn();
        float noise_imag = sigma * randn();

        signal[i] += noise_real + I * noise_imag;
    }
}





// FFT Fuctions
void fft_complex(float complex *x, int N)
{
    // Bit reversal
    int j = 0;
    for (int i = 1; i < N; i++) // start from 1 since bit-reversal of 0 is 0
    {
        int bit = N >> 1;  // this is N/2, used to find the highest bit set in i because N is a power of 2
        while (j & bit) // while the current bit is set in j, we need to clear it and move to the next lower bit
        {
            j ^= bit;  // clear the bit ^= meaning toggle bit
            bit >>= 1;
        }
        j |= bit;

        if (i < j)
        {
            float complex temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
    }

    // FFT stages
    for (int len = 2; len <= N; len <<= 1)
    {
        float angle = -2.0f * PI / len;
        float complex wlen = cosf(angle) + I * sinf(angle);

        for (int i = 0; i < N; i += len)
        {
            float complex w = 1.0f + 0.0f * I;
            for (int j = 0; j < len/2; j++)
            {
                float complex u = x[i+j];
                float complex v = x[i+j+len/2] * w;

                x[i+j] = u + v;
                x[i+j+len/2] = u - v;

                w *= wlen;
            }
        }
    }
}

void ifft_complex(float complex *x, int N)
{
    // Conjugate input
    for (int i = 0; i < N; i++)
        x[i] = conjf(x[i]);

    fft_complex(x, N);

    // Conjugate again and scale
    for (int i = 0; i < N; i++)
        x[i] = conjf(x[i]) / N;
}




void correlate_fft(const float complex *a, int len_a,
                          const float complex *b, int len_b,
                          float complex *out)
{
    int out_len = len_a + len_b - 1;

    int N = 1;
    while (N < out_len)
        N <<= 1;

    static float complex A[4096];
    static float complex B[4096];

    // Zero padding
    for (int i = 0; i < N; i++)
    {
        A[i] = (i < len_a) ? a[i] : 0.0f + 0.0f*I;
        B[i] = (i < len_b) ? b[i] : 0.0f + 0.0f*I;
    }

    fft_complex(A, N);
    fft_complex(B, N);

    // Multiply with conjugate (matched filter)
    for (int i = 0; i < N; i++)
        A[i] *= conjf(B[i]);

    ifft_complex(A, N);

    // Now reorder to match MATLAB xcorr lag ordering
    int shift = len_b - 1;

    for (int i = 0; i < out_len; i++)
        out[i] = A[(i - shift + N) % N];
}

// Simple 1D FFT shift (swap halves)
void fftshift_1d(float complex *x, int N)
{
    int half = N / 2;

    for(int i = 0; i < half; i++)
    {
        float complex temp = x[i];
        x[i] = x[i + half];
        x[i + half] = temp;
    }
}



// Range-Doppler Map computation
void compute_range_doppler(float complex cube[][MAX_PULSES],
                           float complex rd[][MAX_PULSES],
                           int num_range_bins,
                           int num_pulses)
{
    float complex temp[MAX_PULSES];

    for(int r = 0; r < num_range_bins; r++)
    {
        // Copy slow-time vector
        for(int p = 0; p < num_pulses; p++)
            temp[p] = cube[r][p];

        // Doppler FFT
        fft_complex(temp, num_pulses);

        // Scaling (match MATLAB)
        for(int p = 0; p < num_pulses; p++)
            temp[p] /= num_pulses;

        // FFT shift (center zero Doppler)
        fftshift_1d(temp, num_pulses);

        // Store back
        for(int p = 0; p < num_pulses; p++)
            rd[r][p] = temp[p];
    }
}



// Peak finding in Range-Doppler map to get range and velocity estimates
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
                max_val = mag;
                *range_idx = r;
                *doppler_idx = d;
            }
        }
    }
}

void correlate_circular_zc(const float complex *rx, const float complex *tx,
                            float complex *out, int N)
{
    // Allocate local buffers
    float complex RX[NSAMP], TX[NSAMP];
    for(int i = 0; i < N; i++) { RX[i] = rx[i]; TX[i] = tx[i]; }

    fft_complex(RX, N);
    fft_complex(TX, N);

    for(int i = 0; i < N; i++)
        RX[i] *= conjf(TX[i]);

    ifft_complex(RX, N);

    // Scale by N to match MATLAB's ifft(...)./nsamp
    for(int i = 0; i < N; i++)
        out[i] = RX[i];  // ifft_complex already divides by N
}


int main()
{
    //init_platform();

    
    
    compute_range_axis();
    /*for (int i =0; i<NSAMP; i++){
        printf("%f\n",r_axis[i]);
    }*/
    generate_golay();
    build_golay_tx();
    printf("Hello World\n\r");
    /*for (int i =0; i<NSAMP; i++){
        printf("%d\n",ga[i]);
    }*/
    
    /*for (int i =0; i<NSAMP; i++){
        printf("%d\n",gb[i]);
    }*/
    int u =1;
    generate_zc(u);
    
    /*for (int i = 0; i < NSAMP; i++) {
        printf("%f + j%f\n", crealf(zc[i]), cimagf(zc[i]));

    }*/
    build_zc_tx();

    // Variables for FIFO / Timing Parameters
    int SAMPLES_PER_PULSE_GOLAY;
    int SAMPLES_PER_PULSE_ZC;

    int PULSES_PER_BATCH_GOLAY;
    int PULSES_PER_BATCH_ZC;

    int TOTAL_PULSES_GOLAY;
    int TOTAL_PULSES_ZC;

    float PRI_GOLAY;
    float PRF_GOLAY;

    float PRI_ZC;
    float PRF_ZC;

    // Samples per pulse
    SAMPLES_PER_PULSE_GOLAY = L_golay;
    SAMPLES_PER_PULSE_ZC    = L_zc;

    // floor division automatically happens in integer division
    PULSES_PER_BATCH_GOLAY = FIFO_SIZE_SAMPLES / SAMPLES_PER_PULSE_GOLAY;
    PULSES_PER_BATCH_ZC    = FIFO_SIZE_SAMPLES / SAMPLES_PER_PULSE_ZC;

    // Total pulses
    TOTAL_PULSES_GOLAY = PULSES_PER_BATCH_GOLAY * NUM_BATCHES;
    TOTAL_PULSES_ZC    = PULSES_PER_BATCH_ZC * NUM_BATCHES;

    // PRI and PRF
    PRI_GOLAY = (float)L_golay / FS;
    PRF_GOLAY = 1.0f / PRI_GOLAY;

    PRI_ZC = (float)L_zc / FS;
    PRF_ZC = 1.0f / PRI_ZC;

    // Compute Doppler & Velocity axes
    compute_velocity_axes(TOTAL_PULSES_GOLAY,
                      TOTAL_PULSES_ZC,
                      PRF_GOLAY,
                      PRF_ZC);

    /*for(int i = 0; i < 5; i++)
    {
        printf("v_axis_g[%d] = %f\n", i, v_axis_g[i]);
    }*/    




    //initializing the rmse arrays
    for(int i = 0; i < NUM_SNR_LEVELS; i++)
    {
        rmse_r_g[i] = 0.0f;
        rmse_v_g[i] = 0.0f;
        rmse_r_z[i] = 0.0f;
        rmse_v_z[i] = 0.0f;
    }

    printf("\n\n ================= MONTE CARLO SIMULATION =================\n\n");
    printf("-------------------------------------------------------------------\n");
    printf("| SNR | Golay R RMSE | Golay V RMSE | ZC R RMSE | ZC V RMSE |\n");
    printf("-------------------------------------------------------------------\n");



    // ================= MONTE CARLO SIMULATION =================
    for(int s = 0; s < NUM_SNR_LEVELS; s++)
    {
        // Current SNR value in dB
        int snr_db = SNR_START + s;

        // Error accumulators for this SNR level for Golay and ZC
        float err_rg = 0.0f;
        float err_vg = 0.0f;
        float err_rz = 0.0f;
        float err_vz = 0.0f;


        for(int mc = 0; mc < NUM_MC; mc++)
        {
            // Generate random target range (50 to 200 meters)
            float true_r = 50.0f + 150.0f * ((float)rand() / RAND_MAX);

            // Generate random velocity (-100 to +100 m/s)
            float true_v = -100.0f + 200.0f * ((float)rand() / RAND_MAX);


            // Clear Golay cube
            for(int i = 0; i < NSAMP; i++)
            {
                for(int j = 0; j < TOTAL_PULSES_GOLAY; j++)
                {
                    cube_g[i][j] = 0.0f + I*0.0f;
                }
            }

            // Clear ZC cube
            for(int i = 0; i < NSAMP; i++)
            {
                for(int j = 0; j < TOTAL_PULSES_ZC; j++)
                {
                    cube_z[i][j] = 0.0f + I*0.0f;
                }
            }

            // To track time for Doppler phase shift
            float trigger_time_g = 0.0f;        //  Golay slow-time counter
            float trigger_time_z = 0.0f;        //  ZC slow-time counter

            // pulse indexs
            int pulse_idx_g = 0;
            int pulse_idx_z = 0;
            


            ///////// FIFO/BATCH LOOP ///////
            // Compute Doppler frequency once per target
            float fd = (2.0f * true_v * FC) / C;

            //float complex rx[2 * NSAMP];
            //float complex corr[2 * L_golay - 1];  // full correlation output

            // Loop over FIFO batches
            for(int b = 0; b < NUM_BATCHES; b++)
            {
                for(int p = 0; p < (PULSES_PER_BATCH_GOLAY / 2); p++)
                {
                    /******************* GOLAY A *******************/
        
                    float t_now = trigger_time_g;

                    // Instantaneous range (moving target)
                    float r_now = true_r + true_v * t_now;

                    // Round trip delay
                    float tau = (2.0f * r_now) / C;
                    int delay_samples = (int)(tau * FS);

                    

                    // Delay signal
                    delay_signal(tx_ga, rx, L_golay, delay_samples);

                    // Apply Doppler
                    float phase = 2.0f * PI * fd * t_now;
                    float complex shift = cosf(phase) + I * sinf(phase);

                    for(int i = 0; i < L_golay; i++)
                        rx[i] *= shift;

                    // Add noise
                    awgn_measured(rx, L_golay, snr_db);

                    // Correlate
                    correlate_fft(rx, L_golay, tx_ga, L_golay, corr);

                    // Store only valid nsamp samples
                    for(int i = 0; i < NSAMP; i++)
                    {
                        cube_g[i][pulse_idx_g] = corr[L_golay - 1 + i];
                    }

                    pulse_idx_g++;
                    trigger_time_g += PRI_GOLAY;


                    /******************* GOLAY B *******************/

                    t_now = trigger_time_g;

                    r_now = true_r + true_v * t_now;
                    tau = (2.0f * r_now) / C;
                    delay_samples = (int)(tau * FS);

                    delay_signal(tx_gb, rx, L_golay, delay_samples);

                    phase = 2.0f * PI * fd * t_now;
                    shift = cosf(phase) + I * sinf(phase);

                    for(int i = 0; i < L_golay; i++)
                        rx[i] *= shift;

                    awgn_measured(rx, L_golay, snr_db);

                    correlate_fft(rx, L_golay, tx_gb, L_golay, corr);

                    for(int i = 0; i < NSAMP; i++)
                    {
                        cube_g[i][pulse_idx_g] = corr[L_golay - 1 + i];
                    }

                    pulse_idx_g++;
                    trigger_time_g += PRI_GOLAY;


                    /******************* ZC *******************/

                    t_now = trigger_time_z;

                    r_now = true_r + true_v * t_now;
                    tau = (2.0f * r_now) / C;
                    delay_samples = (int)(tau * FS);

                    delay_signal(tx_zc, rx, L_zc, delay_samples);

                    phase = 2.0f * PI * fd * t_now;
                    shift = cosf(phase) + I * sinf(phase);

                    for(int i = 0; i < L_zc; i++)
                        rx[i] *= shift;

                    awgn_measured(rx, L_zc, snr_db);

                    // Remove cyclic prefix (second half only)
                    //float complex rx_win[NSAMP];
                    for(int i = 0; i < NSAMP; i++)
                    {
                        rx_win[i] = rx[NSAMP + i];
                    }
                    // Perform correlation using FFT-based method later
                    // For now, use direct correlation
                    //float complex circ_out[NSAMP];
                    correlate_circular_zc(rx_win, zc, circ_out, NSAMP);
                    for(int i = 0; i < NSAMP; i++)
                    {
                        cube_z[i][pulse_idx_z] = circ_out[i];
                    }

                    pulse_idx_z++;
                    trigger_time_z += PRI_ZC;
                }
            }

             
            // using 2-D FFT
            // After all batches, compute Range-Doppler map for Golay and ZC
            // for golay
            compute_range_doppler(cube_g, rd_map_g, NSAMP, TOTAL_PULSES_GOLAY);
            // for zc
            compute_range_doppler(cube_z, rd_map_z, NSAMP, TOTAL_PULSES_ZC);


            // find peaks 
            int r_g_idx, v_g_idx;
            int r_z_idx, v_z_idx;
            // for golay
            find_peak(rd_map_g, NSAMP, TOTAL_PULSES_GOLAY, &r_g_idx, &v_g_idx);
            // for zc
            find_peak(rd_map_z, NSAMP, TOTAL_PULSES_ZC, &r_z_idx, &v_z_idx);


            // Compute errors
            // Range error is (estimated_range - true_range)^2
            // Velocity error is (estimated_velocity - true_velocity)^2
            // We accumulate these errors over all Monte Carlo runs for each SNR level, and then take the square root of the average to get RMSE.
            err_rg += powf(r_axis[r_g_idx] - true_r, 2);
            err_vg += powf(v_axis_g[v_g_idx] - true_v, 2);

            // ZC errors
            err_rz += powf(r_axis[r_z_idx] - true_r, 2);
            err_vz += powf(v_axis_z[v_z_idx] - true_v, 2);

        }

        // After all Monte Carlo runs for this SNR level, compute RMSE
        rmse_r_g[s] = sqrtf(err_rg / NUM_MC); // RMSE for Golay range
        rmse_v_g[s] = sqrtf(err_vg / NUM_MC); // RMSE for Golay velocity
        rmse_r_z[s] = sqrtf(err_rz / NUM_MC); // RMSE for ZC range
        rmse_v_z[s] = sqrtf(err_vz / NUM_MC); // RMSE for ZC velocity

        // Print results for this SNR level
        printf("| %3d | %12.4f | %12.4f | %10.4f | %10.4f |\n",
       snr_db,
       rmse_r_g[s], rmse_v_g[s],
       rmse_r_z[s], rmse_v_z[s]);
        
    }

    printf("-------------------------------------------------------------------\n");


    //cleanup_platform();

    return 0;
}