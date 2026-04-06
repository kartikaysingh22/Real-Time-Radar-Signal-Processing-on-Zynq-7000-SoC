# Real-Time-Radar-Signal-Processing-on-Zynq-7000-SoC
This project implements a real-time radar signal processing pipeline on the Xilinx Zynq-7000 SoC platform. The system compares the performance of Golay Complementary Pairs (GCP) and Zadoff-Chu (ZC) sequences for target range and radial velocity estimation under varying noise conditions.
The implementation features a complete bare-metal signal processing chain in Embedded C, optimized for high-performance execution on the ARM Cortex-A9 core.

# Key Features
1. Waveform Design: Generation of length-512 Golay pairs and Zadoff-Chu sequences.
2. Custom Signal Processing: Bare-metal implementation of Radix-2 FFT, IFFT, and matched filtering.
3. Hardware Acceleration: Extensive use of compiler optimizations and algorithmic improvements (Twiddle Factor LUTs) to reduce processing time.
4. Simulation Framework: FIFO-based pulse scheduling and Monte Carlo simulations for performance evaluation across SNR levels from -45 dB to -30 dB.

## Optimization Impact

The project explores several optimization tiers to minimize the computational bottleneck in the Receiver (Rx) stage:

| Optimization Level | Golay Total Time | ZC Total Time |
| :--- | :--- | :--- |
| **Baseline (No Optimization)** | [cite_start]5.45 s [cite: 411, 1269, 1270] | [cite_start]1.17 s [cite: 411, 1269, 1270] |
| **Twiddle Factors + Compiler Flags** | [cite_start]1.51 s [cite: 428, 445, 1271] | [cite_start]0.402 s [cite: 428, 446, 1272] |
| **Total Time Reduction** | [cite_start]**~72.3%** [cite: 445] | [cite_start]**~65.6%** [cite: 446] |


# Results
The system demonstrates that Golay processing generally achieves lower Root Mean Square Error (RMSE) at lower SNR levels compared to Zadoff-Chu, showing higher robustness in extremely noisy environments. A sharp drop in estimation error is observed around -40 dB to -38 dB, where the target signal becomes distinguishable from noise. 

more details can be found in the report.
