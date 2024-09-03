# Audio Signal Processing with MATLAB

This repository contains a MATLAB script that demonstrates various audio signal processing techniques. The code performs tasks such as signal extraction, cross-correlation, resampling, Fourier Transform analysis, custom DFT implementation, noise addition, modulation, demodulation, and filtering.

## Key Features
- **Time-Domain Analysis**: Visualization of the original audio signal and its segments in the time domain.
- **Frequency-Domain Analysis**: Fourier Transform and magnitude spectrum visualization for original, downsampled, and upsampled signals.
- **Cross-Correlation**: Analysis of the correlation between different segments of the audio signal.
- **Resampling**: Downsampling and upsampling techniques and their effects on the signal.
- **Noise Addition**: Application of Additive White Gaussian Noise (AWGN) to the signal.
- **Modulation and Demodulation**: Signal modulation and subsequent demodulation with frequency analysis.
- **Filtering**: Application of custom filters to process the demodulated signal.

## Files
- `ghayegh.wav`: Original audio file used in the processing.
- `ghayegh2.wav`: Extracted 10-second segment.
- `ghayegh2_1.wav`: Downsampled version of the extracted segment.
- `ghayegh1_2.wav`: Upsampled version of the extracted segment.
- `ghayegh-n.wav`: Noisy version of the signal with AWGN added.
- `y_modulated.wav`: Modulated signal.
- `y_demodulated.wav`: Demodulated signal.
- `y-demodulated2.wav`: Filtered demodulated signal.

## Usage
To run the script, simply load the `audio_signal_processing.m` file in MATLAB and ensure that the required audio files (`ghayegh.wav`) are available in the same directory.

## Requirements
- MATLAB with Signal Processing Toolbox.

## Author
Ali Pazani

## License
This project is licensed under the MIT License.
