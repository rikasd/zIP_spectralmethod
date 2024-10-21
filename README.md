# Supplementary material for Sugimoto-Dimitrova et al. (2024)

## Quick Start Guide:
1. Download the public data set from Santos and Duarte (2016) and add it to the "Santos2016Data" folder.
2. Ensure that you have MATLAB R2022a or later installed, along with the following toolboxes:
    > Control System Toolbox, Signal Processing Toolbox, Statistics and Machine Learning Toolbox.
3. Navigate to the zIP_spectralmethod folder in MATLAB and run run_main.m to generate figures from Sugimoto-Dimitrova et al. (2024).


## Organization of the repository:

\zIP\             contains functions for computing zIP from data 
                  and from models
                  
\model\           contains functions related to the balance model and
                  simulation
                  
\data_processing\ contains a custom cross-spectral density function that
                  detrends each data segment


## Acknowledgements:
The simulation code is largely based on Shiozawa et al. (2021)

## References:

Santos, D. A., Duarte, M. (2016). A public data set of human balance evaluations. Figshare. https://dx.doi.org/10.6084/m9.figshare.3394432.v2.

Shiozawa, K., Lee, J., Russo, M. et al. (2021). Frequency-dependent force direction elucidates neural control of balance. J NeuroEngineering Rehabil 18, 145. https://doi.org/10.1186/s12984-021-00907-2

Sugimoto-Dimitrova, R., Shiozawa, K., Gruben, K. G., & Hogan, N. (2024). Frequency-domain patterns in foot-force line-of-action: an emergent property of standing balance control. Journal of Neurophysiology. https://doi.org/10.1152/jn.00084.2024
