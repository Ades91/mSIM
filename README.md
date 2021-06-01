# mSIM

A collection of Matlab script and function to process 2D SIM data.

main_SIM.m is the most generic file. It implements, in a succession of cell structure, the whole processing pipeline of SIM reconstruction, namely 
loading, pre-processing, rough parameter estimation, automatic reconstruction (parameter estimation is refined during the reconstruction), rendering and saving.


main_MP_SIM.m is used to process data generated using a unique microscope equiped with a Multi-plane prism. 
In addition to the SIM reconstruction, the script handles the data formating and plane coregistration.


main_SOFI_SIM.m is used to process data of blinking emitters illuminated with a sinusoidal pattern.
After loading, the data is reordered and nth order SOFI is computed for each phase and angle. 
The SOFI images are then recombined using a non-linear SIM reconstruction.