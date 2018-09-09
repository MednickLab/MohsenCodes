# MohsenCodes
ace_analysis: Matlab function, calculates the heart rate burs related changes in delta (SWA), sigma, and HF powers.
Find a GUI for the classic HRV analysis in https://github.com/mohsenbme/sleepHRV
mySOstats : Matlab function, detects the sleep slow oscillations and saves stats such as SO density in csv and mat files.
mySpindleStats: Python function, detects the sleep spindles (Mohsen's method) and EEG power spectra and saves the outputs in mat format.
spndl.ipynb: Jupyther notebook for using the mySpindleStats function.
There is another spindle detection code based on Wamsely et al paper in https://github.com/mohsenbme/EEGspindles. You can use it only of you have saved separate edf files for Stage 2 and SWS. It also gives you other stats such as spindle frequency, duration, and amplitude.
myeegpwr_v2: calculates EEG powers in frequency bands for artifact-rejected epochs. Output format: csv 
myeegpwr_v3: calculates EEG powers in frequency bands, no artifact rejection.

mytopoplot: for EEG headplots. First, install eeglab in Matlab. Based on the number of electrodes use appropriate electrde location files: myloc_file.mat (56 electrodes), myloc_file_18ch (no midline), or myloc_file_24ch. For other caps you need to make your own files. example: figure; mytopoplot(spDensitiesStg2,'myloc_file_24ch.mat');xlim([-.537 0.537]);ylim([-0.525 0.565]);caxis([0 4]) 

morlet_IR: wavelet-based time-frequency representation 
