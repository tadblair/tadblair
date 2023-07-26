# Hippocampal place cell remapping occurs with memory storage of aversive experiences. https://doi.org/10.7554/eLife.80661
---
Within the this repository there are two main folder, one for data and one for scripts that access the data.
#### Directory structure:
#### Data files ####
tadblair-main\Blair_et_al_DATA\
                    |
                    |-> 'cellmaps\'
                    |-> 'prepost\'
                    |-> 'presizecc\'
                    |-> 'pretrn\'
                    |-> 'sessiondata\'
                    --> 'shocktimes\'
#### Scripts for recreating figures ####
tadblair-main\Blair_et_al_MATLAB
                    |
                    |-> 'Behav_effects.mat                              
                    |-> 'Figure1_behavior.m                  
                    |-> 'Figure2_analysis.m                  
                    |-> 'Figure2_plotting.m                  
                    |-> 'Figure3_analysis.m                  
                    |-> 'Figure3_plotting.m  
                    |-> 'Figure4A_analysis.m                      
                    |-> 'Figure4A_plotting.m             
                    |-> 'analyze_shock_responses_drugfree.m  
                    |-> 'analyze_shock_responses_scopo.m     
                    |-> 'between_session_analysis_df.m       
                    |-> 'get_beelines.m                      
                    |-> 'pv_heatmap_cormatrix.m              
                    |-> 'pv_heatmap_tuningcorr.m             
                    |-> 'shadedErrorBar.m                    
                    |-> 'simple_mixed_anova.m                
                    --> 'subplot_tight.m    

##### Script notes ####
# Scripts to run
These scripts should recreate the main figures from our paper. Generally you should run Figure1 first, followed by Figure2_analysis, then Figure2_plotting, and so on (there could be slight dependencies of later figure scripts on earlier ones if you encounter an issue). Below is the recommended order of running. Scripts not listed are dependency scripts called by the figure scripts. Ploing script depends on corresponding analysis script only

allfigs_RUNME.m - will generate all figures sequentially, just change the parent directory name and make sure it's in your matlab path

Or run the other scripts in this order:
--> 'Figure1_behavior.m: generates figure 1 from BehavEffects.mat                  
--> 'Figure2_analysis.m: prepares the data for recreating figure 2
 |-> 'Figure2_plotting.m: generates figure 2
--> 'Figure3_analysis.m: prepares the data for recreating figure 3
 |-> 'Figure3_plotting.m: generates figure 3
--> 'Figure4A_analysis.m: prepares the data for recreating figure 4
 |-> 'Figure4A_plotting.m: generates figure 4A


# provided dependency scripts
'between_session_analysis_df.m': (used in Figure3_plotting) to quantify population vector effect across sessions
'get_beelines.m: extracts beelines from behavior data (no stopping runs along short path)
'pv_heatmap_cormatrix.m: generates the population vector correlation matrix for firing similarity along short path, all pairwise correlation
'pv_heatmap_tuningcorr.m:  generates the population vector tuning correlation for field similarity along short path, only along the diagonal self similarity
'shadedErrorBar.m: plotting utility to plot mean and error range 
'simple_mixed_anova.m: for performing the mixed anova analysis
'subplot_tight.m: for generating tighter plots
'analyze_shock_responses_drugfree.m: plots figure 4 data for drug free shock responses
'analyze_shock_responses_scopo.m:  plots figure 4 data for scopolamine shock responses
## BehavEffects.mat: contains behavior data and statistical effects for generating Figure 1. Each field in 'Effects' is a dependent variable of matrices corresponding to the experimental condition ('shock' first/only shock session, 'barrier', 'scopo' scopolamine only, 'scoposhock' scopolamine and shock, 'shock2' second shock session)
    Fields in 'Effects' variable:
        - rewrate: rewards per minute. NaN's fill animal rows that did not recieve the manip. Each field follows:
            *manipulation*: (animal by session) where sessions (cols) 1:3 are the 3 sessions preceeding the manip, 4 is the first 10 minutes of the session (before the manip is given), 5 is the 5 minutes after the manip is introduced, and cols 6:8 are the 3 sessions immediately following the manip
            *manipulation*_baseline: pre manipulation value
            *manipulation*_mean: by session, mean value across animals
            *manipulation*_p: paired t-test effect of manip, (cols 4 vs 5)
            pre_manip: value before their first manipulation, if multiple were recieved
        - prob_short: proportion of beelines taken along the short (short / short+long)  (for subfield explanation see rewrate)
        - path_length: total distance ran, centimeters  (for subfield explanation see rewrate)
        - num_short: number of short paths  (for subfield explanation see rewrate)
        - num_long: number of long paths (for subfield explanation see rewrate)
        - groups: fields list which animal numbers are in which experimental manipulation. All 15 recieved a shock manipulation, and subsets recieved some combination of others
        - keys: key values are as follows
            -subdays: for *manipulation* (above), the cols corresponding to the pre and post training session
            -animals: per animal index, their experimental number
            -manip_day: per animal number, which session day they recieved a manip
            -manip_key: per animal number, which manip was recieved (order follows manip_day). manip keys:
                's' shock, 'b' barrier, 'scs' scopolamine shock, 'sc' scopolamine only, 'spont' spontaneous avoidance because the track settled into the track and scared the rat and made them avoid but not due to shock (excluded from analysis), 's2' second shock session

##### End script notes ####

## Sharing/Access information
You are free to use this data for educational purposes. Publications using this data should cite this repository dataset and the original paper: https://doi.org/10.7554/eLife.80661. If you have any questions or problems contact either Garrett J. Blair or Hugh Tad Blair (we're not married or related just fun coincidence). See bottom for contact details.

Links to other publicly accessible locations of the data:
  * https://github.com/tadblair/tadblair

Data was derived from the following sources:
  * rat hippocampal brain activity using 1-photon miniscopes (miniscope.org)
  More information on the methodology can be found in our open-access methods paper: https://www.science.org/doi/10.1126/sciadv.adg3918

## Code/Software

# Blair_et_al_DATA: (all matlab files)

#### note - 'LR' and 'RL' notation in all data files corresponds to the 'left to right' or 'right to left' running direction for beelines on the short track. Repetition for LR and RL is denoted with * where applicable

## cellmaps: matching matrices for each animal and condition, ex) Hipp#_condition_cmap.mat. N x 3 matrix (cell# by session) where each entry is the cell id in the session. '0' in cmap indicates the cell was not matched for that session. (adapted from Sheintuch et al. 2017's 'cell_to_index_map')
    VARIABLES:
        - cmap: cell_id by session matrix, used for indexing matched cells in the session data files (non-zero intersection of the given columns)
        - sessionNums: session number of corresponding cmap cols

## prepost: matched data for comparing the effect of training, 2 session per rat (pre and post) per condition. 'pre' is before training, 'post' is after extinction
    'predata' structure VARIABLES: pre-training data
        - *int: [beeline# by 2] corresponds to the start and stop index of each beeline along the track for the LR or RL direction
        - dcurve_*: [cell by position bin] firing field along the short track, events/sec
        - place_score*: [cell] spatial information score (prob of beating chance, -log10(prob))
        - spikes_per_*bee: [cell] average spikes per beeline
        - dcurve_*_bps: bits per spike for the corresponding cels dcurve
        - subspd_median: median speed (cm/sec) along all beelines after subsampling to match behavior between scopo and non-scopo conditions
        - subspd_mean: mean speed (cm/sec) along all beelines after subsampling to match behavior between scopo and non-scopo conditions
    'postdata' structure VARIABLES: post-extinction data
        - shock_code: shock modulation of cell (1 for significantly excited, -1 for surpressed, 0 for neither)
        - all other fields have the same interpretation as in 'predata'

## presizecc: data files for evaluating cells eccentricity (ecc, radial distance from center) and cell size (csize, in pixels) effect of recurrence or place score in Figure 2—figure supplement 3. Each file corresponds to an animals recording session, as in 'sessiondata/', and contains the variables:
    'cecc': cell eccentricity from center (cpos) in pixels
    'cpos': center of recoding, in pixels
    'csize': cell size, in pixels

## pretrn: matched data for comparing before training, 2 session per rat (pre and post) per condition. 'pre' is before training, 'trn' is during training (before shock)
    'predata' structure VARIABLES: pre-training data
        - all other fields have the same meaning as the 'predata' in ##prepost data (above)
    'trndata' structure VARIABLES: training data
        - vmap_*: occupancy map in seconds for either direction in the split session halves
        - fields with *part1 and *part2: correspond to the split half of the session for doing within session analysis
        - all other fields have the same interpretation as in 'predata'

## sessiondata: extracted calcium and position data for each session used, ex) Hipp#_session#_sess.mat. When loaded, a frame# structure variable where # corresponds to the session number. These come in sets of 3 per animal per condition they were included in. Each of the three session files corresponds to the sessions that compose the 'predata', 'trndata', and 'postdata' previously described 
    frame# VARIABLES
        - S: [cell by frame#] deconvolved spike times, output from CaImAn analysis (1== spike inferred)
        - x & y: x and y position per frame along the rectangular maze in centimeteres
        - time: session time per frame since recording start in milliseconds
        - deconv: [cell by frame#] value of deconvolved spike inference
        - posbin: binned x-position along the short track (23 bins, size = 250/23)
        - spd: movement speed (cm/sec) at the given frame
        - vmap_*: occupancy (seconds) in a bin along the short track by direction

## shocktimes: matrices for all animals by condition 
    (barrier, shock, scopolamine shock, or second shock session, 'shock2') giving the time in seconds where each manipulation was encountered (either entering the shock zone or approaching the barrier). Animals not in a given condition have all NaNs in their row. The last approach is duplicated in a row to fill out the matrix. Sessions lasted a total of 900 seconds


# Blair_et_al_MATLAB: (all matlab files, compatible with MATLAB 2016 or later)
    - scripts for generating figures, along with dependencies. ex) Figure#_type.mat (type '_analysis' should be run before '_plotting'). See recommended order at the top of the readme

How to reach the authors:

    Garrett Blair - gblair@g.ucla.edu | garrettjblair.com

    Hugh (Tad) Blair - Professor of Behavioral Neuroscience, UCLA Psychology Department

-Garrett J. Blair, 7/21/2022