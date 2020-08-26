function master_phytospikeclusters
%% HW: 01.08.2020
%script to get spiking info by cluster and trials
%combines multiple analysis scripts in one 

%% EDIT PATH FILES AND OPTIONS HERE
%settings for spiketimes by cluster. gets spike times from phy and
%organizes into cluster and trials as a 3D matrix.
ops.rectype= 'openephys'; %which recording system used
ops.opto= 'yes'; %if optogenetic experiment

ops.byexpt='no'; %if want to seperate trials by time
ops.savenames={'ramp','train','square'};

%for RF/STA experiments, YES if ran, NO is did not.
ops.flash_expt= 'yes'; %if ran sparse noise/STA expt

ops.flashprepost= 'no'; %if prepost included in flash expt (usually no for sprase noise)
ops.flash.blank_pulses=[];
ops.flash.nflashexpts=1;
ops.flash.ishartley=0;
ops.flash.sparsenoise=1;

%path files for create stim_times file %0 1 2 3 5 6

path_files{1}='/Volumes/Seagate_hw/data/ephys_data/20200226/HRV119_2020-02-26_11-14-25_001/';
path_files{2}='/Volumes/Seagate_hw/data/ephys_data/20200226/HRV119_2020-02-26_12-22-44_002/';
path_files{3}='/Volumes/Seagate_hw/data/ephys_data/20200226/HRV119_2020-02-26_13-30-42_003/';
path_files{4}='/Volumes/Seagate_hw/data/ephys_data/20200226/HRV119_2020-02-26_14-39-19_004/';
%path_files{5}='/Volumes/windows_bootcamp/Helen/20200401/HR24_2020-04-01_14-09-33_007/';
%path_files{6}='/Volumes/windows_bootcamp/Helen/20200403/HR25_2020-04-03_15-08-10_007/';
%where to save stim_times file (i.e. on an external drive)
stim_save_path='/Volumes/Seagate_hw/data/ephys_data/20200226/';

%select the analyzer files used
disp('select analyzer files (EXCLUDING ANY STA EXPERIMENTS)');
[file, path]=uigetfile('*.analyzer','MultiSelect','on');
for i= 1: size(file,2)
   analyzer_files{i}=[path file{i}]; 
end

%% get spike times by cluster
[stimulus2, stimulus2_onoff, cluster_id, header2]=HW_spiketimes2_bycluster(path_files, analyzer_files,stim_save_path, ops);

%% if don't run STAs w/ spiketimes
 [psth_norm,rfMap]=sparsenoise_STA(flash_spiketimes, flash_seq);

%% get waveforms
%TO DO:
%figure out wtf peakdur is (2nd peak width)-> bigger = FS???
[templatewf] = getTemplateWaveforms;

%to do:
%classify as FS/RS based on existing pool of wf features collected?

%% take spiketimes and calculate firing rate
%get which units are "significant"
%get mean firing for different sitmulus conditions

%for STF trials
[stimulus_fr_stat,fr_pertrial,keep,header_frlist,stimtypes,pvalues]=...
    stimulus_frmat(stimulus2(1:2196,:,:),stimulus2_onoff(1:2196,:),cluster_id, 1, [], 'yes','yes');


stimulus_frmat(stimulus2(2197:2547,:,:),stimulus2_onoff(2197:2547,:),cluster_id, 1, [], 'yes','yes');

%for como trials
norm=1; %nor to 0 deg/s dots vs. grey screen
[fr_pertrial,keep_best,header_frlist] = ...
    como_frmat(stimulus2(2197:2635,:,:),stimulus2_onoff(2197:2635,:), cluster_id, 1,'', 'yes', 'yes',norm);
%% run CSD
%TO DO: RE-DO HOW YOU DEFINE AL/PM HVA BORDERS (no obvious granular layer,
%so use L2/3 sink instead???)
z_score=0;
combine_columns = 0;
probe='128AN_bottom';
nchans=128;
probe_n=1;
CSD_openEphys(probe, probe_n, nchans, z_score, combine_columns);

%% run OPTO CSD
%To DO:
%create CSD code for parsing out CSDs when light is turned on to
%activate/inactivate cells
%% get cluster dpeth info

%MUST have CSD layer info
%region= 'cortex'
%tosave = 1, save matrix
[norm_cluster_depths]=kilosort_clusterdepth('cortex',1);

%% TO DO:
%plot by shank and depth (STAs) to determine spatial RFs relative to probe
%positioning.
fn_STAs='20200528_STAs'; %rfMaps = need to be transposed by ' to correctly match screen locations if used from MAK's sprase noise STA!!
fn_depth = '20200528_cluster_depths';
probe='128AN_bottom';
cluster_id = load('20200528_spiketimes.mat', 'cluster_id');
cluster_id=cluster_id.cluster_id;
ncols=8; %ncols for figures
%1 figure per shank
plotRFvsPosition(fn_STAs, fn_depth, probe, cluster_id, ncols)

%% calculate OMI
fn_frmat='20200528_fr_nobaselinesub_ORI';
fn_depths = '20200528_cluster_depths';
ops.stimcols=[1:3]; %[1:3]; %which columns contain the stimulus info
ops.lbrcols=[4:5]; %which columns contain the lb on/off and run info (lb first, then run)
ops.run_trials = 0; %0 stat 1 run, which trial conditions to calc OMI for
ops.ledtitle='Blue LED'; %what LED color/title you want to use
%ops.resp_pvalues = pvalues.p_stat; %which p-vals to use when calculating OMI for responsive stim only
main_omi(fn_frmat, fn_depths, ops)
%%

end