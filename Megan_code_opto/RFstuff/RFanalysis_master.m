function RFanalysis_master(exp_path)

% for openEphys data


binsize = 10;      % in ms

%% get experiment name and load raw data
% get name
out = regexp(exp_path,'\\','split');
an_name = out{end-1};       % animal name
if strcmpi(an_name,'LP') || strcmpi(an_name,'V1')
    an_name = strcat(out{end-2},'_',an_name);
elseif contains(an_name,'day','ignorecase',1) 
    an_name = out{end-2};
end
inds = regexp(out{end},'_\d','start');  % looks for what's before an underscore followed immediately by a number; uses that as exp_name
exp_name = out{end}(1:inds(1)-1);
if sum(isstrprop(exp_name,'digit'))>4       % this means it's an OpenEphys file (because it names files starting with the date)
    inds = regexp(out{end},'_\D','start');
    exp_name = strcat(out{end-1},out{end}(inds(1):end));  
    if strfind(out{end}(inds(1):end),an_name)
        exp_name = out{end}(inds(1)+1:end);
    else
        exp_name = strcat(an_name,out{end}(inds(1):end));
    end
    exp_name = exp_name(1:strfind(exp_name,'sparsenoise')-2); % remove "sparsenoise" from end of exp_name
    
end

disp(strcat('Loading raw data from experiment ',exp_name))
cd(exp_path)
if exist(sprintf('%s/data.mat',exp_path),'file')
    load(sprintf('%s/data.mat',exp_path))      
else
    openEphys2matlab_sparsenoise(exp_path);
    load(sprintf('%s/data.mat',exp_path))      % data from intan2matlab.m
end

% data.mat includes the following variables:
% trials - num_trials x 2 matrix of trial start and end times (in seconds)
% field_trials - trial num x 2 matrix of trial start and end samples (in 
    % 1000 Hz sampling rate) - starts from 1!
% time_index - vector of time stamps of each sample (in seconds), with 1000
    % samples per second - starts from 0!
% amp_sr - original amplifier sampling rate
% epoc - vector of 1s and 0s indicating whether each sample was part of a 
    % trial (1) or not (0). Still in original sampling rate (amp_sr)
% mmvt -  vector of 1s and 0s indicating whether mouse was running (1) or  
    % not (0) at time of sample. Still in original sampling rate (amp_sr)
% encd - vector of analog output from movement encoder. Still in original 
    % sampling rate (amp_sr)
% light - is 0 if not an optogenetics experiment, but otherwise, vectors of
    % whether LED was on (1) or off (0)
% re - vector for the timestamps (rising phase of photodiode signal) in 
    % sec. you need to remove the first 2 and the last 2 timestamps. 
    % Timestamps signify the end of a four second black to gray period.
% photo - vector of analog output from photodiode. In original sampling rate
% stim_times - 
% start_times - 

% set up where to store results
if contains(exp_name,'M') || contains(exp_name,'D') || contains(exp_name,'VH')
    full_dir = 'H:\LPproject\LPresults';        % where you want to save data
elseif contains(exp_name,'TH')
    full_dir = 'H:\Tlx3project\Haloresults';
else
    full_dir = 'H:\Tlx3project\ChR2results';
end
main_dir = strcat(full_dir, '\',an_name);
all_fig_dir = sprintf('%s\\RFfigures',main_dir);
if ~exist(all_fig_dir)
    mkdir(all_fig_dir)
end
% if plots
%     if exist(all_fig_dir,'dir')
%         delete(sprintf('%s/*',all_fig_dir));        % delete and remake all figures
%     end
% end

%% load clustering info
disp('Loading clustering data...')
% load kilosort data
spike_times = readNPY('spike_times.npy');           % here, spike times from ALL clusters
clusters = readNPY('spike_clusters.npy');
if exist('..\cluster_group.tsv','file')
    fid=fopen('..\cluster_group.tsv');     % new phy2 output
else
    fid = fopen('..\cluster_groups.csv');
end
% if exist('cluster_groups.csv','file')
%     fid = fopen('cluster_groups.csv');
% elseif exist('..\cluster_groups.csv','file')
%     % get cluster group
%     fid = fopen('..\cluster_groups.csv');
% end

% read column headers
C_text = textscan(fid, '%s', 1, 'delimiter', ',');
% read group info
grp_dat = textscan(fid, '%f %s', 'delimiter', ',');
fclose(fid);
units = grp_dat{1};             % cluster numbers
cluster_group = grp_dat{2};     % 'good', 'noise', etc.
good_units = units(strcmp(cluster_group,'good'));

%% get stimulus info
stim_file = dir(fullfile(exp_path,'_*_*_*.mat'));       % output format for stimulus info .mat file
fprintf(sprintf('Loading stimulus file %s\n',stim_file.name'))
stimdat = load(fullfile(exp_path,stim_file.name),'-mat');
analyze_file = dir(fullfile(exp_path,'*.analyzer'));
fprintf(sprintf('Loading analyzer file %s\n',analyze_file.name'))
load(fullfile(exp_path,analyze_file.name),'-mat');  % under variable name 'Analyzer'
nX = Analyzer.P.param{cellfun(@(x) strcmp(x{1},'Nx'), Analyzer.P.param)}{3};
nY = Analyzer.P.param{cellfun(@(x) strcmp(x{1},'Ny'), Analyzer.P.param)}{3};
nBW = Analyzer.P.param{cellfun(@(x) strcmp(x{1},'bw_bit'), Analyzer.P.param)}{3};
num_reps = length(fieldnames(stimdat))-1;       % assumes one of the fields is for monitor refresh rate ('frate')
num_stim = length(stimdat.randlog_T1.seqs.xseq);

conds = cellfun(@(x) x.val{strcmp(Analyzer.loops.conds{1}.symbol,'light_bit')}, Analyzer.loops.conds);
diff_conds = unique(conds);
trialcond = cellfun(@(x) x.repeats{1}.trialno,Analyzer.loops.conds);
lightcond = zeros(1,num_reps);
for i = 1:length(diff_conds)
    lightcond(ismember(trialcond,find(conds==diff_conds(i)))) = diff_conds(i);
end

% downsample LED to 1000Hz
div             = amp_sr/1000;
zx              = 1:div:length(LED);
izx             = floor(zx);
light = LED(izx);
% pho = photo(izx);

T = 1000/(60/Analyzer.P.param{14}{3});    % length of each stim presentation (in ms). 14th entry is h_per
light_T = Analyzer.LP.pulse_dur(1);       % length of light pulse duration (in ms)
% make matrix of "trials" for each sparse noise stimulus presentation
stim_trials = zeros(num_reps*num_stim,4);       % 4 columns for x pos, y pos, b/w, and lightcond
for n = 1:num_reps
    xx = eval(sprintf('stimdat.randlog_T%d.seqs.xseq',n));
    yy = eval(sprintf('stimdat.randlog_T%d.seqs.yseq',n));
    bw = eval(sprintf('stimdat.randlog_T%d.seqs.bwseq',n'));
    bw(bw==1) = -1; % 1 means BLACK
    bw(bw==2) = 1;  % 2 means WHITE
    for nn = 1:num_stim
        stim_trials((n-1)*num_stim+nn,1:3) = [xx(nn) yy(nn) bw(nn)]; 
        if sum(light(stim_times((n-1)*num_stim+nn):stim_times((n-1)*num_stim+nn)+T)) > round(max(light)-std(light),2)*min(T,light_T)     
%         if find(light(stim_times((n-1)*num_stim+1):stim_times((n)*num_stim))>4) 	% TEMP
            stim_trials((n-1)*num_stim+nn,4) = 1;
        end
    end
end

if sum(stim_trials(:,4)) ~= size(stim_trials,1)/2
    warning('Unequal numbers of light and no-light trials detected')
end

% make "key" matrix for all possible sparse noise stimuli
lconds = unique(stim_trials(:,4));
bwconds = unique(stim_trials(:,3));
count = 1;
stim_key = nan(num_stim,4);
for lc = 1:length(lconds)
    for b = 1:nBW
        for y = 1:nY
            for x = 1:nX
                stim_key(count,:) = [x y bwconds(b) lconds(lc)];
                count = count+1;
            end
        end
    end
end

%% run RF analysis for each unit in this experiment

num_units = length(good_units);
% get units' shanks
full_exp_path = fileparts(exp_path);
load(fullfile(full_exp_path,'rez.mat'));          % load rez.mat output from kilo (for waveform extraction)

for i = 1:num_units
    fprintf(sprintf('Performing RF analysis on %s cluster %s \n',exp_name,num2str(good_units(i))))
    unit_times = spike_times(clusters==good_units(i));
    unit_times_ds = floor(unit_times./(amp_sr/1000));   % change sample #s to account for downsampling
    unit_times_ds = unit_times_ds + 1; % has to be +1 because spike_times starts at 0, but the min possible field_trials value could be 1
    [~,~,shank(i)] = readWaveformsFromRez_K2(good_units(i),full_exp_path,rez);      % default to read waveforms from Rez
    [psth_norm{i}, rfMap{i}, stats(i)] = sparsenoiseRF(unit_times_ds,stim_times,binsize,T,stim_trials,stim_key);
    % save figures
    save_name = sprintf('%s\\%s_%s_Cluster%d',all_fig_dir,exp_name,strcat('shank',num2str(shank(i))),good_units(i));
    set(gcf, 'PaperPosition', [0 0 24 12]);
    print(gcf,'-dpng',save_name)
    close all
end

%% save results
save(fullfile(main_dir,'RFresults'),'good_units','shank','psth_norm','rfMap','stats')

end
