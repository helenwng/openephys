%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast_CSD_64D.m
% Script for extracting LFP data from 64D probes on Intan
%
% Hard coded version for fast analysis during recordings by Ethan McBride
% Adapted from IntanLFP_CSD.m by Megan A Kirchgessner
% Which was adapted from LFPv2_latest from Bryan J. Hansen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function CSD_thalamus(exp_path,probe)
%INPUTS: path, and '64D', '128DN_bottom', etc. 
z_score=0;
combine_columns = 0;

% get necessary data
cd(exp_path)
if exist(sprintf('%s/data.mat',exp_path),'file')
    load(sprintf('%s/data.mat',exp_path))      % data from intanphy2matlab.m
else
    openEphys2matlab(exp_path);      % data from intanphy2matlab.m
    load(sprintf('%s/data.mat',exp_path))
end

%% memory map the amplifier file, filter, & downsample

first_half = exist(fullfile(exp_path,'100_CH1.continuous'),'file');
if first_half
    contfile = fullfile(exp_path,'100_CH1.continuous');
else
    contfile = fullfile(exp_path,'100_CH65.continuous');
end
[data, dataTime, dataInfo] = load_open_ephys_data_faster(contfile);
dataTime = dataTime./dataInfo(1).header.sampleRate;       % uncomment if using load_open_ephys_data_faster
nsamps = length(dataTime);
amp_sr = dataInfo.header.sampleRate;

s=dir;
ch_names = {s(:).name};
nchans = length(cell2mat(cellfun(@(x) strfind(x,'_CH'),ch_names,'UniformOutput',false))); % count number of continuous channel files
short = round(nsamps/4);
v = zeros(short,nchans);   % preallocate matrix for raw data
v(:,1) = data(1:short);
clear data
for i = 2:nchans
    contfile = sprintf('%s/100_CH%d.continuous',exp_path,i);
    [tmp,~,~] = load_open_ephys_data_faster(contfile);
    v(:,i) = tmp(1:short);
end

eventfile = fullfile(exp_path,'all_channels.events')
[events,eventTime,info] = load_open_ephys_data_faster(eventfile);
amp_sr = info.header.sampleRate;

LFPlow = 1000; %low pass forLFP freqs
[bf, af] = butter(2,(LFPlow/(amp_sr/2)), 'low'); 

% LFP=zeros(length(v(:,1))/20,64);

totaltime=0;

var_name = cell(1,nchans);
for ch = 1:nchans
    if ch < 10
        var_name{ch} = sprintf('LFP0%d',ch);
    else
        var_name{ch} = sprintf('LFP%d',ch);
    end
end
cont=struct;

if strcmp(probe,'A')
    field_file_name = 'fields_A';
elseif strcmp(probe,'B')
    field_file_name = 'fields_B';
else
    field_file_name = 'fields';
end
field_output_dir = fullfile(exp_path, field_file_name);
[~,message]=mkdir(field_output_dir);
cd(field_output_dir);
fileID2 = fopen('lfps.txt','w');

LN = size(v,1)/20; %downsample length
div = (amp_sr/20)/1000;
zx = 1:div:LN;
izx = floor(zx);

for i=1:nchans
    tic;
%     if i>1
%         fprintf(fileID2,'\r\n');
%     end
%     ch=sprintf('%s%d', 'Processing channel: ', i);
%     disp(ch);
    if i<10
        save_lfp=sprintf('%s%d','lfp0', i);
    else
        save_lfp=sprintf('%s%d','lfp', i);
    end
    if z_score == 1
        lfp = LineFilterRmFFT(downsample(filtfilt(bf, af, zscore(v(:,i)),0),20),[59 62],0,0.05,0,amp_sr);
    else
        lfp = LineFilterRmFFT(downsample(filtfilt(bf, af, v(:,i)),20),[59 62],0,0.05,0,amp_sr);
    end

%     LFP(1:length(lfp),i) = lfp;
    
    if i>1
        fprintf(fileID2,'\r\n');
    end
    ch=sprintf('%s%d', 'Processing channel: ', i);
    disp(ch);
    if i<10
        save_lfp=sprintf('%s%d','lfp0', i);
    else
        save_lfp=sprintf('%s%d','lfp', i);
    end
    
    cont(i).sig = lfp(izx);
    %     cont(i).sig = lfp;
    cd(field_output_dir);
    %     save(save_lfp,'cont(i).sig');
    fprintf(fileID2,save_lfp);
    
    elapsedtime = toc;
    totaltime = totaltime+elapsedtime;
    timeleft = (64-i)*elapsedtime;
    fprintf('%d channels complete, %.2f sec remaining',i,timeleft)
    fprintf('\n')
    clear lfp
end

fclose(fileID2);



% get probemap
p = eval(sprintf('probemap_%s_func',probe));
[~,inds] = sort(p.z);  % sorts from bottom to top
chan_order = p.channels(inds);
chan_order = flipud(chan_order);    % change from top to bottom
xx = flipud(p.x(inds));
nshanks = length(unique(p.shaft));
% diff_cols = unique(xx);
% x_dif = diff(diff_cols);
% num_cols = length(diff_cols) - sum(diff(diff_cols)<5); % 64D probe has x vals at -20, -16, 0 , 16, and 20, but consider 16 and 20 the same
num_cols = nshanks;
count = 1;
% for i=1:length(x_dif)
%     if x_dif(i) < 5
%         if length(diff_cols)==1
%             xcolvals{count-1} = [xcolvals{count-1}; diff_cols];
%         else
%             xcolvals{count} = diff_cols(1:2);
%             diff_cols(1:2) = [];
%         end
%     else
%         xcolvals{count} = diff_cols(1);
%         diff_cols(1) = [];
%     end
%     count = count+1;
% end
probemap = nan(ceil(nchans/num_cols),num_cols);
shank = flipud(p.shaft(inds));
for i = 1:num_cols
%     col_chans = chan_order(ismember(xx,xcolvals{i}));
    probemap(:,i) = chan_order(shank==i);
end
probemap = probemap+1;          % change so that it stats from 1!!
if combine_columns == 1
    spacing = max(abs(diff(p.z(inds))));
else
    spacing = max(abs(diff(p.z(probemap(:,2)))));
end
% use only first two columns of each shank
twocols = p.x(probemap(:,1))==median(p.x(probemap(:,1)));
probemap = probemap(twocols,:);
%     
pos = probemap(:);
pos(isnan(pos)) = [];
zvals = p.z(pos);       % get corresponding vector of z values
shank = p.shaft(pos);
    
if combine_columns == 1
%     for i=1:numel(probemap(:,1:2))-sum(isnan(probemap(:,1:2)))  % only using channels in first two columns
    for i=1:length(var_name)
        sig{i}=sprintf('%s%d%s','=cont(',pos(i),').sig');
        channel_reorder = evalc([var_name{i}, sig{i}]);   % LFPs are ordered from top to bottom!
    end
else
    for i=1:length(var_name)
        sig{i}=sprintf('%s%d%s','=cont(',pos(i),').sig');
        channel_reorder = evalc([var_name{i}, sig{i}]);   % LFPs are ordered from top to bottom!
    end
end

elapsedtime = toc;
totaltime = totaltime+elapsedtime;
fprintf('Complete! total elapsed time was %.2f seconds',totaltime)
fprintf('\n')

%% find high-intensity opto pulses during blank trials
disp('Find timing using analog the Evoked Response Potential (ERP)')
tic

% new attempt
[~,~,~,trial_type,IVs] = get_exp_params(exp_path,'step');
tris = find(trial_type(:,1)==0 & trial_type(:,3)>0);
bl_field_trials = field_trials(tris,:);
% bl_field_trials = bl_field_trials(1:find(bl_field_trials(:,2)<short./(amp_sr/1000),1,'last'),:);
bl_trials = [trials(tris,1)+.5 trials(tris,2)-1.5];

% new_led=zeros(1,length(stim_times));
% % for i=1:length(stim_times)-1
% %     ldiff(i)=((stim_times(i+1))-(stim_times(i)))./60;
% %     if ldiff(i)>3.5 && ldiff(i)<4.1
% %         new_led(1,i)=(stim_times(i+1))/60;
% %     else
% %         new_led(1,i)=0;
% %     end
% % end
% new_led = stim_times/60;
% disp('Identify timestamps of photodiode')
% new_led(new_led==0)=[];% timestamps for the end of the photo
% for i=1:length(new_led)% the first two and last two pulses are removed
%     %%% this sections is important because it est the timing for the
%     %%% rest of the experiment 
%     %%% 100msec before flip, 400msec after
%     
% %     st_time(i,:)=(new_re(i)-3.00); end_time(i,:)=new_re(i)-1.00; %dark->light
% %     st_time_down(i,:)=(new_re(i)-1.00); end_time_down(i,:)=new_re(i)+1; %light->dark
% %     st_time(i,:)=(new_led(i)-2.50); end_time(i,:)=new_led(i)-1.50; %dark->light
%     st_time(i,:)=(new_led(i)-0.50); end_time(i,:)=new_led(i)+0.50; %light->dark
% end
st_time = trials(:,1) - .5;
end_time = trials(:,1) + .5;
st_time_down = trials(:,2)-.5;
end_time_down = trials(:,2)+.5;
% for i=1:length(st_time)
%     trials (i,1) = find(time_index>=st_time(i)&time_index<=end_time(i),1,'first');  % in samples, starting from 1
%     trials (i,2) = find(time_index>=(time_index(trials(i,1)))&time_index<=(time_index(trials(i,1))+2.000),1,'last');
%     trial_time(i,:)= trials(i,1):1:trials(i,1)+1000;
%     
%     trialsd(i,1) = find(time_index>=st_time_down(i)&time_index<=end_time_down(i),1,'first');  % in samples, starting from 1
%     trialsd(i,2) = find(time_index>=(time_index(trialsd(i,1)))&time_index<=(time_index(trialsd(i,1))+2.000),1,'last');
%     triald_time(i,:)= trialsd(i,1):1:trialsd(i,1)+1000;
% end
% 
% timing=(time_index(trial_time(1,:))-time_index(trial_time(1,1))-0.5);
% timingd=(time_index(triald_time(1,:))-time_index(triald_time(1,1))-0.5);

%% Largely redundant section for movement trials, may remove
% disp('Find the movement data and determine if the mouse was running or stationary')
ts = [ ];

% assumes there will never be more than 1000 events in an interval
maxevents = 1000;
% assumes a block will never be longer than 10000 seconds
maxtime = 10000;
% steps through block in 100 second intervals
steps = maxtime / 100;

mmvt=0;

for i=1:length(tris)
    if mmvt == 0
        mmvt_count(i) = 0;
    else
        mmvt_trials{i,:}    = mmvt(1,bl_trials(i,1):bl_trials(i,2));
        mmvt_count(i,:)     = length(find(mmvt_trials{i}==1));      % added by RK
    end
end
counter1=0;
counter2=0;
for i=1:length(mmvt_count)
    if mmvt_count(i)>=1 %>=2
        counter1=counter1+1;
        trial_running(counter1)=i;
    elseif mmvt_count(i)<1 %<=2
        counter2=counter2+1;
        trial_stationary(counter2)=i;
    end
end
for i=1:length(trial_stationary)
    trials2(i,:)= bl_field_trials(trial_stationary(i),1)+500:1:bl_field_trials(trial_stationary(i),1)+1500;
%     trials2d(i,:)= trialsd(trial_stationary(i),1):1:trialsd(trial_stationary(i),1)+1000;
end
% lastgood = find(trials2(:,end)<round(short/20),1,'last');
% trials2(lastgood+1:end,:) = [];
% trials2d(lastgood+1:end,:)=[];

trials2 = trials2(1:find(bl_field_trials(:,2)<short./(amp_sr/1000),1,'last'),:);

%% Remove trials with high diff from median!
disp('Removing noisy trials')
multiplier = 2; %stdev cutoff for trial removal
% DIFF = bsxfun(@minus,LFP01,median(LFP01,1));
% ERP01(std(bsxfun(@minus,ERP01,median(ERP01,1)),0,2)>multiplier*median(std(bsxfun(@minus,ERP01,median(ERP01,1)),0,2)),:)=[];
% ERP01=mean(ERP01,1);

for i=1:64
% ERP01 = LFP01(trials2);
% ERP01 = mean(ERP01(std(bsxfun(@minus,LFP01(trials2),median(LFP01(trials2),1)),0,2)<=multiplier*median(std(bsxfun(@minus,LFP01(trials2),median(LFP01(trials2),1)),0,2)),:),1);
    if i<10
        eval(['ERP0' num2str(i) '= LFP0' num2str(i) '(trials2);']);
        eval(['ERP0' num2str(i) '= mean(ERP0' num2str(i) '(std(bsxfun(@minus,LFP0' num2str(i) '(trials2),median(LFP0'...
            num2str(i) '(trials2),1)),0,2)<=multiplier*median(std(bsxfun(@minus,LFP0'...
            num2str(i) '(trials2),median(LFP0' num2str(i) '(trials2),1)),0,2)),:),1);']);
%         eval(['ERPd0' num2str(i) '= LFP0' num2str(i) '(trials2d);']);
%         eval(['ERPd0' num2str(i) '= mean(ERPd0' num2str(i) '(std(bsxfun(@minus,LFP0' num2str(i) '(trials2d),median(LFP0'...
%             num2str(i) '(trials2d),1)),0,2)<=multiplier*median(std(bsxfun(@minus,LFP0'...
%             num2str(i) '(trials2d),median(LFP0' num2str(i) '(trials2d),1)),0,2)),:),1);']);
    else
        eval(['ERP' num2str(i) '= LFP' num2str(i) '(trials2);']);
        eval(['ERP' num2str(i) '= mean(ERP' num2str(i) '(std(bsxfun(@minus,LFP' num2str(i) '(trials2),median(LFP'...
            num2str(i) '(trials2),1)),0,2)<=multiplier*median(std(bsxfun(@minus,LFP'...
            num2str(i) '(trials2),median(LFP' num2str(i) '(trials2),1)),0,2)),:),1);']);
%         eval(['ERPd' num2str(i) '= LFP' num2str(i) '(trials2d);']);
%         eval(['ERPd' num2str(i) '= mean(ERPd' num2str(i) '(std(bsxfun(@minus,LFP' num2str(i) '(trials2d),median(LFP'...
%             num2str(i) '(trials2d),1)),0,2)<=multiplier*median(std(bsxfun(@minus,LFP'...
%             num2str(i) '(trials2d),median(LFP' num2str(i) '(trials2d),1)),0,2)),:),1);']);
    end
end

if strfind(probe,'64D')
    ERP_shk1=vertcat(ERP01,ERP02,ERP03,ERP04,ERP05,ERP06,ERP07,ERP08,ERP09,ERP10,ERP11,ERP12,ERP13,ERP14,ERP15,ERP16,ERP17,ERP18,ERP19,ERP20,ERP21);
    ERP_shk2= vertcat(ERP22,ERP23,ERP24,ERP25,ERP26,ERP27,ERP28,ERP29,ERP30,ERP31,ERP32,ERP33,ERP34,ERP35,ERP36,ERP37,ERP38,ERP39,ERP40,ERP41,ERP42,ERP43);
    ERP_shk3 = vertcat(ERP44,ERP45,ERP46,ERP47,ERP48,ERP49,ERP50,ERP51,ERP52,ERP53,ERP54,ERP55,ERP56,ERP57,ERP58,ERP59,ERP60,ERP61,ERP62,ERP63,ERP64);

    % ERP_shk1=(aux_shk1/(1*10^8))*(1*10^6);% Scale factor for  ERP
    % ERP_shk2=(aux_shk2/(1*10^8))*(1*10^6);% Scale factor for  ERP
    % ERP_shk3=(aux_shk3/(1*10^8))*(1*10^6);% Scale factor for  ERP

    ERPd_shk1=vertcat(ERPd01,ERPd02,ERPd03,ERPd04,ERPd05,ERPd06,ERPd07,ERPd08,ERPd09,ERPd10,ERPd11,ERPd12,ERPd13,ERPd14,ERPd15,ERPd16,ERPd17,ERPd18,ERPd19,ERPd20,ERPd21);
    ERPd_shk2=vertcat(ERPd22,ERPd23,ERPd24,ERPd25,ERPd26,ERPd27,ERPd28,ERPd29,ERPd30,ERPd31,ERPd32,ERPd33,ERPd34,ERPd35,ERPd36,ERPd37,ERPd38,ERPd39,ERPd40,ERPd41,ERPd42,ERPd43);
    ERPd_shk3= vertcat(ERPd44,ERPd45,ERPd46,ERPd47,ERPd48,ERPd49,ERPd50,ERPd51,ERPd52,ERPd53,ERPd54,ERPd55,ERPd56,ERPd57,ERPd58,ERPd59,ERPd60,ERPd61,ERPd62,ERPd63,ERPd64);

    % ERPd_shk1=(aux_shk1d/(1*10^8))*(1*10^6);% Scale factor for  ERP      <<<< MAYBE CHANGE THIS
    % ERPd_shk2=(aux_shk2d/(1*10^8))*(1*10^6);% Scale factor for  ERP
    % ERPd_shk3=(aux_shk3d/(1*10^8))*(1*10^6);% Scale factor for  ERP
    
elseif strfind(probe,'64G')
    ERP_shk1=vertcat(ERP01,ERP02,ERP03,ERP04,ERP05,ERP06,ERP07,ERP08,ERP09,ERP10,ERP11,ERP12,ERP13,ERP14,ERP15,ERP16,ERP17,ERP18,ERP19,ERP20,ERP21,ERP22,ERP23,ERP24,ERP25,ERP26,ERP27,ERP28,ERP29,ERP30,ERP31,ERP32);
    ERP_shk2= vertcat(ERP33,ERP34,ERP35,ERP36,ERP37,ERP38,ERP39,ERP40,ERP41,ERP42,ERP43,ERP44,ERP45,ERP46,ERP47,ERP48,ERP49,ERP50,ERP51,ERP52,ERP53,ERP54,ERP55,ERP56,ERP57,ERP58,ERP59,ERP60,ERP61,ERP62,ERP63,ERP64);

%     ERPd_shk1=vertcat(ERPd01,ERPd02,ERPd03,ERPd04,ERPd05,ERPd06,ERPd07,ERPd08,ERPd09,ERPd10,ERPd11,ERPd12,ERPd13,ERPd14,ERPd15,ERPd16,ERPd17,ERPd18,ERPd19,ERPd20,ERPd21,ERPd22,ERPd23,ERPd24,ERPd25,ERPd26,ERPd27,ERPd28,ERPd29,ERPd30,ERPd31,ERPd32);
%     ERPd_shk2=vertcat(ERPd33,ERPd34,ERPd35,ERPd36,ERPd37,ERPd38,ERPd39,ERPd40,ERPd41,ERPd42,ERPd43,ERPd44,ERPd45,ERPd46,ERPd47,ERPd48,ERPd49,ERPd50,ERPd51,ERPd52,ERPd53,ERPd54,ERPd55,ERPd56,ERPd57,ERPd58,ERPd59,ERPd60,ERPd61,ERPd62,ERPd63,ERPd64);

end

%only first two columns
ERP_shk1 = ERP_shk1(twocols,:);
ERP_shk2 = ERP_shk2(twocols,:);


% trials2 = floor(trials2/(amp_sr/1000)); % downsample


% %% Get the mean of ERP and ERPd
% 
% for i=1:length(ERP_shk1(:,1))
%     ERP_shk1_testmean(i,:) = mean([ERP_shk1(i,:);ERPd_shk1(i,:)]);
% end
% % ERPd_shk1=(aux_shk1d/(1*10^8))*(1*10^6);
% ERP_shk1 = ERP_shk1_testmean;
% 
% % if combine_columns==0
%     for i=1:length(ERP_shk2(:,1))
%         ERP_shk2_testmean(i,:) = mean([ERP_shk2(i,:);ERPd_shk2(i,:)]);
%     end
% %     for i=1:length(ERP_shk3(:,1))
% %         ERP_shk3_testmean(i,:) = mean([ERP_shk3(i,:);ERPd_shk3(i,:)]);
% %     end
%     ERP_shk2 = ERP_shk2_testmean;
% %     ERP_shk3 = ERP_shk3_testmean;
% % end

%%
ERP_shk2_clean = ERP_shk2;
ERP_shk2_clean(:,500:510)=0;
ERP_shk2_clean(:,511:end)=ERP_shk2_clean(:,511:end)-repmat(ERP_shk2(:,510),1,size(ERP_shk2_clean(:,511:end),2));
figure;plot(ERP_shk2_clean')

ERP_shk1_clean = ERP_shk1;
ERP_shk1_clean(:,500:510)=0;
ERP_shk1_clean(:,511:end)=ERP_shk1_clean(:,511:end)-repmat(ERP_shk1(:,510),1,size(ERP_shk1_clean(:,511:end),2));
ERP_shk1 = ERP_shk1_clean;
ERP_shk2 = ERP_shk2_clean;
%% Plot and clean up ERP figures %%%MEAN NORMALIZATION

probemap_chrem = probemap; %map of removed channels = nan
% probemapd_chrem = probemap;

%TMP
timing = time_index(1:1001)-.5;

[ERP_shk1 probemap_chrem]=clean_ERP_MAK(ERP_shk1,timing,probemap_chrem,1);
ERP_shk1_nonan=ERP_shk1;

[ERP_shk2 probemap_chrem]=clean_ERP_MAK(ERP_shk2,timing,probemap_chrem,2);
% [ERP_shk3 probemap_chrem]=clean_ERP_MAK(ERP_shk3,timing,probemap_chrem,3);

% remove nan's by replacing with the average of the nearest neighbors
%check if there are non-nan channels adjacent
nanidx1 = find(isnan(ERP_shk1(:,1)));
nonnanidx1 = find(~isnan(ERP_shk1(:,1)));
adjidx1 = zeros(length(nanidx1),2);
for i=1:length(nanidx1)
    temp = nonnanidx1-nanidx1(i);%distance from nan chan to other chans
%     if nanidx1(i) == 1 %if it's channel 1
%         adjidx1(i,1) = nanidx1(i); 
%     else
        [tempmin tempidx] = min(abs(temp(temp<0)));%find first min distance
        if ~isempty(tempmin)
            adjidx1(i,1) = nonnanidx1(tempidx);
            temp(1:tempidx) = nan; %make first min nan 
        else
            adjidx1(i,1) = nan;
        end
        
%     end
%     if nanidx1(i) == length(ERP_shk1(:,1)) %if it's the last channel
%         adjidx1(i,2) = nanidx1(i);
%     else
        [tempmin tempidx] = min(abs(temp));%find second min distance
        if ~isempty(tempmin) && ~isnan(tempmin)
            adjidx1(i,2) = nonnanidx1(tempidx);
        else
            adjidx1(i,2) = nan;
        end
%     end
end
%find the distance to the nearest non-nan channels and then do weighted averages!
%shank1
for i=1:length(nanidx1)
    if isnan(adjidx1(i,1)) %if channels at end of probe are removed
        ERP_shk1_nonan(nanidx1(i),:) = ERP_shk1(adjidx1(i,2),:); %just duplicate nearest neighbor
    elseif isnan(adjidx1(i,2))
        ERP_shk1_nonan(nanidx1(i),:) = ERP_shk1(adjidx1(i,1),:);
    elseif adjidx1(i,2)-adjidx1(i,1)>=2 %if nan chan is between channels
        d1 = nanidx1(i)-adjidx1(i,1); %distances to nearest non-nan channels
        d2 = adjidx1(i,2)-nanidx1(i);
        if d1 == d2
            ERP_shk1_nonan(nanidx1(i),:) = mean([ERP_shk1(adjidx1(i,1),:);ERP_shk1(adjidx1(i,2),:)],1);
        else
            %WEIGHTED AVERAGE....NOT SURE I DID THIS RIGHT!?!
            ERP_shk1_nonan(nanidx1(i),:) = (ERP_shk1(adjidx1(i,1),:)*(1/d1) + ERP_shk1(adjidx1(i,2),:)*(1/d2))/(1/d1+1/d2);
        end
    elseif adjidx1(i,2)-adjidx1(i,1)==1 %if nan chan is at the end of the probe
        if nanidx1(i) == 1
            ERP_shk1_nonan(nanidx1(i),:) = ERP_shk1(adjidx1(i,2),:); %just duplicate nearest neighbor
        elseif nanidx1(i) == length(ERP_shk1(:,1))
            ERP_shk1_nonan(nanidx1(i),:) = ERP_shk1(adjidx1(i,1),:);
        end
    elseif adjidx1(i,2)-adjidx1(i,1)==0
        %%%
    end
end
% normalize ERPS i.e. subtract common noise (from monitor?)<< but already
% did this!! (MAK)
% if combine_columns == 0
    ERP_shk1_mean = nanmean(ERP_shk1_nonan,1);
    
%     for i = 1:length(ERP_shk1_nonan(:,1))
%         if isnan(ERP_shk1_nonan(i,1))
%             ERP_shk1_norm(i,1:length(ERP_shk1(1,:))) = 0;
%         else
%             ERP_shk1_norm(i,1:length(ERP_shk1(1,:))) = ERP_shk1_nonan(i,:) - ERP_shk1_mean;
%         end
%     end

%try normalizing each column before combining    
% else
%     ERP_shk1_mean1 = nanmedian(ERP_shk1_nonan(pos_col==1),1);
%     ERP_shk1_mean2 = nanmedian(ERP_shk1_nonan(pos_col==2),1);
%     
%     for i = 1:length(ERP_shk1_nonan(:,1))
%         if isnan(ERP_shk1_nonan(i,1))
%             ERP_shk1_norm(i,1:length(ERP_shk1(1,:))) = 0;  
%         else
%             if pos_col(i)==1
%                 ERP_shk1_norm(i,1:length(ERP_shk1(1,:))) = ERP_shk1_nonan(i,:) - ERP_shk1_mean1;
%             elseif pos_col(i)==0
%                 ERP_shk1_norm(i,1:length(ERP_shk1(1,:))) = ERP_shk1_nonan(i,:) - ERP_shk1_mean2;
%             end
%         end
%     end
%     
% end


% if combine_columns==0
    
% [ERPd_shk1 probemapd_chrem]=clean_ERP_EGM(ERPd_shk1,timing,probemapd_chrem,1);
% [ERPd_shk2 probemapd_chrem]=clean_ERP_EGM(ERPd_shk2,timing,probemapd_chrem,2);
% [ERPd_shk3 probemapd_chrem]=clean_ERP_EGM(ERPd_shk3,timing,probemapd_chrem,3);

% for i=1:length(ERP_shk1(:,1))
%     ERP_shk1_both(i,:) = nanmedian([ERP_shk1(i,:);ERPd_shk1(i,:)]);
% end
% for i=1:length(ERP_shk2(:,1))
%     ERP_shk2_both(i,:) = nanmedian([ERP_shk2(i,:);ERPd_shk2(i,:)]);
% end
% for i=1:length(ERP_shk3(:,1))
%     ERP_shk3_both(i,:) = nanmedian([ERP_shk3(i,:);ERPd_shk3(i,:)]);
% end
% 
% ERP_shk1 = ERP_shk1_both;
% ERP_shk2 = ERP_shk2_both;
% ERP_shk3 = ERP_shk3_both;

% % normalize ERPS i.e. subtract common noise (from monitor?)
% ERP_shk1_mean = nanmean(ERP_shk1,1);
% ERP_shk2_mean = nanmean(ERP_shk2,1);
% ERP_shk3_mean = nanmean(ERP_shk3,1);
% ERP_shk_all = [ERP_shk1; ERP_shk2; ERP_shk3];
% ERP_shk_all_mean = nanmean(ERP_shk_all,1);

% remove nan's by replacing with the average of the nearest neighbors

%check if there are non-nan channels adjacent
ERP_shk2_nonan=ERP_shk2;
% ERP_shk3_nonan=ERP_shk3;

nanidx2 = find(isnan(ERP_shk2(:,1)));
nonnanidx2 = find(~isnan(ERP_shk2(:,1)));
adjidx2 = zeros(length(nanidx2),2);
for i=1:length(nanidx2)
    temp = nonnanidx2-nanidx2(i);%distance from nan chan to other chans
    [tempmin tempidx] = min(abs(temp(temp<0)));%find first min distance
    if ~isempty(tempmin)
        adjidx2(i,1) = nonnanidx2(tempidx);
        temp(1:tempidx) = nan; %make first min nan
    else
        adjidx2(i,1) = nan;
    end
    [tempmin tempidx] = min(abs(temp));%find second min distance
    if ~isempty(tempmin) && ~isnan(tempmin)
        adjidx2(i,2) = nonnanidx2(tempidx);
    else
        adjidx2(i,2) = nan;
    end
end

% nanidx3 = find(isnan(ERP_shk3(:,1)));
% nonnanidx3 = find(~isnan(ERP_shk3(:,1)));
% adjidx3 = zeros(length(nanidx3),2);
% for i=1:length(nanidx3)
%     temp = nonnanidx3-nanidx3(i);%distance from nan chan to other chans
%     [tempmin tempidx] = min(abs(temp(temp<0)));%find first min distance
%     if ~isempty(tempmin)
%         adjidx3(i,1) = nonnanidx3(tempidx);
%         temp(1:tempidx) = nan; %make first min nan
%     else
%         adjidx3(i,1) = nan;
%     end
%     [tempmin tempidx] = min(abs(temp));%find second min distance
%     if ~isempty(tempmin) && ~isnan(tempmin)
%         adjidx3(i,2) = nonnanidx3(tempidx);
%     else
%         adjidx3(i,2) = nan;
%     end
% end

%find the distance to the nearest non-nan channels and then do weighted averages!


% %if there are many consecutive nan channels just delete them!
% consec1 = nanidx1(find(diff(nanidx1)==1));
% if ~isempty(consec1)
%     ERP_shk1_nonan([consec1; consec1(end)+1],:)=[];
% end

%shank2
for i=1:length(nanidx2)
    if isnan(adjidx2(i,1)) %if channels at end of probe are removed
        ERP_shk2_nonan(nanidx2(i),:) = ERP_shk2(adjidx2(i,2),:); %just duplicate nearest neighbor
    elseif isnan(adjidx2(i,2))
        ERP_shk2_nonan(nanidx2(i),:) = ERP_shk2(adjidx2(i,1),:);
    elseif adjidx2(i,2)-adjidx2(i,1)>=2 %if nan chan is between channels
        d1 = nanidx2(i)-adjidx2(i,1); %distances to nearest non-nan channels
        d2 = adjidx2(i,2)-nanidx2(i);
        if d1 == d2
            ERP_shk2_nonan(nanidx2(i),:) = mean([ERP_shk2(adjidx2(i,1),:);ERP_shk2(adjidx2(i,2),:)],1);
        else
            %WEIGHTED AVERAGE....NOT SURE I DID THIS RIGHT!?!
            ERP_shk2_nonan(nanidx2(i),:) = (ERP_shk2(adjidx2(i,1),:)*(1/d1) + ERP_shk2(adjidx2(i,2),:)*(1/d2))/(1/d1+1/d2);
        end
    else %if nan chan is at the end of the probe
        if nanidx2(i) == 1
            ERP_shk2_nonan(nanidx2(i),:) = ERP_shk2(adjidx2(i,2),:); %just duplicate nearest neighbor
        elseif nanidx2(i) == length(ERP_shk2(:,1))
            ERP_shk2_nonan(nanidx2(i),:) = ERP_shk2(adjidx2(i,1),:);
        end
    end
end
% %if there are many consecutive nan channels just delete them!
% consec2 = nanidx2(find(diff(nanidx2)==1));
% if ~isempty(consec2)
%     ERP_shk2_nonan([consec2; consec2(end)+1],:)=[];
% end

% %shank3
% for i=1:length(nanidx3)
%     if isnan(adjidx3(i,1)) %if channels at end of probe are removed
%         ERP_shk3_nonan(nanidx3(i),:) = ERP_shk3(adjidx3(i,2),:); %just duplicate nearest neighbor
%     elseif isnan(adjidx3(i,2))
%         ERP_shk3_nonan(nanidx3(i),:) = ERP_shk3(adjidx3(i,1),:);
%     elseif adjidx3(i,2)-adjidx3(i,1)>=2 %if nan chan is between channels
%         d1 = nanidx3(i)-adjidx3(i,1); %distances to nearest non-nan channels
%         d2 = adjidx3(i,2)-nanidx3(i);
%         if d1 == d2
%             ERP_shk3_nonan(nanidx3(i),:) = mean([ERP_shk3(adjidx3(i,1),:);ERP_shk3(adjidx3(i,2),:)],1);
%         else
%             %WEIGHTED AVERAGE....NOT SURE I DID THIS RIGHT!?!
%             ERP_shk3_nonan(nanidx3(i),:) = (ERP_shk3(adjidx3(i,1),:)*(1/d1) + ERP_shk3(adjidx3(i,2),:)*(1/d2))/(1/d1+1/d2);
%         end
%     else %if nan chan is at the end of the probe
%         if nanidx3(i) == 1
%             ERP_shk3_nonan(nanidx3(i),:) = ERP_shk3(adjidx3(i,2),:); %just duplicate nearest neighbor
%         elseif nanidx3(i) == length(ERP_shk3(:,1))
%             ERP_shk3_nonan(nanidx3(i),:) = ERP_shk3(adjidx3(i,1),:);
%         end
%     end
% end
% % consec3 = nanidx3(find(diff(nanidx3)==1));
% % if ~isempty(consec3)
% %     ERP_shk3_nonan([consec3; consec3(end)+1],:)=[];
% % end

% normalize ERPS i.e. subtract common noise (from monitor?)

ERP_shk2_mean = nanmean(ERP_shk2_nonan,1);
% ERP_shk3_mean = nanmean(ERP_shk3_nonan,1);
ERP_shk_all = [ERP_shk1_nonan; ERP_shk2_nonan];
% ERP_shk_all = [ERP_shk1_nonan; ERP_shk2_nonan; ERP_shk3_nonan];
ERP_shk_all_mean = nanmean(ERP_shk_all,1);

% for i = 1:length(ERP_shk2_nonan(:,1))
%     if isnan(ERP_shk2_nonan(i,1))
%         ERP_shk2_norm(i,1:length(ERP_shk2(1,:))) = 0;
%     else
%         ERP_shk2_norm(i,1:length(ERP_shk2(1,:))) = ERP_shk2_nonan(i,:) - ERP_shk2_mean;
%     end
% end
% for i = 1:length(ERP_shk3_nonan(:,1))
%     if isnan(ERP_shk3_nonan(i,1))
%         ERP_shk3_norm(i,1:length(ERP_shk3(1,:))) = 0;
%     else
%         ERP_shk3_norm(i,1:length(ERP_shk3(1,:))) = ERP_shk3_nonan(i,:) - ERP_shk3_mean;
%     end
% end

% ERP_columns12([1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41],:) = ERP_shk1_norm;
% ERP_columns12([2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,43],:) = ERP_shk2_norm;

% end
%%%FOR NOW JUST MAKING REMOVED CHANNELS == 0 ... THERE'S PROBABLY A BETTER WAY!!



%% The final ERP figure is sent to the screen to view and now the next section deals
%  with the CSD plotter and average analysis
% pause(5);
file_name = 'CSD_LFPdata';
save (file_name,'ERP_shk1_nonan','ERP_shk2_nonan');
disp ('saving mat')
HC=CSDplotter;
disp ('Use CSD plotter')
clc
pause(20)
disp ('ERP: normal LFP time series w/o moving trials')
reply = input('Ready to continue? Y/N: ','s');
if reply == 'Y';
    close all;
    disp ('Open CSD matrix file for Shank 1')
    filename = uigetfile;
    load (filename);
    CSD_matrix1=CSD_matrix;
    fprintf('%s%s','Filename:   ',filename);
    close all
    clc
    disp ('Plot final ERP and CSD side by side')
    X2=timing.*1000;
    clear CSD_matrix
    
    if combine_columns == 0
    disp ('Open CSD matrix file for Shank 2')
    filename = uigetfile;
    load (filename);
    CSD_matrix2=CSD_matrix;
    fprintf('%s%s','Filename:   ',filename);
    close all
    clc
    
%     disp ('Open CSD matrix file for Shank 3')
%     filename = uigetfile;
%     load (filename);
%     CSD_matrix3=CSD_matrix;
%     fprintf('%s%s','Filename:   ',filename);
%     close all
%     clc
    end
    
    disp ('Plot final ERP and CSD side by side')
    %% Plot ERPs
    FontName = 'MyriadPro-Regular'; % or choose any other font
    FontSize = 14;
    figure_width = 14;
    figure_height = 10;
    figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
    ERP_fig1 = figure('Name','ERP1');
    set(ERP_fig1,'Visible', figuresVisible)
    set(ERP_fig1, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
    set(ERP_fig1, 'PaperPositionMode', 'auto');
    set(ERP_fig1, 'Color', [1 1 1]); % Sets figure background
    set(ERP_fig1, 'Color', [1 1 1]); % Sets axes background
    hsp = subplot(1,1,1, 'Parent', ERP_fig1);
    set(hsp,'Position',[0.15 0.17 0.75 0.80]);
    plot (X2, ERP_shk1_nonan');
    axis on;      % display axis
    axis tight;   % no white borders
    set(gca, ...
        'Box'         , 'off'      , ...
        'TickDir'     , 'out'      , ...
        'TickLength'  , [.015 .015] , ...
        'XMinorTick'  , 'off'      , ...
        'YMinorTick'  , 'off'     , ...
        'XGrid'       , 'off'     , ...
        'YGrid'       , 'off'     , ...
        'XColor'      , [.0 .0 .0], ...
        'YColor'      , [.0 .0 .0], ...
        'LineWidth'   , 0.6        );
    set(gca,'Xlim',[-25 225]);
    yaxis=ylim;
    set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
    prestim_offset_y            = yaxis(1):1:yaxis(2);
    prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
    hold on;plot(prestim_offset_t, prestim_offset_y, 'k','linewidth', 1);
    xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
    yLabelText = 'Amplitude (uV)';
    % save handles to set up label properties
    hXLabel = xlabel(xLabelText);
    hYLabel = ylabel(yLabelText);
    set([gca, hXLabel, hYLabel], ...
        'FontSize'   , FontSize    , ...
        'FontName'   , FontName);
    fig_title=sprintf('%s','ERP 1 ');
    set(gca,'Layer', 'top');
    drawnow
    %export_fig (fig_title, '-pdf')
    %         export_fig(ERP_fig1, '-png','-r600','-zbuffer');
    
    if combine_columns==0
    FontName = 'MyriadPro-Regular'; % or choose any other font
    FontSize = 14;
    figure_width = 14;
    figure_height = 10;
    figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
    ERP_fig2 = figure('Name','ERP2');
    set(ERP_fig2,'Visible', figuresVisible)
    set(ERP_fig2, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
    set(ERP_fig2, 'PaperPositionMode', 'auto');
    set(ERP_fig2, 'Color', [1 1 1]); % Sets figure background
    set(ERP_fig2, 'Color', [1 1 1]); % Sets axes background
    hsp = subplot(1,1,1, 'Parent', ERP_fig2);
    set(hsp,'Position',[0.19 0.19 0.75 0.80]);
    plot (X2, ERP_shk2_nonan');
    axis on;      % display axis
    axis tight;   % no white borders
    set(gca, ...
        'Box'         , 'off'      , ...
        'TickDir'     , 'out'      , ...
        'TickLength'  , [.015 .015] , ...
        'XMinorTick'  , 'off'      , ...
        'YMinorTick'  , 'off'     , ...
        'XGrid'       , 'off'     , ...
        'YGrid'       , 'off'     , ...
        'XColor'      , [.0 .0 .0], ...
        'YColor'      , [.0 .0 .0], ...
        'LineWidth'   , 0.6        );
    set(gca,'Xlim',[-25 225]);
    yaxis=ylim;
    set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
    prestim_offset_y            = yaxis(1):1:yaxis(2);
    prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
    hold on;plot(prestim_offset_t, prestim_offset_y, 'k','linewidth', 1);
    xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
    yLabelText = 'Amplitude (uV)';
    % save handles to set up label properties
    hXLabel = xlabel(xLabelText);
    hYLabel = ylabel(yLabelText);
    set([gca, hXLabel, hYLabel], ...
        'FontSize'   , FontSize    , ...
        'FontName'   , FontName);
    fig_title=sprintf('%s','ERP 2 ');
    set(gca,'Layer', 'top');
    drawnow
    %     export_fig (fig_title, '-pdf')
    %         export_fig (ERP_fig2, '-png','-r600','-opengl')
    
    FontName = 'MyriadPro-Regular'; % or choose any other font
    FontSize = 14;
    figure_width = 14;
    figure_height = 10;
    figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
    ERP_fig3 = figure('Name','ERP3');
    set(ERP_fig3,'Visible', figuresVisible)
    set(ERP_fig3, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
    set(ERP_fig3, 'PaperPositionMode', 'auto');
    set(ERP_fig3, 'Color', [1 1 1]); % Sets figure background
    set(ERP_fig3, 'Color', [1 1 1]); % Sets axes background
    hsp = subplot(1,1,1, 'Parent', ERP_fig3);
    set(hsp,'Position',[0.19 0.19 0.75 0.80]);
    plot (X2, ERP_shk3_nonan');
    axis on;      % display axis
    axis tight;   % no white borders
    set(gca, ...
        'Box'         , 'off'      , ...
        'TickDir'     , 'out'      , ...
        'TickLength'  , [.015 .015] , ...
        'XMinorTick'  , 'off'      , ...
        'YMinorTick'  , 'off'     , ...
        'XGrid'       , 'off'     , ...
        'YGrid'       , 'off'     , ...
        'XColor'      , [.0 .0 .0], ...
        'YColor'      , [.0 .0 .0], ...
        'LineWidth'   , 0.6        );
    set(gca,'Xlim',[-25 225]);
    yaxis=ylim;
    set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
    prestim_offset_y            = yaxis(1):1:yaxis(2);
    prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
    hold on;plot(prestim_offset_t, prestim_offset_y, 'k','linewidth', 1);
    xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
    yLabelText = 'Amplitude (uV)';
    % save handles to set up label properties
    hXLabel = xlabel(xLabelText);
    hYLabel = ylabel(yLabelText);
    set([gca, hXLabel, hYLabel], ...
        'FontSize'   , FontSize    , ...
        'FontName'   , FontName);
    fig_title=sprintf('%s','ERP 3 ');
    set(gca,'Layer', 'top');
    drawnow
    end
    %     export_fig (fig_title, '-pdf')
    %         export_fig (ERP_fig3, '-png','-r600','-opengl')
    
    %% Plot ERPs stacked
    
    if combine_columns==1
        channels= fliplr([1:1:size(ERP_shk1,1)]);
        FontName = 'MyriadPro-Regular'; % or choose any other font
        FontSize = 14;
        figure_width = 14;
        figure_height = 14;
        figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
        ERP_stacked=figure('Name','ERP_stacked');% figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        set(ERP_stacked, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
        set(ERP_stacked, 'PaperPositionMode', 'auto');
        set(ERP_stacked, 'Color', [1 1 1]); % Sets figure background
        set(ERP_stacked, 'Color', [1 1 1]); % Sets axes background
        hsp = subplot(1,1,1, 'Parent', ERP_stacked);
        set(hsp,'Position',[0.15 0.17 0.75 0.80]);
        depth           = 0;
        depth_spacing   = 25;
        max_y           = 0;
        hold all;
        % for each channel
        chan_legends    = {};
        for j=1:length(channels);
            for i = channels(j)
                averaged_ERP= ERP_shk1_nonan(i,:);
                plot(X2, averaged_ERP*depth_spacing + depth,'LineWidth',2);
                max_y= max_y + max(averaged_ERP);
                depth= depth + depth_spacing;
                chan_legends= [chan_legends, num2str(i)];
            end
        end
        yaxis=ylim;
        set(gca,'Ylim',[-depth_spacing depth_spacing*length(channels)])
        
        set(gca, 'ytick', [0:depth_spacing:(depth_spacing*length(channels))-depth_spacing],'tickdir','out','yticklabel',[channels]);
        axis on;      % display axis
        set(gca, ...
            'Box'         , 'off'      , ...
            'TickDir'     , 'out'      , ...
            'TickLength'  , [0 0] , ...
            'XMinorTick'  , 'off'      , ...
            'YMinorTick'  , 'off'     , ...
            'XGrid'       , 'off'     , ...
            'YGrid'       , 'off'     , ...
            'XColor'      , [.0 .0 .0], ...
            'YColor'      , [.0 .0 .0], ...
            'LineWidth'   , 0.6        );
        set(gca,'Xlim',[-25 225]);
        set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
        xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
        yLabelText = 'Electrode number (Sup-->Deep)';
        % save handles to set up label properties
        hXLabel = xlabel(xLabelText);
        hYLabel = ylabel(yLabelText);
        set([gca, hXLabel, hYLabel], ...
            'FontSize'   , FontSize    , ...
            'FontName'   , FontName);
        prestim_offset_y            = yaxis(1):1:yaxis(2);
        prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
        plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',2);
        % poststimulus onset
    else
        left_channels= fliplr([1:1:size(ERP_shk1,1)]);
        right_channels= fliplr([size(ERP_shk1,1)+1:1:size(ERP_shk2,1)*2]);
        channels={left_channels';right_channels'}';
        
        FontName = 'MyriadPro-Regular'; % or choose any other font
        FontSize = 14;
        figure_width = 28;
        figure_height = 14;
        figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
        ERP_stacked=figure('Name','ERP_stacked');% figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        set(ERP_stacked, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
        set(ERP_stacked, 'PaperPositionMode', 'auto');
        set(ERP_stacked, 'Color', [1 1 1]); % Sets figure background
        set(ERP_stacked, 'Color', [1 1 1]); % Sets axes background
        hsp = subplot(1,1,1, 'Parent', ERP_stacked);
        set(hsp,'Position',[0.15 0.17 0.75 0.80]);
        for ii = 1:size(channels,2)
            subplot(1, size(channels,2), ii);
            depth           = 0;
            depth_spacing   = 50;
            max_y           = 0;
            hold all;
            % for each channel
            chan_legends    = {};
            for j=1:length(channels{ii});
                for i = channels{ii}(j)
                    if ii==1
                        averaged_ERP= (ERP_shk1_nonan(i,:)/(1*10^8))*(1*10^6);      % scale factor for plotting(not totally sure where this came from?)
                    elseif ii==2
                        i=i-size(ERP_shk1,1);
                        averaged_ERP= (ERP_shk2_nonan(i,:)/(1*10^8))*(1*10^6);
                    end
                    
                    plot(X2, averaged_ERP*depth_spacing + depth,'LineWidth',2);
                    max_y= max_y + max(averaged_ERP);
                    depth= depth + depth_spacing;
                    chan_legends= [chan_legends, num2str(i)];
                end
            end
            yaxis=ylim;
            set(gca,'Ylim',[-depth_spacing depth_spacing*length(channels{ii})])
            if ii==1
                set(gca, 'ytick', [0:depth_spacing:(depth_spacing*length(channels{ii}))-depth_spacing],'tickdir','out','yticklabel',[left_channels]);
                axis on;      % display axis
                set(gca, ...
                    'Box'         , 'off'      , ...
                    'TickDir'     , 'out'      , ...
                    'TickLength'  , [0 0] , ...
                    'XMinorTick'  , 'off'      , ...
                    'YMinorTick'  , 'off'     , ...
                    'XGrid'       , 'off'     , ...
                    'YGrid'       , 'off'     , ...
                    'XColor'      , [.0 .0 .0], ...
                    'YColor'      , [.0 .0 .0], ...
                    'LineWidth'   , 0.6        );
                set(gca,'Xlim',[-25 225]);
                set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
                xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
                yLabelText = 'Electrode number (Sup-->Deep)';
                % save handles to set up label properties
                hXLabel = xlabel(xLabelText);
                hYLabel = ylabel(yLabelText);
                set([gca, hXLabel, hYLabel], ...
                    'FontSize'   , FontSize    , ...
                    'FontName'   , FontName);
                prestim_offset_y            = yaxis(1):1:yaxis(2);
                prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
                plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',2);
                % poststimulus onset
            elseif ii==2;
                set(gca, 'ytick', [0:depth_spacing:(depth_spacing*length(channels{ii}))-depth_spacing],'tickdir','out','yticklabel',[right_channels]);
                axis on;      % display axis
                set(gca, ...
                    'Box'         , 'off'      , ...
                    'TickDir'     , 'out'      , ...
                    'TickLength'  , [0 0] , ...
                    'XMinorTick'  , 'off'      , ...
                    'YMinorTick'  , 'off'     , ...
                    'XGrid'       , 'off'     , ...
                    'YGrid'       , 'off'     , ...
                    'XColor'      , [.0 .0 .0], ...
                    'YColor'      , [.0 .0 .0], ...
                    'LineWidth'   , 0.6        );
                set(gca,'Xlim',[-25 225]);
                set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
                xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
                yLabelText = 'Electrode number (Sup-->Deep)';
                % save handles to set up label properties
                hXLabel = xlabel(xLabelText);
                hYLabel = ylabel(yLabelText);
                set([gca, hXLabel, hYLabel], ...
                    'FontSize'   , FontSize    , ...
                    'FontName'   , FontName);
                prestim_offset_y            = yaxis(1):1:yaxis(2);
                prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
                plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',2);
            elseif ii==3;
                set(gca, 'ytick', [0:depth_spacing:(depth_spacing*length(channels{ii}))-depth_spacing],'tickdir','out','yticklabel',[center_channels]);
                axis on;      % display axis
                set(gca, ...
                    'Box'         , 'off'      , ...
                    'TickDir'     , 'out'      , ...
                    'TickLength'  , [0 0] , ...
                    'XMinorTick'  , 'off'      , ...
                    'YMinorTick'  , 'off'     , ...
                    'XGrid'       , 'off'     , ...
                    'YGrid'       , 'off'     , ...
                    'XColor'      , [.0 .0 .0], ...
                    'YColor'      , [.0 .0 .0], ...
                    'LineWidth'   , 0.6        );
                set(gca,'Xlim',[-25 225]);
                set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
                xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
                yLabelText = 'Electrode number (Sup-->Deep)';
                % save handles to set up label properties
                hXLabel = xlabel(xLabelText);
                hYLabel = ylabel(yLabelText);
                set([gca, hXLabel, hYLabel], ...
                    'FontSize'   , FontSize    , ...
                    'FontName'   , FontName);
                prestim_offset_y            = yaxis(1):1:yaxis(2);
                prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
                plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',2);
            end
        end
    end
    fig_title=sprintf('%s','ERP Stacked ');
    set(gca,'Layer', 'top');
    drawnow
    %export_fig (fig_title, '-pdf')
        export_fig (fig_title, '-png');
    
    %%     Plot CSDs
    
    figure_width = 30;
    figure_height = 14;
    FontSize = 12;
    FontName = 'MyriadPro-Regular'; % or choose any other font
    % --- setup plot windows
    figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
    CSD_fig1 = figure('Name','CSD1');
    set(CSD_fig1,'Visible', figuresVisible)
    set(CSD_fig1, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
    set(CSD_fig1, 'PaperPositionMode', 'auto');
    set(CSD_fig1, 'Renderer','Zbuffer');
    set(CSD_fig1, 'Color', [1 1 1]); % Sets figure background
    set(CSD_fig1, 'Color', [1 1 1]); % Sets axes background
    % --- dimensions and position of plot
%     hsp = subplot(1,1,1, 'Parent', CSD_fig1);
%     set(hsp,'Position',[0.15 0.15 0.60 0.80]);
    colorDepth = 1000;
    colormap(flipud(jet(colorDepth)));
    %     pcolor(X2, 1:1:size(CSD_matrix1,1), CSD_matrix1);
    %     imagesc(X2,[],CSD_matrix1)
    subplot(1, 2, 1);
    title('Column 1')
    hold all;
    imagesc(X2(476:576),[],CSD_matrix1(:,476:576))
    shading interp; % do not interpolate pixels
    axis on; % display axis
    axis tight;% no white borders
    set(gca, ...
        'Box'         , 'off'      , ...
        'TickDir'     , 'in'      , ...
        'Ydir'        , 'reverse', ...
        'TickLength'  , [.01 .01] , ...
        'XMinorTick'  , 'off'      , ...
        'YMinorTick'  , 'off'     , ...
        'XGrid'       , 'off'     , ...
        'YGrid'       , 'off'     , ...
        'XColor'      , [.0 .0 .0], ...
        'YColor'      , [.0 .0 .0], ...
        'LineWidth'   , 0.6        );
    set(gca,'Xlim',[-25 75]);
    set(gca,'XTickLabel',[0 50], 'Xtick', [0 50])
    set(gca, 'yticklabel',1:1:size(CSD_matrix1,1), 'Ytick', 1:1:size(CSD_matrix1,1));
    xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
    yLabelText = 'Electrode number';
    hXLabel = xlabel(xLabelText);
    hYLabel = ylabel(yLabelText);
    fig_title=sprintf('%s','Shank 1 ');
    yaxis=ylim;
    prestim_offset_y            = yaxis(1):1:yaxis(2);
    prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
    hold on;plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',1);
    %     caxis([min(min(CSD_matrix1)) max(max(CSD_matrix1))]);
    zLabelText = 'nA / mm^3';  % greek letters in LaTeX Syntax
%     hcb = colorbar('eastoutside');
%     h_bar = findobj(gcf,'Tag','Colorbar');
%     initpos = get(h_bar,'Position');
%     set(h_bar, ...
%         'Position',[initpos(1)+initpos(3)*2.5 initpos(2)+initpos(4)*0.3 ...
%         initpos(3)*0.4 initpos(4)*0.4]);
%     hcLabel = ylabel(hcb,zLabelText);
%     set(hcb,'YTickLabel',{'Sink','Source'}, 'Ytick', [min(min(CSD_matrix1)) max(max(CSD_matrix1))])
%     set(hcb, ...
%         'Box'         , 'on'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.010 .010] , ...
%         'LineWidth'   , 0.6);
%     set([gca, hcb, hXLabel, hYLabel, hcLabel], ...
%         'FontSize'   , FontSize    , ...
%         'FontName'   , FontName);
%     ylabh=get(hcb,'Ylabel');
%     set(ylabh,'Position',get(ylabh,'Position')-[8 0 0]);
%     set(gca,'Layer', 'top');
%     drawnow
    
%     nrm = input('Do you want to normalize the CSD plot? Y/N:','s');
%     if nrm == 'Y'
%         clear nrm_CSD_matrix
%         prestim_csd =  mean(CSD_matrix1(:,976:1025),2); % prestim = 200ms before flip
%         for c = 1: size(CSD_matrix1,1)
%             nrm_CSD_matrix(c,:) = CSD_matrix1(c,976:end)-prestim_csd(c);
%         end
%         CSD_fig1 = figure('Name','CSD1_norm');
%         set(CSD_fig1,'Visible', figuresVisible)
%         set(CSD_fig1, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
%         set(CSD_fig1, 'PaperPositionMode', 'auto');
%         set(CSD_fig1, 'Renderer','Zbuffer');
%         set(CSD_fig1, 'Color', [1 1 1]); % Sets figure background
%         set(CSD_fig1, 'Color', [1 1 1]); % Sets axes background
%         % --- dimensions and position of plot
%         hsp = subplot(1,1,1, 'Parent', CSD_fig1);
%         set(hsp,'Position',[0.15 0.15 0.60 0.80]);
%         colorDepth = 1000;
%         colormap(flipud(jet(colorDepth)));
%         %             pcolor(X2(1001:end), 1:1:size(CSD_matrix1,1), nrm_CSD_matrix);
%         imagesc(X2(976:end), [], nrm_CSD_matrix); hold on
%         shading interp; % do not interpolate pixels
%         axis on; % display axis
%         axis tight;% no white borders
%         set(gca, ...
%             'Box'         , 'off'      , ...
%             'TickDir'     , 'in'      , ...
%             'Ydir'        , 'reverse', ...
%             'TickLength'  , [.01 .01] , ...
%             'XMinorTick'  , 'off'      , ...
%             'YMinorTick'  , 'off'     , ...
%             'XGrid'       , 'off'     , ...
%             'YGrid'       , 'off'     , ...
%             'XColor'      , [.0 .0 .0], ...
%             'YColor'      , [.0 .0 .0], ...
%             'LineWidth'   , 0.6        );
%         set(gca,'Xlim',[-25 225]);
%         set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
%         set(gca, 'yticklabel',1:1:size(CSD_matrix1,1), 'Ytick', 1:1:size(CSD_matrix1,1));
%         xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
%         yLabelText = 'Electrode number';
%         hXLabel = xlabel(xLabelText);
%         hYLabel = ylabel(yLabelText);
%         fig_title=sprintf('%s','Shank 1 - Normalized');
%         yaxis=ylim;
%         prestim_offset_y            = yaxis(1):1:yaxis(2);
%         prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
%         hold on;plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',1);
%         caxis([min(min(CSD_matrix1)) max(max(CSD_matrix1))]);
%         zLabelText = 'nA / mm^3';  % greek letters in LaTeX Syntax
%         hcb = colorbar('eastoutside');
%         h_bar = findobj(gcf,'Tag','Colorbar');
%         initpos = get(h_bar,'Position');
%         set(h_bar, ...
%             'Position',[initpos(1)+initpos(3)*2.5 initpos(2)+initpos(4)*0.3 ...
%             initpos(3)*0.4 initpos(4)*0.4]);
%         hcLabel = ylabel(hcb,zLabelText);
%         set(hcb,'YTickLabel',{'Sink','Source'}, 'Ytick', [min(min(CSD_matrix1)) max(max(CSD_matrix1))])
%         set(hcb, ...
%             'Box'         , 'on'     , ...
%             'TickDir'     , 'in'     , ...
%             'TickLength'  , [.010 .010] , ...
%             'LineWidth'   , 0.6);
%         set([gca, hcb, hXLabel, hYLabel, hcLabel], ...
%             'FontSize'   , FontSize    , ...
%             'FontName'   , FontName);
%         ylabh=get(hcb,'Ylabel');
%         set(ylabh,'Position',get(ylabh,'Position')-[8 0 0]);
%         set(gca,'Layer', 'top');
%         plot([25 25],get(gca,'ylim'),'-k','Linewidth',.5);
%         drawnow
%         hold on
%     end
    %         export_fig (fig_title, '-png','-r600','-zbuffer');
    
    if combine_columns==0
%     figure_width = 12;
%     figure_height = 10;
%     FontSize = 12;
%     FontName = 'MyriadPro-Regular'; % or choose any other font
%     % --- setup plot windows
%     CSD_fig2 = figure('Name','CSD2');
%     figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
%     set(CSD_fig2,'Visible', figuresVisible)
%     set(CSD_fig2, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
%     set(CSD_fig2, 'PaperPositionMode', 'auto');
%     set(CSD_fig2, 'Renderer','Zbuffer');
%     set(CSD_fig2, 'Color', [1 1 1]); % Sets figure background
%     set(CSD_fig2, 'Color', [1 1 1]); % Sets axes background
%     % --- dimensions and position of plot
%     hsp = subplot(1,1,1, 'Parent', CSD_fig2);
%     set(hsp,'Position',[0.15 0.15 0.60 0.80]);
    colorDepth = 1000;
    colormap(flipud(jet(colorDepth)));
    %         pcolor(X2, 1:1:size(CSD_matrix2,1), CSD_matrix2);
    subplot(1, 2, 2);
    imagesc(X2(476:576),[],CSD_matrix2(:,476:576))
    shading interp; % do not interpolate pixels
    axis on; % display axis
    axis tight;% no white borders
    set(gca, ...
        'Box'         , 'off'      , ...
        'TickDir'     , 'in'      , ...
        'Ydir'        , 'reverse', ...
        'TickLength'  , [.01 .01] , ...
        'XMinorTick'  , 'off'      , ...
        'YMinorTick'  , 'off'     , ...
        'XGrid'       , 'off'     , ...
        'YGrid'       , 'off'     , ...
        'XColor'      , [.0 .0 .0], ...
        'YColor'      , [.0 .0 .0], ...
        'LineWidth'   , 0.6        );
    set(gca,'Xlim',[-25 75]);
    set(gca,'XTickLabel',[0 50], 'Xtick', [0 50])
    set(gca, 'yticklabel',1:1:size(CSD_matrix2,1), 'Ytick', 1:1:size(CSD_matrix2,1));
    xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
    yLabelText = 'Electrode number';
    hXLabel = xlabel(xLabelText);
    hYLabel = ylabel(yLabelText);
%     fig_title=sprintf('%s','Shank 2 ');
    yaxis=ylim;
    prestim_offset_y            = yaxis(1):1:yaxis(2);
    prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
    hold on;plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',1);
    %         caxis([min(min(CSD_matrix2)) max(max(CSD_matrix2))]);
    zLabelText = 'nA / mm^3';  % greek letters in LaTeX Syntax
%     hcb = colorbar('eastoutside');
%     h_bar = findobj(gcf,'Tag','Colorbar');
%     initpos = get(h_bar,'Position');
%     set(h_bar, ...
%         'Position',[initpos(1)+initpos(3)*2.5 initpos(2)+initpos(4)*0.3 ...
%         initpos(3)*0.4 initpos(4)*0.4]);
%     hcLabel = ylabel(hcb,zLabelText);
%     set(hcb,'YTickLabel',{'Sink','Source'}, 'Ytick', [min(min(CSD_matrix2)) max(max(CSD_matrix2))])
%     set(hcb, ...
%         'Box'         , 'on'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.010 .010] , ...
%         'LineWidth'   , 0.6);
%     set([gca, hXLabel, hYLabel], ...
%         'FontSize'   , FontSize    , ...
%         'FontName'   , FontName);
%     set([hcb, hcLabel], ...
%         'FontSize'   , 10    , ...
%         'FontName'   , FontName);
%     ylabh=get(hcb,'Ylabel');
%     set(ylabh,'Position',get(ylabh,'Position')-[8 0 0]);
%     set(gca,'Layer', 'top');
%     drawnow
%     
%     nrm = input('Do you want to normalize the CSD plot? Y/N:','s');
%     clear nrm_CSD_matrix
%     if nrm == 'Y'
%         prestim_csd =  mean(CSD_matrix2(:,900:1000),2); % prestim = 200ms before flip
%         for c = 1: size(CSD_matrix2,1)
%             nrm_CSD_matrix(c,:) = CSD_matrix2(c,1001:end)-prestim_csd(c);
%         end
%         CSD_fig2 = figure('Name','CSD2_norm');
%         set(CSD_fig2,'Visible', figuresVisible)
%         set(CSD_fig2, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
%         set(CSD_fig2, 'PaperPositionMode', 'auto');
%         set(CSD_fig2, 'Renderer','Zbuffer');
%         set(CSD_fig2, 'Color', [1 1 1]); % Sets figure background
%         set(CSD_fig2, 'Color', [1 1 1]); % Sets axes background
%         % --- dimensions and position of plot
%         hsp = subplot(1,1,1, 'Parent', CSD_fig1);
%         set(hsp,'Position',[0.15 0.15 0.60 0.80]);
%         colorDepth = 1000;
%         colormap(flipud(jet(colorDepth)));
%         %             pcolor(X2(1001:end), 1:1:size(CSD_matrix1,1), nrm_CSD_matrix);
%         imagesc(X2(1001:end), [], nrm_CSD_matrix);
%         shading interp; % do not interpolate pixels
%         axis on; % display axis
%         axis tight;% no white borders
%         set(gca, ...
%             'Box'         , 'off'      , ...
%             'TickDir'     , 'in'      , ...
%             'Ydir'        , 'reverse', ...
%             'TickLength'  , [.01 .01] , ...
%             'XMinorTick'  , 'off'      , ...
%             'YMinorTick'  , 'off'     , ...
%             'XGrid'       , 'off'     , ...
%             'YGrid'       , 'off'     , ...
%             'XColor'      , [.0 .0 .0], ...
%             'YColor'      , [.0 .0 .0], ...
%             'LineWidth'   , 0.6        );
%         set(gca,'Xlim',[-25 225]);
%         set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
%         set(gca, 'yticklabel',1:1:size(CSD_matrix2,1), 'Ytick', 1:1:size(CSD_matrix2,1));
%         xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
%         yLabelText = 'Electrode number';
%         hXLabel = xlabel(xLabelText);
%         hYLabel = ylabel(yLabelText);
%         fig_title=sprintf('%s','Shank 1 - Normalized');
%         yaxis=ylim;
%         prestim_offset_y            = yaxis(1):1:yaxis(2);
%         prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
%         hold on;plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',1);
%         caxis([min(min(CSD_matrix2)) max(max(CSD_matrix2))]);
%         zLabelText = 'nA / mm^3';  % greek letters in LaTeX Syntax
%         hcb = colorbar('eastoutside');
%         h_bar = findobj(gcf,'Tag','Colorbar');
%         initpos = get(h_bar,'Position');
%         set(h_bar, ...
%             'Position',[initpos(1)+initpos(3)*2.5 initpos(2)+initpos(4)*0.3 ...
%             initpos(3)*0.4 initpos(4)*0.4]);
%         hcLabel = ylabel(hcb,zLabelText);
%         set(hcb,'YTickLabel',{'Sink','Source'}, 'Ytick', [min(min(CSD_matrix2)) max(max(CSD_matrix2))])
%         set(hcb, ...
%             'Box'         , 'on'     , ...
%             'TickDir'     , 'in'     , ...
%             'TickLength'  , [.010 .010] , ...
%             'LineWidth'   , 0.6);
%         set([gca, hcb, hXLabel, hYLabel, hcLabel], ...
%             'FontSize'   , FontSize    , ...
%             'FontName'   , FontName);
%         ylabh=get(hcb,'Ylabel');
%         set(ylabh,'Position',get(ylabh,'Position')-[8 0 0]);
%         set(gca,'Layer', 'top');
%         drawnow
%     end
%             export_fig (fig_title, '-png','-r600','-zbuffer');
    
    %third column
%     figure_width = 12;
%     figure_height = 10;
%     FontSize = 12;
%     FontName = 'MyriadPro-Regular'; % or choose any other font
%     % --- setup plot windows
%     figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
%     CSD_fig3 = figure('Name','CSD3');
%     set(CSD_fig3,'Visible', figuresVisible)
%     set(CSD_fig3, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
%     set(CSD_fig3, 'PaperPositionMode', 'auto');
%     set(CSD_fig3, 'Renderer','Zbuffer');
%     set(CSD_fig3, 'Color', [1 1 1]); % Sets figure background
%     set(CSD_fig3, 'Color', [1 1 1]); % Sets axes background
%     % --- dimensions and position of plot
%     hsp = subplot(1,1,1, 'Parent', CSD_fig3);
%     set(hsp,'Position',[0.15 0.15 0.60 0.80]);
%     subplot(1, 3, 3);
%     colorDepth = 1000;
%     colormap(flipud(jet(colorDepth)));
%     %         pcolor(X2, 1:1:size(CSD_matrix2,1), CSD_matrix2);
%     imagesc(X2(476:726),[],CSD_matrix3(:,476:726))
%     shading interp; % do not interpolate pixels
%     axis on; % display axis
%     axis tight;% no white borders
%     set(gca, ...
%         'Box'         , 'off'      , ...
%         'TickDir'     , 'in'      , ...
%         'Ydir'        , 'reverse', ...
%         'TickLength'  , [.01 .01] , ...
%         'XMinorTick'  , 'off'      , ...
%         'YMinorTick'  , 'off'     , ...
%         'XGrid'       , 'off'     , ...
%         'YGrid'       , 'off'     , ...
%         'XColor'      , [.0 .0 .0], ...
%         'YColor'      , [.0 .0 .0], ...
%         'LineWidth'   , 0.6        );
%     set(gca,'Xlim',[-25 225]);
%     set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
%     set(gca, 'yticklabel',1:1:size(CSD_matrix3,1), 'Ytick', 1:1:size(CSD_matrix3,1));
%     xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
%     yLabelText = 'Electrode number';
%     hXLabel = xlabel(xLabelText);
%     hYLabel = ylabel(yLabelText);
% %     fig_title=sprintf('%s','Shank 3 ');
%     yaxis=ylim;
%     prestim_offset_y            = yaxis(1):1:yaxis(2);
%     prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
%     hold on;plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',1);
%     %         caxis([min(min(CSD_matrix2)) max(max(CSD_matrix2))]);
%     zLabelText = 'nA / mm^3';  % greek letters in LaTeX Syntax
%     hcb = colorbar('eastoutside');
%     h_bar = findobj(gcf,'Tag','Colorbar');
%     initpos = get(h_bar,'Position');
%     set(h_bar, ...
%         'Position',[initpos(1)+initpos(3)*2.5 initpos(2)+initpos(4)*0.3 ...
%         initpos(3)*0.4 initpos(4)*0.4]);
%     hcLabel = ylabel(hcb,zLabelText);
%     set(hcb,'YTickLabel',{'Sink','Source'}, 'Ytick', [min(min(CSD_matrix3)) max(max(CSD_matrix3))])
%     set(hcb, ...
%         'Box'         , 'on'     , ...
%         'TickDir'     , 'in'     , ...
%         'TickLength'  , [.010 .010] , ...
%         'LineWidth'   , 0.6);
%     set([gca, hXLabel, hYLabel], ...
%         'FontSize'   , FontSize    , ...
%         'FontName'   , FontName);
%     set([hcb, hcLabel], ...
%         'FontSize'   , 10    , ...
%         'FontName'   , FontName);
%     ylabh=get(hcb,'Ylabel');
%     set(ylabh,'Position',get(ylabh,'Position')-[8 0 0]);
%     set(gca,'Layer', 'top');
%     drawnow
%     
%     nrm = input('Do you want to normalize the CSD plot? Y/N:','s');
%     clear nrm_CSD_matrix
%     if nrm == 'Y'
%         prestim_csd =  mean(CSD_matrix3(:,900:1000),2); % prestim = 200ms before flip
%         for c = 1: size(CSD_matrix3,1)
%             nrm_CSD_matrix(c,:) = CSD_matrix3(c,1001:end)-prestim_csd(c);
%         end
%         CSD_fig3 = figure('Name','CSD3_norm');
%         set(CSD_fig3,'Visible', figuresVisible)
%         set(CSD_fig3, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
%         set(CSD_fig3, 'PaperPositionMode', 'auto');
%         set(CSD_fig3, 'Renderer','Zbuffer');
%         set(CSD_fig3, 'Color', [1 1 1]); % Sets figure background
%         set(CSD_fig3, 'Color', [1 1 1]); % Sets axes background
%         % --- dimensions and position of plot
%         hsp = subplot(1,1,1, 'Parent', CSD_fig1);
%         set(hsp,'Position',[0.15 0.15 0.60 0.80]);
%         colorDepth = 1000;
%         colormap(flipud(jet(colorDepth)));
%         %             pcolor(X2(1001:end), 1:1:size(CSD_matrix1,1), nrm_CSD_matrix);
%         imagesc(X2(1001:end), [], nrm_CSD_matrix);
%         shading interp; % do not interpolate pixels
%         axis on; % display axis
%         axis tight;% no white borders
%         set(gca, ...
%             'Box'         , 'off'      , ...
%             'TickDir'     , 'in'      , ...
%             'Ydir'        , 'reverse', ...
%             'TickLength'  , [.01 .01] , ...
%             'XMinorTick'  , 'off'      , ...
%             'YMinorTick'  , 'off'     , ...
%             'XGrid'       , 'off'     , ...
%             'YGrid'       , 'off'     , ...
%             'XColor'      , [.0 .0 .0], ...
%             'YColor'      , [.0 .0 .0], ...
%             'LineWidth'   , 0.6        );
%         set(gca,'Xlim',[-25 225]);
%         set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
%         set(gca, 'yticklabel',1:1:size(CSD_matrix3,1), 'Ytick', 1:1:size(CSD_matrix3,1));
%         xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
%         yLabelText = 'Electrode number';
%         hXLabel = xlabel(xLabelText);
%         hYLabel = ylabel(yLabelText);
%         fig_title=sprintf('%s','Shank 1 - Normalized');
%         yaxis=ylim;
%         prestim_offset_y            = yaxis(1):1:yaxis(2);
%         prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
%         hold on;plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',1);
%         caxis([min(min(CSD_matrix3)) max(max(CSD_matrix3))]);
%         zLabelText = 'nA / mm^3';  % greek letters in LaTeX Syntax
%         hcb = colorbar('eastoutside');
%         h_bar = findobj(gcf,'Tag','Colorbar');
%         initpos = get(h_bar,'Position');
%         set(h_bar, ...
%             'Position',[initpos(1)+initpos(3)*2.5 initpos(2)+initpos(4)*0.3 ...
%             initpos(3)*0.4 initpos(4)*0.4]);
%         hcLabel = ylabel(hcb,zLabelText);
%         set(hcb,'YTickLabel',{'Sink','Source'}, 'Ytick', [min(min(CSD_matrix3)) max(max(CSD_matrix3))])
%         set(hcb, ...
%             'Box'         , 'on'     , ...
%             'TickDir'     , 'in'     , ...
%             'TickLength'  , [.010 .010] , ...
%             'LineWidth'   , 0.6);
%         set([gca, hcb, hXLabel, hYLabel, hcLabel], ...
%             'FontSize'   , FontSize    , ...
%             'FontName'   , FontName);
%         ylabh=get(hcb,'Ylabel');
%         set(ylabh,'Position',get(ylabh,'Position')-[8 0 0]);
%         set(gca,'Layer', 'top');
%         drawnow
%     end
    end
            export_fig (fig_title, '-png');
    
end

%% Determine the center of mass
t1=find(X2>=-200 & X2<=500,1,'first');
t2=find(X2>=-200 & X2<=500,1,'last');
[mean_x1 mean_y1]=centroid_me(CSD_matrix1(:,t1:t2));
mean_x1=X2((t1-1)+mean_x1);
fprintf ('%s%d%s%d','The centroid value of shank1 is: ',round(mean_x1), ' ms at contact number ',mean_y1)
fprintf('\n')

if combine_columns==0
t1=find(X2>=-200 & X2<=500,1,'first');
t2=find(X2>=-200 & X2<=500,1,'last');
[mean_x2 mean_y2]=centroid_me(CSD_matrix2(:,t1:t2));
mean_x2=X2((t1-1)+mean_x2);
fprintf ('%s%d%s%d','The centroid value of shank2 is: ',round(mean_x2), ' ms at contact number ',mean_y2)
fprintf('\n')

t1=find(X2>=-200 & X2<=500,1,'first');
t2=find(X2>=-200 & X2<=500,1,'last');
[mean_x3 mean_y3]=centroid_me(CSD_matrix3(:,t1:t2));
mean_x3=X2((t1-1)+mean_x3);
fprintf ('%s%d%s%d','The centroid value of shank3 is: ',round(mean_x3), ' ms at contact number ',mean_y3)
fprintf('\n')
end

%% identify layers
disp ('Look at the CSD and find the layers which contain the sink(red)')
st_sink1 = input('Left column begining of the sink: ','s');
end_sink1 = input('Left column Ending of the sink: ','s');
st_5B1 = input('Left column begining of L5B sink: ','s');
end_5B1 = input('Left column Ending of L5B sink: ','s');
end_61 = input('Left column Ending of L6: ', 's');
gran1=str2double(st_sink1):str2double(end_sink1);
L5B1 = str2double(st_5B1):str2double(end_5B1);
L6end1 = str2double(end_61);

if combine_columns==0
st_sink2 = input('Mid column begining of the sink: ','s');
end_sink2 = input('Mid column Ending of the sink: ','s');
st_5B2 = input('Mid column begining of L5B sink: ','s');
end_5B2 = input('Mid column Ending of L5B sink: ','s');
end_62 = input('Mid column Ending of L6: ', 's');
gran2=str2double(st_sink2):str2double(end_sink2);
L5B2 = str2double(st_5B2):str2double(end_5B2);
L6end2 = str2double(end_62);

st_sink3 = input('Right column begining of the sink: ','s');
end_sink3 = input('Right column Ending of the sink: ','s');
st_5B3 = input('Right column begining of L5B sink: ','s');
end_5B3 = input('Right column Ending of L5B sink: ','s');
end_63 = input('Right column Ending of L6: ', 's');
gran3=str2double(st_sink3):str2double(end_sink3);
L5B3 = str2double(st_5B3):str2double(end_5B3);
L6end3 = str2double(end_63);

% disp ('Look at the CSD and find the layers which contain the sink(red)')
% st_sink3 = input('Shank 3 begining of the sink: ','s');
% end_sink3 = input('Shank 3 Ending of the sink: ','s');
% gran3=str2double(st_sink3):str2double(end_sink3);

probemap_gran = zeros(size(probemap));
probemap_gran(1:gran1(1)-1,1) = 2.5;
probemap_gran(gran1,1) = 4;
probemap_gran(gran1(end)+1:L5B1(1)-1,1) = 5;
probemap_gran(L5B1,1) = 5.5;
probemap_gran(L5B1(end)+1:L6end1,1) = 6;
probemap_gran(1:gran2(1)-1,2) = 2.5;
probemap_gran(gran2,2) = 4;
probemap_gran(gran2(end)+1:L5B2(1)-1,2) = 5;
probemap_gran(L5B2,2) = 5.5;
probemap_gran(L5B2(end)+1:L6end2,2) = 6;
probemap_gran(1:gran3(1)-1,3) = 2.5;
probemap_gran(gran3,3) = 4;
probemap_gran(gran3(end)+1:L5B2(1)-1,3) = 5;
probemap_gran(L5B3,3) = 5.5;
probemap_gran(L5B3(end)+1:L6end3,3) = 6;
% probemap_gran(1:gran3(1)-1,3) = -1;
% probemap_gran(gran3(end)+1:21,3) = 1;
else
    granidx = zeros(43,1);
    gran=gran1;
    granidx(gran)=1;
    gran1=find(granidx(pos_col==1));
    gran2=find(granidx(pos_col==2));
    gran3=gran1;
    probemap_gran = nan(size(probemap));
    probemap_gran(22,:) = nan;
    probemap_gran(1:gran1(1)-1,1) = -1;
    probemap_gran(gran1(end)+1:21,1) = 1;
    probemap_gran(1:gran2(1)-1,2) = -1;
    probemap_gran(gran2(end)+1:22,2) = 1;
    probemap_gran(1:gran3(1)-1,3) = -1;
    probemap_gran(gran3(end)+1:21,3) = 1;
end

% %% plot layers
% %%% Now this and final section creates the average CSD across contacts
% %%% given
% %%% the user defined top and bottom of the sink(red) region in the CSD
% %%% plots. The mean and the STD are calaucated for each layer. and the
% %%% envelope of the std is ploted using the jbfill fx
% axis_time=t1:t2;
% if gran1(1)==2;
%     SG_c=(CSD_matrix1(1,axis_time));
%     SG_std=std(CSD_matrix1(1,axis_time))./1;
%     SG_p=SG_c+SG_std;  SG_m=SG_c-SG_std;
% else
%     SG_c=mean(CSD_matrix1(1:(gran1(1)-1),axis_time),1);
%     SG_std=(std(CSD_matrix1((1:gran1(1)-1),axis_time),0,1))./size(CSD_matrix1(1:gran1(1)-1),2);
%     SG_p=SG_c+SG_std;  SG_m=SG_c-SG_std;
% end
% G_c=mean(CSD_matrix1(gran1,axis_time),1);
% G_std=(std(CSD_matrix1(gran1,axis_time),0,1))./size(CSD_matrix1(gran1),2);
% G_p=G_c+G_std;     G_m=G_c-G_std;
% 
% IG_c=mean(CSD_matrix1((gran1(length(gran1))+1):size(CSD_matrix1,1),axis_time),1);
% IG_std=(std(CSD_matrix1((gran1(length(gran1))+1):size(CSD_matrix1,1),axis_time),0,1))...
%     ./(length((gran1(length(gran1))+1):size(CSD_matrix1)));
% IG_p=IG_c+IG_std;  IG_m=IG_c-IG_std;
% 
% figure_width = 12;
% figure_height = 12;
% FontSize = 12;
% FontName = 'MyriadPro-Regular'; % or choose any other font
% % --- setup plot windows
% figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
% Avg_CSD1 = figure('Name','CSD_layers1');
% set(Avg_CSD1,'Visible', figuresVisible)
% set(Avg_CSD1, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
% set(Avg_CSD1, 'PaperPositionMode', 'auto');
% set(Avg_CSD1, 'Color', [1 1 1]); % Sets figure background
% set(Avg_CSD1, 'Color', [1 1 1]); % Sets axes background
% % --- dimensions and position of plot
% hsp = subplot(1,1,1, 'Parent', Avg_CSD1);
% set(hsp,'Position',[0.10 0.15 0.80 0.80]);
% x=X2(axis_time)';
% hold on;handle_vector(:,1) = plot(x,SG_c,'r','LineWidth',2);
% hold on; handle_vector(:,2) = jbfill(x,SG_p,SG_m,'r','r',1,.4);
% hold on; handle_vector(:,3) = plot(x,G_c,'b','LineWidth',2);
% hold on; handle_vector(:,4) = jbfill(x,G_p,G_m,'b','b',1,.4);
% hold on; handle_vector(:,5) = plot(x,IG_c,'color',[0 .5 0],'LineWidth',2);
% hold on; handle_vector(:,6) = jbfill(x,IG_p,IG_m,[0 .5 0],[0 .5 0],1,.4);
% % will remove the 3rd legend entry.
% hasbehavior(handle_vector(2),'legend',false);
% hasbehavior(handle_vector(4),'legend',false);
% hasbehavior(handle_vector(6),'legend',false);
% % will remove the 3rd legend entry.
% set(gca, ...
%     'Box'         , 'off'      , ...
%     'TickDir'     , 'out'      , ...
%     'TickLength'  , [0 0] , ...
%     'XMinorTick'  , 'off'      , ...
%     'YMinorTick'  , 'off'     , ...
%     'XGrid'       , 'off'     , ...
%     'YGrid'       , 'off'     , ...
%     'XColor'      , [.0 .0 .0], ...
%     'YColor'      , [.0 .0 .0], ...
%     'LineWidth'   , 0.6        );
% set(gca,'Xlim',[-25 225]);
% axis off
% %export_fig ('Avg1', '-png','-r600','-opengl')
% axis on
% set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
% set(gca, 'ytick', [],'tickdir','out');
% set(gca, 'YTickLabel', num2str(get(gca, 'YTick')'))
% yaxis=ylim;
% prestim_offset_y            = yaxis(1):1:yaxis(2);
% prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
% hold on;plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',1);
% xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
% yLabelText = 'nA / mm^3';
% h=legend({'Supragranular','Granular','Infragranular'});
% set(h, 'Box', 'off','location', 'Best')
% hXLabel = xlabel(xLabelText);
% hYLabel = ylabel(yLabelText);
% fig_title=sprintf('%s','Avg1 ');
% set([gca, hXLabel, hYLabel, h], ...
%     'FontSize'   , FontSize    , ...
%     'FontName'   , FontName);
% set(h, ...
%     'FontSize'   , 10    , ...
%     'FontName'   , FontName);
% set(gca,'Layer', 'top');
% drawnow
% %     export_fig (fig_title, '-png','-r600','-opengl')
% 
% if combine_columns==0
% if gran2(1)==2;
%     SG_c=(CSD_matrix2(1,axis_time));
%     SG_std=std(CSD_matrix2(1,axis_time))./1;
%     SG_p=SG_c+SG_std;  SG_m=SG_c-SG_std;
% else
%     SG_c=mean(CSD_matrix2(1:(gran2(1)-1),axis_time),1);
%     SG_std=(std(CSD_matrix2((1:gran2(1)-1),axis_time),0,1))./sqrt(size(CSD_matrix2(1:gran2(1)-1),2));       % MK added square roots!
%     SG_p=SG_c+SG_std;  SG_m=SG_c-SG_std;
% end
% G_c=mean(CSD_matrix2(gran2,axis_time),1);
% G_std=(std(CSD_matrix2(gran2,axis_time),0,1))./sqrt(size(CSD_matrix2(gran2),2));
% G_p=G_c+G_std;     G_m=G_c-G_std;
% 
% IG_c=mean(CSD_matrix2((gran2(length(gran2))+1):size(CSD_matrix2,1),axis_time),1);
% IG_std=(std(CSD_matrix2((gran2(length(gran2))+1):size(CSD_matrix2,1),axis_time),0,1))...
%     ./sqrt(length((gran1(length(gran1))+1):size(CSD_matrix2)));
% IG_p=IG_c+IG_std;  IG_m=IG_c-IG_std;
% 
% figure_width = 14;
% figure_height = 12;
% FontSize = 12;
% FontName = 'MyriadPro-Regular'; % or choose any other font
% % --- setup plot windows
% figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
% Avg_CSD2 = figure('Name','CSD_layers2');
% set(Avg_CSD2,'Visible', figuresVisible)
% set(Avg_CSD2, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
% set(Avg_CSD2, 'PaperPositionMode', 'auto');
% set(Avg_CSD2, 'Color', [1 1 1]); % Sets figure background
% set(Avg_CSD2, 'Color', [1 1 1]); % Sets axes background
% % --- dimensions and position of plot
% hsp = subplot(1,1,1, 'Parent', Avg_CSD2);
% set(hsp,'Position',[0.10 0.15 0.80 0.80]);
% x=X2(axis_time)';
% hold on;handle_vector(:,1) = plot(x,SG_c,'r','LineWidth',2);
% hold on; handle_vector(:,2) = jbfill(x,SG_p,SG_m,'r','r',1,.4);
% hold on; handle_vector(:,3) = plot(x,G_c,'b','LineWidth',2);
% hold on; handle_vector(:,4) = jbfill(x,G_p,G_m,'b','b',1,.4);
% hold on; handle_vector(:,5) = plot(x,IG_c,'color',[0 .5 0],'LineWidth',2);
% hold on; handle_vector(:,6) = jbfill(x,IG_p,IG_m,[0 .5 0],[0 .5 0],1,.4);
% % will remove the 3rd legend entry.
% hasbehavior(handle_vector(2),'legend',false);
% hasbehavior(handle_vector(4),'legend',false);
% hasbehavior(handle_vector(6),'legend',false);
% % will remove the 3rd legend entry.
% set(gca, ...
%     'Box'         , 'off'      , ...
%     'TickDir'     , 'out'      , ...
%     'TickLength'  , [0 0] , ...
%     'XMinorTick'  , 'off'      , ...
%     'YMinorTick'  , 'off'     , ...
%     'XGrid'       , 'off'     , ...
%     'YGrid'       , 'off'     , ...
%     'XColor'      , [.0 .0 .0], ...
%     'YColor'      , [.0 .0 .0], ...
%     'LineWidth'   , 0.6        );
% set(gca,'Xlim',[-25 225]);
% axis off
% %export_fig ('Avg2', '-png','-r600','-opengl')
% axis on
% set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
% set(gca, 'ytick', [],'tickdir','out');
% set(gca, 'YTickLabel', num2str(get(gca, 'YTick')'))
% yaxis=ylim;
% prestim_offset_y            = yaxis(1):1:yaxis(2);
% prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
% hold on; plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',1);
% xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
% yLabelText = 'nA / mm^3';
% h=legend({'Supragranular','Granular','Infragranular'});
% set(h, 'Box', 'off','location', 'Best')
% hXLabel = xlabel(xLabelText);
% hYLabel = ylabel(yLabelText);
% fig_title=sprintf('%s','Avg2 ');
% set([gca, hXLabel, hYLabel, h], ...
%     'FontSize'   , FontSize    , ...
%     'FontName'   , FontName);
% set(h, ...
%     'FontSize'   , 10    , ...
%     'FontName'   , FontName);
% set(gca,'Layer', 'top');
% drawnow
% 
% 
% if gran3(1)==2;
%     SG_c=(CSD_matrix3(1,axis_time));
%     SG_std=std(CSD_matrix3(1,axis_time))./1;
%     SG_p=SG_c+SG_std;  SG_m=SG_c-SG_std;
% else
%     SG_c=mean(CSD_matrix3(1:(gran3(1)-1),axis_time),1);
%     SG_std=(std(CSD_matrix3((1:gran3(1)-1),axis_time),0,1))./sqrt(size(CSD_matrix3(1:gran3(1)-1),2));       % MK added square roots!
%     SG_p=SG_c+SG_std;  SG_m=SG_c-SG_std;
% end
% G_c=mean(CSD_matrix3(gran3,axis_time),1);
% G_std=(std(CSD_matrix3(gran3,axis_time),0,1))./sqrt(size(CSD_matrix3(gran3),2));
% G_p=G_c+G_std;     G_m=G_c-G_std;
% 
% IG_c=mean(CSD_matrix3((gran3(length(gran3))+1):size(CSD_matrix3,1),axis_time),1);
% IG_std=(std(CSD_matrix3((gran3(length(gran3))+1):size(CSD_matrix3,1),axis_time),0,1))...
%     ./sqrt(length((gran1(length(gran1))+1):size(CSD_matrix3)));
% IG_p=IG_c+IG_std;  IG_m=IG_c-IG_std;
% 
% figure_width = 14;
% figure_height = 12;
% FontSize = 12;
% FontName = 'MyriadPro-Regular'; % or choose any other font
% % --- setup plot windows
% figuresVisible = 'on'; % 'off' for non displayed plots (will still be exported)
% Avg_CSD3 = figure('Name','CSD_layers3');
% set(Avg_CSD3,'Visible', figuresVisible)
% set(Avg_CSD3, 'units', 'centimeters', 'pos', [5 5 figure_width figure_height])
% set(Avg_CSD3, 'PaperPositionMode', 'auto');
% set(Avg_CSD3, 'Color', [1 1 1]); % Sets figure background
% set(Avg_CSD3, 'Color', [1 1 1]); % Sets axes background
% % --- dimensions and position of plot
% hsp = subplot(1,1,1, 'Parent', Avg_CSD3);
% set(hsp,'Position',[0.10 0.15 0.80 0.80]);
% x=X2(axis_time)';
% hold on;handle_vector(:,1) = plot(x,SG_c,'r','LineWidth',2);
% hold on; handle_vector(:,2) = jbfill(x,SG_p,SG_m,'r','r',1,.4);
% hold on; handle_vector(:,3) = plot(x,G_c,'b','LineWidth',2);
% hold on; handle_vector(:,4) = jbfill(x,G_p,G_m,'b','b',1,.4);
% hold on; handle_vector(:,5) = plot(x,IG_c,'color',[0 .5 0],'LineWidth',2);
% hold on; handle_vector(:,6) = jbfill(x,IG_p,IG_m,[0 .5 0],[0 .5 0],1,.4);
% % will remove the 3rd legend entry.
% hasbehavior(handle_vector(2),'legend',false);
% hasbehavior(handle_vector(4),'legend',false);
% hasbehavior(handle_vector(6),'legend',false);
% % will remove the 3rd legend entry.
% set(gca, ...
%     'Box'         , 'off'      , ...
%     'TickDir'     , 'out'      , ...
%     'TickLength'  , [0 0] , ...
%     'XMinorTick'  , 'off'      , ...
%     'YMinorTick'  , 'off'     , ...
%     'XGrid'       , 'off'     , ...
%     'YGrid'       , 'off'     , ...
%     'XColor'      , [.0 .0 .0], ...
%     'YColor'      , [.0 .0 .0], ...
%     'LineWidth'   , 0.6        );
% set(gca,'Xlim',[-25 225]);
% axis off
% %export_fig ('Avg2', '-png','-r600','-opengl')
% axis on
% set(gca,'XTickLabel',[0 200], 'Xtick', [0 200])
% set(gca, 'ytick', [],'tickdir','out');
% set(gca, 'YTickLabel', num2str(get(gca, 'YTick')'))
% yaxis=ylim;
% prestim_offset_y            = yaxis(1):1:yaxis(2);
% prestim_offset_t            = ones(1, length(prestim_offset_y))*0;
% hold on; plot(prestim_offset_t, prestim_offset_y, 'k', 'linewidth',1);
% xLabelText = 'Time from stimulus onset (ms)';  % greek letters in LaTeX Syntax
% yLabelText = 'nA / mm^3';
% h=legend({'Supragranular','Granular','Infragranular'});
% set(h, 'Box', 'off','location', 'Best')
% hXLabel = xlabel(xLabelText);
% hYLabel = ylabel(yLabelText);
% fig_title=sprintf('%s','Avg3 ');
% set([gca, hXLabel, hYLabel, h], ...
%     'FontSize'   , FontSize    , ...
%     'FontName'   , FontName);
% set(h, ...
%     'FontSize'   , 10    , ...
%     'FontName'   , FontName);
% set(gca,'Layer', 'top');
% drawnow
% end
% 
% savefile=sprintf('%s%s','FinalCSD','.mat');
% save (savefile);

%% Save layers!!!

layers = nan(1,sum(~isnan(probemap(:))));
[~,inds] = sort(zvals);
inds = flipud(inds);
layers_bycol = probemap_gran(~isnan(probemap(:)));
layers = layers_bycol(inds);
discreps = find(diff(layers_bycol(inds))<0) +1; % discrepancies in layers
layers(discreps) = layers(discreps-1);

cd ..
    expdir = cd;
    dirname= ...
        uigetdir(expdir,'In which experiment folders do you want to save the layer information?');
        save(sprintf('%s\\layers.mat',dirname),'layers')
        save(sprintf('%s\\layer_info.mat',dirname),'probemap','probemap_chrem','probemap_gran')
    moredirs = input('Do you want to save the layers in more places? Y or N: ','s');
    if strcmp(moredirs,'Y')
        howmany = input('How many more directories do you want to save the layer info in? ','s');
        for n = 1:str2double(howmany)
            dirname = uigetdir(expdir,'In which experiment folders do you want to save the layer information?');
            save(sprintf('%s\\layers.mat',dirname),'layers')
            save(sprintf('%s\\layer_info.mat',dirname),'probemap','probemap_chrem','probemap_gran')
        end
    end