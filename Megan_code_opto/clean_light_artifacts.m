 
function [clean_spike_times,clean_clusters] = clean_light_artifacts(field_trials,LED,spike_times,clusters,amp_sr)
% drop spikes that occur within 1ms of light start, in case it's
% artifactual
% edited 8/2/17 to make faster (MAK)

% load('data.mat','field_trials','LED') 
% spike_times = readNPY('spike_times.npy');           % here, spike times from ALL clusters
% clusters = readNPY('spike_clusters.npy');           
khz = amp_sr/1000;

tic
dropped_spikes = [];
samps_per_t = max(diff(field_trials,[],2))*khz;                     % do this in case each trial doens't have exact same number of samples
% start_samp = zeros(1,size(field_trials,1));
% light_start = start_samp;
for t = 1:size(field_trials,1)      % for each trial
    start_samp = field_trials(t,1)*khz;  % multiply by samples per ms cuz using ORIGINAL sampling rate (e.g. 20kHz)
    light_out = LED(start_samp:start_samp+samps_per_t);      % get actual values of LED input to Intan
    thresh = mean(light_out(1:1000))+10*std(light_out(1:1000));      % threshold for determining light on or off
    light_on = find(light_out>thresh);
    light_diff = diff(light_on);
    light_start = [find(light_out>thresh,1,'first') light_on(find(light_diff>1)+1)];     % in case of trains experiment
    spks2drop = cell(1,length(light_start));
    for ii = 1:length(light_start)
        spks2drop{ii} = find(ismember(spike_times,start_samp+light_start(ii)-1:start_samp+light_start(ii)+19))';    % drop if its within 20 samples (1ms)
    end
    badspikes = [spks2drop{cellfun(@(x) ~isempty(x),spks2drop)}];
    if ~isempty(badspikes)
        dropped_spikes = [dropped_spikes badspikes];   
        disp(strcat('Dropping spike # ',num2str(badspikes)));
%     dropped_spikes = [dropped_spikes spike_times(ismember(spike_times,start_samp+light_start-1:start_samp+light_start+19))'];   % drop if its within 20 samples (1ms)
    end
end
toc

clusters(dropped_spikes)=[];
spike_times(dropped_spikes)=[];    
clean_spike_times = spike_times;
clean_clusters = clusters;

save('dropped_spikes.mat', 'dropped_spikes')
end
