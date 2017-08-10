function [psth,psth_FR] = make_psth_v2(binsize,edges,which_trials,spike_raster,all_light_trials)

% binsize = size of bin in seconds
% edges = vector 0:binsize:totaltrialtime
% which_trials = 1xnum_trials vector of 1s and 0s, 1s indicating trials to look
% at
% spike_raster = num_trials x time (totaltime*1000, so in 1000Hz samp rate)
    % matrix of 1s and 0s indicating when spikes occurred
% all_light_trials = 1xnum_trials vector indicating which light condition
    % each trial belongs to
% output: psth = numconds x length(edges)-1 matrix of spike counts (psth) or firing rate(psth_FR, in spikes/s)

conds = unique(all_light_trials);
num_conds = length(conds);
psth = zeros(num_conds,length(edges));

light_trials = all_light_trials(which_trials);
rast = spike_raster(which_trials,:);

for i = 1:num_conds
    trials = rast(light_trials==conds(i),:);
    ntrials(i) = size(trials,1);
    timevec = repmat((edges(1):size(spike_raster,2)-1)/1000,ntrials(i),1);           % in ms
    spikes = timevec.*trials;
    spikes = spikes(:);
    spikes = spikes(spikes~=0);
    psth(i,:) = histc(spikes,edges);
end

psth = psth(:,1:end-1);     % because last column from histc is values that equal edges(end)
psth_FR = (psth./repmat(ntrials',1,size(psth,2)))./binsize;


end

