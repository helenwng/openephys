function [waveforms,max_ch,shank] = readWaveformsFromRez_K2(unit,exp_path,rez)
% unit = cluster identifier (e.g. 47)
% exp_path = experiment directory
% rez = rez.mat output from kilosort
    % load(sprintf('%s\\rez.mat',exp_path))

clusters = readNPY(fullfile(exp_path,'spike_clusters.npy'));
templates = readNPY(fullfile(exp_path,'spike_templates.npy'));
chans = readNPY(fullfile(exp_path,'channel_positions.npy'));
templates_w = readNPY(fullfile(exp_path,'templates.npy'));    % whitened templates from kilosort2 - MAK 6/13/19
% templates_unw = readNPY(fullfile(exp_path,'templates_unw.npy'));    % whitened templates from kilosort2 - MAK 6/13/19
% templates_ind = readNPY(fullfile(exp_path,'templates_ind.npy'));    % whitened templates from kilosort2 - MAK 6/13/19
% load(fullfile(exp_path,'rez2.mat'))

n = hist(templates(clusters==unit),0:double(max(templates)));
[~,m] = max(n);

% % TEMPPPP NEED TO FIX
% if m>size(rez.Wraw,3)
%     m=size(rez.Wraw,3);
% end

% waveforms(:,:) = rez.Wraw(:,:,m)';
waveforms(:,:) = squeeze(templates_w(m,:,:));   %MAK 6/13/19
waveforms = waveforms(22:61,:).*.195;     % NEED TO FIGURE OUT UNITS
nchs = size(chans,1);
if range(chans(:,1))<50
    shankchs = nchs;
else
    shankchs = nchs/ceil(range(chans(:,1))/max(diff(unique(chans(:,1)))));     % chs per shank 
end

% figure out max_ch   
[~,max_chN] = max(range(waveforms));
shank = ceil(max_chN/shankchs)-1;       
[~,pos] = sort(chans(shank*shankchs+1:(shank+1)*shankchs,2));
max_ch = find(flipud(pos)==(max_chN-(shank*shankchs)));

% if shankchs == nchs 
%     shank = [];
% end

end