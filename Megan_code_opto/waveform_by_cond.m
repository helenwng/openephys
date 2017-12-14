function waveform_by_cond(unit,chansInDat,exp_system,exp_path)
cd(exp_path)
% unit = cluster identifier (from kilosort)
spike_times = readNPY('spike_times.npy');
clusters = readNPY('spike_clusters.npy');
load('data.mat','field_trials','LED','amp_sr')
load('light_params.mat','all_light')
if strcmpi(exp_system,'openephys')
    clusters = readNPY(sprintf('%s\\spike_clusters.npy',exp_path));
    templates = readNPY(sprintf('%s\\spike_templates.npy',exp_path));
    load('rez.mat');
end
[spike_times,clusters] = clean_light_artifacts(field_trials,LED,spike_times,clusters,amp_sr);      % NEW MAK addition 2/26/17 - get rid of light artifacts mistaken for spikes - modified and confirmed works 8/9/17
if strcmpi(exp_system,'openephys')
    load('dropped_spikes.mat')
    templates(dropped_spikes')=[];
end
sampleTimes = spike_times(clusters==unit);
window = [-12 16];
nToRead = [];

s = dir(exp_path); 
if strcmpi(exp_system,'intan')
    key = 'amplifier';
else
    key = '.dat';
end
for i=1:length(s)
    if strfind(s(i).name,key) 
        datFilename = s(i).name;
    end
end
datFilename = (sprintf('%s\\%s',exp_path,datFilename));

FileInf = dir(datFilename);
nSampsInDat = (FileInf.bytes/chansInDat/2);
rawData = memmapfile(datFilename, 'Format', {'int16', [chansInDat, nSampsInDat], 'x'});

if isempty(nToRead)
    theseTimes = sampleTimes;
    nToRead = length(sampleTimes);
else
    if nToRead>length(sampleTimes)
        nToRead = length(sampleTimes);
    end
    q = randperm(nToRead);
    theseTimes = sampleTimes(sort(q(1:nToRead)));
end

theseTimes = theseTimes(theseTimes>-window(1) & theseTimes<nSampsInDat-window(2)-1);
nToRead = numel(theseTimes);

allWF = zeros(chansInDat, diff(window)+1, nToRead, 'int16');

for i=1:nToRead
    allWF(:,:,i) = (double(rawData.Data.x(1:chansInDat,theseTimes(i)+window(1):theseTimes(i)+window(2))));  % scale appropriately for Intan (get units in microvolts)
end

if strcmpi(exp_system,'openephys')
    n = hist(templates(clusters==unit),0:double(max(templates)));
    [~,m] = max(n);
    W(:,:) = allWF(rez.ypos(m),:,:);
else
    [~,m] = max(range(mean(allWF,3),2));
    W(:,:) = allWF(m,:,:);
end
W = W';
W = double(W)-repmat(mean(W(:,1:6),2),1,size(W,2));    % normalize so that each waveform starts from 0
meanW = mean(W,1);

samplesbyLight = nan(size(sampleTimes));
for i = 1:size(field_trials,1)
    samplesbyLight(ismember(sampleTimes,[field_trials(i,1):field_trials(i,2)])) = all_light(i);
end

lightconds = unique(all_light);
    color_matA = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4; 0 .7 .2]; % for graphing purposes (first is black, last is green)
    color_matB = [.7 .7 .7; .6 .8 1; .8 .8 1; .4 0.5 .4; .7 .7 .2]; % for graphing purposes (first is black, last is green)
fig=figure;
hold on
for ii = 1:length(lightconds)
    plot(W(samplesbyLight==lightconds(ii),:,:)','Color',color_matB(ii,:))
end
for ii = 1:length(lightconds)
    plot(mean(W(samplesbyLight==lightconds(ii),:,:),1)','Color',color_matA(ii,:),'linewidth',8)
end
xlim([1 size(W,2)])
print2eps(strcat('waveforms_',num2str(unit)),fig)
end

