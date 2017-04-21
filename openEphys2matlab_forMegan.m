function openEphys2matlab(exp_path)
%extracts info from open Ephys files
%hard coded for 32 channel for now - to do: write code to detect all
%continuous files in exp_path, separate aux and electrodes and load all
%channels

% files = dir(exp_path);
% nfiles = length(files)-2;
% eventfile = {};
% allFiles = cell(1,nfiles);
% for i = 1:nfiles
%     allFiles{i} = files(i+2).name;
%     if strcomp(allFiles{i},'all_channels.events')
%         eventfile = fullfile(exp_path,'all_channels.events');
%     elseif i==nfiles & isempty(eventfile)
%         error('all_channels.events is missing')
%     end
% end
% 
% 
% nADC = 0;

    

%todo - initialize channels w/ size info

%load electrode channels
contfile = fullfile(exp_path,'100_CH1.continuous');
[~, dataTime, dataInfo] = load_open_ephys_data(contfile);
% dataTime = dataTime./dataInfo(1).header.sampleRate;       % uncomment if using load_open_ephys_data_faster
nsamples = length(dataTime);

%load ADC inputs
for i = 1:2
    [ADCin(i,:),~,ADCinfo(i,:)] = load_open_ephys_data(fullfile(exp_path,sprintf('100_ADC%d.continuous',i)));
end

%load all_channels.events
eventfile = fullfile(exp_path,'all_channels.events')
[events,eventTime,info] = load_open_ephys_data_faster(eventfile);
 amp_sr = info.header.sampleRate;
 eventIdx = floor((eventTime-dataTime(1))*amp_sr+1);

%define digital channels (subtract 1 from label on i/o boards)
epocCH = 0;
encdACH = 1;
encdBCH = 2;

%define analog channels
photo = ADCin(1,:);
LED = ADCin(2,:);

%digital events
epoc = zeros(size(dataTime));
encdA = zeros(size(dataTime));
encdB = zeros(size(dataTime));

%find epoc on and off times and fill out binary vector same size as data
epocOn = eventIdx(events==epocCH&info.eventId&info.eventType==3);
epocOff = eventIdx(events==epocCH&~info.eventId&info.eventType==3);

%in case of corrupted files:
epocOn = epocOn(ismember(epocOn,1:length(dataTime)));
epocOff = epocOff(ismember(epocOff,1:length(dataTime)));

% downsample to 1000Hz to save memory
LN              = length(dataTime);
div             = dataInfo(1).header.sampleRate/1000;
zx              = 1:div:LN;
izx             = floor(zx);
time      = dataTime(izx);    % downsample from 20000 hz to 1000 hz
% epocOn = epocOn(izx);
% epocOff = epocOff(izx);

if epocOn(1)<epocOff(1)
    if size(epocOn,1 )> size(epocOff,1)   % if experiment got cut off in middle of trial...
        epocOn = epocOn(1:size(epocOff,1));         % ...drop the last trial
    end
    if size(epocOn) == size(epocOff)
        field_trials = [epocOn epocOff];
        trials = [dataTime(epocOn) dataTime(epocOff)];    % find and ismember in case of corrupted files
    end
    for i = 1:length(epocOn)
        if i<=length(epocOff)
            epoc(epocOn(i):epocOff(i)) = ones(size(epocOn(i):epocOff(i)));
        else
            epoc(epocOn(i):end) = ones(1,nsamples-epocOn(i));
        end
    end
end

%find encdA on and off times and fill out binary vector same size as data
encdAOn = eventIdx(events==encdACH&info.eventId&info.eventType==3);
encdAOff = eventIdx(events==encdACH&~info.eventId&info.eventType==3);

if encdAOn(1)>encdAOff(1)
    encdA(1:encdAOff(1)) = ones(encdAOff(1),1);
    for i = 1:length(encdAOn)
        if i+1<=length(encdAOff)
            encdA(encdAOn(i):encdAOff(i+1)) = ones(size(encdAOn(i):encdAOff(i+1)));
        else
            encdA(encdAOn(i):end) = ones(nsamples-encdAOn(i)+1,1);
        end
    end
else
    for i = 1:length(encdAOn)
        if i<=length(encdAOff)
            encdA(encdAOn(i):encdAOff(i)) = ones(size(encdAOn(i):encdAOff(i)));
        else
            encdA(encdAOn(i):end) = ones(nsamples-encdAOn(i)+1,1);
        end
    end
end


%find encdB on and off times and fill out binary vector same size as data
encdBOn = eventIdx(events==encdBCH&info.eventId&info.eventType==3);
encdBOff = eventIdx(events==encdBCH&~info.eventId&info.eventType==3);

if encdBOn(1)>encdBOff(1)
    encdB(1:encdBOff(1)) = ones(encdBOff(1),1);
    for i = 1:length(encdBOn)
        if i+1<=length(encdBOff)
            encdB(encdBOn(i):encdBOff(i+1)) = ones(size(encdBOn(i):encdBOff(i+1)));
        else
            encdB(encdBOn(i):end) = ones(nsamples-encdBOn(i)+1,1);
        end
    end
else
    for i = 1:length(encdBOn)
        if i<=length(encdBOff)
            encdB(encdBOn(i):encdBOff(i)) = ones(size(encdBOn(i):encdBOff(i)));
        else
            encdB(encdBOn(i):end) = ones(nsamples-encdBOn(i)+1,1);
        end
    end
end

%getSyncTimes for "re"
[re,x]= getSyncTimesRevCorr_AC(photo',(1/amp_sr));

%todo - write code for when recording starts during a trial ->
%epocOff(1)<epocOn(1)

%save variables in data.mat
cd(exp_path)
save('data.mat', 'trials','field_trials','amp_sr','photo','LED','epoc','encdA','encdB','re')

end