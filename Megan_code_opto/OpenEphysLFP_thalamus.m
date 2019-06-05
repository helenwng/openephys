function openEphysLFP_thalamus(exp_path,exp_type,probe)
% modified from IntanLFP_fastfilt (MAK 2/4/19)

% get necessary data
cd(exp_path)
if exist(sprintf('%s/data.mat',exp_path),'file')
    load(sprintf('%s/data.mat',exp_path))      % data from intanphy2matlab.m
else
    openEphys2matlab(exp_path);       % data from openEphys2matlab.m
    load(sprintf('%s/data.mat',exp_path)) 
end

% % Determine which probe was used
% tmp = dir(exp_path);
% possible_probes = {'A2x16_2.prb','A1x32.prb','NeuroNexus_8FCS.prb','A1x32_A48.prb'};
% for i = 1:length(tmp)
%     isprobename = strcmp(tmp(i).name,possible_probes);
%     if sum(isprobename)
%         probetype = find(isprobename);
%     end
% end
% num_channels = 32;              % currently hard-coded for 32 channel probe and 25um spacing
% spacing = 25;
num_channels = str2double(probe(regexp(probe,'\d')));

%% Get data for LFPs

if ~exist(sprintf('%s/LFP_all.mat',exp_path),'file')
    first_half = exist(fullfile(exp_path,'100_CH1.continuous'),'file');
    if first_half
        ch_offset = 0;
        [data, ~, dataInfo] = load_open_ephys_data_faster(sprintf('%s\\100_CH1.continuous',exp_path));
    else
        ch_offset = 64;
        [data, ~, dataInfo] = load_open_ephys_data_faster(sprintf('%s\\100_CH65.continuous',exp_path));
    end
    % downsample
    amp_sr = dataInfo.header.sampleRate;
    div = amp_sr/1000;
    zx = 1:div:length(data);
    izx = floor(zx);
    lfp_data = zeros(num_channels,length(izx));
    for n = 1:num_channels
        [data, ~, ~] = load_open_ephys_data(sprintf('%s\\100_CH%d.continuous',exp_path,n+ch_offset));
        lfp_data(n,:) = data(izx);
    end

% filter LFPs
    fRawo_ln = newfilter(lfp_data,1000,0);
    save('LFP_all.mat','fRawo_ln','-nocompression','-v7.3')
else
    load(sprintf('%s/LFP_all.mat',exp_path))
end

%% split LFPs and average according to shank
p = eval(sprintf('probemap_%s_func',probe));
shks = unique(p.shaft);
mean_lfp = zeros(length(shks),size(fRawo_ln,2));
for nn=1:length(shks)
    shk_inds(nn,:) = find(p.shaft==nn);
    mean_lfp(nn,:) = mean(fRawo_ln(p.channels(shk_inds(nn,:))+1,:),1);
end

if exist(sprintf('%s\\light_params.mat',exp_path),'file')
    load(sprintf('%s\\light_params.mat',exp_path))
else
    [all_light, pulse_dur, ~,av_light_start] = get_lightstim_v2(exp_path,exp_type);
end

% startT = av_light_start(1)*1000;
% pulseT = round(pulse_dur(2),2)*1000;
startT = 1001;
pulseT = 1000;
for l = 1:size(mean_lfp,1)
    for t = 1:size(field_trials,1)
        lfp_all{l}(:,t) = mean_lfp(l,field_trials(t,1)+startT:field_trials(t,1)+startT+pulseT)';
    end

    T = size(lfp_all{l},1);
    df = 1/(T/1000);
    f = [0:df:1000-df];
    t_chunk = 500;
    light_trials = find(all_light);
    nolight_trials = find(all_light == 0);
    for t = 1:length(all_light)
        X(t,:) = fft(lfp_all{l}(:,t));
        P(t,:) = (abs(X(t,:))).^2;
    end


    mean_Plight = mean(P(light_trials,:),1);
    mean_Pnolight = mean(P(nolight_trials,:),1);
    SE_Plight = std(P(light_trials,:),1)/sqrt(length(light_trials));
    SE_Pnolight = std(P(nolight_trials,:),1)/sqrt(length(nolight_trials));
    mean_P = [mean_Pnolight; mean_Plight];
    SE_P = [SE_Pnolight; SE_Plight];

    figure;
    shadedErrorBar(f,mean_Pnolight,SE_Pnolight,'b')
    hold on;
    shadedErrorBar(f,mean_Plight,SE_Plight,'g')
    xlim([0 80])
    legend('no light','light')
    
    figure;plot(f(1:100),log10(mean_P(:,1:100)))
    legend('no light','light')
end

%% using chronux

[prestim,poststim,stimtime,trial_type,IVs] = get_exp_params(exp_path,exp_type);
vis_trials = find(trial_type(:,1));
blank_trials = find(trial_type(:,1)==0);

% make appropriate figure legends
if strcmp(exp_type,'trains')
    for i = 2:length(lightconds)        % for graphing purposes
        legend_labels{i-1} = sprintf('%dHz',lightconds(i));
    end
elseif strcmp(exp_type,'intensities')
    legend_labels = {'Low light','Medium light','High light'};
    full_legend_labels = {'Light OFF','Low light','Medium light','High light'};
else
    legend_labels = {'Light ON - short','Light ON - long'};
end

params.tapers = [3 5];
params.fpass = [0 200];
params.err = [2 .05];
params.Fs = 1000;
params.trialave = 1;

conds = unique(all_light);
for ii = 1:length(conds)
    trial_conds{ii,:} = find(all_light == conds(ii));
end

layer_titles = {'Shank 1' 'Shank 2' 'Shank 3' 'Shank4'};
clear f
colors = {'k' 'g' 'r' 'b' 'm'}; 
    color_mat = [0 0 0; 0 .8 1; 0 0 1; 0 0.5 .4; 0 .7 .2]; % for graphing purposes (first is black, last is green)
spec_fig = figure;
for l = 1:length(lfp_all)
    for ii = 1:length(conds)
        [S(:,ii),f(ii,:),Serr{ii}] = mtspectrumc(lfp_all{l}(:,trial_conds{ii,:}),params);
%     [Slight_inf,flight_inf,errorlight_inf] = mtspectrumc(lfp_all{l}(:,light_trials),params); 
    end
    subplot(1,length(lfp_all),l)
%     for ii = 1:length(conds)
% %        plot(S(:,ii),f(ii,:),'l',Serr{ii},colors{ii})
%         h(ii) = plot(f(ii,:),10*log10(S(:,ii)),colors{ii});
%         hold on
%     end
    for ii = 1:length(conds)
        h = shadedErrorBar(f(ii,:),log10(S(:,ii))', [log10(S(:,ii))'-log10(Serr{ii}(1,:));-(log10(S(:,ii))'-log10(Serr{ii}(2,:)))],{'Color',color_mat(ii,:),'linewidth',2});    % shadedErrorBar requires error from mean; Serr gives actual value of errorbars
        leg_handles(ii) = h.mainLine;       % or h.patch
%         plot(f(ii,:),10*log10(Serr{ii}(2,:)),sprintf('%s:',colors{ii}))
        hold on
    end
    hold on
%     plot_vector(S(:,ii),f(ii,:),'n',Serr{l},'g')
    xlim([0 100])
    xlabel('Frequency (Hz)','fontsize',24)
    ylabel('Log power','fontsize',24)
    title(layer_titles(l),'fontsize',24)
    hold off
    %     legend(leg_handles,'No light','1Hz','10Hz','20Hz','40Hz');
    set(gca,'fontsize',18)
    l=legend(leg_handles,'Light OFF',legend_labels{:});
    set(l,'fontsize',18)
end
    

% save chronux plot
xSize = 24; ySize = 11;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0 0 xSize*50 ySize*50])
save_name= 'Power_by_shank';
save('power_by_shank.mat','S','f','Serr')
print(gcf,'-dpng',save_name)

%% running
move_trials = trial_type(:,4);
run_conds = unique(move_trials);
for l = 1:length(lfp_all)
    for ii = 1:length(run_conds)
        [S(:,ii),f(ii,:),Serr{ii}] = mtspectrumc(lfp_all{l}(:,intersect(trial_conds{1,:},find(move_trials==run_conds(ii)))),params);
%     [Slight_inf,flight_inf,errorlight_inf] = mtspectrumc(lfp_all{l}(:,light_trials),params); 
    end
    subplot(1,length(lfp_all),l)
%     for ii = 1:length(conds)
% %        plot(S(:,ii),f(ii,:),'l',Serr{ii},colors{ii})
%         h(ii) = plot(f(ii,:),10*log10(S(:,ii)),colors{ii});
%         hold on
%     end
    for ii = 1:length(run_conds)
        h = shadedErrorBar(f(ii,:),log10(S(:,ii))', [log10(S(:,ii))'-log10(Serr{ii}(1,:));-(log10(S(:,ii))'-log10(Serr{ii}(2,:)))],{'Color',color_mat(ii,:),'linewidth',2});
        leg_handles(ii) = h.mainLine;       % or h.patch
%         plot(f(ii,:),10*log10(Serr{ii}(2,:)),sprintf('%s:',colors{ii}))
        hold on
    end
    hold on
%     plot_vector(S(:,ii),f(ii,:),'n',Serr{l},'g')
    xlim([0 80])
    xlabel('Frequency (Hz)')
    ylabel('Log power')
    title(layer_titles(l))
    hold off
%     legend(leg_handles,'No light','1Hz','10Hz','20Hz','40Hz');
    legend(leg_handles(ii),'Stationary','Running');
%     [a,b,c,d] = legend('No light','Low intensity','Medium intensity','High intensity');
end
% save chronux plot
xSize = 24; ySize = 11;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0 0 xSize*50 ySize*50])
save_name= 'Power_by_layer_runningvstationary';
save('power_by_layer_runningvstationary.mat','S','f','Serr')
print(gcf,'-dpng',save_name)

figure;
title('Running trials')
for l = 1:length(lfp_all)
    for ii = 1:length(conds)
        [Srun(:,ii),f(ii,:),Srunerr{ii}] = mtspectrumc(lfp_all{l}(:,intersect(trial_conds{ii,:},find(move_trials))),params);
%     [Slight_inf,flight_inf,errorlight_inf] = mtspectrumc(lfp_all{l}(:,light_trials),params); 
    end
    subplot(1,length(lfp_all),l)
%     for ii = 1:length(conds)
% %        plot(S(:,ii),f(ii,:),'l',Serr{ii},colors{ii})
%         h(ii) = plot(f(ii,:),10*log10(S(:,ii)),colors{ii});
%         hold on
%     end
    for ii = 1:length(conds)
        h = shadedErrorBar(f(ii,:),log10(Srun(:,ii))', [log10(Srun(:,ii))'-log10(Srunerr{ii}(1,:));-(log10(Srun(:,ii))'-log10(Srunerr{ii}(2,:)))],{'Color',color_mat(ii,:),'linewidth',2});
        leg_handles(ii) = h.mainLine;       % or h.patch
%         plot(f(ii,:),10*log10(Serr{ii}(2,:)),sprintf('%s:',colors{ii}))
        hold on
    end
    hold on
%     plot_vector(S(:,ii),f(ii,:),'n',Serr{l},'g')
    xlim([0 80])
    xlabel('Frequency (Hz)')
    ylabel('Log power')
    title(layer_titles(l))
    hold off
%     legend(leg_handles,'No light','1Hz','10Hz','20Hz','40Hz');
    legend(leg_handles,'No light','Low intensity','Medium intensity','High intensity');
%     [a,b,c,d] = legend('No light','Low intensity','Medium intensity','High intensity');
end
% save chronux plot
xSize = 24; ySize = 11;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0 0 xSize*50 ySize*50])
save_name= 'Power_by_layer_runningtrials';
save('power_by_layer_running.mat','S','f','Serr')
print(gcf,'-dpng',save_name)

figure;
title('Stationary trials')
for l = 1:length(lfp_all)
    for ii = 1:length(conds)
        [Sstat(:,ii),f(ii,:),Sstaterr{ii}] = mtspectrumc(lfp_all{l}(:,intersect(trial_conds{ii,:},intersect(vis_trials,find(~move_trials)))),params);
%     [Slight_inf,flight_inf,errorlight_inf] = mtspectrumc(lfp_all{l}(:,light_trials),params); 
    end
    subplot(1,length(lfp_all),l)
%     for ii = 1:length(conds)
% %        plot(S(:,ii),f(ii,:),'l',Serr{ii},colors{ii})
%         h(ii) = plot(f(ii,:),10*log10(S(:,ii)),colors{ii});
%         hold on
%     end
    for ii = 1:length(conds)
        h = shadedErrorBar(f(ii,:),log10(Sstat(:,ii))', [log10(Sstat(:,ii))'-log10(Sstaterr{ii}(1,:));-(log10(Sstat(:,ii))'-log10(Sstaterr{ii}(2,:)))],{'Color',color_mat(ii,:),'linewidth',2});
        leg_handles(ii) = h.mainLine;       % or h.patch
%         plot(f(ii,:),10*log10(Serr{ii}(2,:)),sprintf('%s:',colors{ii}))
        hold on
    end
    hold on
%     plot_vector(S(:,ii),f(ii,:),'n',Serr{l},'g')
    xlim([0 100])
    xlabel('Frequency (Hz)','fontsize',24)
    ylabel('Log power','fontsize',24)
            title(layer_titles(l),'fontsize',24)
        hold off
    %     legend(leg_handles,'No light','1Hz','10Hz','20Hz','40Hz');
        set(gca,'fontsize',18)
        l=legend(leg_handles,'Light OFF',legend_labels{:});
        set(l,'fontsize',18)
end
% save chronux plot
xSize = 24; ySize = 11;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0 0 xSize*50 ySize*50])
save_name= 'Power_by_layer_stationarytrials';
save('power_by_layer_stationary.mat','S','f','Serr')
print(gcf,'-dpng',save_name)

return