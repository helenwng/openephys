function [tuning_curve tuning_curve_SE] = makeTCfromRast(spike_mat,preMS,stimMS,sumPeriod,oridom,trial_num,trial_subset,normalizing,figtitle)
% this function will generate a tuning curve from a raster plot

% [tuning_curve tuning_curve_SE] = makeTCfromRast(spike_mat,preMS,stimMS,sumPeriod,oridom,trials,normalizing)
%   spike_mat: your matrix of spikes / raster plot
%   preMS: the pre stimulus period, in MS
%   stimMS: length of the stimulus, in MS
%   sumPeriod: how much of your stimulus period would you like to use to generate the tuning period (e.g. first 200 ms, all/preMS)
%   oridom: your orientation domain (e.g. 0:30:330)
%   trial_num: your matrix of all trial information. if this has two
%   columns, both will be plotted
%   trial_subset: which trials to use. enter 'all' if this is irrelevant
%    normalizing: normalize to 'blank','prestim', or 'none'


 % make pre and stimulus rasters   
 prestim_rast=spike_mat(1:preMS,:);
 stim_rast=spike_mat(preMS:preMS+sumPeriod,:); %tuning curve for designated period

 % if you want to use all trials, just create a vector of 1:length(trial_num)
 if strcmp(trial_subset,'all')
     trial_subset = 1:length(trial_num);
 end

assert(ndims(trial_num) ==2,'This code currently only works for 2 dimensions (e.g. ori and LED)'); % if you have two variables (e.g. orientation/LED)
 
    for i=1:length(oridom) %get average firing rate for each orientation
        [TC1rast{i},c]=find(trial_num(trial_subset,1)==oridom(i)&trial_num(trial_subset,2)==0); %for the first condition
        TC1(i)=sum(sum(stim_rast(:,TC1rast{i}))'/(size(stim_rast,1)/1000))/length(TC1rast{i}); %get mean firing rate (/s)          
        TC1_SE(i)= std(sum(stim_rast(:,TC1rast{i}))'/(size(stim_rast,1)/1000))/sqrt(length(TC1rast{i})); %get standard error of mean FR
            if strcmp(normalizing,'prestim');
                TC1spont(i)= sum(sum(prestim_rast(:,TC1rast{i}))'/(size(prestim_rast,1)/1000))/length(TC1rast{i});
            end
                
        [TC2rast{i},c]=find(trial_num(trial_subset,1)==oridom(i)&trial_num(trial_subset,2)==1); %for the second condition
        TC2(i)=sum(sum(stim_rast(:,TC2rast{i}))'/(size(stim_rast,1)/1000))/length(TC2rast{i}); %get mean firing rate            
        TC2_SE(i)= std(sum(stim_rast(:,TC1rast{i}))'/(size(stim_rast,1)/1000))/sqrt(length(TC2rast{i})); %get standard error of mean FR
            if strcmp(normalizing,'prestim');
                TC2spont(i)= sum(sum(prestim_rast(:,TC2rast{i}))'/(size(prestim_rast,1)/1000))/length(TC2rast{i});
            end

    end
    
    if strcmp(normalizing,'blank'); %currently just using ALL trials because there aren't LED off blanks
        %get average of blank trials
        blank_trials=find(trial_num(trial_subset,1)==256); %256 is a random number determined by looper script
        %blank_check = exist('blank_trials','var');
        blank_check = length(blank_trials);
        assert(blank_check>0,'You specified to use blank trials. There are no blank trials.');
        TC1spont=ones(1,12)*sum(sum(stim_rast(:,blank_trials))'/(size(stim_rast,1)/1000))/length(blank_trials);
        blank_SE= std(sum(stim_rast(:,blank_trials))')/sqrt(length(blank_trials));
        
        blank_trials2=find(trial_num(trial_subset,1)==256); %256 is a random number determined by looper script
        blank_check2 = exist('blank_trials','var');
        assert(blank_check2==1,'You specified to use blank trials. There are no blank trials.');
        TC2spont=ones(1,12)*sum(sum(stim_rast(:,blank_trials2))'/(size(stim_rast,1)/1000))/length(blank_trials2);
        blank_SE2= std(sum(stim_rast(:,blank_trials2))')/sqrt(length(blank_trials2));                

    elseif strcmp(normalizing,'none');
        TC1spont(i)= 0; %set TCspont to 0 so it has no effect on tuning curve (no normalization)
        TC2spont(i)= 0; %set TCspont to 0 so it has no effect on tuning curve (no normalization)

    end
    
    TC1_corrected=TC1-TC1spont;
    TC2_corrected=TC2-TC2spont;

    tuning_curve(1,:) = TC1_corrected;
    tuning_curve(2,:) = TC2_corrected;
    tuning_curve_SE(1,:) = TC1_SE;
    tuning_curve_SE(2,:) = TC2_SE;
    
    
    if sum(isnan(TC1_corrected))>0 %make sure there aren't NaN trials
        TC1_corrected2=TC1_corrected(~isnan(TC1_corrected));
        deg=oridom(~isnan(TC1_corrected));
        plot(deg,TC1_corrected2,'k','linewidth',2);
    else
        shadedErrorBar(oridom,TC1_corrected,TC1_SE,{'k','linewidth',2},1);       
    end
    
    if sum(isnan(TC2_corrected))>0
        TC2_corrected2=TC2_corrected(~isnan(TC2_corrected));
        deg2=oridom(~isnan(TC2_corrected));
        hold on;plot(deg2,TC2_corrected2,'g','linewidth',2);
    else
        %    hold on;plot(oridom,ori_run,'b','linewidth',2);
        hold on;shadedErrorBar(oridom,TC2_corrected,TC2_SE,{'color',[0 .75 .75],'linewidth',2},1);
    end
    
    set(gca,'Fontsize', 10);
    set(gca,'Xlim',[0 oridom(end)]);
    set(gca, 'xtick', [0:60:330],'tickdir','out');        
    xlabel({'Orientation (degrees)'},'Fontsize',10);
    set(gca,'tickdir','out');     
    ylabel({'Mean FR (spks/s)'},'Fontsize',10);  
    title(strrep(figtitle,'_','\_'),'fontsize',10);  


