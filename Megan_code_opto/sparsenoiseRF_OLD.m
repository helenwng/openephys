function sparsenoiseRF(unit_times_ds,stim_times,binsize,T,LED,amp_sr,exp_path)
  

%% extract spike times 
% spike_raster = make_raster_v2(unit_times_ds,field_trials,1000,4)';      % time x trials matrix of spike times (0s and 1s)
% spike_raster = make_raster_v2(unit_times_ds,stim_times(1:240:end),1000,60)';
spike_raster = make_raster_v2(unit_times_ds,stim_times,1000,T/1000)';
% [psth,~] = make_psth_v2(binsize/1000,[0:binsize:T]./1000,ones(1,size(spike_raster,2)),spike_raster',ones(1,size(spike_raster,2)));
% psth = psth';
% big_psth = psth(:);
% big_rast = spike_raster(:);     % turn into single time vector of spikes

% refresh = 4;        % vis stimuli per sec (hz)
% T = 60;        % time duration
% dt = 1/refresh;
% t = linspace(0,T,refresh*T+1);

stim_file = dir(fullfile(exp_path,'_*_*_*.mat'));
fprintf(sprintf('loading stimulus file %s\n',stim_file.name'))
stimdat = load(fullfile(exp_path,stim_file.name),'-mat');
analyze_file = dir(fullfile(exp_path,'*.analyzer'));
fprintf(sprintf('loading analyzer file %s\n',analyze_file.name'))
load(fullfile(exp_path,analyze_file.name),'-mat');  % under variable name 'Analyzer'
nX = Analyzer.P.param{find(cellfun(@(x) strcmp(x{1},'Nx'), Analyzer.P.param))}{3};
nY = Analyzer.P.param{find(cellfun(@(x) strcmp(x{1},'Ny'), Analyzer.P.param))}{3};
nBW = Analyzer.P.param{find(cellfun(@(x) strcmp(x{1},'bw_bit'), Analyzer.P.param))}{3};

num_reps = length(fieldnames(stimdat))-1;       % assumes one of the fields is for monitor refresh rate ('frate')
num_stim = length(stimdat.randlog_T1.seqs.xseq);

conds = cellfun(@(x) x.val{find(strcmp(Analyzer.loops.conds{1}.symbol,'light_bit'))}, Analyzer.loops.conds);
diff_conds = unique(conds);
trialcond = cellfun(@(x) x.repeats{1}.trialno,Analyzer.loops.conds);
lightcond = zeros(1,num_reps);
for i = 1:length(diff_conds)
    lightcond(ismember(trialcond,find(conds==diff_conds(i)))) = diff_conds(i);
end

div             = amp_sr/1000;
zx              = 1:div:length(LED);
izx             = floor(zx);
light = LED(izx);
% pho = photo(izx);

% stim_mat = zeros(size(big_rast,1),nX,nY);  % initialize 3d matrix for stimulus (ms in all trials x xposition x yposition)
% trial_stiminds = zeros(num_stim+1,num_reps);
% light_mat = zeros(size(big_rast,1),1);
% for n = 1:num_reps
%     xx = eval(sprintf('stimdat.randlog_T%d.seqs.xseq',n));
%     yy = eval(sprintf('stimdat.randlog_T%d.seqs.yseq',n));
%     bw = eval(sprintf('stimdat.randlog_T%d.seqs.bwseq',n'));
%     bw(bw==2) = -1; 
%     inds = stim_times((n-1)*num_stim+1:(n-1)*num_stim+num_stim)-stim_times((n-1)*num_stim+1)+1; 
%     trial_stiminds(:,n) = [inds + (n-1)*60000; n*60000];
%     for nn = 1:num_stim
%         stim_mat(trial_stiminds(nn,n):trial_stiminds(nn+1,n)-1,xx(nn),yy(nn)) = bw(nn);
%         if sum(light(stim_times((n-1)*num_stim+nn):stim_times((n-1)*num_stim+nn)+250)) > floor(max(light))*250
%             light_mat(trial_stiminds(nn,n):trial_stiminds(nn+1,n)-1) = 1;
%         end
%     end
% end
% stim_mat_ext = [zeros(T,nX,nY); stim_mat];

stim_trials = zeros(size(spike_raster,2),4);       %
spks_bystim = zeros(num_stim*length(diff_conds),size(spike_raster,1));
for n = 1:num_reps
    xx = eval(sprintf('stimdat.randlog_T%d.seqs.xseq',n));
    yy = eval(sprintf('stimdat.randlog_T%d.seqs.yseq',n));
    bw = eval(sprintf('stimdat.randlog_T%d.seqs.bwseq',n'));
    bw(bw==1) = -1; % 1 means BLACK
    bw(bw==2) = 1;  % 2 means WHITE
    for nn = 1:num_stim
        stim_trials((n-1)*num_stim+nn,1:3) = [xx(nn) yy(nn) bw(nn)]; 
        if sum(light(stim_times((n-1)*num_stim+nn):stim_times((n-1)*num_stim+nn)+size(spike_raster,1))) > floor(max(light))*size(spike_raster,1)
%         if find(light(stim_times((n-1)*num_stim+1):stim_times((n)*num_stim))>4) 	% TEMP
            stim_trials((n-1)*num_stim+nn,4) = 1;
        end
    end
end

lconds = unique(stim_trials(:,4));
bwconds = unique(stim_trials(:,3));
count = 1;
N = sum(sum(spike_raster));
% psth = nan(nX*nY*nBW*length(diff_conds),(T/binsize)/1000);
% psth = nan(nX,nY,(T/binsize)/1000,nBW,length(diff_conds));
psth = nan(T/binsize,num_stim*length(lconds));
stim_mat = zeros(nX,nY,nBW,T/binsize,num_stim,length(lconds));
% for lc = 1:length(diff_conds)
%     for x = 1:nX
%         for y = 1:nY
%             for b = 1:nBW
%                 which_trials = (stim_trials(:,1) == x)&(stim_trials(:,2) == y)&(stim_trials(:,3)==b)&(stim_trials(:,4)==lconds(lc));
% %                 [psth(count,:),~] = make_psth_v2(binsize,0:binsize:T/1000,which_trials',spike_raster',ones(1,size(spike_raster,2)));
%                 [psth(x,y,:,b,lc),~] = make_psth_v2(binsize,0:binsize:T/1000,which_trials',spike_raster',ones(1,size(spike_raster,2)));
%                 stim_mat(x,y,count) = b;
%                 count = count+1;
%             end
%         end
%     end
% end
stim_key = nan(num_stim,4);
for lc = 1:length(lconds)
    for b = 1:nBW
        for y = 1:nY
            for x = 1:nX
                which_trials = (stim_trials(:,1) == x)&(stim_trials(:,2) == y)&(stim_trials(:,3)==bwconds(b))&(stim_trials(:,4)==lconds(lc));
                [psth(:,count),~] = make_psth_v2(binsize/1000,[0:binsize:T]./1000,which_trials',spike_raster',ones(1,length(which_trials)));
                stim_mat(x,y,b,:,:,lc) = bwconds(b);
                stim_key(count,:) = [x y bwconds(b) lconds(lc)];
                count = count+1;
            end
        end
    end
end
% find peak response for each subfield and each light condition
% psth_bls = psth - repmat(psth(1,:),size(psth,1),1);
for lc = 1:length(lconds)
    for b = 1:nBW
        K = psth(:,stim_key(:,3)==bwconds(b) & stim_key(:,4)==lconds(lc));
        [max_bin{b,lc}, max_stim{b,lc}] = find(ismember(K, max(K(:))));
        max_r(b,lc) = K(max_bin{b,lc}(1),max_stim{b,lc}(1));
    end
end
max_r = max(max_r);        % get max for each light cond across subfields
max_r_mat = zeros(size(psth));    
for lc = 1:length(lconds)
    max_r_mat(:,stim_key(:,4)==lconds(lc)) = max_r(lc);
end
psth_norm = psth./max_r_mat;
% all_max = max(max_r(:));        % get cell's max firing rate
% psth_norm = psth./all_max;
[lowx,lowy] = find(psth_norm < .35);
for i=1:length(lowx)
    psth_norm(lowx(i),lowy(i)) = 0;
end
RFmap_on = zeros(nX,nY,T/binsize,length(lconds));
RPmap_off = RFmap_on;
for lc = 1:length(lconds)
    for x=1:nX
        for y = 1:nY
            for t=1:T/binsize
                RFmap_on(x,y,t,lc) = psth_norm(t,stim_key(:,1)==x & stim_key(:,2)==y & stim_key(:,3)==1 & stim_key(:,4)==lconds(lc));
                RFmap_off(x,y,t,lc) = psth_norm(t,stim_key(:,1)==x & stim_key(:,2)==y & stim_key(:,3)==-1 & stim_key(:,4)==lconds(lc));
            end
        end
    end
end
% figure;
% title('Light off, off-subfields')
% for t=1:T/binsize
%     subplot(5,6,t)
%     imagesc(1:nX,1:nY,squeeze((RFmap_off(:,:,t,1)))',[0 1]);
% end
% figure;
% title('Light off, on-subfields')
% for t=1:T/binsize
%     subplot(5,6,t)
%     imagesc(1:nX,1:nY, squeeze((RFmap_on(:,:,t,1)))',[0 1]);
% end
% figure;
% title('Light on, off-subfields')
% for t=1:T/binsize
%     subplot(5,6,t)
%     imagesc(1:nX,1:nY,squeeze((RFmap_off(:,:,t,2)))',[0 1]);
% end
% figure;
% title('Light on, on-subfields')
% for t=1:T/binsize
%     subplot(5,6,t)
%     imagesc(1:nX,1:nY,squeeze((RFmap_on(:,:,t,2)))',[0 1]);
% end

% Xlight = reshape(stim_mat(:,:,:,:,:,1),nX*nY*nBW,T/binsize*nX*nY*nBW);
% Xnolight = reshape(stim_mat(:,:,:,:,:,2),nX*nY*nBW,T/binsize*nX*nY*nBW);
% Ylight = reshape(psth(:,(find(lconds>0)-1)*num_stim+1:(find(lconds>0))*num_stim),1,num_stim*T/binsize);
% Ynolight = reshape(psth(:,(find(lconds==0)-1)*num_stim+1:(find(lconds==0))*num_stim),1,num_stim*T/binsize);
% STAlight = (1/sum(Ylight)).*(Xlight*Ylight');
% STAnolight = (1/sum(Ynolight)).*(Xnolight*Ynolight');
% RFlight = reshape(STAlight,nX,nY,nBW);


% % X = reshape(stim_mat(:,:,1:num_stim),nX*nY,T);
% X = reshape(psth(:,:,:,1,1),nX*nY,(T/binsize)/1000);
% baseline = mean(X(:,1)); % take the first bin as the baseline
% % X(X==2) = -1;
% [U,S,V] = svd(X - baseline,'econ');
% baseline = mean(psth(1,:));

% psth_norm = zeros(size(psth));
psth_norm = zeros(size(psth,1)-4,size(psth,2));
Model = zeros(size(psth_norm'));
baseline = mean(psth(1,stim_key(:,4)==0));  % use when only did short pulses?
for n = 1:size(psth_norm,2)/(nX*nY)
%     baseline = mean(psth(1,(n-1)*(nX*nY)+1:n*nX*nY));     % use when doing sustained light stim?
%     psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY) = psth(:,(n-1)*(nX*nY)+1:n*nX*nY)-baseline;
    psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY) = psth(5:end,(n-1)*(nX*nY)+1:n*nX*nY)-baseline;    % use when only did short pulses?
    [U,S,V] = svd(psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY)','econ');    
    rfMapVec{n} = U(:,1);
    timeCourse{n} = V(:,1)';
    peakTimePoint(n) = find(abs(timeCourse{n})==max(abs(timeCourse{n})),1);    
    ts = sign(timeCourse{n}(peakTimePoint(n)));
    if ts<0
        rfMapVec{n} = -rfMapVec{n};
        timeCourse{n} = -timeCourse{n};
    end
    rfMap{n} = reshape(rfMapVec{n},nX, nY);
    Scalar(n) = S(1,1);
    Model((n-1)*(nX*nY)+1:n*nX*nY,:) = rfMapVec{n}*timeCourse{n}*Scalar(n) + baseline;
%     Model((n-1)*(nX*nY)+1:n*nX*nY,:) = rfMapVec{n}*timeCourse{n}*Scalar(n);
    stats.timeCourse{n} = timeCourse{n};
    maxZ(n) = (max(rfMap{n}(:))-mean(rfMap{n}(:)))./std(rfMap{n}(:));
    minZ(n) = (min(rfMap{n}(:))-mean(rfMap{n}(:)))./std(rfMap{n}(:));
    stats.peakZscore(n) = max(abs([minZ(n) maxZ(n)]));
end
Residual = psth_norm' - Model;

figure;
timeBins = 0:binsize:T;
timeBins = timeBins(1:end-1)+binsize/2;
nCol = 4;
timeBins=timeBins(5:end);       %TEMP

key_leg = {'OFF subunit, light OFF', 'ON subunit, light OFF', 'OFF subunit, light ON', 'ON subunit, light ON'};
for n = 1:size(psth_norm,2)/(nX*nY)
    
    subplot(4,nCol,n);
    plot(timeBins, timeCourse{n});
    xlim([timeBins(1) timeBins(end)]);
    title('response time course')
    xlabel('time')
    
    ax = subplot(4,nCol,4+n);
    imagesc(1:size(psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY),2), timeBins,  psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY));
    %         cax = caxis();
    %         caxis(max(abs(cax))*[-1 1]);
    % colormap(colorcet('L3'));
    %         axis image
    colorbar
    title(sprintf('all PSTHs - %s',key_leg{n}));
    xlabel('space')
    ylabel('time')
    
    subplot(4,nCol,8+n); 
    imagesc(Residual((n-1)*(nX*nY)+1:n*nX*nY,:)');
    %         cax = caxis();
    %         caxis(max(abs(cax))*[-1 1]);     
    %         axis image
    colorbar
    title(sprintf('residual - %s',key_leg{n}))
    xlabel('space')
    ylabel('time')
    
    subplot(4,nCol,12+n);
    imagesc(1:nX, 1:nY, rfMap{n}(:,:)');
    axis image
    colorbar
    title(sprintf('%s map, peakZ = %.2f', key_leg{n}, stats.peakZscore(n)));
end





% % 
% % psth_long = reshape(psth',1,size(psth,2)*size(psth,1));
% % STA = (1/N).*(X*psth_long);
% 
% 
% 
% N = sum(big_rast);          % total number of spikes
% spks = find(big_rast);
% spks_light = find(big_rast&light_mat);
% spks_nolight = find(big_rast&~light_mat);
% N_light = length(spks_light);
% N_nolight = length(spks_nolight);
% X_light = zeros(nX*nY*(T+1),N_light);
% for i=1:N_light
%     X_light(:,i) = reshape(permute(stim_mat_ext(spks_light(i):spks_light(i)+T,:,:),[2 3 1]),nX*nY*(T+1),1);
% end
% STAlight = (1/N_light).*(X_light*big_rast(spks_light));
% STAlight = reshape(STAlight,nX,nY,T+1);
% figure;imagesc(mean(STAlight(:,:,1:150),3))
% colorbar
% figure
% for n = 1:20
%     subplot(4,5,n)
%     imagesc(mean(STAlight(:,:,(n-1)*10+1:n*10),3))
% end
% 
% clear X_light
% X_nolight = zeros(nX*nY*(T+1),N_nolight);
% for i=1:N_nolight
%     X_nolight(:,i) = reshape(permute(stim_mat_ext(spks_nolight(i):spks_nolight(i)+T,:,:),[2 3 1]),nX*nY*(T+1),1);
% end
% STAnolight = (1/N_nolight).*(X_nolight*big_rast(spks_nolight));
% STAnolight = reshape(STAnolight,nX,nY,T+1);
% figure;imagesc(mean(STAnolight(:,:,1:150),3))
% colorbar
% figure
% for n = 1:20
%     subplot(4,5,n)
%     imagesc(mean(STAnolight(:,:,(n-1)*10+1:n*10),3))
% end
% 
% clear X_nolight
% X = zeros(nX*nY*(T+1),N);
% for i=1:N
%     X(:,i) = reshape(permute(stim_mat_ext(spks(i):spks(i)+T,:,:),[2 3 1]),nX*nY*(T+1),1);
% end
% STA = (1/N).*(X*big_rast(spks));
% STA = reshape(STA,nX,nY,T+1);
% figure;imagesc(mean(STA(:,:,1:150),3))
% colorbar
% figure
% for n = 1:20
%     subplot(4,5,n)
%     imagesc(mean(STA(:,:,(n-1)*10+1:n*10),3))
% end
% 
% % sta = ifft(fft(big_rast).*conj(fft(stim_mat,[],1)),[],1);
% % sta_nolight = ifft(fft(big_rast(light_mat==0)).*conj(fft(stim_mat(light_mat==0,:,:),[],1)),[],1);
% % N_nolight = sum(big_rast(light_mat==0));
% % sta_light = ifft(fft(big_rast(light_mat==1)).*conj(fft(stim_mat(light_mat==1,:,:),[],1)),[],1);
% % N_light = sum(big_rast(light_mat==1));
% 
% 
%     
return