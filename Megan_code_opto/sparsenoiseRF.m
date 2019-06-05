function [psth_norm, rfMap, stats] = sparsenoiseRF(unit_times_ds,stim_times,binsize,T,stim_trials,stim_key)
% description here:

% unit_times_ds = 
% stim_times = 
% binsize = 
% T = stimulus duration


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

N = sum(sum(spike_raster));
% psth = nan(nX*nY*nBW*length(diff_conds),(T/binsize)/1000);
% psth = nan(nX,nY,(T/binsize)/1000,nBW,length(diff_conds));
num_stim = size(stim_key,1);
X = unique(stim_trials(:,1));
nX = length(X);
Y = unique(stim_trials(:,2));
nY = length(Y);
BW = unique(stim_trials(:,3));
nBW = length(BW);
lconds = unique(stim_trials(:,4));
psth = nan(T/binsize,num_stim);
count = 1;
for lc = 1:length(lconds)
    for b = 1:nBW
        for y = 1:nY
            for x = 1:nX
                which_trials = (stim_trials(:,1) == x)&(stim_trials(:,2) == y)&(stim_trials(:,3)==BW(b))&(stim_trials(:,4)==lconds(lc));
                [psth(:,count),~] = make_psth_v2(binsize/1000,[0:binsize:T]./1000,which_trials',spike_raster',ones(1,length(which_trials)));
                count = count+1;
            end
        end
    end
end

% %% trying simple STA (11/30/18)
% test = reshape(psth(:,1:120)',nX,nY,size(psth,1));
% % test_norm = test./max(test(:));
% % test_bs = test_norm(:,:,1);
% test_bl = test(:,:,1);
% test_bs = test - repmat(test_bl,1,1,size(test,3));
% test_norm = test_bs./max(abs(test_bs(:)));
% figure;
% for i=1:size(test,3)
%     subplot(5,5,i)
% %     imagesc((squeeze(test_norm(:,:,i)-test_bs))')
%     imagesc(squeeze(test_norm(:,:,i)'))
%     colorbar
% end
% % [a,b] = max((test_norm-test_bs),[],3);
% % [~,maxind] = max(a(:));
% % maxT = b(maxind);
% % peakVar = var(reshape(test_norm(:,:,maxT),1,nX*nY));
% % noiseVar = mean(var(reshape(test_norm(:,:,21:end),5,nX*nY),0,2));
% % SNR(1) = peakVar/noiseVar;
% [a,b] = max(abs(test_norm),[],3);
% [~,maxind] = max(abs(a(:)));
% maxT = b(maxind);
% peakVar = var(reshape(test_norm(:,:,maxT),1,nX*nY));
% noiseVar = mean(var(reshape(test_norm(:,:,21:end),5,nX*nY),0,2));
% 
% test = reshape(psth(:,121:240)',nX,nY,size(psth,1));
% test_norm = test./max(test(:));
% test_bs = test_norm(:,:,1);
% figure;
% for i=1:size(test,3)
%     subplot(5,5,i)
%     imagesc((squeeze(test_norm(:,:,i)-test_bs))')
%     colorbar
% end
% [a,b] = max((test_norm-test_bs),[],3);
% [~,maxind] = max(a(:));
% maxT = b(maxind);
% peakVar = var(reshape(test_norm(:,:,maxT),1,nX*nY));
% noiseVar = mean(var(reshape(test_norm(:,:,21:end),5,nX*nY),0,2));
% SNR(2) = peakVar/noiseVar;
% 
% test = reshape(psth(:,241:360)',nX,nY,size(psth,1));
% test_norm = test./max(test(:));
% test_bs = test_norm(:,:,1);
% figure;
% for i=1:size(test,3)
%     subplot(5,5,i)
%     imagesc((squeeze(test_norm(:,:,i)-test_bs))')
%     colorbar
% end
% [a,b] = max((test_norm-test_bs),[],3);
% [~,maxind] = max(a(:));
% maxT = b(maxind);
% peakVar = var(reshape(test_norm(:,:,maxT),1,nX*nY));
% noiseVar = mean(var(reshape(test_norm(:,:,21:end),5,nX*nY),0,2));
% SNR(3) = peakVar/noiseVar;
% 
% 
% test = reshape(psth(:,361:end)',nX,nY,size(psth,1));
% test_norm = test./max(test(:));
% test_bs = test_norm(:,:,1);
% figure;
% for i=1:size(test,3)
%     subplot(5,5,i)
%     imagesc((squeeze(test_norm(:,:,i)-test_bs))')
%     colorbar
% end
% [a,b] = max((test_norm-test_bs),[],3);
% [~,maxind] = max(a(:));
% maxT = b(maxind);
% peakVar = var(reshape(test_norm(:,:,maxT),1,nX*nY));
% noiseVar = mean(var(reshape(test_norm(:,:,21:end),5,nX*nY),0,2));
% SNR(4) = peakVar/noiseVar;


% % find peak response for each subfield and each light condition
% % psth_bls = psth - repmat(psth(1,:),size(psth,1),1);
% for lc = 1:length(lconds)
%     for b = 1:nBW
%         K = psth(:,stim_key(:,3)==BW(b) & stim_key(:,4)==lconds(lc));
%         [max_bin{b,lc}, max_stim{b,lc}] = find(ismember(K, max(K(:))));
%         max_r(b,lc) = K(max_bin{b,lc}(1),max_stim{b,lc}(1));
%     end
% end
% max_r = max(max_r);        % get max for each light cond across subfields
% max_r_mat = zeros(size(psth));    
% for lc = 1:length(lconds)
%     max_r_mat(:,stim_key(:,4)==lconds(lc)) = max_r(lc);
% end
% psth_norm = psth./max_r_mat;
% % all_max = max(max_r(:));        % get cell's max firing rate
% % psth_norm = psth./all_max;
% [lowx,lowy] = find(psth_norm < .35);
% for i=1:length(lowx)
%     psth_norm(lowx(i),lowy(i)) = 0;
% end
% RFmap_on = zeros(nX,nY,T/binsize,length(lconds));
% RPmap_off = RFmap_on;
% for lc = 1:length(lconds)
%     for x=1:nX
%         for y = 1:nY
%             for t=1:T/binsize
%                 RFmap_on(x,y,t,lc) = psth_norm(t,stim_key(:,1)==x & stim_key(:,2)==y & stim_key(:,3)==1 & stim_key(:,4)==lconds(lc));
%                 RFmap_off(x,y,t,lc) = psth_norm(t,stim_key(:,1)==x & stim_key(:,2)==y & stim_key(:,3)==-1 & stim_key(:,4)==lconds(lc));
%             end
%         end
%     end
% end

% psth_norm = zeros(size(psth));
psth_norm = zeros(size(psth,1)-2,size(psth,2));
Model = zeros(size(psth_norm'));
baseline = zeros(1,size(psth_norm,2)/(nX*nY));
% for nn = 1:nBW
%     baseline([nn nn+2]) = mean(psth(1,(nn-1)*(nX*nY)+1:nn*nX*nY));     % baseline = start of first time bin for each bw condition in light OFF - use for exps with single pulses?? trials only
% end
for n = 1:size(psth_norm,2)/(nX*nY)
    baseline(n) = mean(psth(1,(n-1)*(nX*nY)+1:n*nX*nY));     % baseline = start of first time bin for each bw and lcond condition -- use for most exps
%     psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY) = psth(:,(n-1)*(nX*nY)+1:n*nX*nY)-baseline;
    psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY) = psth(3:end,(n-1)*(nX*nY)+1:n*nX*nY)-baseline(n);    % start from third time bin to avoid light artifacts
    [U,S,V] = svd(psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY)','econ');    
    rfMapVec{n} = U(:,1);
    timeCourse{n} = V(:,1)';
    peakTimePoint(n) = find(abs(timeCourse{n})==max(abs(timeCourse{n})),1);    
    ts(n) = sign(timeCourse{n}(peakTimePoint(n)));
    if n>2 
        if sign(timeCourse{n}(peakTimePoint(n-2))) < 0  || (ts(n) < 0 && ts(n-2) > 0) % matches sign flip with no light condition
            rfMapVec{n} = -rfMapVec{n};
            timeCourse{n} = -timeCourse{n};
        end            
    elseif ts(n)<0
        rfMapVec{n} = -rfMapVec{n};
        timeCourse{n} = -timeCourse{n};
    end
    rfMap{n} = reshape(rfMapVec{n},nX, nY);
    Scalar(n) = S(1,1);
    Model((n-1)*(nX*nY)+1:n*nX*nY,:) = rfMapVec{n}*timeCourse{n}*Scalar(n) + baseline(n);
%     Model((n-1)*(nX*nY)+1:n*nX*nY,:) = rfMapVec{n}*timeCourse{n}*Scalar(n);
    stats.timeCourse{n} = timeCourse{n};
    maxZ(n) = (max(rfMap{n}(:))-mean(rfMap{n}(:)))./std(rfMap{n}(:));
    minZ(n) = (min(rfMap{n}(:))-mean(rfMap{n}(:)))./std(rfMap{n}(:));
    stats.peakZscore(n) = max(abs([minZ(n) maxZ(n)]));
    
end
Residual = psth_norm' - Model;



% % % test if RF is significant (chi-squared test of independence across
% % % neighboring pixels of peak)
% for i=1:size(psth_norm,2)/(nX*nY)
%     tmp_psth = psth_norm(:,(i-1)*(nX*nY)+1:i*nX*nY);
%     small_stim_key = stim_key((i-1)*(nX*nY)+1:i*nX*nY,:);
% %     [~,b] = max(abs(tmp_psth(:)));
% %     peak = small_stim_key(ceil(b/size(tmp_psth,1)),1:2);
%     [~,max_which] =sort(abs(rfMap{i}(:)));
%     peak = small_stim_key(max_which(end),1:2);
% %     [~,max_inds] = sort(max(abs(tmp_psth)));
% %     peaks = small_stim_key(max_inds(end-2:end),1:2)
%     if peak(1)==1
%         first_xneighb = 1;
%         last_xneighb = peak(1)+1;
%     elseif peak(1) == nX
%         first_xneighb = peak(1)-1;
%         last_xneighb = nX;
%     else
%         first_xneighb = peak(1)-1;
%         last_xneighb = peak(1)+1;
%     end
%     if peak(2) == 1
%         first_yneighb = 1;
%         last_yneighb = peak(2) +1;
%     elseif peak(2) == nY
%         first_yneighb = peak(2)-1;
%         last_yneighb = nY;
%     else
%         first_yneighb = peak(2)-1;
%         last_yneighb = peak(2) +1;
%     end
%     [p,q] = meshgrid(first_xneighb:last_xneighb,first_yneighb:last_yneighb);
%     pairs = [p(:) q(:)];
%     which = zeros(1,length(pairs));
%     for ii = 1:length(pairs)
%         which(ii) = find((small_stim_key(:,1)==pairs(ii,1))&(small_stim_key(:,2)==pairs(ii,2)));
%     end
% %     [~,max_which] =sort(abs(rfMap{i}(:)));
% % [~,max_which] =sort(max(abs(tmp_psth(:,which))));
% %     neighb_psth = num2cell(tmp_psth(:,which(max_which(end-6:end))),1);      % int8 to save memory
%     try
%         neighb_psth = num2cell(int8(tmp_psth(:,which)),1);      % int8 to save memory
%         [~, ~,pval,~] = crosstab(neighb_psth{:});
%     catch
%         fprintf('9 is too many\n')
%         [~,max_which] =sort(max(abs(tmp_psth(:,which))));
%         neighb_psth = num2cell(int8(tmp_psth(:,which(max_which(end-7:end)))),1);      % int8 to save memory
%         try
%             [~, ~,pval,~] = crosstab(neighb_psth{:});
%         catch
%             fprintf('8 is still too many\n')
%             neighb_psth = num2cell(int8(tmp_psth(:,which(max_which(end-6:end)))),1); 
%             [~, ~,pval,~] = crosstab(neighb_psth{:});
%         end
%     end
%     if pval < .05
%         stats.issig(i) = 1;
%     else
%         stats.issig(i) = 0;
%         stats.pvals{i} = nan(size(rfMap{i}));
%     end
% end
% 
% %% make null distributions for statistical comparison
% % if find(stats.issig)      % only do it if any statistically sig RFs were identified (because this takes time)
%     tic
%     rfMap_shuf = nan(nX,nY,100);
%     reps = 1000;
%     nshufstim = num_stim/(length(lconds)*nBW);
%     trial_mat = repmat(1:nshufstim,1,length(which_trials)/nshufstim);
%     for ii = 1:reps
%         shuf_trials = trial_mat(randperm(length(which_trials)));    % randomly assign each stim presentation to one of nX*nY stimulus conditions
%         [psth_shuf,~] = make_psth_v2(binsize/1000,[0:binsize:T]./1000,1:length(which_trials),spike_raster',shuf_trials);
%         bs = mean(psth_shuf(:,1));
%         psth_shuf = psth_shuf(:,3:end) - bs;
%         [Ushuf,~,~] = svd(psth_shuf,'econ');        % get spatial filter for baseline-subtracted shuffled psth
%         rfMapVec_shuf = Ushuf(:,1);
%         rfMap_shuf(:,:,ii) = reshape(rfMapVec_shuf,nX, nY);
%     end
%     toc
% 
%     for i=1:size(psth_norm,2)/(nX*nY)
%         pvals{i} = zeros(size(rfMap{i}));
%         for xx = 1:nX
%             for yy = 1:nY
%                 if rfMap{i}(xx,yy)<0
%                     pvals{i}(xx,yy) = -sum(rfMap_shuf(xx,yy,:)<=rfMap{i}(xx,yy))/reps;
%                 else
%                     pvals{i}(xx,yy) = sum(rfMap_shuf(xx,yy,:)>=rfMap{i}(xx,yy))/reps;
%                 end
%             end
%         end
%     %     pvals{i}(abs(pvals{i})>=.05/(nX*nY)) = nan;   % because two-sided test. do I need a multiple comparisons correction?
%         [sort_pvals,ord] = sort(abs(pvals{i}(:)));      % treat all pvals as +
%         testp = zeros(1,length(sort_pvals));
%         for k=1:length(sort_pvals)
%             testp(k) =  sort_pvals(k)<=(k/(nX*nY))*.05;      % Benjamini-Hochberg FDR correction for multiple comparisons
%         end
%         pvals_ord = pvals{i}(ord);
%         max_p = max(abs(pvals_ord(find(testp))));
%         if isempty(max_p)
%             max_p = 0;
%         end
%         pvals{i}(abs(pvals{i})>max_p) = nan; 
%         stats.pvals{i} = pvals{i};
% 
%         if sum(~isnan(pvals{i}(:)))
%             stats.issig(i) = 1;
%         else
%             stats.issig(i) = 0;
%         end
%     end
% % else
% %     for i=1:4
% %         stats.pvals{i} = nan(size(rfMap{1}));
% %     end
% % end

%%
figure;
timeBins = 0:binsize:T;
timeBins = timeBins(1:end-1)+binsize/2;
nCol = 4;
timeBins=timeBins(3:end);       %TEMP

key_leg = {'OFF subunit, light OFF', 'ON subunit, light OFF', 'OFF subunit, light ON', 'ON subunit, light ON'};
for n = 1:size(psth_norm,2)/(nX*nY)
    
    subplot(4,nCol,n);
    plot(timeBins, timeCourse{n});
    xlim([timeBins(1) timeBins(end)]);
    title('response time course')
    xlabel('time')
    
    ax = subplot(4,nCol,nCol+n);
    imagesc(1:size(psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY),2), timeBins,  psth_norm(:,(n-1)*(nX*nY)+1:n*nX*nY));
    %         cax = caxis();
    %         caxis(max(abs(cax))*[-1 1]);
    % colormap(colorcet('L3'));
    %         axis image
    colorbar
    title(sprintf('all PSTHs - %s',key_leg{n}));
    xlabel('space')
    ylabel('time')
    
    subplot(4,nCol,2*nCol+n); 
    imagesc(Residual((n-1)*(nX*nY)+1:n*nX*nY,:)');
    %         cax = caxis();
    %         caxis(max(abs(cax))*[-1 1]);     
    %         axis image
    colorbar
    title(sprintf('residual - %s',key_leg{n}))
    xlabel('space')
    ylabel('time')
    
    subplot(4,nCol,3*nCol+n);
    imagesc(1:nX, 1:nY, rfMap{n}(:,:)');
    axis image
    colorbar
    title(sprintf('%s map, peakZ = %.2f', key_leg{n}, stats.peakZscore(n)));
    
%     subplot(5,nCol,16+n)
%     imagesc(1:nX,1:nY,stats.pvals{n}');
%     axis image
%     caxis([-.05 .05])
%     colorbar
%     title(sprintf('%s RF sig map (issig = %.2f)',key_leg{n},stats.issig(n)));
end

%% Trying new thing - look at variance across X and Y dimensions of calculated spatial filter and compare that against "shuffled" spatial filters
for i=1:length(rfMap)
    varY = var(rfMap{i},0,2); % get the variance across Y values for every X
    varX = var(rfMap{i},0,1); % get the variance across X values for every Y
    [peakVarY,whereVarY] = max(varY); % get X coordinate and value with maximum variance in Y
    [peakVarX,whereVarX] = max(varX); % get Y coordinate and value with maximum variance in X
    shufmat = zeros(nX,nY,1000);
    for n=1:1000
        tmp =randperm(numel(rfMap{i}));
        shufmat(:,:,n) = reshape(tmp,size(rfMap{i}));   % shuffle the spatial filter values across positions in nX by nY grid
    end
    shufvarY = squeeze(var(shufmat,0,2)); % get the variance across Y values for every X in every shuffle
    shufvarX = squeeze(var(shufmat,0,1)); % get the variance across X values for every Y in every shuffle
    shufSNR_Y = zeros(1,1000);
    if whereVarY>nX/2   % if peak in Y variances is in right half of grid
        spSNR_Y = peakVarY/mean(varY(1:2));   % get ratio of the peak Y variance to the mean of Y variances in first two columns (aka "spatial SNR")    
    else
        spSNR_Y = peakVarY/mean(varY(nX-1:nX));  % if the peak in Y variances is in the left half of the grid, get ratio of peak Y variance to the mean of Y variances in last two columns
    end
    [shufpeaksY,wherepeaksY] = max(shufvarY);   % get peak of X variance for each shuffle
    shufSNR_Y(wherepeaksY>nX/2) = shufpeaksY(wherepeaksY>nX/2)./mean(shufvarY(1:2,wherepeaksY>nX/2)); % get same SNR val from every shuffle to get a null distribution of SNRs (using max variance from each shuffle)
    shufSNR_Y(wherepeaksY<=nX/2) = shufpeaksY(wherepeaksY<=nX/2)./mean(shufvarY(nX-1:nX,wherepeaksY<=nX/2));
    
    shufSNR_X = zeros(1,1000);
    if whereVarX>nY/2   % if peak in X variances is in lower half of grid
        spSNR_X = peakVarX/mean(varX(1:2));   % get ratio of the peak X variance to the mean of X variances in first two rows (aka "spatial SNR")    
    else
        spSNR_X = peakVarX/mean(varX(nY-1:nY));  % if the peak in X variances is in the upper half of the grid, get ratio of peak X variance to the mean of X variances in last two rows
    end
    [shufpeaksX,wherepeaksX] = max(shufvarX);   % get peak of X variance for each shuffle
    shufSNR_X(wherepeaksX>nY/2) = shufpeaksX(wherepeaksX>nY/2)./mean(shufvarX(1:2,wherepeaksX>nY/2)); % get same SNR val from every shuffle to get a null distribution of SNRs (using max variance from each shuffle)
    shufSNR_X(wherepeaksX<=nY/2) = shufpeaksX(wherepeaksX<=nY/2)./mean(shufvarX(nY-1:nY,wherepeaksX<=nY/2));
    
    probSigSNR_Y(i) = sum(shufSNR_Y>=spSNR_Y)/1000;
    probSigSNR_X(i) = sum(shufSNR_X>spSNR_X)/1000;
    
    if probSigSNR_Y(i) <= .01 && probSigSNR_X(i) <= .01 % consider it a significant receptive field if there's a significant (<=.01) spatial SNR in BOTH x and y dimensions
        stats.issig(i) = 1;
    else
        stats.issig(i) = 0;
    end
    subplot(4,nCol,3*nCol+i);
    title(sprintf('%s map, peakZ = %.2f, sig = %d', key_leg{i}, stats.peakZscore(i),stats.issig(i)));

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