% load('H:\LPproject\LPresults\VH3_LP_PVChR2\RFresults')
% load('H:\LPproject\LPresults\VH3_LP_PVChR2\VH3_LP_PVChR2_cluster_quality.mat')
cd('H:\LPproject\LPresults\MH24_LP')
load('RFresults.mat')
% load('H:\LPproject\LPresults\MH18_LP\MH18_LP_cluster_quality.mat')
% good_units = intersect(find(SNR>=1.5&refr_idx<1),find(shank==0));
% allsig = reshape([stats.issig],4,length(stats))';
% sig_mat = [allsig(:,1)+allsig(:,3) allsig(:,2)+allsig(:,4)];
% good_units(sum(sig_mat(good_units,:),2)<1) = [];

if ~exist('RFmap_analysis','dir')
    mkdir('RFmap_analysis');
end
cd('RFmap_analysis')

sig = reshape([stats.issig],length(stats(1).issig),length(good_units))';
sig2 = reshape([stats.issig2],length(stats(1).issig),length(good_units))';
good_units(sum(sig2,2)<1) = [];
stats(sum(sig2,2)<1) = [];
rfMap(sum(sig2,2)<1) = [];
sig2(sum(sig2,2)<1,:) = [];
sig(sum(sig2,2)<1,:) = [];
nX = size(rfMap{1}{1},1);
nY = size(rfMap{1}{1},2);
sig_byfield = [sum(sig2(:,[1 3]),2) sum(sig2(:,[2 4]),2)];

area = nan(length(good_units),4);
overlap = nan(length(good_units),2);

for n = 1:length(good_units)
    unit = good_units(n);
%     reps = 10000;
%     test = zeros(reps,numel(psth_norm{unit}(:,1:120)));
%     test_psth = zeros(size(psth_norm{unit}(:,1:120),1),size(psth_norm{unit}(:,1:120),2),reps);
%     rfShuf = zeros(12,10,reps);
%     rfShuf_int = zeros(12*8,10*8,reps);
%     for n=1:reps
%         test(n,:) = randperm(numel(psth_norm{unit}(:,1:120)));
%         test_psth(:,:,n) = reshape(psth_norm{unit}(test(n,:)),size(psth_norm{unit},1),120);
%         [Ushuf,~,~] = svd(test_psth(:,:,n)','econ');
%         rfShuf(:,:,n) = reshape(Ushuf(:,1),12,10);
%         rfShuf_int(:,:,n) = imresize(rfShuf(:,:,n),[12*8 10*8],'bilinear');
%     end
%     pvals = zeros(size(RF));
%     for xx = 1:size(RF,1)
%         for yy = 1:size(RF,2)
%             if RF(xx,yy)<0
%                 pvals(xx,yy) = -sum(rfShuf(xx,yy,:)<=RF(xx,yy))/reps;
%             else
%                 pvals(xx,yy) = sum(rfShuf(xx,yy,:)>=RF(xx,yy))/reps;
%             end
%         end
%     end
%     pvals(abs(pvals)>=.05) = nan;   % do I need a multiple comparisons correction?
% 
%     % shuf_psth = mean(test_psth,3);
%     figure;imagesc(1:120,20:10:240,shuf_psth)
%     [U,S,V] = svd(shuf_psth','econ');
%     rfTest = U(:,1);
%     rfTest = reshape(rfTest,12,10);
%     figure;imagesc(1:12,1:10,rfTest')

    f1 = figure;
    for i = 1:length(rfMap{n})
%         if sig_byfield(n,rem(i+1,2)+1)
            RF = rfMap{n}{i}(:,:);
            ZRF = (RF-mean(RF(:)))./std(RF(:));
            pmap = stats(n).pvals{i};
            if sig2(n,i) && sig_byfield(n,rem(i+1,2)+1)==2     % only compare maps for which there was a significant RF in both light on and light off conds
                zmap = abs(ZRF)>=2;
            else
                zmap = zeros(size(RF));
            end
            subplot(length(rfMap{n}),5,1+(i-1)*5)
            imagesc(1:nX,1:nY,ZRF')
            RFnorm = RF./max(abs(RF(:)));
%             RFnorm(isnan(pmap)) = 0;
            RFnorm(~zmap) = 0;
            subplot(length(rfMap{n}),5,2+(i-1)*5)
            imagesc(1:nX,1:nY,RFnorm')
            testRF = imresize(RFnorm,[nX*8 nY*8],'bilinear'); % bilinear interpolation to 1deg resolution
            subplot(length(rfMap{n}),5,3+(i-1)*5)
            imagesc(1:nX,1:nY,testRF')
            filtRF = imgaussfilt(testRF, 5); % smooth with 2d gaussian filter (5 stdev = 5*(1deg resolution?) = 5degs??)
            subplot(length(rfMap{n}),5,4+(i-1)*5)
            imagesc(1:nX,1:nY,filtRF')
    %         filtRF(abs(filtRF(:))<.5) = 0;
            subplot(length(rfMap{n}),5,5+(i-1)*5)
    %         imagesc(1:12,1:10,filtRF')
            colorbar
            I{i} = imbinarize(abs(filtRF)).*sign(filtRF);
            imagesc(1:nX,1:nY,I{i}')
            bw{i} = bwboundaries(I{i});
            % get rid of "rfs" not part of the main one
            if ~isempty(bw{i})
                [~,which_bord(i)] = max(cellfun(@(x) length(x),bw{i},'uniformoutput',1));      % if multiple borders, only use biggest one
                [p,q] = meshgrid(min(bw{i}{which_bord(i)}(:,1)):max(bw{i}{which_bord(i)}(:,1)),min(bw{i}{which_bord(i)}(:,2)):max(bw{i}{which_bord(i)}(:,2)));
                filt_grid = zeros(size(I{i}));
                filt_grid(p(:),q(:)) = 1;
                I{i} = I{i}.*filt_grid;
                perim(i) = size(bw{i}{which_bord(i)},1);
            end
            area(n,i)= bwarea(I{i});
            if i>2
                overlap(n,i-2) = (sum(abs(I{i-2}(:)).*abs(I{i}(:))))/sum(abs(I{i-2}(:))); % as % of no-light condition
            end
%         end
    end


    f2 = figure; 
    hold on
    for i=1:length(rfMap{n})

        if rem(i,2)
            c = 'b-';
            subplot(2,2,ceil(i/2))
            hold on
        else
            c = 'r-';
        end
        if ~isempty(bw{i})
            plot(bw{i}{which_bord(i)}(:,1),bw{i}{which_bord(i)}(:,2),c,'linewidth',4)
        end
        view(0,270)
        xlim([0 size(filtRF,1)])
        ylim([0 size(filtRF,2)])

    end

    subplot(223)
    imagesc(1:nX,1:nY,[I{3}-I{1}]')
    title('OFF subfields: Light ON - Light OFF')
    subplot(224)
    imagesc(1:nX,1:nY,[I{4}-I{2}]')
    title('ON subfields: Light ON - Light OFF')
    save_name1 = sprintf('Fields_Cluster%d',good_units(n));
    save_name2 = sprintf('Optocomp_Cluster%d',good_units(n));
    print(f1,'-dpng',save_name1)
    print(f2,'-dpng',save_name2)
    close all
   
end

