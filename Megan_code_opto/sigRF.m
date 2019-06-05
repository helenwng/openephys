function [pvals, issig] = sigRF(psth_mat,rfMap)

nX = size(rfMap{1},1);
nY = size(rfMap{1},2);
tbins = size(psth_mat,1);
nstim = size(psth_mat,2)/length(rfMap);
reps = 1000;
for i = 1:length(rfMap)
    RF = rfMap{i}(:,:);
    %
    rfShuf = zeros(nX,nY,reps);
%     rfShuf_int = zeros(nX*8,nY*8,reps);
    for n=1:reps
        shuf = randperm(tbins*nstim);
        tmp_mat = psth_mat(:,(i-1)*nstim+1:i*nstim);
        shuf_psth = reshape(tmp_mat(shuf),tbins,nstim);
        [Ushuf,~,~] = svd(shuf_psth','econ');
        rfShuf(:,:,n) = reshape(Ushuf(:,1),12,10);
%         rfShuf_int(:,:,n) = imresize(rfShuf(:,:,n),[12*8 10*8],'bilinear');
    end
    pvals{i} = zeros(size(RF));
    for xx = 1:nX
        for yy = 1:nY
            if RF(xx,yy)<0
                pvals{i}(xx,yy) = -sum(rfShuf(xx,yy,:)<=RF(xx,yy))/reps;
            else
                pvals{i}(xx,yy) = sum(rfShuf(xx,yy,:)>=RF(xx,yy))/reps;
            end
        end
    end
    pvals{i}(abs(pvals{i})>=.025) = nan;   % because two-sided test. do I need a multiple comparisons correction?

    if sum(~isnan(pvals{i}(:)))
        issig(i) = 1;
    else
        issig(i) = 0;
    end

end