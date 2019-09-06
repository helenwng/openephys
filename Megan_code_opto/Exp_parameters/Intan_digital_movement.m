function move_trials = Intan_digital_movement(field_trials,mmvtA,mmvtB,show_plot)

%% Movement decoding for the digital encoder w/ intan
% Ethan McBride - 5/5/2016

num_trials = size(field_trials,1);
distVect = zeros(2,2,num_trials);

for t=1:num_trials

    tempA = mmvtA(field_trials(t,1)*20:field_trials(t,2)*20);
    tempB = mmvtB(field_trials(t,1)*20:field_trials(t,2)*20);
    k=2;
    for i=2:length(tempA)
        if tempA(i)~=tempA(i-1) %if there is a change in mmvtA
            if tempA(i)==1 %Aup
                if tempB(i)==0 %if B is down, forward motion
                    %forward
                    distVect(k,1,t) = distVect(k-1,1,t)+0.25;
                    distVect(k,2,t) = i/20000;
                    k=k+1;
                elseif tempB(i)==1 %if B is up, backward motion
                    %backward
                    distVect(k,1,t) = distVect(k-1,1,t)-0.25;
                    distVect(k,2,t) = i/20000;
                    k=k+1;
                end
            elseif tempA(i)==0 %Adown
                if tempB(i)==1 %if B is up, forward motion
                    %forward
                    distVect(k,1,t) = distVect(k-1,1,t)+0.25;
                    distVect(k,2,t) = i/20000;
                    k=k+1;
                elseif tempB(i)==0 %if B is down, backward motion
                    %backward
                    distVect(k,1,t) = distVect(k-1,1,t)-0.25;
                    distVect(k,2,t) = i/20000;
                    k=k+1;
                end
            end
        elseif tempB(i)~=tempB(i-1) %if there is a change in mmvtB
            if tempB(i)==1 %Bup
                if tempA(i)==1 %if A is up, forward motion
                    %forward
                    distVect(k,1,t) = distVect(k-1,1,t)+0.25;
                    distVect(k,2,t) = i/20000;
                    k=k+1;
                elseif tempA(i)==0 %if A is down, backward motion
                    %backward
                    distVect(k,1,t) = distVect(k-1,1,t)-0.25;
                    distVect(k,2,t) = i/20000;
                    k=k+1;
                end
            elseif tempB(i)==0 %Bdown
                if tempA(i)==0 %if A is down, forward motion
                    %forward
                    distVect(k,1,t) = distVect(k-1,1,t)+0.25;
                    distVect(k,2,t) = i/20000;
                    k=k+1;
                elseif tempA(i)==1 %if A is up, backward motion
                    %backward
                    distVect(k,1,t) = distVect(k-1,1,t)-0.25;
                    distVect(k,2,t) = i/20000;
                    k=k+1;
                end
            end
        end
    end
        
    
end


    distVect(:,1,:) = distVect(:,1,:)./256 .* 47.75; %if the circumference is 47.75cm
    
    speed = diff(distVect(:,1,:))./diff(distVect(:,2,:));

    
    maxdistances = reshape(max(abs(distVect(:,1,:))),num_trials,1);
    maxspeeds = reshape(max(abs(speed(:,1,:))),num_trials,1); % TEMP
    meanspeeds = reshape(nanmean(abs(speed(:,1,:))),num_trials,1);

    meanspeeds(isnan(meanspeeds))=0;
    maxspeeds(isnan(maxspeeds))=0;
    
    move_trials = zeros(num_trials,1);
    move_trials(find(meanspeeds>=2)) = 1;
  
 if show_plot
    figure;
    subplot(211)
    hist(maxdistances,20)
    ylabel('#trials')
    xlabel('distance traveled (cm)')
    subplot(212)
    hist(meanspeeds,20)
    ylabel('#trials')
    xlabel('mean speed (cm/sec)')
 end

return