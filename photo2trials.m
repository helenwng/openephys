function [phototrial, trial_idx] = photo2trials(photo)

    %find threshold crossings of photo signal
    photoOn = [];
    photoOff = [];
    thresh = min(photo)+0.5*(max(photo)-min(photo));
    photo_zero = photo-thresh;
    eps = 0.005;
    photoXs = find(abs(photo_zero)<eps);
    
    %^that will give redundant crossings for on and off times, so this part
    %picks the closest one to the threshold for each on/off and categorizes
    %them as on or off
    for x = 1:length(photoXs)

        %if current crossing is ON (positive slope)
        if photo_zero(photoXs(x)-3)<0&photo_zero(photoXs(x)+3)>0,
            photoOn = [photoOn photoXs(x)]; %Add it to the On list
            
            %if current crossing is within 10 data points away from the
            %last one, it's probably redundant, so pick the closest and
            %throw the other one out
            if length(photoOn)>1&(photoXs(x)-photoOn(end-1))<10,
                if abs(photo_zero(photoXs(x)))<abs(photo_zero(photoOn(end-1)))
                    photoOn(end-1) = [];
                else
                    photoOn(end) = [];
                end
            end
        
            %else if current crossing is OFF (negative slope)
        elseif photo_zero(photoXs(x)-3)>0&photo_zero(photoXs(x)+3)<0,
            photoOff = [photoOff photoXs(x)]; %Add it to the off list
            
            %if current crossing is within 10 data points away from the
            %last one, it's probably redundant, so pick the closest and
            %throw the other one out
            if length(photoOff)>1&(photoXs(x)-photoOff(end-1))<10,
                if abs(photo_zero(photoXs(x)))<abs(photo_zero(photoOff(end-1)))
                    photoOff(end-1) = [];
                else
                    photoOff(end) = [];
                end
            end
        else
            photoXs(x) = NaN; %neither positive or negative slope, this would be a local extrema so just jitter - get rid of it
        end

    end
    
    %populate a vector, size(photo), with trial on/off signal
    phototrial = zeros(size(photo));
    IPI = diff(photoOn); %find the intervals between photo onsets - mode will be the small interval within a trial
    trialOn = photoOn([1 1+find(IPI>2*mode(IPI))]); %picks the first onset and all that are more than twice the short interval after the previous
    trialOff = [(floor((photoOn(find(IPI>2*mode(IPI)))+photoOff(find(IPI>2*mode(IPI))))./2)) floor((photoOn(end)+photoOff(end))/2)]; %Epoc off signal is halfway between the on and offsets of the last photo pulse, so this averages the on and offsets of the last pulse in each trial 
    
    %make ones
    for i = 1:length(trialOn)
        phototrial(trialOn(i):trialOff(i)) = ones(size(trialOn(i):trialOff(i)));
    end
    
    trial_idx = [trialOn' trialOff'];
end