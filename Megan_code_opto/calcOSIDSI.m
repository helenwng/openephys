function [OSI,OSI_CV,DSI,DSI_CV] = calcOSIDSI(tuning_curve,deg)

%tuning_curve = vector of firing rates at each orientation 
%deg = vector of different orientations used, in degrees

tuning_curve_rep = repmat(tuning_curve,1,2);
oridom_rad = deg.*pi/180';      % convert to radians
dsiDistance = length(deg)/2;    % assumes even number of orientations
osiDistance = dsiDistance/2;

peaks(1) = find(tuning_curve_rep(1:length(tuning_curve_rep)/2)==max(tuning_curve_rep),1,'first');  % get indice of orientation with peak response
peaks(2) = find(tuning_curve_rep(length(tuning_curve_rep)/2+1:end)==max(tuning_curve_rep),1,'first')+length(tuning_curve_rep)/2;
if peaks(1) <= dsiDistance
    peak = peaks(2);
else
    peak = peaks(1);
end
OSI = (tuning_curve_rep(peak)-(mean([tuning_curve_rep(peak+osiDistance) tuning_curve_rep(peak-osiDistance)])))/(tuning_curve_rep(peak)+(mean([tuning_curve_rep(peak+osiDistance) tuning_curve_rep(peak-osiDistance)]))); %Determine OSI, max - response at orthogonal ori (90 degrees away)
DSI = (tuning_curve_rep(peak)-(mean([tuning_curve_rep(peak+dsiDistance) tuning_curve_rep(peak-dsiDistance)])))/(tuning_curve_rep(peak)+(mean([tuning_curve_rep(peak+dsiDistance) tuning_curve_rep(peak-dsiDistance)]))); %Determine DSI, max - response at ori 180 degrees away

if isnan(OSI)
    OSI=0;
end
if isnan(DSI)
    DSI=0;
end

OSI_CV = abs((sum(tuning_curve.*exp(1i*2*oridom_rad)))/sum(tuning_curve)); % Determine OSI, 1-CV
DSI_CV = abs((sum(tuning_curve.*exp(1i*oridom_rad)))/sum(tuning_curve)); % Determine OSI, 1-CV
return