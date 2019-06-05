function [refr_idx, SNR] = eval_clusters(good_units,spike_times,clusters,nchs,amp_sr,exp_system,exp_path)
% evaluate cluster quality

for n = 1:length(good_units)
    % get autocorrelogram and refractory period violations
    disp(strcat('Cluster quality analysis: ', num2str(good_units(n))))
    unit_times = spike_times(clusters==good_units(n));
    unit_times_ms = round(unit_times./(amp_sr/1000));
    unit_spiketrain = zeros(max(unit_times_ms),1);
    unit_spiketrain(unit_times_ms) = 1;
    [~, lags] = xcorr(unit_spiketrain);
    refr_idx(n) = sum(lags<=2&lags>=-2)/length(unit_times_ms)*100;     % %refractory period violations

%     % get SNR for unit
%     SNR(n) = waveformSNR(good_units(n),unit_times,nchs,exp_system,exp_path);
    SNR(n) = nan;
end

end