function [wfMean_mean_shift,splineField_shift,spiketimes,wf_shift,single_spike_map,spike_time_tab,shuffled_map_tab] = get_WFs_NO_CSD(sp,nrn_ID, nrn, info, plotParams, gwfparams, NP_probe,drift_nrn,batch_id,templateYpos,hippocampal_units)

	% This function first averages field potential across all 385 channels.
	% Then it deals with parallel channels: meanField is averaged signal
	% for parallel channels and splineField is taking two columns of
	% channels with highest amplitudes and doing cubic spline interpolation
	% along the missing distances in each, since they have
	% a 40 um intersite interval. Spline also interpolates over channels that are reference or
	% defined as being broken.  It then averages two values that correspond
	% to the same depth (see PowerPoint).
	% CSD is then done the same for both of these methods by just taking
	% the second spatial derivative (to get meanCSD and splineCSD).
	
	% even though meanField and meanCSD aren't used later on they are still
	% calculated and the same functions can be used on them as on
	% splineField ('getSoma', 'getFeats', 'plot_meanWFs', 'plot_trace',
	% 'plot_dist_amp' etc.)
	
	% TODO: use xColumns to check that the same two columns are used for
	% spline under different conditions!!! + include number of spikes used
	% for plotting!!!
	

    fprintf(['Number of spikes for ' num2str(nrn_ID) ' is: ' num2str(numel(gwfparams.spikeTimes)) '\n'])
    
	
%     if numel(gwfparams.spikeTimes) < gwfparams.nWf - 3 % the minus 3 operation here is to prevent a bug that it skips a unit because gwfparams.spikeTimes was truncated so you don't index negatively in the getWaveForms function when the unit has very early spike times 
% %   
% %         fprintf(2, ['not enough spikes for: ' num2str(nrn_ID), ' index: ', num2str(nrn) '\n'])
% 
%     %   return
%     end

    %% FOR SELECTED SPIKES

    fprintf(['Getting spikes for neuron: ' num2str(nrn_ID)  ' (index: ', num2str(nrn) ')\n'])
	
        %% could literally just feed spike times into the next function (aka ISOLATE based on clu first, then on WHEEL VELOCITY)
%         gwfparams.spikeTimes = ceil(sp.st(sp.clu==nrn_ID)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
%         gwfparams.spikeClusters = sp.clu(sp.clu==nrn_ID);
%         
%         gwfparams.spikeTimes = ceil(selected_st(selected_clu==nrn_ID)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
%         gwfparams.spikeClusters = selected_clu(selected_clu==nrn_ID);
if isequal(NP_probe, '1.0')
[drift_stamps,spiketimes,nrn_idx]  = get_drift(nrn,nrn_ID,sp,gwfparams,drift_nrn,batch_id,templateYpos,NP_probe);

else 
drift_stamps = [];
spiketimes = gwfparams.spikeTimes;
end
   
% OG getWaveForms (no drift correction)
%   		[wf] = getWaveForms(gwfparams);

% new getWaveForms2 = with drift correction, added by GM
% tic
% [wf_shift] = getWaveForms2(gwfparams,drift_stamps,spiketimes,templateYpos,nrn,hippocampal_units,NP_probe,info);    
% toc

tic 
[wf_shift] = getWaveForms2_parallelcomp(gwfparams,drift_stamps,spiketimes,templateYpos,nrn,hippocampal_units,NP_probe,info, plotParams);    
toc

% 	if isequal(gwfparams.analyse_data, 'field')
% fprintf(['Taking ',num2str(gwfparams.nWf),' spikes for bAP analysis \n'])
if info.single_spike_maps 

if info.take_all_spikes
n_of_spikes = gwfparams.use_nWfs;

elseif ~(info.take_all_spikes)
%  only used to get concatenated voltage_maps of n spikes
if numel(spiketimes) < 300
n_of_spikes = numel(spiketimes); % minimum to compare differences between single spikes
else
n_of_spikes = gwfparams.use_nWfs;
end
end

try
[single_spike_map,spike_time_tab,shuffled_map_tab] = get_single_spike_maps(gwfparams,drift_stamps,spiketimes,templateYpos,nrn, plotParams, info,n_of_spikes,hippocampal_units, NP_probe);
	
catch ME
	rethrow(ME)

end

else
	
single_spike_map = []; spike_time_tab = []; shuffled_map_tab = [];
end

% elseif isequal(gwfparams.analyse_data, 'LFP')
% 		wf = getWaveForms_LFP_medianSubt(gwfparams); % this is getWaveForms_LFP2 with template-, median- and mean subtraction
% 	end
	
 
  %% 'preprocessing'

if size(wf_shift.waveFormsMean,1) > 1   % running parallel computing getWaveforms will give a 385x401 matrix that doesn't need to be squeezed
	wfMean_shift =  wf_shift.waveFormsMean;
else
	wfMean_shift = squeeze(wf_shift.waveFormsMean(1,:,:));
end
	
	if isequal(gwfparams.analyse_data, 'both')
	wfMean_LFP = squeeze(wf_shift.waveFormsMeanLFP(1,:,:));	
	
	sub_mean_shift = repmat(mean(wfMean_LFP, 2), 1, size(wfMean_LFP, 2));
    wfMean_mean_LFP = wfMean_LFP - sub_mean_shift;  
	end
  % subtracting channel mean
 
sub_mean_shift = repmat(mean(wfMean_shift, 2), 1, size(wfMean_shift, 2));
wfMean_mean_shift = wfMean_shift - sub_mean_shift;  

if ~isequal(NP_probe, 'UHD_4')
% median subtraction
medianwv = median(wfMean_mean_shift, 1, 'omitnan');
sub_median_mat_shift = repmat(medianwv, size(wfMean_shift, 1), 1);
wfMean_median_shift = wfMean_mean_shift - sub_median_mat_shift; 
end

    %% getting spline interpolation
	
    win_size = gwfparams.wfWin(2) - gwfparams.wfWin(1);  % window around spike time

	 if isequal(NP_probe, '2.0')
  [splineField_shift, xColumns] = computeSplineLFP_NoHistology_2(wfMean_mean_shift, wfMean_mean_shift, plotParams.equalChans, info, win_size);
	 elseif isequal(NP_probe, '1.0')
 [splineField_shift, xColumns] = computeSplineLFP_NoHistology(wfMean_mean_shift, wfMean_mean_shift, plotParams.equalChans, info, win_size);
		if isequal(gwfparams.analyse_data, 'both')
    win_size = gwfparams.wfWinLFP(2) - gwfparams.wfWinLFP(1);  % window around spike time			
 [splineFieldLFP, xColumns] = computeSplineLFP_NoHistology(wfMean_mean_LFP, wfMean_mean_LFP, plotParams.equalChans, info, win_size);
		end			
	 elseif isequal(NP_probe, 'UHD') 
	xColumns = [18 12];
	splineField_shift = wfMean_mean_shift;

	 elseif isequal(NP_probe, 'UHD_4')
 [splineField_shift, xColumns] = computeSpline_UHD(wfMean_mean_shift, wfMean_mean_shift, plotParams.equalChans, info, win_size);
		 

	 end
	 
  
			
end
