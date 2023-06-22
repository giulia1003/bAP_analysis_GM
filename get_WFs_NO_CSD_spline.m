function [splineField] = get_WFs_NO_CSD_spline(info, plotParams, gwfparams, NP_probe,spiketimes)

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
	    
	
%     if numel(gwfparams.spikeTimes) < gwfparams.nWf - 3 % the minus 3 operation here is to prevent a bug that it skips a unit because gwfparams.spikeTimes was truncated so you don't index negatively in the getWaveForms function when the unit has very early spike times 
% %   
% %         fprintf(2, ['not enough spikes for: ' num2str(nrn_ID), ' index: ', num2str(nrn) '\n'])
% 
%     %   return
%     end

    %% FOR SELECTED SPIKES


	gwfparams.spikeTimes = spiketimes;

 try


[wf_shuf] = getWaveForms_shuf(gwfparams);
	

	
 catch ME
	 rethrow(ME)

 end
 
			     %% 'preprocessing'

	if size(wf_shuf.waveFormsMean,1) > 1
        
	wfMean = wf_shuf.waveFormsMean;
	else
    wfMean = squeeze(wf_shuf.waveFormsMean(1,:,:));
    end
    
    % subtracting channel mean

    sub_mean_mat = repmat(mean(wfMean, 2), 1, size(wfMean, 2));
    wfMean_mean = wfMean - sub_mean_mat;    
	 
	
    win_size = gwfparams.wfWin(2) - gwfparams.wfWin(1);
	

    %% getting spline interpolation
	
	
     if isequal(NP_probe, '2.0')
    [splineField, xColumns] = computeSplineLFP_NoHistology_2(wfMean_mean, wfMean_mean, plotParams.equalChans, info, win_size);
	 elseif isequal(NP_probe, '1.0')
 [splineField, xColumns] = computeSplineLFP_NoHistology(wfMean_mean, wfMean_mean, plotParams.equalChans, info, win_size);
 
	 elseif isequal(NP_probe, 'UHD') 
		 
	xColumns = [18 12];
	splineField = wfMean_mean;

	 elseif isequal(NP_probe, 'UHD_4')
	
	 [splineField, xColumns] = computeSpline_UHD(wfMean_mean, wfMean_mean, plotParams.equalChans, info, win_size);

	 end
	
	 
    
			
end
