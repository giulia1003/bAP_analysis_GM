function [norm_amplitudes_vect,full_map_bAP_times_and_channels] = get_amplitudes(voltage_map, bAP_times_and_channels, nrn, plotParams, NP_probe) 

 % inputs: voltage_map, nrn, plotParams
   % getting amplitude of signal at soma and above around spike time
  try
 if isempty(bAP_times_and_channels)
	 norm_amplitudes_vect = [];
	 full_map_bAP_times_and_channels = [];
	 
 else
	 
     min_value = min(voltage_map(:)); 
    [soma_depth,soma_time] = find(voltage_map == min_value);
	
	% maps have been centered around the soma and somatic spike time so you need to get the right times and channels for bAP extent
	diff_ch = soma_depth-(bAP_times_and_channels(1,2));  % first value is soma ch
	diff_t = soma_time-(bAP_times_and_channels(1,1));    % first value is soma time
	
	bAP_times = bsxfun(@plus, bAP_times_and_channels(:,1), diff_t);
	bAP_chans = bsxfun(@plus, bAP_times_and_channels(:,2), diff_ch);

	full_map_bAP_times_and_channels = [bAP_times bAP_chans];
	
	amplitude_at_soma = voltage_map(soma_depth,soma_time);
	
	bAP_amplitudes_vect = [];

	% get vector of amplitudes acorss significant channels and time points
	
	if isequal(NP_probe, '1.0')
	for i = 1:size(bAP_chans,1)
	bAP_amplitudes = voltage_map(bAP_chans(i),bAP_times(i));
	bAP_amplitudes_vect = [bAP_amplitudes_vect bAP_amplitudes];
	end
	
	elseif isequal(NP_probe, 'UHD') || isequal(NP_probe, 'UHD_4')
    bAP_chans = bAP_chans-50;
	bAP_times = bAP_times-141;
	for i = 1:size(bAP_chans,1)
	bAP_amplitudes = voltage_map(bAP_chans(i),bAP_times(i));
	bAP_amplitudes_vect = [bAP_amplitudes_vect bAP_amplitudes];
	end		
	
	end
	% normalise to somatic value
	norm_amplitudes_vect = bAP_amplitudes_vect./amplitude_at_soma;

 end
	
  catch ME 
	  rethrow(ME)
	  
full_map_bAP_times_and_channels = [];	  
norm_amplitudes_vect = [];	  
  end
% OLD METHOD 

%     % taking window around spike time
%     window = round([size(voltage_map, 2)/2 + plotParams.soma_window(1), size(voltage_map, 2)/2 + plotParams.soma_window(2)]);
%     
% 	% get minimum values around spike time (map is cut around soma spike)
% 	% (column index = ch_min_times_sp)
%     [waveForms_min_val, ch_min_times_sp] = min(voltage_map(:, window(1):window(2))');
% 	
% 	% same for maximum values
%     [waveForms_max_val, ch_max_times_sp] = max(voltage_map(:, window(1):window(2))');
%     
% 	% getting time point where minima and maxima where recorded (rows are channels)
%     ch_min_times_sp = floor(ch_min_times_sp + (size(voltage_map, 2)/2 + plotParams.soma_window(1)));
%     ch_max_times_sp = floor(ch_max_times_sp + (size(voltage_map, 2)/2 + plotParams.soma_window(1)));
%     
%     meanWFs_diff_val_window = waveForms_max_val - waveForms_min_val;
%    
% 	% normalised amplitudes vector
% 	% max(meanWFs_diff_val_window) should correspond to amplitude at soma at time of spike
% 	
%     norm_amp_val = meanWFs_diff_val_window / max(meanWFs_diff_val_window);
%     
%     soma_depth = find(norm_amp_val == max(norm_amp_val));   % max val should be 1
%     soma_time = ch_min_times_sp (soma_depth);  
%    
    
    %% splitting normalised vector into 'below soma' and 'above soma' - getting AUC
    
%     soma_depth = find(norm_amplitudes_vect == max(norm_amplitudes_vect));
%     
%     below_soma = norm_amplitudes_vect(1:soma_depth);
%     above_soma = norm_amplitudes_vect(soma_depth:end);   
%     % above soma could also start 150um from soma (signal amplitude differences will be less evident close to the soma)
%     AUC_ratio = trapz(above_soma)/trapz(below_soma); 
%     AUC_diff = trapz(above_soma)-trapz(below_soma); 
% 
% 
%     AUC_ratio_neurons(nrn) = AUC_ratio; % saving value into vector
%     AUC_diff_neurons(nrn) = AUC_diff;


	  	
end

