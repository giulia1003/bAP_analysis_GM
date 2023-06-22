function [bAP_map_info] = get_bap_map3(info,analyse_nrn,min_map,bAP_times_and_channels,voltage_map,...
	soma_depth_sp_Field,bAP_extent,signal_axon_only,soma_trace,above_trace,ind_analyse_nrn,bAP_map_info,high_trace,wf_shift,velocity,single_spike_map,spike_time_tab,shuffled_map_tab,amplitudes_vect, correlation_matrix)
% bAP_map_info = [];
if size(analyse_nrn,2) == 1
	ind_analyse_nrn = 1;
end
	
	
bAP_map_info(ind_analyse_nrn).nrn = analyse_nrn(ind_analyse_nrn);
	bAP_map_info(ind_analyse_nrn).voltage_map = voltage_map;
    bAP_map_info(ind_analyse_nrn).single_spike_map = single_spike_map;
	bAP_map_info(ind_analyse_nrn).single_spike_times = spike_time_tab;
	bAP_map_info(ind_analyse_nrn).shuffled_spike_map = shuffled_map_tab;
    bAP_map_info(ind_analyse_nrn).bAP_times_and_channels = bAP_times_and_channels;
	bAP_map_info(ind_analyse_nrn).min_map = min_map;
% 	bAP_map_info(ind_analyse_nrn).spline = spline;
	bAP_map_info(ind_analyse_nrn).soma_depth_sp_Field = soma_depth_sp_Field;
	bAP_map_info(ind_analyse_nrn).soma_amplitude = wf_shift.min_value_soma;
    bAP_map_info(ind_analyse_nrn).correlation_coeff = correlation_matrix;
	if isempty(voltage_map)   % setting fields for units thar didn't have enough spikes
    bAP_map_info(ind_analyse_nrn).single_spike_map = [];
	bAP_map_info(ind_analyse_nrn).frequency_sp = [];
    bAP_map_info(ind_analyse_nrn).frequency_bursts = [];
	bAP_map_info(ind_analyse_nrn).n_of_single_spikes = [];
	
	if ~(info.analyse_whisk_behaviour)
	if isequal(info.analyse_running,'running')
	bAP_map_info(ind_analyse_nrn).n_all_spikes_run = [];
	elseif isequal(info.analyse_running,'stationary')
	bAP_map_info(ind_analyse_nrn).n_all_spikes_stat = [];
	end
	end
	
	bAP_map_info(ind_analyse_nrn).bAP_extent = [];
	bAP_map_info(ind_analyse_nrn).axonal_extent = [];
	bAP_map_info(ind_analyse_nrn).soma_trace =  [];
	bAP_map_info(ind_analyse_nrn).above_soma_trace = [];
	bAP_map_info(ind_analyse_nrn).high_trace = [];
	bAP_map_info(ind_analyse_nrn).bAP_conduction_velocity = [];
	
	else % for neurons with enough spikes
	
	bAP_map_info(ind_analyse_nrn).bAP_extent = bAP_extent;
	bAP_map_info(ind_analyse_nrn).axonal_extent = signal_axon_only;
	bAP_map_info(ind_analyse_nrn).amplitudes = amplitudes_vect;
	bAP_map_info(ind_analyse_nrn).soma_trace =  soma_trace;
	bAP_map_info(ind_analyse_nrn).above_soma_trace = above_trace;
	bAP_map_info(ind_analyse_nrn).high_trace = high_trace;
% 	bAP_map_info(ind_analyse_nrn).soma_amplitude_across_time = wf_shift.soma_amplitude_over_time;
	bAP_map_info(ind_analyse_nrn).bAP_conduction_velocity = velocity;
   
% 	return

	end
	
end