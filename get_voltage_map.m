
function [bAP_map_info,plot_splineField_shift] = get_voltage_map(nrn,analyse_nrn,ind_analyse_nrn, info, plotParams, gwfparams, NP_probe,...
bAP_map_info,w_drive,rec_date,isi_type,isi,sp,drift_nrn,batch_id,templateYpos,hippocampal_units,analyse_bursts,take_burst_type, high_network_analysis)

 soma_depth_sp_Field_tab = []; 
 single_spike_map = []; 
spike_time_tab = [];
shuffled_map_tab = [];  
% 	idx = find(sp.cids == nrn);
	nrn_ID = double(sp.cids(nrn));
%     nrn_ID = double(sp.cids(idx));
	
	neuron_list = extractfield(isi_type,'nrn');
	
	index_nrn = find(nrn == neuron_list);
	

    gwfparams.spikeTimes = isi_type(index_nrn).spiketimes;
	gwfparams.spikeClusters = ones(1,size(gwfparams.spikeTimes,1));
	

warning('off');

gwfparams.nWf = gwfparams.use_nWfs;

	
if ~isempty(gwfparams.spikeTimes) & numel(gwfparams.spikeTimes) >= gwfparams.use_nWfs  
	try
% 		nrn = nrn(ind_analyse_nrn);
		
	[wfMean_mean_shift,voltage_map_shift,...
	spiketimes,wf_shift,single_spike_map,spike_time_tab,...
		shuffled_map_tab] = get_WFs_NO_CSD(sp,nrn_ID, nrn, info, plotParams, gwfparams, NP_probe,drift_nrn,batch_id,templateYpos,hippocampal_units);
		
	catch ME
% 	
      			 rethrow(ME)
 		fprintf(2, 'Error while extracting spikes.\n')
		voltage_map = []; bAP_extent_shift = []; voltage_map_shift = []; single_spike_map = []; spike_time_tab = []; shuffled_map_tab = []; 
				soma_time_sp_Field = 0; wf_shift = []; wf = []; amplitudes_vect = []; soma_depth_sp_Field_tab = [soma_depth_sp_Field_tab; soma_time_sp_Field]; 
%  		    		 continue
	end
else 
		voltage_map = []; voltage_map_shift = [];
		soma_time_sp_Field = 0; bAP_times_and_channels = []; 
		soma_depth_sp_Field_tab = [soma_depth_sp_Field_tab; soma_time_sp_Field];
	
end
	
%% Get spline and bAP extent, plot bAP 

if isempty(gwfparams.spikeTimes)
	bAP_times_and_channels = []; bAP_extent = []; soma_depth = []; soma_time = []; min_map = []; plot_splineField = [];
	wf_shift = []; wf = []; single_spike_map = []; bAP_times_and_channels = [];
end

if ~isempty(voltage_map_shift) &&  plotParams.plot_fig
	
try


[bAP_times_and_channels,plot_splineField_shift,plot_CSD,bAP_extent_shift,...
    signal_axon_only,soma_depth,soma_time,min_map,velocity] = get_voltage_map_shift(voltage_map_shift,nrn_ID, nrn, info, plotParams, gwfparams, NP_probe,spiketimes,isi,w_drive,rec_date,hippocampal_units,analyse_bursts,take_burst_type, high_network_analysis);

[amplitudes_vect,bAP_times_and_channels] = get_amplitudes(voltage_map_shift, bAP_times_and_channels, nrn, plotParams,NP_probe); 

catch ME
     		rethrow(ME)
         
bAP_times_and_channels = []; bAP_channels = []; bAP_extent = [];bAP_extent_shift = []; soma_depth = []; soma_time = []; min_map = []; plot_splineField = []; plot_splineField_shift = []; plot_CSD = [];
soma_depth_shift = []; soma_time = []; min_map_shift = []; min_value = []; min_value_shift = []; velocity =[];signal_axon_only = [];
amplitudes_vect = []; 
end


end

	
if ~(plotParams.plot_fig)
plot_splineField_shift = []; ch_p = [];
end

%%  Look at bAP extent for each unit 
% using 2 metrics : amplitude ~160um above the soma
% spline extent from soma to channels above
if ~(exist('plot_splineField_shift'))
	plot_splineField_shift = [];
end

if ~isempty(plot_splineField_shift) 
get_min_traces = true;

min_spline_shift = min(voltage_map_shift(:));
[soma_depth,soma_time] = find(voltage_map_shift == min_spline_shift);
select_height = soma_depth;
[soma_trace, above_trace,high_trace, correlation_matrix] = get_amplitude_above_soma(info,sp,select_height,wfMean_mean_shift,isi,w_drive,rec_date,nrn,NP_probe);
wf_shift.min_value_soma = min(soma_trace);

% convert to uV (NP 1.0)
% Gain  = 500;
% UVPerBit = 1.2/2^10/Gain*1e6;
% soma_trace = soma_trace.*UVPerBit;

elseif isempty(plot_splineField_shift) 
	bAP_times_and_channels = []; soma_depth = []; soma_time = []; min_map = []; plot_splineField = []; wf = []; wf_shift.min_value_soma = []; velocity =[]; signal_axon_only = [];
	voltage_map = []; spline = []; gwfparams.spikeTimes = []; firing_freq_change = []; above_trace = []; plot_CSD = []; amplitudes_vect = [];  correlation_matrix = [];

	spikes_bhv_unit = [];  bAP_extent= []; bAP_extent_shift = []; amps_vect_sp_Field = []; soma_trace = []; high_trace = []; plot_splineField_shift = []; 
end

[bAP_map_info] = get_bap_map3(info,analyse_nrn,min_map,bAP_times_and_channels,voltage_map_shift,...
	soma_depth,bAP_extent_shift,signal_axon_only,soma_trace,above_trace,ind_analyse_nrn,bAP_map_info,high_trace,wf_shift,velocity,single_spike_map,spike_time_tab,shuffled_map_tab,amplitudes_vect, correlation_matrix);




end