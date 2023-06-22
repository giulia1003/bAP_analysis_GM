function [spike_nums_5type1,spike_nums_5type2,spike_nums_20type1,spike_nums_20type2,spike_nums_100type2,spike_nums_100type1,spike_nums_50type2,spike_nums_50type1,spike_nums_500type2,spike_nums_500type1,gwfparams,analyse_bursts,take_burst_type]  = get_spiketimes(info,load_spikes,infodir,analyse_nrn,gwfparams,...
	running_runningTimes,stat_runningTimes,sesh_ON_periods,sp,LED_trig_data,...
	all_sesh_ON_periods,plotParams,LED_data_all,time_motion,not_whisking_times,whisking_times)
% if analyse_LED_spikes = false and baseline_rec = true
% only spikes fired during initial baseline recordings will be considered
% as true LED OFF spikes (this is necessary for long-lasting opsins like eOPN3)
% if analyse_LED_spikes = false and baseline_rec = false
% no restrictions, all spikes during led off periods will be used 
% if analyse_LED_spikes = true and baseline_rec = false
% same but with led on conditions

if load_spikes
    
 analyse_bursts = false;
 take_burst_type = [];
 
file_name = 'spike_nums.mat'; file_name_led = 'spike_nums_LED.mat';
file_name_hippo = 'spike_nums_hippo.mat'; file_name_whisk = 'spike_nums_whisking.mat';

if ~(info.analyse_LED_spikes) && ~(info.analyse_whisk_behaviour) && ~(info.analyse_network_activity) &&  exist([infodir,'\',file_name])
	load(strcat(infodir,'\',file_name)) 
	
elseif info.analyse_LED_spikes && ~(info.analyse_whisk_behaviour) && ~(info.analyse_network_activity) && exist([infodir,'\',file_name_led])
	load(strcat(infodir,'\',file_name_led)) 
	
elseif ~(info.analyse_LED_spikes) && exist([infodir,'\',file_name_hippo]) && ~(info.analyse_network_activity) && hippocampal_units 
	load(strcat(infodir,'\',file_name_hippo)) 

elseif info.analyse_whisk_behaviour && ~(info.analyse_network_activity) && exist([infodir,'\',file_name_whisk])
	load(strcat(infodir,'\',file_name_whisk))
else
end



else


hippocampal_units = false;

prompt1 = 'Burst analysis: true/false '
info.analyse_bursts = input(prompt1);

if info.analyse_bursts
prompt2 = 'take isolated bursts? true/false '
info.isolated_bursts = input(prompt2);
end

% get spikes during high network and low network states
if info.analyse_network_activity

[firing_rate_network, high_net_activity,low_net_activity,thresh_firing] = get_network_activity(sp, analyse_nrn);

else
high_net_activity = []; low_net_activity = [];
end

if info.analyse_bursts
prompt2 = 'Choose burst type to analyse: "firsts" / "seconds" / "thirds" / "all" '   
% firsts: first spike of a burst
% seconds: second spikes 
% others: all spikes within a burst apart from the first one
take_burst_type = input(prompt2);

[spike_nums_5type1,spike_nums_5type2,spike_nums_20type1,spike_nums_20type2,spike_nums_100type2,spike_nums_100type1,spike_nums_50type2,spike_nums_50type1,spike_nums_500type2,spike_nums_500type1,gwfparams,spikes_bhv_unit,spikes_bhv_b_unit,time_motion,LED_data,...
	T,bAP_map_info,ISI_range,take_burst_type] = get_spike_nums_ISI(info,analyse_nrn,gwfparams,...
	running_runningTimes,stat_runningTimes,sesh_ON_periods,sp,LED_trig_data,...
	all_sesh_ON_periods,plotParams,LED_data_all,time_motion,not_whisking_times,whisking_times,take_burst_type,high_net_activity,low_net_activity);


else % not looking at bursts
	
take_burst_type = []; 


if info.analyse_LED_spikes && info.baseline_rec

[all_sesh_ON_periods] = extend_LED_on_periods(LED_trig_data,all_sesh_ON_periods);

end


[spike_nums_5type1,spike_nums_5type2,spike_nums_20type1,spike_nums_20type2,spike_nums_100type2,spike_nums_100type1,spike_nums_50type2,spike_nums_50type1,spike_nums_500type2,spike_nums_500type1,gwfparams,spikes_bhv_unit,spikes_bhv_b_unit,time_motion,LED_data,...
	T,bAP_map_info,ISI_range,take_burst_type] = get_spike_nums_ISI(info,analyse_nrn,gwfparams,...
	running_runningTimes,stat_runningTimes,sesh_ON_periods,sp,LED_trig_data,...
	all_sesh_ON_periods,plotParams,LED_data_all,time_motion,not_whisking_times,whisking_times,take_burst_type, high_net_activity,low_net_activity);

if ~(info.analyse_LED_spikes) && info.baseline_rec
	% takes only spikes fired before the LED was turned on for the first time
[spike_nums_5type1,spike_nums_5type2,spike_nums_20type1,spike_nums_20type2,spike_nums_50type2,spike_nums_50type1,spike_nums_100type2,spike_nums_100type1,spike_nums_500type2,spike_nums_500type1] = keep_only_baseline_spikes(spike_nums_5type1,spike_nums_5type2,spike_nums_20type1,spike_nums_20type2,spike_nums_50type2,spike_nums_50type1,...
	spike_nums_100type2,spike_nums_100type1,spike_nums_500type2,spike_nums_500type1,LED_trig_data);
end


end


cd(infodir)
if info.analyse_LED_spikes
save('spike_nums_LED.mat','spike_nums_5type1','spike_nums_5type2','spike_nums_20type1','spike_nums_20type2','spike_nums_100type1','spike_nums_100type2','spike_nums_50type1','spike_nums_50type2','spike_nums_500type1','spike_nums_500type2');
elseif ~(info.analyse_LED_spikes) && ~(info.analyse_whisk_behaviour)
save('spike_nums.mat','spike_nums_5type1','spike_nums_5type2','spike_nums_20type1','spike_nums_20type2','spike_nums_100type1','spike_nums_100type2','spike_nums_50type1','spike_nums_50type2','spike_nums_500type1','spike_nums_500type2');
elseif ~(info.analyse_LED_spikes) && info.analyse_whisk_behaviour
save('spike_nums_whisking.mat','spike_nums_5type1','spike_nums_5type2','spike_nums_20type1','spike_nums_20type2','spike_nums_100type1','spike_nums_100type2','spike_nums_50type1','spike_nums_50type2','spike_nums_500type1','spike_nums_500type2');
elseif ~(info.analyse_LED_spikes) && info.analyse_network_activity
save('spike_nums_network.mat','spike_nums_5type1','spike_nums_5type2','spike_nums_20type1','spike_nums_20type2','spike_nums_100type1','spike_nums_100type2','spike_nums_50type1','spike_nums_50type2','spike_nums_500type1','spike_nums_500type2');

end

analyse_bursts = info.analyse_bursts;
end


end % end function

