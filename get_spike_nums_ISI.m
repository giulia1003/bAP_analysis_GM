function [spike_nums_5type1,spike_nums_5type2,spike_nums_20type1,spike_nums_20type2,spike_nums_100type2,spike_nums_100type1,spike_nums_50type2,spike_nums_50type1,spike_nums_500type2,spike_nums_500type1,gwfparams,spikes_bhv_unit,spikes_bhv_b_unit,time_motion,LED_data,...
	T,bAP_map_info,ISI_range,take_burst_type] = get_spike_nums_ISI(info,analyse_nrn,gwfparams,...
	running_runningTimes,stat_runningTimes,sesh_ON_periods,sp,LED_trig_data,...
	all_sesh_ON_periods,plotParams,LED_data_all,time_motion,not_whisking_times,whisking_times,take_burst_type,high_net_activity,low_net_activity)

 spike_nums_5type1 = []; spike_nums_5type2 = []; spike_nums_50type2 = []; spike_nums_50type1 = []; spike_nums_500type2 = []; spike_nums_500type1 = [];
spike_nums_20type1 = []; spike_nums_20type2 = []; spike_nums_100type2 = []; spike_nums_100type1 = []; T = []; bAP_map_info = []; spike_nums = [];
 

conditions = [1:10];

if ~(info.analyse_network_activity)
high_net_activity = [];
low_net_activity = [];

end
info.analyse_licking = false;

% to focus on bursts only (first of a doublet or triplet, second spikes etc)
if info.analyse_bursts
	if info.analyse_locomotion
		info.analyse_whisking = false;
		conditions_bursts = [1:2];
		high_network_analysis = false;
		
		for j = 1:length(conditions_bursts)
		if conditions_bursts(j) == 1
		info.analyse_running = 'stationary';
		prespike_restriction = [1:10]; %250:18000000	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
        postspike_restriction = prespike_restriction;	% same for the time of silence after the burst or after single spikes
		elseif conditions_bursts(j) == 2
		info.analyse_running = 'running';
		prespike_restriction = [1:10]; %250:18000000	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
        postspike_restriction = prespike_restriction;	% same for the time of silence after the burst or after single spikes
		end
		
[gwfparams,spikes_bhv_unit,spikes_bhv_b_unit,time_motion,LED_data,T,bAP_map_info,ISI_range,take_burst_type] = get_spike_num_nrn(info,analyse_nrn,gwfparams,prespike_restriction, postspike_restriction,...
	running_runningTimes,stat_runningTimes,sesh_ON_periods,sp,LED_trig_data,...
	all_sesh_ON_periods,plotParams,LED_data_all,T,time_motion,bAP_map_info,not_whisking_times,whisking_times,take_burst_type,high_network_analysis, high_net_activity, low_net_activity);
       
if conditions(j) == 1
	     spike_nums_5type1 = T;
        elseif conditions(j) == 2
	     spike_nums_5type2 = T;
		end
		
		end
	

	
	elseif info.analyse_network_activity
		
		info.analyse_whisking = false;
		info.analyse_locomotion= false;
		conditions_bursts = [1:2];

		
		for j = 1:length(conditions_bursts)
		if conditions_bursts(j) == 1
	high_network_analysis = false;
		prespike_restriction = [1:10]; %250:18000000	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
        postspike_restriction = prespike_restriction;	% same for the time of silence after the burst or after single spikes
		elseif conditions_bursts(j) == 2
			high_network_analysis = true;
		prespike_restriction = [1:10]; %250:18000000	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
        postspike_restriction = prespike_restriction;	% same for the time of silence after the burst or after single spikes
		end
		
[gwfparams,spikes_bhv_unit,spikes_bhv_b_unit,time_motion,LED_data,T,bAP_map_info,ISI_range,take_burst_type] = get_spike_num_nrn(info,analyse_nrn,gwfparams,prespike_restriction, postspike_restriction,...
	running_runningTimes,stat_runningTimes,sesh_ON_periods,sp,LED_trig_data,...
	all_sesh_ON_periods,plotParams,LED_data_all,T,time_motion,bAP_map_info,not_whisking_times,whisking_times,take_burst_type,high_network_analysis, high_net_activity, low_net_activity);
       
if conditions(j) == 1
	     spike_nums_5type1 = T;
        elseif conditions(j) == 2
	     spike_nums_5type2 = T;
end
		
		end
		
	end
	
else  % to look at spikes across all ISI (not just bursts)
		
if info.analyse_locomotion && ~(info.analyse_network_activity)

info.analyse_whisking = false; 
high_network_analysis = false;


for i = 1:length(conditions)
	
if conditions(i) == 1
info.analyse_running = 'stationary';         % 'running','stationary','no';
prespike_restriction = [250:18000000]; %250:18000000	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;	% same for the time of silence after the burst or after single spikes
elseif conditions(i) == 2
info.analyse_running = 'running';         % 'running','stationary','no';
prespike_restriction = [250:18000000];	
postspike_restriction = prespike_restriction;
elseif conditions(i) == 3
info.analyse_running = 'stationary';         % 'running','stationary','no';
info.analyse_licking = false;
prespike_restriction = [100:250];	% 100:250
postspike_restriction = prespike_restriction;
elseif conditions(i) == 4
info.analyse_running = 'running';         % 'running','stationary','no';
prespike_restriction = [100:250];	% 100:250
postspike_restriction = prespike_restriction;
elseif conditions(i) == 5
info.analyse_running = 'stationary';         % 'running','stationary','no';
prespike_restriction = [50:100];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
elseif conditions(i) == 6
info.analyse_running = 'running';         % 'running','stationary','no';
prespike_restriction = [50:100];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
elseif conditions(i) == 7   
info.analyse_running = 'stationary';         % 'running','stationary','no';
prespike_restriction = [10:30];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
elseif conditions(i) == 8    
info.analyse_running = 'running';         % 'running','stationary','no';
prespike_restriction = [10:30];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;

elseif conditions(i) == 9    
info.analyse_running = 'stationary';         % 'running','stationary','no';
prespike_restriction = [1:10];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
elseif conditions(i) == 10    
info.analyse_running = 'running';         % 'running','stationary','no';
prespike_restriction = [1:10];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
else 
end

[gwfparams,spikes_bhv_unit,spikes_bhv_b_unit,time_motion,LED_data,T,bAP_map_info,ISI_range,take_burst_type] = get_spike_num_nrn(info,analyse_nrn,gwfparams,prespike_restriction, postspike_restriction,...
	running_runningTimes,stat_runningTimes,sesh_ON_periods,sp,LED_trig_data,...
	all_sesh_ON_periods,plotParams,LED_data_all,T,time_motion,bAP_map_info,not_whisking_times,whisking_times,take_burst_type,high_network_analysis, high_net_activity, low_net_activity);

if conditions(i) == 1
	spike_nums_500type1 = T;
elseif conditions(i) == 2
	spike_nums_500type2 = T;
elseif conditions(i) == 3
	spike_nums_100type1 = T;
elseif conditions(i) == 4	
spike_nums_100type2 = T;
elseif conditions(i) == 5
spike_nums_50type1 = T;
elseif conditions(i) == 6
spike_nums_50type2 = T;
elseif conditions(i) == 7	
spike_nums_20type1 = T;
elseif conditions(i) == 8	
spike_nums_20type2 = T;
elseif conditions(i) == 9	
spike_nums_5type1 = T;
elseif conditions(i) == 10	
spike_nums_5type2 = T;
else
end
end


elseif info.analyse_whisk_behaviour && ~(info.analyse_network_activity)

info.analyse_running = false;
high_network_analysis = false;

for i = 1:length(conditions)
	
if conditions(i) == 1
info.analyse_whisking = true;      
prespike_restriction = [250:18000000]; %250:600	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;	% same for the time of silence after the burst or after single spikes
elseif conditions(i) == 2
info.analyse_whisking = false;     
prespike_restriction = [250:18000000];	
postspike_restriction = prespike_restriction;
elseif conditions(i) == 3
info.analyse_whisking = true;        
prespike_restriction = [100:250];	% 100:250
postspike_restriction = prespike_restriction;
elseif conditions(i) == 4
info.analyse_whisking = false;          
prespike_restriction = [100:250];	% 100:250
postspike_restriction = prespike_restriction;
elseif conditions(i) == 5
info.analyse_whisking = true;          
prespike_restriction = [50:100];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
elseif conditions(i) == 6
info.analyse_whisking = false;        
prespike_restriction = [50:100];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
elseif conditions(i) == 7   
info.analyse_whisking = true;          % 'running','stationary','no';
prespike_restriction = [10:50];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
elseif conditions(i) == 8    
info.analyse_whisking = false;        % 'running','stationary','no';
prespike_restriction = [10:50];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;

elseif conditions(i) == 9    
info.analyse_whisking = true;         % 'running','stationary','no';
prespike_restriction = [1:10];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
elseif conditions(i) == 10    
analyse_whisking = false;         % 'running','stationary','no';
prespike_restriction = [1:10];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
else 
end

[gwfparams,spikes_bhv_unit,spikes_bhv_b_unit,time_motion,LED_data,T,bAP_map_info,ISI_range,take_burst_type] = get_spike_num_nrn(info,analyse_nrn,gwfparams,prespike_restriction, postspike_restriction,...
	running_runningTimes,stat_runningTimes,sesh_ON_periods,sp,LED_trig_data,...
	all_sesh_ON_periods,plotParams,LED_data_all,T,time_motion,bAP_map_info,not_whisking_times,whisking_times,take_burst_type,high_network_analysis, high_net_activity, low_net_activity);

if conditions(i) == 1
spike_nums_500type1 = T;
elseif conditions(i) == 2
spike_nums_500type2 = T;
elseif conditions(i) == 3
spike_nums_100type1 = T;
elseif conditions(i) == 4	
spike_nums_100type2 = T;
elseif conditions(i) == 5
spike_nums_50type1 = T;
elseif conditions(i) == 6
spike_nums_50type2 = T;
elseif conditions(i) == 7	
spike_nums_20type1 = T;
elseif conditions(i) == 8	
spike_nums_20type2 = T;
elseif conditions(i) == 9	
spike_nums_5type1 = T;
elseif conditions(i) == 10	
spike_nums_5type2 = T;
else
end
end

elseif info.analyse_network_activity && ~(info.analyse_locomotion) && ~(info.analyse_whisk_behaviour)
	
info.analyse_whisking = false;  
info.analyse_running = false;

for i = 1:length(conditions)
	
if conditions(i) == 1
high_network_analysis = false;      
prespike_restriction = [250:18000000]; %250:600	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;	% same for the time of silence after the burst or after single spikes
elseif conditions(i) == 2
high_network_analysis = true;     
prespike_restriction = [250:18000000];	
postspike_restriction = prespike_restriction;
elseif conditions(i) == 3
high_network_analysis = false;        
prespike_restriction = [100:250];	% 100:250
postspike_restriction = prespike_restriction;
elseif conditions(i) == 4
high_network_analysis = true;          
prespike_restriction = [100:250];	% 100:250
postspike_restriction = prespike_restriction;
elseif conditions(i) == 5
high_network_analysis = false;          
prespike_restriction = [50:100];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
elseif conditions(i) == 6
high_network_analysis = true;        
prespike_restriction = [50:100];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
elseif conditions(i) == 7   
high_network_analysis = false;          % 'running','stationary','no';
prespike_restriction = [10:30];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
elseif conditions(i) == 8    
high_network_analysis = true;        % 'running','stationary','no';
prespike_restriction = [10:30];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;

elseif conditions(i) == 9    
high_network_analysis = false;         % 'running','stationary','no';
prespike_restriction = [1:10];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
elseif conditions(i) == 10    
high_network_analysis = true;         % 'running','stationary','no';
prespike_restriction = [1:10];	% specify the interval in ms of silence before the first spike of a burst or between single spikes (!)
postspike_restriction = prespike_restriction;
else 
end

[gwfparams,spikes_bhv_unit,spikes_bhv_b_unit,time_motion,LED_data,T,bAP_map_info,ISI_range,take_burst_type] = get_spike_num_nrn(info,analyse_nrn,gwfparams,prespike_restriction, postspike_restriction,...
	running_runningTimes,stat_runningTimes,sesh_ON_periods,sp,LED_trig_data,...
	all_sesh_ON_periods,plotParams,LED_data_all,T,time_motion,bAP_map_info,not_whisking_times,whisking_times,take_burst_type,high_network_analysis, high_net_activity, low_net_activity);

if conditions(i) == 1
spike_nums_500type1 = T;
elseif conditions(i) == 2
spike_nums_500type2 = T;
elseif conditions(i) == 3
spike_nums_100type1 = T;
elseif conditions(i) == 4	
spike_nums_100type2 = T;
elseif conditions(i) == 5
spike_nums_50type1 = T;
elseif conditions(i) == 6
spike_nums_50type2 = T;
elseif conditions(i) == 7	
spike_nums_20type1 = T;
elseif conditions(i) == 8	
spike_nums_20type2 = T;
elseif conditions(i) == 9	
spike_nums_5type1 = T;
elseif conditions(i) == 10	
spike_nums_5type2 = T;
else
end
end


	
	

end
end