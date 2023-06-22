function [sesh_ON_periods,all_sesh_ON_periods,sesh_ON_starts,templates_table,gwfparams,running_running,running_runningTimes,stat_runningTimes,time_motion,whisking_times,not_whisking_times,whisking_vect] = load_bhv_params(info,gwfparams,LED_trig_data,...
	procDataDir_General,data_path,sp,running_threshold, stationary_threshold,bAP_nrn_ind,analyse_multiple_sessions,session)

sesh_ON_starts = []; all_sesh_ON_periods = []; sesh_ON_periods = []; not_whisking_times = []; whisking_times = [];
running_running = []; running_runningTimes = []; templates_table = []; 
stat_runningTimes = []; licking_lickingTimes = []; licking = []; time_motion = []; whisking_vect = []; 


[interpolated_firing_rate] = get_network_activity_interp(sp);

if ~isempty(LED_trig_data)
	sesh_ON_periods = [];
if analyse_multiple_sessions
% sesh_ON_periods = LED_trig_data.session_g6.raw.LED_on_periods;             % LED session to be analysed
%  all_sesh_ON_periods = LED_trig_data.LED_on_periods_full.cat;                % all LED session, concatenated files
all_sesh_ON_periods = LED_trig_data.all.offset;
else 
all_sesh_ON_periods = [LED_trig_data.(session).raw];
end

if info.analyse_LED_starts 
 sesh_ON_starts = LED_trig_data.session_g6.raw.LED_trig_starts;            % set this to look at LED start triggers
end
end

if exist([procDataDir_General, '/touch_spikes_table.mat'], 'file')
	load([procDataDir_General, '/touch_spikes_table.mat'])
end
	
if info.analyse_locomotion && ~(info.analyse_whisk_behaviour)
[running_running,running_runningTimes,stat_runningTimes, time_motion] = load_running_gm(data_path,sp,running_threshold,stationary_threshold);
disp('Analysing spike-subset for locomotion')
end

if info.analyse_whisk_behaviour && ~(info.analyse_locomotion)
[whisking_times,not_whisking_times] = load_whisking_gm(data_path);
disp('whisking times loaded to the workspace')
end

if info.analyse_whisk_behaviour && info.analyse_locomotion
[running_running,running_runningTimes,stat_runningTimes, time_motion] = load_running_gm(data_path,sp,running_threshold,stationary_threshold);
[whisking_times,not_whisking_times,whisking_vect] = load_whisking_gm(data_path);
disp('Locomotion and whisking times loaded to the workspace')
end



end
