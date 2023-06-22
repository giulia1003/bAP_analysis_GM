function [bAP_map_info_type1_all,bAP_map_info_type1_500,bAP_map_info_type2_500,bAP_map_info_type1_100,bAP_map_info_type2_100,bAP_map_info_type2_50,bAP_map_info_type1_50,bAP_map_info_type2_20,bAP_map_info_type1_20,...
	bAP_map_info_type2_5,bAP_map_info_type1_5] = get_voltage_maps_isi(analyse_nrn,spike_nums_5type1,spike_nums_5type2,spike_nums_500type1,spike_nums_500type2,spike_nums_100type1,spike_nums_100type2,...
	spike_nums_20type1,spike_nums_20type2,spike_nums_50type1,spike_nums_50type2,gwfparams,w_drive,rec_date,info, plotParams, NP_probe,sp,sp_threshold,drift_nrn,batch_id,templateYpos,hippocampal_units,analyse_bursts,take_burst_type, use_fixed_n_spikes)
wfMean_mean_shuffle = [];

bAP_map_info_type1_500 = []; bAP_map_info_type2_500 = []; bAP_map_info_type1_100 = []; bAP_map_info_type2_100 = [];
bAP_map_info_type2_50 = []; bAP_map_info_type1_50 = []; bAP_map_info_type2_20 = []; bAP_map_info_type1_20 = [];
bAP_map_info_type1_5 = []; bAP_map_info_type2_5 = [];  bAP_map_info_type1_all = [];

%% taking single spike voltage maps of bAPs using spikes across all ISI

if info.single_spike_maps && info.take_all_spikes  % save voltage maps for neurons using single spikes across different ISI

gwfparams.use_nWfs = sp_threshold;  % make voltage maps with 2000 spikes 

		info.analyse_whisking = false;
		high_network_analysis = false;
        info.analyse_locomotion = false;
        info.analyse_running = false;
		isi = 'all';
		bAP_map_info = [];
        nrn = analyse_nrn;   
        isi_type = [];
        
		% no distinction between running and stationary spikes
		
       for ind_spike_nrn = 1:size(analyse_nrn,2)

        isi_type(ind_spike_nrn).nrn = analyse_nrn(ind_spike_nrn);
% 	neuron_list = extractfield(isi_type,'nrn');
%        spike_times_all = extractfield(
		spike_times_all = []; spike_times_long_isi = []; spike_times_short_isi = [];
 
		% get spikes for every neuron (long ISI only)
		spike_times_long_isi = cat(1, spike_times_long_isi, spike_nums_500type1(ind_spike_nrn).spiketimes, spike_nums_500type2(ind_spike_nrn).spiketimes, spike_nums_100type1(ind_spike_nrn).spiketimes, spike_nums_100type2(ind_spike_nrn).spiketimes);
		spike_times_long_isi = spike_times_long_isi(randperm(gwfparams.use_nWfs/2));
		
		% get spikes for every neuron (short ISI only)
		spike_times_short_isi = cat(1, spike_times_short_isi,spike_nums_20type1(ind_spike_nrn).spiketimes, spike_nums_20type2(ind_spike_nrn).spiketimes, spike_nums_5type1(ind_spike_nrn).spiketimes, spike_nums_5type2(ind_spike_nrn).spiketimes);
		spike_times_short_isi = spike_times_short_isi(randperm(gwfparams.use_nWfs/2));
		
		% make vector with both long and short ISI spikes
		spike_times_all = cat(1, spike_times_all, spike_times_short_isi, spike_times_long_isi);
	
        isi_type(ind_spike_nrn).spiketimes = spike_times_all;
%         isi_type(ind_spike_nrn).isi = diff(spike_times_all)./30000; % in ms

       end

        plotParams.plot_fig = false;



for ind_analyse_nrn = 1:size(analyse_nrn,2)  
	
nrn = analyse_nrn(ind_analyse_nrn);

[bAP_map_info,plot_splineField_shift] = get_voltage_map(nrn,analyse_nrn,ind_analyse_nrn, info, plotParams, gwfparams, NP_probe,...
bAP_map_info,w_drive,rec_date,isi_type,isi,sp,drift_nrn,batch_id,templateYpos,hippocampal_units,analyse_bursts,take_burst_type, high_network_analysis);
pause(2)
close all

end


bAP_map_info_type1_all = bAP_map_info;


end





%% looking at bAPs across different ISI, or during bursts

if ~(info.take_all_spikes)
       

if analyse_bursts   % concentrate on bursts
	
	if ~(info.single_spike_maps)
% 	if analyse_locomotion
		info.analyse_whisking = false;
		high_network_analysis = false;
		
		conditions_bursts = [1:2];
		
		for j = 1:length(conditions_bursts)
		if j == 1
		bAP_map_info = [];
        nrn = analyse_nrn;

        isi_type = spike_nums_5type1;
        isi = 5;

		plotParams.plot_fig = true;
		info.analyse_running = 'stationary';  % stationary running
		
		elseif j == 2
		bAP_map_info = [];
        nrn = analyse_nrn;

       isi_type = spike_nums_5type2;
        isi = 5;

		plotParams.plot_fig = true;
		info.analyse_running = 'running';  % stationary running		
		
		end
		
		
for ind_analyse_nrn = 1:size(analyse_nrn,2)  

	if size(analyse_nrn,2) == 1
		nrn = analyse_nrn;
	else
nrn = analyse_nrn(ind_analyse_nrn);
	end
% spike_n = [spike_nums_5s(ind_analyse_nrn),spike_nums_5r(ind_analyse_nrn)];
% 
% spike_n = [spike_n(:).n_spikes];  min_spike = min(spike_n(spike_n >sp_threshold));  % take minimum number of spikes at least > 100 spikes

min_spike = sp_threshold;
gwfparams.use_nWfs = min_spike;

[bAP_map_info,plot_splineField_shift] = get_voltage_map(nrn,analyse_nrn,ind_analyse_nrn, info, plotParams, gwfparams, NP_probe,...
bAP_map_info,w_drive,rec_date,isi_type,isi,sp,drift_nrn,batch_id,templateYpos,hippocampal_units,analyse_bursts,take_burst_type, high_network_analysis);
 
close all

end

if j == 1
bAP_map_info_type1_5 = bAP_map_info;
elseif j == 2
bAP_map_info_type2_5 = bAP_map_info;
end

		end


		
elseif single_spike_maps
	
		% getting single maps for bursts
		
info.analyse_whisking = false;
info.analyse_locomotion = false;
info.analyse_running = false;
high_network_analysis = false; 
isi = 5;
bAP_map_info = [];
nrn = analyse_nrn;   
isi_type = [];	
  

for ind_spike_nrn = 1:size(analyse_nrn,2)
spike_times_short_isi = [];

isi_type(ind_spike_nrn).nrn = analyse_nrn(ind_spike_nrn);
spike_times_short_isi = cat(1, spike_times_short_isi, spike_nums_5type1(ind_spike_nrn).spiketimes, spike_nums_5type2(ind_spike_nrn).spiketimes);

if length(spike_times_short_isi) >= gwfparams.use_nWfs
spike_times_short_isi = spike_times_short_isi(randperm(gwfparams.use_nWfs));

else
	gwfparams.use_nWfs = length(spike_times_short_isi);
end

isi_type(ind_spike_nrn).spiketimes = spike_times_short_isi;

plotParams.plot_fig = false;
end

for ind_analyse_nrn = 1:size(analyse_nrn,2)  
	
nrn = analyse_nrn(ind_analyse_nrn);

[bAP_map_info,plot_splineField_shift] = get_voltage_map(nrn,analyse_nrn,ind_analyse_nrn, info, plotParams, gwfparams, NP_probe,...
bAP_map_info,w_drive,rec_date,isi_type,isi,sp,drift_nrn,batch_id,templateYpos,hippocampal_units,analyse_bursts,take_burst_type, high_network_analysis);
pause(2)
close all

end


bAP_map_info_type1_all = bAP_map_info;


end		
		
	
else % consider spikes during all ISI (not just bursts)
	
	
	
if info.analyse_locomotion

info.analyse_whisking = [];
high_network_analysis = false;

conditions = [1:10];

for i = 1:length(conditions)

if i == 1
% %  %% ISI 500 ms
bAP_map_info = [];

isi_type = spike_nums_500type1;
isi = 500;
plotParams.plot_fig = true;
info.analyse_running = 'stationary';  % stationary running

elseif i == 2
bAP_map_info = [];
isi_type = spike_nums_500type2;
isi = 500;
plotParams.plot_fig = true;
info.analyse_running = 'running'; % stationary running

elseif i == 3

bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_100type1;
isi = 100;
plotParams.plot_fig = true;
info.analyse_running = 'stationary';  % stationary running

elseif i == 4

bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_100type2;
isi = 100;
plotParams.plot_fig = true;
info.analyse_running = 'running';  % stationary running

elseif i == 5

bAP_map_info = [];
isi_type = spike_nums_50type1;
isi = 50;
plotParams.plot_fig = true;
info.analyse_running = 'stationary';  % stationary running

elseif i == 6

bAP_map_info = [];
isi_type = spike_nums_50type2;
isi = 50;
plotParams.plot_fig = true;
info.analyse_running = 'running';  % stationary running

elseif i == 7
	
bAP_map_info = [];
isi_type = spike_nums_20type1;
isi = 30;

plotParams.plot_fig = true;
info.analyse_running = 'stationary';   % stationary running

elseif i == 8 

bAP_map_info = [];
isi_type = spike_nums_20type2;
isi = 30;
plotParams.plot_fig = true;
info.analyse_running = 'running';  % stationary running

elseif i == 9
	
bAP_map_info = [];
isi_type = spike_nums_5type1;
isi = 5;
plotParams.plot_fig = true;
info.analyse_running = 'stationary';  % stationary running

elseif i == 10
	
bAP_map_info = [];
isi_type = spike_nums_5type2;
isi = 5;
plotParams.plot_fig = true;
info.analyse_running = 'running';  % stationary running

end

for ind_analyse_nrn = 1:size(analyse_nrn,2)  

nrn = analyse_nrn(ind_analyse_nrn);

spike_n = [spike_nums_500type1(ind_analyse_nrn),spike_nums_500type2(ind_analyse_nrn),spike_nums_100type1(ind_analyse_nrn),spike_nums_100type2(ind_analyse_nrn),...
	spike_nums_20type1(ind_analyse_nrn),spike_nums_20type2(ind_analyse_nrn),spike_nums_50type1(ind_analyse_nrn),spike_nums_50type2(ind_analyse_nrn),spike_nums_5type1(ind_analyse_nrn),spike_nums_5type2(ind_analyse_nrn)];

if use_fixed_n_spikes
min_spike =  sp_threshold;
else
spike_n = [spike_n(:).n_spikes];  min_spike = min(spike_n(spike_n >sp_threshold));  % take minimum number of spikes at least > 100 spikes
end

gwfparams.use_nWfs = min_spike;
try
[bAP_map_info,plot_splineField_shift] = get_voltage_map(nrn,analyse_nrn,ind_analyse_nrn, info, plotParams, gwfparams, NP_probe,...
bAP_map_info,w_drive,rec_date,isi_type,isi,sp,drift_nrn,batch_id,templateYpos,hippocampal_units,analyse_bursts,take_burst_type, high_network_analysis);
catch ME
	rethrow(ME)
end
pause(2)
close all

end

if i == 1
bAP_map_info_type1_500 = bAP_map_info;
elseif i == 2
bAP_map_info_type2_500 = bAP_map_info;
elseif i == 3
bAP_map_info_type1_100 = bAP_map_info;
elseif i == 4
bAP_map_info_type2_100 = bAP_map_info;
elseif i == 5
bAP_map_info_type1_50 = bAP_map_info;
elseif i == 6
bAP_map_info_type2_50 = bAP_map_info;
elseif i == 7
bAP_map_info_type1_20 = bAP_map_info;
elseif i == 8
bAP_map_info_type2_20 = bAP_map_info;
elseif i == 9
bAP_map_info_type1_5 = bAP_map_info;
elseif i == 10
bAP_map_info_type2_5 = bAP_map_info;
end

end

elseif info.analyse_whisk_behaviour

analyse_running = [];
high_network_analysis = false;

conditions = [1:10];

for i = 1:length(conditions)

if i == 1
% %  %% ISI 500 ms
bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_500type1;
isi = 500;

plotParams.plot_fig = true;
info.analyse_whisking = true; 

elseif i == 2
bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_500type2;
isi = 500;

plotParams.plot_fig = true;
info.analyse_whisking = false; 

elseif i == 3

bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_100type1;
isi = 100;

plotParams.plot_fig = true;
info.analyse_whisking = true; 

elseif i == 4

bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_100type2;
isi = 100;

plotParams.plot_fig = true;
info.analyse_whisking = false; 

elseif i == 5

bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_50type1;
isi = 50;

plotParams.plot_fig = true;
info.analyse_whisking = true; 

elseif i == 6

bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_50type2;
isi = 50;

plotParams.plot_fig = true;
info.analyse_whisking = false; 

elseif i == 7
	
bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_20type1;
isi = 30;

plotParams.plot_fig = true;
info.analyse_whisking = true; 

elseif i == 8 

bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_20type2;
isi = 30;

plotParams.plot_fig = true;
info.analyse_whisking = false; 

elseif i == 9
	
bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_5type1;
isi = 5;

plotParams.plot_fig = true;
info.analyse_whisking = true; 

elseif i == 10
	
bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_5type2;
isi = 5;

plotParams.plot_fig = true;
info.analyse_whisking = false; 

end

for ind_analyse_nrn = 1:size(analyse_nrn,2)  
	
nrn = analyse_nrn(ind_analyse_nrn);

spike_n = [spike_nums_500type1(ind_analyse_nrn),spike_nums_500type2(ind_analyse_nrn),spike_nums_100type1(ind_analyse_nrn),spike_nums_100type2(ind_analyse_nrn),...
	spike_nums_20type1(ind_analyse_nrn),spike_nums_20type2(ind_analyse_nrn),spike_nums_50type1(ind_analyse_nrn),spike_nums_50type2(ind_analyse_nrn),spike_nums_5type1(ind_analyse_nrn),spike_nums_5type2(ind_analyse_nrn)];

if use_fixed_n_spikes
	min_spike =  sp_threshold;
	
else
spike_n = [spike_n(:).n_spikes];  min_spike = min(spike_n(spike_n >sp_threshold));  % take minimum number of spikes at least > 100 spikes
end

gwfparams.use_nWfs = min_spike;

[bAP_map_info,plot_splineField_shift] = get_voltage_map(nrn,analyse_nrn,ind_analyse_nrn, info, plotParams, gwfparams, NP_probe,...
bAP_map_info,w_drive,rec_date,isi_type,isi,sp,drift_nrn,batch_id,templateYpos,hippocampal_units,analyse_bursts,take_burst_type, high_network_analysis);

close all

end

if i == 1
bAP_map_info_type1_500 = bAP_map_info;
elseif i == 2
bAP_map_info_type2_500 = bAP_map_info;
elseif i == 3
bAP_map_info_type1_100 = bAP_map_info;
elseif i == 4
bAP_map_info_type2_100 = bAP_map_info;
elseif i == 5
bAP_map_info_type1_50 = bAP_map_info;
elseif i == 6
bAP_map_info_type2_50 = bAP_map_info;
elseif i == 7
bAP_map_info_type1_20 = bAP_map_info;
elseif i == 8
bAP_map_info_type2_20 = bAP_map_info;
elseif i == 9
bAP_map_info_type1_5 = bAP_map_info;
elseif i == 10
bAP_map_info_type2_5 = bAP_map_info;
end

end
	
	
		
elseif info.analyse_network_activity

info.analyse_running = [];
info.analyse_whisk_behaviour = false;
info.analyse_whisking = [];


conditions = [1:10];

for i = 1:length(conditions)

if i == 1
% %  %% ISI 500 ms
bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_500type1;
isi = 500;

plotParams.plot_fig = true;
high_network_analysis = false; 

elseif i == 2
bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_500type2;
isi = 500;

plotParams.plot_fig = true;
high_network_analysis = true; 

elseif i == 3

bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_100type1;
isi = 100;

plotParams.plot_fig = true;
high_network_analysis = false; 

elseif i == 4

bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_100type2;
isi = 100;

plotParams.plot_fig = true;
high_network_analysis = true; 

elseif i == 5

bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_50type1;
isi = 50;

plotParams.plot_fig = true;
high_network_analysis = false; 

elseif i == 6

bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_50type2;
isi = 50;

plotParams.plot_fig = true;
high_network_analysis = true; 

elseif i == 7
	
bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_20type1;
isi = 30;

plotParams.plot_fig = true;
high_network_analysis = false; 

elseif i == 8 

bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_20type2;
isi = 30;

plotParams.plot_fig = true;
high_network_analysis = true; 

elseif i == 9
	
bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_5type1;
isi = 5;

plotParams.plot_fig = true;
high_network_analysis = false; 

elseif i == 10
	
bAP_map_info = [];
nrn = analyse_nrn;

isi_type = spike_nums_5type2;
isi = 5;

plotParams.plot_fig = true;
high_network_analysis = true; 

end


for ind_analyse_nrn = 1:size(analyse_nrn,2)  
	
nrn = analyse_nrn(ind_analyse_nrn);

spike_n = [spike_nums_500type1(ind_analyse_nrn),spike_nums_500type2(ind_analyse_nrn),spike_nums_100type1(ind_analyse_nrn),spike_nums_100type2(ind_analyse_nrn),...
	spike_nums_20type1(ind_analyse_nrn),spike_nums_20type2(ind_analyse_nrn),spike_nums_50type1(ind_analyse_nrn),spike_nums_50type2(ind_analyse_nrn),spike_nums_5type1(ind_analyse_nrn),spike_nums_5type2(ind_analyse_nrn)];

if use_fixed_n_spikes
	min_spike =  sp_threshold;
	
else
spike_n = [spike_n(:).n_spikes];  min_spike = min(spike_n(spike_n >sp_threshold));  % take minimum number of spikes at least > 100 spikes
end
gwfparams.use_nWfs = min_spike;
 
[bAP_map_info,plot_splineField_shift] = get_voltage_map(nrn,analyse_nrn,ind_analyse_nrn, info, plotParams, gwfparams, NP_probe,...
bAP_map_info,w_drive,rec_date,isi_type,isi,sp,drift_nrn,batch_id,templateYpos,hippocampal_units,analyse_bursts,take_burst_type, high_network_analysis);

close all

end
if i == 1
bAP_map_info_type1_500 = bAP_map_info;
elseif i == 2
bAP_map_info_type2_500 = bAP_map_info;
elseif i == 3
bAP_map_info_type1_100 = bAP_map_info;
elseif i == 4
bAP_map_info_type2_100 = bAP_map_info;
elseif i == 5
bAP_map_info_type1_50 = bAP_map_info;
elseif i == 6
bAP_map_info_type2_50 = bAP_map_info;
elseif i == 7
bAP_map_info_type1_20 = bAP_map_info;
elseif i == 8
bAP_map_info_type2_20 = bAP_map_info;
elseif i == 9
bAP_map_info_type1_5 = bAP_map_info;
elseif i == 10
bAP_map_info_type2_5 = bAP_map_info;
end

	
end 




end





end

end

























