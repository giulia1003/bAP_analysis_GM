function [bAP_times_and_channels,plot_splineField,plot_CSD,bAP_extent_shift,signal_axon_only,soma_depth,soma_time,min_map,velocity] = get_voltage_map_shift(voltage_map,nrn_ID, nrn, info, plotParams, gwfparams, NP_probe,spikeTimes,isi,w_drive,rec_date,hippocampal_units,analyse_bursts,take_burst_type, high_network_analysis)

% this function determines the extent of backpropagating action potentials from the soma to the apical dendrites
% using a linear shift test to compute where the minima of every channel above the soma
% is significantly above noise levels (doesn't fall within a null distribution of channel probabilities from voltage maps
% where spike times got shifted by hundreds of s)
% the outputs of the function consist of: an average voltage map centered around the spike times at the soma,
% its CSD, a vector of the real map's channel probabilities (ch_p) for every channel above soma, 
% the extent of the bAP in um, its conduction velocity (in cm/s)


% compute CSD 
plotParams.FieldSmooth_win = 5;     % window (number of depths) for gaussian smoothing of map before doing CSD
 gaussFilter = gausswin(plotParams.FieldSmooth_win);   
 voltage_map_smooth = conv2(gaussFilter,1, voltage_map, 'same');
 
csd_map = computeCSD(voltage_map_smooth);
 
addt = 120;  % pick how many ms after spike you want to include in the plot (120 = 4 ms)
addm = 40; % same thing but for channels above soma

%% get vector of shuffled spikes and shuffled voltage map

% for loop that makes 1) vector of spike times + each element of rand num vectors
% 2) a voltage map with shuffled spiketimes
% voltage_map_shuf_tab = [];voltage_map_shuf_tab.map = []; voltage_map_shuf_tab.csd = [];

tic

parpool;

p = gcp('nocreate');    addAttachedFiles(gcp,["get_WFs_NO_CSD_spline.m"]); 
addAttachedFiles(gcp,["getWaveForms_shuf.m"]);  addAttachedFiles(gcp,["computeSplineLFP_NoHistology.m"]);   

n = 300;

v_map = cell(1, n);   csd_map_shuf = cell(1, n);
spiketimes_shuf_all = cell(1, n);
shuffle_vec = [1:1:n];

% get vector of random nums 
shuffle_nums = round((rand(n,1)-0.5)*10000);  % 10000
shuffle_nums(shuffle_nums<0) = shuffle_nums(shuffle_nums<0) - 250;
shuffle_nums(shuffle_nums>=0) = shuffle_nums(shuffle_nums>=0) + 250;

% Load .dat and KiloSort/Phy output, needed to make shuffled maps
fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
gwfparams.mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'});
gwfparams.chMap = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab

spiketimesSub = spikeTimes(1:200);

parfor idx = 1 : numel(shuffle_vec)
    
   	spiketimes_shuf = spiketimesSub + shuffle_nums(idx);
    spiketimes_shuf = round(spiketimes_shuf);
     
	[voltage_map_shuf] = get_WFs_NO_CSD_spline(info, plotParams, gwfparams, NP_probe,spiketimes_shuf);
	       
 	 csd_map_shuf{idx} = computeCSD(voltage_map_shuf); % get CSD of shuffled map
	 
     v_map{idx} = voltage_map_shuf;
     
	
end
 delete(gcp('nocreate'));

toc
% get mean map of shuffled voltage maps
reshapedArrays = reshape(v_map, [1, 1, size(v_map,2)]);
mean_spline_shuf = cat(3, reshapedArrays{:});
shuf_spline = mean(mean_spline_shuf,3, 'omitnan');

% get mean of csd of shuffled maps
reshapedArrays = reshape(csd_map_shuf, [1, 1, size(csd_map_shuf,2)]);
mean_csd_shuf = cat(3, reshapedArrays{:});
shuf_csd = mean(mean_csd_shuf,3, 'omitnan');


%% non parametric test starts here
% done for voltage map and CSD

 run_min_shuf = [];
 
[nCh, nT, nSh] = size(shuf_spline); 
 min_value = min(voltage_map(:)); 
%  soma_depth = wf_shift.mean_soma_depth;
%  some_time = wf_shift.mean_soma_time;
 [soma_depth_whole_map,soma_time] = find(voltage_map == min_value);

 min_value = abs(min(voltage_map(:))); % normalise map
voltage_map = voltage_map./min_value;

% some maps have noise around tip of the probe (first few channels) which have to removed

if soma_depth_whole_map == 1
voltage_map = voltage_map(20:180,:);
 min_value = min(voltage_map(:)); 
 [soma_depth_whole_map,soma_time] = find(voltage_map == min_value);
end

try
 centered_spline = voltage_map(soma_depth_whole_map-30:soma_depth_whole_map+addm,soma_time-60:soma_time+addt);
 centered_csd = csd_map(soma_depth_whole_map-30:soma_depth_whole_map+addm,soma_time-60:soma_time+addt);
catch 

centered_spline = voltage_map;
centered_csd = csd_map;
end


 min_value = min(centered_spline(:)); % minimum of voltage map

[soma_depth,soma_time] = find(centered_spline == min_value); 


bf_soma = 10; % factor that decides how many timestamps before the somatic spike
% does your minima window start
try
min_map = centered_spline(:, soma_time-bf_soma:soma_time+70); % takes the timestamps from somatic spike onwards
max_csd = centered_csd(:, soma_time-30:soma_time+60);
min_map_full_voltage_map = voltage_map(:,soma_time-bf_soma:soma_time+70);

catch
min_map = centered_spline(:, soma_time-10:soma_time+10);
max_csd = centered_csd(:, soma_time-10:soma_time+10);
min_map_full_voltage_map = voltage_map(:,soma_time-bf_soma:soma_time+10);
end

for iSh = 1:nSh 
	% TO DO: get minima after spike time!!!!!!!!!!
	run_min_shuf(:, :, iSh) = movmin(shuf_spline(:, :, iSh), 60, 2);  % gets the running minimum of each time point
	% in the shuffled maps (row by row)
	
	run_max_csd_shuf(:,:, iSh) = movmax(shuf_csd(:, :, iSh), 60, 2); % same for CSD
end

nCh = size(centered_spline,1);

% mean subtraction
for channel = 1:nCh

mean_channel = mean(centered_spline(channel));
	
centered_spline(channel,1:soma_time-1)-mean_channel;
end



%% getting channel probability vector, computing stats test
 for iCh = 1:nCh
	 
	[ch_min(iCh), min_t(iCh)] = min(min_map(iCh,:));  % gets the minimum for each time point in 
	 % actual voltage map

	 ch_shuffle = squeeze(run_min_shuf(iCh, :, :)); 
	 ch_shuffle = ch_shuffle(:);
	 
	 ch_p(iCh) = 1 - sum(ch_shuffle>ch_min(iCh))/numel(ch_shuffle); % what's the probability to find that
	 % minimum in the real data?
	 % *1-  getting the shuffled n smaller than the actual minimum
 end

 
 % getting channel probability in full map (not centered around soma)
 nCh_all = size(voltage_map,1);
 
  for iCh = 1:nCh_all
	 
	[ch_min_all(iCh), min_t_all(iCh)] = min(min_map_full_voltage_map(iCh,:));  % gets the minimum for each time point in 
	 % actual voltage map

	 ch_shuffle = squeeze(run_min_shuf(iCh, :, :)); 
	 ch_shuffle = ch_shuffle(:);
	 
	 ch_p_all(iCh) = 1 - sum(ch_shuffle>ch_min_all(iCh))/numel(ch_shuffle); % what's the probability to find that
	 % minimum in the real data?
	 % *1-  getting the shuffled n smaller than the actual minimum
 end
 
 
 %%
 % repeat the same for CSD - sometimes gives issues with nCh
 try	 
for iCh = 1:nCh	 
[ch_max_csd(iCh), max_t(iCh)] = max(max_csd(iCh,:));
ch_shuffle_csd = squeeze(run_max_csd_shuf(iCh, :, :)); 
 ch_shuffle_csd = ch_shuffle_csd(:);
ch_p_csd(iCh) = 1 - sum(ch_shuffle_csd<ch_max_csd(iCh))/numel(ch_shuffle_csd);
  end
 catch
ch_p_csd = [];	 
 end

 % P-VALUE
p = 0.01; % divide per n of channels you consider 

 if ~isempty(ch_p_csd)
 signif_bAP_ch_csd = find(ch_p_csd<p); 
signif_bAP_t_csd = max_t(signif_bAP_ch_csd)+soma_time-1; 
% get rid of channels before somatic channel CSD (for plotting, can be removed to look at axon)
soma_pos = find(signif_bAP_ch_csd == soma_depth);
low_ch = find(signif_bAP_ch_csd < soma_depth);
signif_bAP_ch_csd(low_ch) = [];
signif_bAP_t_csd = max_t(signif_bAP_ch_csd)+soma_time-1; 
signif_bAP_ch_csd = signif_bAP_ch_csd(1:size(signif_bAP_t_csd,2));
 else
 end
 
 figure; plot(ch_p_all,'LineWidth',1.2), hold on
  plot(soma_depth_whole_map,0,'r*', 'MarkerSize',14), hold on % plots probability distribution for each channel
% yline(0.05,'--','p-value = 0.05','FontSize',10,'LineWidth',0.9)
yline(0.01,'--','LineWidth',0.4); ax = gca;
yticks([0:0.05:0.9]); ylabel('channel probability','FontSize',12); xlim([0 200])
xlabel('channel index','FontSize',12); title('Mean bAP extent across channels')
set(gca, 'TickDir', 'out'); set (gca,'Xdir','reverse'); set(gca,'view',[90 90]);  box off; axis square


signif_bAP_ch = find(ch_p< p); % ch_p<0.01 set p value to include channels
signif_bAP_t = min_t(signif_bAP_ch)+ (soma_time-1); % gets coordinates of significant data points
signif_bAP_t_axon = min_t(signif_bAP_ch)+ (soma_time-1);
signif_bAP_ch_axon = signif_bAP_ch;

% get rid of channels before somatic channel (for plotting, can be removed to look at axon)
soma_pos = find(signif_bAP_ch == soma_depth); low_ch = find(signif_bAP_ch < soma_depth);
signif_bAP_ch(low_ch) = []; signif_bAP_t = min_t(signif_bAP_ch)+soma_time-1; 

signif_bAP_ch = signif_bAP_ch(1:size(signif_bAP_t,2)); % in case for some reasons they are two different sizes

%% plot VOLTAGE MAP

plotParams.plot_fig = 1;

if plotParams.plot_fig == true 
	
plot_splineField = figure; hold on


min_spline = min(centered_spline(:));
[soma_depth,soma_time] = find(centered_spline == min_spline);

 min_voltagemap = abs(min(centered_spline(:))); 


centered_spline = centered_spline./min_voltagemap;

imagesc(centered_spline)

 set(gcf, 'Position', [80, 100, 630, 850]);

 % if map doesn't have the standard centering around soma (because of location in the brain)
 if size(centered_spline,2) == size(voltage_map,2)
	 
 ylim([soma_depth-50 soma_depth+30])
xlim([soma_time-60 soma_time+80])

 else  % every other map will have the soma in the same location within the map 
 
 if size(centered_spline,1) == 71
 ylim([0 66])
 xlim([0 180])
	 else
 ylim([0 82])
 xlim([0 150]) 
	 end
 end
 
% set(gca,'YTickLabel',[]);
  set(gca,'YDir','normal') %   set(gca,'YDir','reverse')

tb.Visible = 'off';
 scatter(soma_time,soma_depth, 'w', '^', 'LineWidth', 2);


 include_chans = zeros(1,size(signif_bAP_t, 2));

for channel = 1:size(signif_bAP_t, 2) % making logical vector for including channels
	% don-t consider next channel if there's a big gap between timepoints 
	 	if channel > 1 && signif_bAP_t(channel) - signif_bAP_t(channel-1) >= 15
			
		include_chans(channel) = 0;
		break
elseif channel > 10 && signif_bAP_t(channel) < signif_bAP_t(channel-1) 
				  % get rid of sequential channels
				  % whose timepoints seem to go back in time 
		include_chans(channel) = 0;
		break
% 		elseif channel > 10 && signif_bAP_t(channel) <= signif_bAP_t(channel-4) 
% 				% dont consider bap points taking place at the same time (network artefacts)
% 		include_chans(channel) = 0;
% 		break
% 		% when different channels are recorded at the same time (artefact)
% 	   elseif channel > 15 && signif_bAP_t(channel) <= signif_bAP_t(channel-5) 
% 		include_chans(channel) = 0;
% 		break
		elseif channel > 10 && signif_bAP_t(channel) <  signif_bAP_t(channel-10) 
% 		% for bAPs that seem to go back in time (network effect)
		include_chans(channel) = 0;
        break
% 		elseif channel > 1 && signif_bAP_t(channel) - signif_bAP_t(channel-1) > 10
% 		% stop if there's a temporal gap between data points 	
% 		include_chans(channel) = 0;
% 		break
		else
		include_chans(channel) = 1;
		end

end

y = find(include_chans == 1);

try
% extrapolate to get smoother signal from spike time to the last significant time point
signif_time = signif_bAP_t(y)-bf_soma; last_timep = signif_time(end);
signif_time = unique(signif_time(1:2:end));
signif_ch = signif_bAP_ch(y); signif_ch = signif_ch(1:2:end);

last_time_point = find(signif_time == last_timep);

query_time_points = [signif_time(1):2:signif_time(last_time_point)];

interp_bap = interp1(signif_time, signif_ch(1:length(signif_time)), query_time_points, 'linear', 'extrap');

if isempty(interp_bap)
signif_time = signif_bAP_t(y)-bf_soma; 
signif_time = sort(signif_time); 
signif_time = unique(signif_time(1:2:end)); last_timep = signif_time(end);
last_time_point = find(signif_time == last_timep);
query_time_points = [signif_time(1):2:signif_time(last_time_point)];
interp_bap = interp1(signif_time, signif_ch(1:length(signif_time)), query_time_points, 'linear', 'extrap');
end

% plot(query_time_points, interp_bap, '-ok')

catch
plot(signif_bAP_t(y)-bf_soma,signif_bAP_ch(y),'ok-')
query_time_points = signif_bAP_t(y);
interp_bap = signif_bAP_ch(y);
end

% the bf_soma is normally -10/-30 is because you start taking the minima from 10/30 timestamps)
% before spike time

% 
plot(signif_bAP_t(y)-bf_soma,signif_bAP_ch(y),'ok-')
% 
query_time_points = signif_bAP_t(y);
interp_bap = signif_bAP_ch(y);
% 
 bAP_times_and_channels = [round(query_time_points)' round(interp_bap)'];
 
% bAP_times_and_channels = [(signif_bAP_t(y)-30)' (signif_bAP_ch(y))']; 
%  dist = abs(signif_bAP_ch_axon - soma_depth);
% min_dist = min(dist(:));
% idx = find(dist == min_dist);

 ch_signif_axon = []; % tab to collect significant channels in the axon
 
%% to plot axon activity

% conditions to plot axon extent:
% 1) distance between consecutive points should not be > 0.5 ms 
% or more than 5 channels 
% 2) not plotting channels above soma (we are only looking at the axon here)

 for idx = size(signif_bAP_t_axon,2):-1:1
	 
time_point = signif_bAP_t_axon(idx);   

channel_time_point = signif_bAP_ch_axon(idx);  

if idx == 1
previous_channel = signif_bAP_ch_axon(1);
previous_time_point = signif_bAP_t_axon(1);  
else
previous_channel = signif_bAP_ch_axon(idx-1);
previous_time_point = signif_bAP_t_axon(idx-1); 
end

% only take points around somatic spike
	
if time_point - (soma_time+30) > 30 % see point 1)
		continue
	end

	% if spatial distance between two consecutive points is more than
	% 5 channels OR if channel is above soma
	% skip to next point
if channel_time_point - previous_channel  > 5 || channel_time_point > soma_depth  % see point 2)
continue
 
   % if temporal distance between two consecutive points
   % is more than half a millisecond
   % skip to next point
elseif time_point - previous_time_point > 10
break

% if the point before spike time AND it's far from the soma 
% interrupt the loop (probably network artefact)

elseif soma_depth - channel_time_point > 5 && time_point < soma_time

	break

	% if there's a big spatial gap between consecutive significant points
	% stop loop
elseif soma_depth - channel_time_point > 15 && channel_time_point - previous_channel >= -4
	break
else
	
	 sig_ch_axon = channel_time_point;
	 plot(time_point-bf_soma, channel_time_point,'-d','Color',[.3 .3 .3]), hold on
	 ch_signif_axon = [ch_signif_axon sig_ch_axon];
end



 end

 
 
%  end % end loop
 
axon_extent_channel = min(ch_signif_axon);

%%

meanWFs_wind = centered_spline(:);
colormap(colormap_BlueWhiteRed);
%colormap(c_map);
caxis(plotParams.caxisLFP); hcb = colorbar;
% colorTitleHandle = get(hcb,'Title');
% titleString = 'Voltage (\muV)';
% set(colorTitleHandle ,'String',titleString);

box off
if size(centered_spline,2) == size(voltage_map,2) 

if isequal(NP_probe,'1.0')
% set(gca,'YTickLabel',[]);
yticks([1:5:size(voltage_map, 1)])
yticklabels([-500:100:2500])
xticks(0:50:400)
xticklabels([-6:1.5:6])
else
	
yticks([1:10:size(voltage_map, 1)])
yticklabels([-450:150:800])  % NP 2.0 have 15um spacing!!
xticks(soma_time-30:30:400)
xticklabels([-1:3])
	
end

else % if maps are centered around somatic spikes

if isequal(NP_probe,'1.0')
% set(gca,'YTickLabel',[]);
yticks([1:5:size(voltage_map, 1)])
yticklabels([-600:100:800])
xticks(soma_time-30:30:400)
xticklabels([-1:4])
else
	
yticks([1:10:size(voltage_map, 1)])
yticklabels([-450:150:800])  % NP 2.0 have 15um spacing!!
xticks(soma_time-30:30:400)
xticklabels([-1:3])
	
end

end

ylabel('distance from soma ( \mum)')
xlabel('time (ms)')

end

%% get bAP extent above soma and of signal travelling to the axon

try
bAP_extent_shift = 	interp_bap(end) - soma_depth;

% bAP_extent_shift = signif_bAP_ch(y(end)) - soma_depth; % extent of significant channels
catch
	bAP_extent_shift = []; 
end

try
	
	signal_axon_only = soma_depth - axon_extent_channel;
catch
	signal_axon_only = []; signal_axon_dendrite = []; 
end

if isequal(NP_probe,'1.0')
% bAP_extent_micro = bAP_extent*20;
if isempty(bAP_extent_shift)
	bAP_extent_micro_shift = 0;
end
bAP_extent_micro_shift = bAP_extent_shift*20;
axon_extent_micro_shift = signal_axon_only*20;
else
% bAP_extent_micro = bAP_extent*15;
bAP_extent_micro_shift = bAP_extent_shift*15;
end

%% title main figure

if analyse_bursts
	
if info.analyse_locomotion
	
if isequal(info.analyse_running,'running')

if isequal(take_burst_type, 'firsts')
title('bAP running, first spikes of bursts, LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
elseif isequal(take_burst_type, 'seconds')
title('bAP running, second spikes within bursts, LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
elseif isequal(take_burst_type, 'thirds')
title('bAP running, third spikes whitin bursts, LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});

elseif isequal(take_burst_type, 'others')
title('bAP running, all spikes of bursts (except for first), LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
end

elseif isequal(info.analyse_running,'stationary')

if isequal(take_burst_type, 'firsts')
title('bAP stationary, first spikes of bursts, LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
elseif isequal(take_burst_type, 'seconds')
title('bAP stationary, second spikes within bursts, LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
elseif isequal(take_burst_type, 'thirds')
title('bAP stationary, third spikes within bursts, LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});

elseif isequal(take_burst_type, 'others')
title('bAP stationary, all spikes of bursts (except for first), LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
end

end

elseif info.analyse_network_activity
	
	if high_network_analysis
if isequal(take_burst_type, 'firsts')
title('bAP high network, first spikes of bursts, LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
elseif isequal(take_burst_type, 'seconds')
title('bAP high network, second spikes within bursts, LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
elseif isequal(take_burst_type, 'thirds')
title('bAP high network, third spikes whitin bursts, LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});

elseif isequal(take_burst_type, 'others')
title('bAP high network, all spikes of bursts (except for first), LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
end
	
	elseif ~(high_network_analysis)
	
if isequal(take_burst_type, 'firsts')
title('bAP low network, first spikes of bursts, LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
elseif isequal(take_burst_type, 'seconds')
title('bAP low network, second spikes within bursts, LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
elseif isequal(take_burst_type, 'thirds')
title('bAP low network, third spikes whitin bursts, LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});

elseif isequal(take_burst_type, 'others')
title('bAP low network, all spikes of bursts (except for first), LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
end	
	
    end
end


else % if not concentrating on bursts

if ~(info.analyse_whisk_behaviour) && info.analyse_locomotion && ~(info.analyse_network_activity)
	
if isequal(info.analyse_running,'running') && ~(info.analyse_LED_spikes)
title('bAP running LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
elseif isequal(info.analyse_running,'stationary') && ~(info.analyse_LED_spikes)
	title('bAP stationary LED OFF')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
elseif isequal(info.analyse_running,'running') && info.analyse_LED_spikes
	title('bAP running LED ON')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
elseif isequal(info.analyse_running,'stationary') && info.analyse_LED_spikes
	title('bAP stationary LED ON')
subtitle({[strcat('nrn',num2str(nrn),',ind:',num2str(nrn_ID))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});

end

elseif info.analyse_whisk_behaviour && ~(info.analyse_locomotion) && ~(info.analyse_network_activity)
	
if info.analyse_whisking && ~(info.analyse_LED_spikes)
title('bAP whisking LED OFF')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]});
elseif ~(analyse_whisking) && ~(analyse_LED_spikes)

title('bAP not-whisking LED OFF')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]});
elseif info.analyse_whisking && info.analyse_LED_spikes
title('bAP whisking LED ON')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]});
elseif ~(info.analyse_whisking) && info.anlyse_LED_spikes
title('bAP not whisking whisking LED ON')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]});

end

elseif ~(info.analyse_locomotion) && ~(info.analyse_whisk_behaviour) && info.analyse_network_activity
	if high_network_analysis
	title('bAP - high network state')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]});
	else
	title('bAP - low network state')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_shift),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]});
	end
		
end

end

hold off

%% get conduction velocity of bAP
% exclude points just above soma to get better estimate of slope (normally very vertical just above soma)

n = size(signif_bAP_ch(y),2);   % n of data points in bAP extent
x1 = signif_bAP_t(y)-30; time = x1(end) - x1(1); time = time/30000;  % time coordinates, in s
y1 = signif_bAP_ch(y); space = y1(end) - y1(1);
space = (space*20)/10000; % channel extent (for NP 1.0 1 channel is 20um), in cm
velocity = space/time; % cm/s


%% plot csd

if plotParams.plot_fig == true && ~isempty(ch_p_csd)
	
plot_CSD = figure; hold on


max_csd = max(centered_csd(:));
[soma_depth,soma_time] = find(centered_csd == max_csd);

 max_csd = abs(max(centered_csd(:))); 
% 
 centered_csd = centered_csd./max_csd;
 imagesc(centered_csd)

 set(gcf, 'Position', [170, 300, 560, 620]);
 
 % if map doesn't have the standard centering around soma (because of location in the brain)
 if size(centered_spline,2) == size(voltage_map,2)
	 
 ylim([soma_depth-50 soma_depth+30])
xlim([soma_time-60 soma_time+80])

 else  % every other map will have the soma in the same location within the map 
 
 if size(centered_spline,1) == 71
 ylim([0 66])
 xlim([0 150])
	 else
 ylim([0 82])
 xlim([0 150]) 
	 end
 end
% set(gca,'YTickLabel',[]);
set(gca,'YDir','normal')

tb.Visible = 'off';
 scatter(soma_time,soma_depth, 'w', '^', 'LineWidth', 2);
 include_chans = zeros(1,size(signif_bAP_t_csd, 2));

for channel = 1:size(signif_bAP_t_csd, 2) % making logical vector for including channels
	% don-t consider next channel if there's a big gap between timepoints (0.5 ms)
	 	if channel > 1 && signif_bAP_t_csd(channel) - signif_bAP_t_csd(channel-1) > 5
			
		include_chans(channel) = 0;
		continue
		elseif channel > 2 && signif_bAP_t_csd(channel) < signif_bAP_t_csd(channel-1) 
				  % this is a convoluted temporary way to get rid of sequential channels
				  % whose timepoints seem to go back in time
		include_chans(channel) = 0;
		continue
		elseif channel > 10 && signif_bAP_t_csd(channel) <  signif_bAP_t_csd(10) 
% 		% for bAPs that seem to go back in time (network effect)
		include_chans(channel) = 0;

		else
		include_chans(channel) = 1;
		end

end
y = find(include_chans == 1);

include_time = zeros(1,size(y,2));
for channel = 1:size(y,2)
	if channel > 1 && signif_bAP_t_csd(channel) < signif_bAP_t_csd(channel-1)
		include_time(channel) = 0;
		continue
	else
		include_time(channel) = 1;
	end
end

y = find(include_time == 1);
  
plot(signif_bAP_t_csd(y)-30,signif_bAP_ch_csd(y),'ok-')
% remove the -10 when you include time points before somatic spot
meancsd_wind = centered_csd(:);
colormap(colormap_BlueWhiteRed);
	%colormap(c_map);
	caxis(plotParams.caxisLFP);
	hcb = colorbar;
	colorTitleHandle = get(hcb,'Title');
	titleString = 'Current density (AU)';
	set(colorTitleHandle ,'String',titleString);
box off

if size(centered_spline,2) == size(voltage_map,2) 

if isequal(NP_probe,'1.0')
% set(gca,'YTickLabel',[]);
yticks([1:5:size(voltage_map, 1)])
yticklabels([-3500:100:2500])
xticks(soma_time-30:30:400)
xticklabels([-1:3])
else
	
yticks([1:10:size(voltage_map, 1)])
yticklabels([-450:150:800])  % NP 2.0 have 15um spacing!!
xticks(soma_time-30:30:400)
xticklabels([-1:3])
	
end

else % if maps are centered around somatic spike
	
if isequal(NP_probe,'1.0')
% set(gca,'YTickLabel',[]);
yticks([1:5:size(voltage_map, 1)])
yticklabels([-600:100:800])
xticks(soma_time-30:30:400)
xticklabels([-1:3])
else
	
yticks([1:10:size(voltage_map, 1)])
yticklabels([-450:150:800])  % NP 2.0 have 15um spacing!!
xticks(soma_time-30:30:400)
xticklabels([-1:3])
	
end

end

ylabel('distance from soma ( \mum)')
xlabel('time (ms)')


try
bAP_extent_csd = signif_bAP_ch_csd(y(end)) - soma_depth;

catch
bAP_extent_csd = [];

end

if isequal(NP_probe,'1.0')
% bAP_extent_micro = bAP_extent*20;
if isempty(bAP_extent_csd)
	bAP_extent_micro_csd = 0;
end

bAP_extent_micro_csd = bAP_extent_csd*20;

else
% bAP_extent_micro = bAP_extent*15;
bAP_extent_micro_csd = bAP_extent_csd*15;
end

%% titles CSD

if ~(info.analyse_whisk_behaviour) 
	
	if ~(info.analyse_network_activity)
if isequal(info.analyse_running,'running') && ~(info.analyse_LED_spikes)
	title('CSD - bAP running LED OFF')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_csd),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_run(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
elseif isequal(info.analyse_running,'stationary') && ~(info.analyse_LED_spikes)
	title('CSD - bAP stationary LED OFF')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_csd),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});

elseif isequal(info.analyse_running,'stationary') && info.analyse_LED_spikes
	title('CSD - bAP stationary LED ON')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_csd),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});

elseif isequal(info.analyse_running,'running') && info.analyse_LED_spikes
	title('CSD - bAP running LED ON')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_csd),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
else
end
 
	else % when looking at network
		
		if high_network_analysis
		title('CSD - bAP high network state')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_csd),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});
		else
		title('CSD - bAP low network state')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_csd),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]}); %['number of spikes:',num2str(size(gwfparams.spikeTimes,1))],['spikes firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_sp)],['bursts firing rate:',num2str(bAP_map_info_stat(ind_analyse_nrn).frequency_bursts)],['bAP extent:',num2str(bAP_extent_micro),' \mum']});		
		end
end
	
elseif info.analyse_whisk_behaviour

if info.analyse_whisking && ~(info.analyse_LED_spikes)
title('CSD - bAP whisking LED OFF')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_csd),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]});
elseif ~(info.analyse_whisking) && ~(info.analyse_LED_spikes)
title('CSD - bAP not-whisking LED OFF')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_csd),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]});
elseif info.analyse_whisking && info.analyse_LED_spikes
title('CSD - bAP whisking LED ON')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_csd),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]});
elseif ~(analyse_whisking) && anlyse_LED_spikes
title('CSD - bAP not whisking whisking LED ON')
subtitle({[strcat('nrn',num2str(nrn))],['bAP extent:',num2str(bAP_extent_micro_csd),' \mum'],['Number of spikes used: ',num2str(gwfparams.nWf)]});

end
end


hold off

else
plot_CSD = [];
end

%% save figures 

figure_dir = strcat(w_drive,'\','nrn_',num2str(nrn),'\','drift_correction','\',date);
if ~exist(figure_dir,'dir')
	mkdir(figure_dir)
end
if analyse_bursts
	if info.analyse_locomotion
if isequal(info.analyse_running,'running')
	if isequal(take_burst_type, 'firsts')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_running_first_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_running_first_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
	elseif isequal(take_burst_type, 'seconds')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_running_second_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_running_second_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
elseif isequal(take_burst_type, 'thirds')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_running_third_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
try
exportgraphics(plot_CSD,[figure_dir,'\','bAP_running_third_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
catch 
end
elseif isequal(take_burst_type, 'others')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_running_other_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_running_other_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
		
	end

elseif isequal(info.analyse_running,'stationary')		
		if isequal(take_burst_type, 'firsts')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_stationary_first_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_stationary_first_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
		
	elseif isequal(take_burst_type, 'seconds')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_stationary_second_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_stationary_second_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
	elseif isequal(take_burst_type, 'thirds')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_stationary_third_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_stationary_third_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
		
	elseif isequal(take_burst_type, 'others')
		
exportgraphics(plot_splineField,[figure_dir,'\','bAP_stationary_other_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_stationary_other_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
		
		end
		
end
		
elseif info.analyse_network_activity
	
	
	
	if high_network_analysis
	if isequal(take_burst_type, 'firsts')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_high_net_first_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_high_net_first_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
	elseif isequal(take_burst_type, 'seconds')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_high_net_second_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_high_net_second_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
elseif isequal(take_burst_type, 'thirds')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_high_net_third_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
try
exportgraphics(plot_CSD,[figure_dir,'\','bAP_high_net_third_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
catch 
end
elseif isequal(take_burst_type, 'others')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_high_net_other_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_high_net_other_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
		
	end

elseif ~(high_network_analysis)
	
		if isequal(take_burst_type, 'firsts')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_low_net_first_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_low_net_first_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
		
	elseif isequal(take_burst_type, 'seconds')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_low_net_second_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_low_net_second_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
	elseif isequal(take_burst_type, 'thirds')
exportgraphics(plot_splineField,[figure_dir,'\','bAP_low_net_third_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_low_net_third_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
		
	elseif isequal(take_burst_type, 'others')
		
exportgraphics(plot_splineField,[figure_dir,'\','bAP_low_nety_other_sp_bursts',num2str(isi),'_ISI.jpg'], 'Resolution',600)
exportgraphics(plot_CSD,[figure_dir,'\','bAP_low_net_other_sp_bursts_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
		
		end
		
	end
	end
else   % if not looking at bursts specifically
	
if info.analyse_locomotion
	
	
if isequal(info.analyse_running,'stationary') && ~(info.analyse_LED_spikes)
exportgraphics(plot_splineField,[figure_dir,'\','bAP_stationary_drift_corrected',num2str(isi),'_ISI.jpg'], 'Resolution',600)
try
exportgraphics(plot_CSD,[figure_dir,'\','bAP_stationary_drift_corrected_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
catch
end
elseif isequal(info.analyse_running,'running') && ~(info.analyse_LED_spikes)
exportgraphics(plot_splineField,[figure_dir,'\','bAP_running_drift_corrected',num2str(isi),'_ISI.jpg'], 'Resolution',600)
try
exportgraphics(plot_CSD,[figure_dir,'\','bAP_running_drift_corrected_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
catch
end
elseif isequal(info.analyse_running,'running') && info.analyse_LED_spikes
exportgraphics(plot_splineField,[figure_dir,'\','bAP_running_drift_corrected_LED',num2str(isi),'_ISI.jpg'], 'Resolution',600)
try
exportgraphics(plot_CSD,[figure_dir,'\','bAP_running_drift_corrected_CSD_LED',num2str(isi),'_ISI.jpg'], 'Resolution',600)
catch 
end
elseif isequal(info.analyse_running,'stationary') && info.analyse_LED_spikes
exportgraphics(plot_splineField,[figure_dir,'\','bAP_stationary_drift_corrected_LED',num2str(isi),'_ISI.jpg'], 'Resolution',600)
try
exportgraphics(plot_CSD,[figure_dir,'\','bAP_stationary_drift_corrected_CSD_LED',num2str(isi),'_ISI.jpg'], 'Resolution',600)
catch
end
end

elseif info.analyse_whisk_behaviour
	
if info.analyse_whisking && ~(info.analyse_LED_spikes)
exportgraphics(plot_splineField,[figure_dir,'\','bAP_whisking_drift_corrected',num2str(isi),'_ISI.jpg'], 'Resolution',600)
try
exportgraphics(plot_CSD,[figure_dir,'\','bAP_whisking_drift_corrected_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
catch
end
elseif ~(info.analyse_whisking) && ~(info.analyse_LED_spikes)
exportgraphics(plot_splineField,[figure_dir,'\','bAP_not_whisking_drift_corrected',num2str(isi),'_ISI.jpg'], 'Resolution',600)	
try
exportgraphics(plot_CSD,[figure_dir,'\','bAP_not_whisking_drift_corrected_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
catch
end
elseif info.analyse_whisking && info.analyse_LED_spikes
exportgraphics(plot_splineField,[figure_dir,'\','bAP_whisking_LED_drift_corrected',num2str(isi),'_ISI.jpg'], 'Resolution',600)
try
exportgraphics(plot_CSD,[figure_dir,'\','bAP_whisking_LED_drift_corrected_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
catch
end
elseif ~(info.analyse_whisking) && info.analyse_LED_spikes
exportgraphics(plot_splineField,[figure_dir,'\','bAP_not_whisking_LED_drift_corrected',num2str(isi),'_ISI.jpg'], 'Resolution',600)		
try
exportgraphics(plot_CSD,[figure_dir,'\','bAP_not_whisking_LED_drift_corrected_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
catch
end
end

elseif info.analyse_network_activity
	if high_network_analysis
	exportgraphics(plot_splineField,[figure_dir,'\','bAP_high_network_drift_corrected',num2str(isi),'_ISI.jpg'], 'Resolution',600)		
	try
	exportgraphics(plot_CSD,[figure_dir,'\','bAP_high_network_drift_corrected_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
	catch
	end
	
	else
	exportgraphics(plot_splineField,[figure_dir,'\','bAP_low_network_drift_corrected',num2str(isi),'_ISI.jpg'], 'Resolution',600)		
	try
	exportgraphics(plot_CSD,[figure_dir,'\','bAP_low_network_drift_corrected_CSD',num2str(isi),'_ISI.jpg'], 'Resolution',600)
	catch
	end
	end
end
end


end
