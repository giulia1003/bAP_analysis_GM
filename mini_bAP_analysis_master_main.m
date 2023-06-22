%% Reduced bAP script - only looking at bAP maps comparisons
% set paths
% .paq files should be in the same folder as the recording you want to analyse

code_folder = 'W:\code';
cd(code_folder)
addpath(genpath('NP2_bAP_analysis_GM_current'));

cortexfolder = 'C:\Users\gmast\OneDrive\Desktop\spikes-master';
addpath(genpath(cortexfolder));

%%
data_path = 'Z:\data\raw_data\cb387\20230509\cb387_20230509_dstripe'; 
% W:\data\raw_data\cb365\20220613\cb365_20220613_all 
% W:\data\raw_data\cb330\20230113\cb330_20230113_all
% W:\data\raw_data\cb341\20210830\cb341_20210830_all
% W:\data\raw_data\cb369\20220823\cb369_20220823_all 
% W:\data\raw_data\cb370\20220831\cb370_20220831_all
% W:\data\raw_data\cb365\20220613\cb365_20220613_all 
% W:\data\raw_data\cb374\20221109\cb374_20221109_all
% W:\data\raw_data\cb375\20221121\cb375_20221121_all
% W:\data\raw_data\cb377\20221127\cb377_20221127_all
% W:\data\raw_data\cb378\20221130\cb378_20221130_all
% W:\data\raw_data\cb381\20230207\cb381_20230207_all
% W:\data\raw_data\cb386\20230425\cb386_20230425_all
% W:\data\raw_data\cb389\20230614\cb389_20230614_all
data_np_split = 'W:\data\raw_data\cb381\20230207\'; % C:\raw_data\cb349\20210804  %'F:\raw_data\cb341\20210830\'; %'F:\raw_data\cb348\20210723\cb348_20210723_all_11' % where non-concatenated files are stored
% F:\raw_data\cb341\20210830   E:\raw_data\cb336\20210708
%W:\data\raw_data\cb297\20200418\cb297_20200418_g6_14  % W:\data\raw_data\cb365\20220613
NP_probe = 'UHD';   analyse_multiple_sessions = true; session = 'session_g6'; % LED analysis
[mouse,session,linux,rec_date,insertion_angle_name,myKsDir,procDataDir_General,infodir] = get_path(data_path,NP_probe);
if isequal(NP_probe,'1.0')
chan_map_path = 'W:\code\chanMap\neuropixPhase3A_kilosortChanMap'; % 'C:\Users\tex_analysis\Desktop\Dendrites_project\NP2_bAP_analysis_GM_v20220104\chanMap_NP2\NP2_kilosortChanMap'; load(chan_map_path)
else
chan_map_path = 'W:\code\chanMap\NPUHD2_inner_vstripe_ref0';
end
w_drive = strcat(data_path,'\','figures');
load(chan_map_path)
running_threshold = 0.42;   	% choose threshold for running assignment in cm/sec (could be 2) 
stationary_threshold = 0.1;      % could be 0.5
infodir = [myKsDir, '/info_nrn/'];
if exist(infodir, 'dir') ~= 7
	mkdir(infodir)
end

%% Loading spiking data from kilosort/phy easily and visualising features for cortex
% params.excludeNoise = true;
sp = loadKSdir(myKsDir)%, params)    % sp.st are spike times in seconds
                                     % sp.clu are cluster identities
if size(sp.cgs,1) == 1
    sp.cgs = sp.cgs';
end
% spikes from clusters labeled "noise" have already been omitted
if size(sp.cids,1) == 1
    sp.cids = sp.cids';
end


% identifying cortical neuron

% templateYpos is in um! 1 = tip of probe

procDataDir_General = [data_path, '\general'];                     % this is where variables such as cortical channels will be stored
previously_run = exist([procDataDir_General, '\crtx_brkn_chans.mat'], 'file');

if previously_run
	load([procDataDir_General, '\crtx_brkn_chans.mat']); % fix so that it goes across days as well
	load([procDataDir_General, '\templateYpos.mat']);
	fprintf(['Cortical depth set at previous script run. Bottom channel set to: ' num2str(crtx_bottCh) ', top channel set to: ' num2str(crtx_topCh) ' \nDisplaying plots from previous script run.\n'])
	connected = true(384,1); connected(bad_chans) = 0;
	
elseif ~previously_run
	mkdir(procDataDir_General)
	fprintf('Plotting spike depths, amplitudes and LFP power spectrum for determining cortical channels...\n\n')
		[bAP_nrn_ind,templateYpos, crtx_bottCh, crtx_topCh,driftEvents] = define_crtx(myKsDir, sp, linux, procDataDir_General, previously_run); % NOTE: this still needs some corrections
	
	% 4. Identify broken channels
	% do not choose reference channels here, this is done separately
	plotLFP(myKsDir,linux, procDataDir_General, crtx_bottCh, crtx_topCh, connected)         % plotting lfp plot again to check for broken channels
	bad_chans = get_broken_chans();
	connected = true(384,1); connected(bad_chans) = 0;
	save([procDataDir_General, '/crtx_brkn_chans'], 'crtx_bottCh', 'crtx_topCh', 'bad_chans', 'bAP_nrn_ind') 	% saving info about cortex extent and broken channel IDs so that the                                                                                                      % section for identification doesn't have to be repeated upon each run
	save([procDataDir_General, '/templateYpos'], 'templateYpos')
end

% making 'info' file for data without histology (our data)
% taking into account additional broken channels (besides the reference)

info.region = 'all';
info.chanMap = chanMap(connected);          % this is including reference and broken channels - to get y positions!
info.chans = connected;                     % NOTE: This struct field is something else in Marta's analysis!
info.connected = connected;
info.connectedPar = getConnectedPar(connected);         % this keeps track of which channels are interpolated over when parallel channels are combined
% FK, this seems to just indicate if in this row there is a ref chan instead of two data chans
info.connectedParSpline = info.connectedPar(4:end);     % for alignment, because spline gets rid of channels at the side
info.chanCoords = [xcoords(connected), ycoords(connected)];
info.neuronID = [unique(sp.clu)];
info.chanMapAll = chanMap;                  % this is including reference and broken channels - to get y positions!
info.chanCoordsAll = [xcoords, ycoords];

% Setting general parameters for extracting raw data
% choose if you want to analyse AP or LFP data

gwfparams.analyse_data = 'field';	% setting can be 'field' or 'LFP' or 'both'
[gwfparams] = get_gwfparams(myKsDir,gwfparams);

% setting specific parameters for extracting multiple spikes for averaging and plotting
close all
[plotParams,gwfparams] = get_plotParams(gwfparams,info,sp,NP_probe);

%% important switches !!

info.analyse_LED = false;
info.analyse_LED_spikes = false;   
info.baseline_rec = false; 

info.analyse_whisk_behaviour = false;  % true to look at whisking
info.analyse_locomotion = false;     % true to look at at locomotion (running vs stationary)
info.analyse_network_activity = true;

info.single_spike_maps = false;  % true if you want to get voltage maps for every spike (takes a long time)
info.take_all_spikes = false; % true to get all spikes across ISI in voltage maps 

%% Get LED trigger times

info.analyse_LED_starts = false;
triggers_per_session = [5 300];
file_name = 'LED_trig_data.mat';

if info.analyse_LED
if exist([infodir,'\',file_name])
	try
	load(strcat(infodir,'\',file_name))
	catch
	load(strcat(infodir,file_name))	
	end
	if isempty(LED_trig_data.all.offset.LED_on_periods)
		disp('LED periods are empty')
	end
else	
prompt = 'was the FIRST session a baseline recording? 1/0: ';
baseline_start = input(prompt); % true first recording is baseline (no LED stims);

if baseline_start
prompt = 'type number of baseline recordings at the beginning of the experiment: ';
baseline_start_n = input(prompt);
else
baseline_start_n = [];
end
prompt = 'was the LAST session a baseline recording? 1/0: ';
baseline_end = input(prompt); % true first recording is baseline (no LED stims);


	% get the LED trigger starts and ends and the periods where the LED was on
	%to run the function you need a folder called 'packIO' with all necessary packIO files in your data directory 
		% could also be made to take packIO files directly from Y -- MAKE SURE THE Y DRIVE IS CONNECTED
 [LED_trig_data] = LEDpower_trigger_gm(data_path, data_np_split, mouse, rec_date, triggers_per_session,NP_probe,baseline_start,baseline_end,baseline_start_n);
 if ~isempty(LED_trig_data)
if exist(infodir, 'dir') ~= 7
	mkdir(infodir)
end
 save ([infodir,'\',file_name], 'LED_trig_data')
 end
end

if info.analyse_LED_starts == true
		plotParams.plotting_LED_starts = true;
		analyse_nrn = 1; 
	end
else
	LED_trig_data = []; all_sesh_ON_periods.LED_on_periods = [];
end

%% initialise some variables %

bAP_extent_table = []; LED_data_sum = []; spikes_bhv = []; spikes_bhv_b = []; LED_data_all = []; bAP_map_info = [];
spikenums = []; field_pool = []; soma_depth_sp_Field_tab = []; amps_vect_sp_collect = [];nrn_info_file_all = []; diff_plot = [];
burst_freq = []; AUC_run = []; AUC_stat = [];


%% Get locomotion, whisking, and LED times

[sesh_ON_periods,all_sesh_ON_periods,sesh_ON_starts,templates_table,gwfparams,running_running,running_runningTimes,stat_runningTimes,time_motion,whisking_times,not_whisking_times,whisking_vect] = load_bhv_params(info,gwfparams,LED_trig_data,...
	procDataDir_General,data_path,sp,running_threshold, stationary_threshold,bAP_nrn_ind,analyse_multiple_sessions,session);

%% Load probe drift from Kilosort 3.0

load([data_path,'\','probe_drift.mat'])
drift_nrn = drift; % drift is a function of some financial matlab packages, creates issues as variable name
clear drift
load([data_path,'\','batch_id_drift.mat'])  
% batch_id: column 1) spike times, column 2) cluster id, column 5) batch (time segment)
% use rez.st3 to extract batch_id

%% Select the units you want to analyse

analyse_nrn = [61];  %bAP_nrn_ind';  to plot for all...
% check sp.cids to get correct nrn (inseert INDEX here)

% cb341 [216,218,219,221]; 
% cb365 [373,386,393,397]  L2/3 [528, 541] 
% cb370 [330,339,350,358,368,376,378,384];
% day 1: cb374 [279,285,290,291,294,306,323];     
% day 2: cb374 [307,309,321,324,327,330,347];
% cb373 [64,78,157] HIPPO
% cb374 [175,181,226]; HIPPO
% cb374 [279,285,290,291,294,306,323];
% cb375 [17,18,33,35,36,40,41,42,48,51,54,58,61];   HIPPO
% cb375 [114,132,148,153,158,160,161,163];
% cb377 day 1 [183,184,185,194,195,196,197,199,201,202,203,215,216,257,284];
% cb377 day 2 [252,254,284,286,289,293,294,295];
% cb378 [81,82,83,84,95,96,99,100,101,103,105,110,115,116,119,124,126,260];
% cb330/380 [236,240,247,248,257,258,261,271,272];   % interneuron: 279
% cb381 [179, 181, 187, 188,192];

% cb386 [252,255];
nrn = analyse_nrn;
close all   

%% get isi properties and bursty neurons IDs 
plot_figs = true;
% isi distribution is in ms
[isi_properties,bursty_neurons,non_bursty_neurons] = get_isi_distribution(analyse_nrn,sp,infodir,plot_figs);

%% plot multiunit activity and LFP trace
% for the whole recording

cd(infodir)

short_stims = true;
[mua_plot] = plot_mua(sp,running_running,LED_trig_data,whisking_vect,short_stims);
saveas(mua_plot,'MUA_running_led.fig');
close(mua_plot);

% get raster plot of multiunit activity

get_raster_activity(data_path,LED_trig_data,sp);

% [lfp_plot] = plot_lfp(myKsDir);
% saveas(lfp_plot,'LFP_plot.fig');

%% Get spiketimes for each neuron across different ISI

load_spikes = false;% true if you want to load spiketimes that were previously computed, false if you want to recompute them

[spike_nums_5type1,spike_nums_5type2,spike_nums_20type1,spike_nums_20type2,spike_nums_100type2,...
spike_nums_100type1,spike_nums_50type2,spike_nums_50type1,spike_nums_500type2,spike_nums_500type1,...
gwfparams,analyse_bursts,take_burst_type]  = get_spiketimes(info,load_spikes,infodir,analyse_nrn,gwfparams,...
	running_runningTimes,stat_runningTimes,sesh_ON_periods,sp,LED_trig_data,...
	all_sesh_ON_periods,plotParams,LED_data_all,time_motion,not_whisking_times,whisking_times);

%% set structures where you can collect the data later on
% type1 = stationary  // whisking
% type2 = running  // not whisking 

bAP_map_info_type1_500 = []; bAP_map_info_type2_20 =[]; bAP_map_info_type1_5 = []; bAP_map_info_type1_5 = []; bAP_map_info_type1_all = []; bAP_map_info_type2_5 =[]; 


%% Plot spike-triggered voltage maps 
% for each unit looking at different ISI and behavioural conditions (running vs stationary, LED ON vs LED OFF peridos)

% for LED OFF periods, make sure you use the same fixed n of spikes
% as in the LED ON periods (not the minimum above a threshold)
delete(gcp('nocreate'));

sp_threshold = 500; % minimum n of  spikes to make voltage maps
use_fixed_n_spikes = true; % use the same number of spikes for each condition (not minimum), indicated by sp_threshold
 
hippocampal_units = false;
if hippocampal_units  % warning in case you select the wrong condition
	disp('analsyng hippocampal neurons - flipping maps')
end

tic

[bAP_map_info_type1_all,bAP_map_info_type1_500,bAP_map_info_type2_500,bAP_map_info_type1_100,bAP_map_info_type2_100,bAP_map_info_type2_50,bAP_map_info_type1_50,bAP_map_info_type2_20,bAP_map_info_type1_20,...
	bAP_map_info_type2_5,bAP_map_info_type1_5] = get_voltage_maps_isi(analyse_nrn,spike_nums_5type1,spike_nums_5type2,spike_nums_500type1,spike_nums_500type2,spike_nums_100type1,spike_nums_100type2,...
	spike_nums_20type1,spike_nums_20type2,spike_nums_50type1,spike_nums_50type2,gwfparams,w_drive,rec_date,info, plotParams, NP_probe,sp,sp_threshold,drift_nrn,batch_id,templateYpos,hippocampal_units,analyse_bursts,take_burst_type, use_fixed_n_spikes);  
toc
%% get bAP amplitudes 
% using the same channels across conditions

info_file = bAP_map_info_type2_20;

[new_info_file] = get_amplitudes_vect_equal(bAP_map_info_type1_500, info_file);

bAP_map_info_type2_20 = new_info_file;

%% plot amplitudes across conditions

nrn = 5;
amp_dir = strcat(w_drive,'\','bAP_amplitudes');
if ~exist(amp_dir,'dir')
	mkdir(amp_dir)
end
cd(amp_dir)

plot_amplitudes_across_conditions(info,nrn, bAP_map_info_type1_500, bAP_map_info_type1_100, bAP_map_info_type1_50, bAP_map_info_type1_20, bAP_map_info_type1_5, ...
	bAP_map_info_type2_500, bAP_map_info_type2_100, bAP_map_info_type2_50, bAP_map_info_type2_20, bAP_map_info_type2_5)


%% save data
close all

bap_dir = ([infodir,'\',date]);
if ~exist(bap_dir)
	mkdir(bap_dir)
end
cd(bap_dir)
if info.take_all_spikes
save ([bap_dir, '\','nrn_bAP_map_all_single_spikes.mat'], 'bAP_map_info_type1_all','-v7.3')

elseif ~(info.take_all_spikes)
	
if analyse_bursts
	if isequal(take_burst_type, 'firsts')
	save ([bap_dir,'\','nrn_bAP_map_type1_first_spikes_bursts.mat'], 'bAP_map_info_type1_5','-v7.3')
	save ([bap_dir,'\','nrn_bAP_map_type2_first_spikes_bursts.mat'], 'bAP_map_info_type2_5','-v7.3')
	elseif isequal(take_burst_type, 'seconds')
	save ([bap_dir,'\','nrn_bAP_map_type1_second_spikes_bursts.mat'], 'bAP_map_info_type1_5','-v7.3')
	save ([bap_dir,'\','nrn_bAP_map_type2_second_spikes_bursts.mat'], 'bAP_map_info_type2_5','-v7.3')
	elseif isequal(take_burst_type, 'thirds')
	save ([bap_dir,'\','nrn_bAP_map_type1_third_spikes_bursts.mat'], 'bAP_map_info_type1_5','-v7.3')
	save ([bap_dir,'\','nrn_bAP_map_type2_third_spikes_bursts.mat'], 'bAP_map_info_type2_5','-v7.3')
	elseif isequal(take_burst_type, 'all')
		save ([bap_dir,'\','nrn_bAP_map_all_spikelets_bursts.mat'], 'bAP_map_info_type1_all','-v7.3')
	end
elseif ~(analyse_bursts)

if info.analyse_LED_spikes && info.analyse_locomotion
save ([bap_dir,'\','nrn_bAP_map_run_100_LED.mat'], 'bAP_map_info_type2_100','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_run_500_LED.mat'], 'bAP_map_info_type2_500','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_run_50_LED.mat'], 'bAP_map_info_type2_50','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_run_20_LED.mat'], 'bAP_map_info_type2_20','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_run_5_LED.mat'], 'bAP_map_info_type2_5','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_stat_100_LED.mat'], 'bAP_map_info_type1_100','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_stat_500_LED.mat'], 'bAP_map_info_type1_500','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_stat_50_LED.mat'], 'bAP_map_info_type1_50','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_stat_20_LED.mat'], 'bAP_map_info_type1_20','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_stat_5_LED.mat'], 'bAP_map_info_type1_5','-v7.3')

elseif ~(info.analyse_LED_spikes) && info.analyse_locomotion
save ([bap_dir,'\','nrn_bAP_map_run_100.mat'], 'bAP_map_info_type2_100','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_run_500.mat'], 'bAP_map_info_type2_500','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_run_50.mat'], 'bAP_map_info_type2_50','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_run_20.mat'], 'bAP_map_info_type2_20','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_stat_100.mat'], 'bAP_map_info_type1_100','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_stat_500.mat'], 'bAP_map_info_type1_500','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_stat_50.mat'], 'bAP_map_info_type1_50','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_stat_20.mat'], 'bAP_map_info_type1_20','-v7.3')	
save ([bap_dir,'\','nrn_bAP_map_stat_5.mat'], 'bAP_map_info_type1_5','-v7.3')		
save ([bap_dir,'\','nrn_bAP_map_run_5.mat'], 'bAP_map_info_type2_5','-v7.3')	

elseif info.analyse_LED_spikes && info.analyse_whisk_behaviour
	
save ([bap_dir,'\','nrn_bAP_map_not_whisk_100_LED.mat'], 'bAP_map_info_type2_100','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_not_whisk_500_LED.mat'], 'bAP_map_info_type2_500','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_not_whisk_50_LED.mat'], 'bAP_map_info_type2_50','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_not_whisk_20_LED.mat'], 'bAP_map_info_type2_20','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_not_whisk_5_LED.mat'], 'bAP_map_info_type2_5','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_whisk_100_LED.mat'], 'bAP_map_info_type1_100','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_whisk_500_LED.mat'], 'bAP_map_info_type1_500','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_whisk_50_LED.mat'], 'bAP_map_info_type1_50','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_whisk_20_LED.mat'], 'bAP_map_info_type1_20','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_whisk_5_LED.mat'], 'bAP_map_info_type1_5','-v7.3')

elseif ~(info.analyse_LED_spikes) && info.analyse_whisk_behaviour

save ([bap_dir,'\','nrn_bAP_map_not_whisk_100.mat'], 'bAP_map_info_type2_100','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_not_whisk_500.mat'], 'bAP_map_info_type2_500','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_not_whisk_50.mat'], 'bAP_map_info_type2_50','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_not_whisk_20.mat'], 'bAP_map_info_type2_20','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_whisk_100.mat'], 'bAP_map_info_type1_100','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_whisk_500.mat'], 'bAP_map_info_type1_500','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_whisk_50.mat'], 'bAP_map_info_type1_50','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_whisk_20.mat'], 'bAP_map_info_type1_20','-v7.3')	
save ([bap_dir,'\','nrn_bAP_map_whisk_5.mat'], 'bAP_map_info_type1_5','-v7.3')		
save ([bap_dir,'\','nrn_bAP_map_whisk_5.mat'], 'bAP_map_info_type2_5','-v7.3')

elseif info.analyse_network_activity

save ([bap_dir,'\','nrn_bAP_map_high_network_100.mat'], 'bAP_map_info_type2_100','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_high_network_500.mat'], 'bAP_map_info_type2_500','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_high_network_50.mat'], 'bAP_map_info_type2_50','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_high_network_20.mat'], 'bAP_map_info_type2_20','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_high_network_5.mat'], 'bAP_map_info_type2_5','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_low_network_100.mat'], 'bAP_map_info_type1_100','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_low_network_500.mat'], 'bAP_map_info_type1_500','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_low_network_50.mat'], 'bAP_map_info_type1_50','-v7.3')
save ([bap_dir,'\','nrn_bAP_map_low_network_20.mat'], 'bAP_map_info_type1_20','-v7.3')	
save ([bap_dir,'\','nrn_bAP_map_low_network_5.mat'], 'bAP_map_info_type1_5','-v7.3')		


	
end
end
end

%% Get PSTHs for LED ON periods
% % looking at neuronal activity at stim onset and offset 

eventTimes = LED_trig_data.session_g8.offset.LED_trig_starts./30000;  % in SECONDS
% 
window = [-0.8 30];                   % look at spike times from 0.3 sec before each event to 1 sec after
 trialGroups = ones(size(eventTimes)); 
 
% 
psthViewer(sp.st, sp.clu, eventTimes, window, trialGroups);	


%% compare bAP amplitudes in relation to soma distance
% subtract dist amp matrices to compare conditions

nrn = analyse_nrn; 

for ind_analyse_nrn = 1:size(bAP_map_info.bAP_map_info_run,2)	
subplot6 = figure;
subplot(1,2,1)
plot_dist_amp2(ind_analyse_nrn,bAP_map_info,plotParams, nrn,NP_probe)
title('Mean spike amplitude, LED OFF')
subplot(1,2,2)
plot_dist_amp3(ind_analyse_nrn,bAP_map_info,plotParams, nrn,NP_probe)
title('Mean spike amplitude, LED ON')
set(gcf, 'Position', [180, 60, 1700, 980]);
exportgraphics(subplot6,[w_drive,'\',rec_date,'\','nrn_',num2str(nrn(ind_analyse_nrn)),'\','bAP_dist_amplitude.jpg'], 'Resolution',600)
close(subplot6)
end


 