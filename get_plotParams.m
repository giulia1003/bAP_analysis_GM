function [plotParams,gwfparams] = get_plotParams(gwfparams,info,sp,NP_probe)
plotParams = [];
plotParams.soma_window = [-10 40]; % window around spike time for determining min/max (in time stamps at 30khz)

 plotParams.trace_spacing = 15;   % how far apart should lines be in trace view (just visualisation)

plotParams.smoothing = 2;           % do not change - THIS IS FOR 'mean' ('Jure' Method).

plotParams.caxisLFP = [-0.5 0.5];   % [-6 6]

% plotParams.caxisLFP = [-0.5 0.5]; % for z-scored data
% sensitivity of color axis for field [-4 4] without std division [-0.5 0.5] with division

plotParams.anim_caxis = [-100 100];   % for animation (logarithmic scale)

plotParams.n_stdevs = 3; %6;        % for defining signal vs. no signal - (Christina: this is what I will be looking at next week for defining signal)
plotParams.prop_amp_lim = 0.06;     % proportion for amplitude clipping when assigning significance of signal at a particular depth - applied to channel minuimum
plotParams.bAP_restr_jitter = 2;    % in units of timepoints (1/30) - how much the depth above can jump as compared to depth below (should be based on physiological velocity of bAP)
% plotParams.n_chans_skip = 2;      % TODO(currently still hard-coded): if there is a gap larger than this ammount of depths, minima will not be assigned in EAP or CSD depth/time plots

plotParams.FieldSmooth = 'n';       % can be 'y' or 'n' for spatial smoothing
plotParams.FieldSmooth_win = 5;     % window (number of depths) for gaussian smoothing of LFP before doing CSD


plotParams.equalChans = unique(info.chanCoordsAll(:,2))'; % this is needed for spline interpolation - All because we don't want to skip any depths where for example two channels are broken!!!

if isequal(gwfparams.analyse_data, 'field')
% 	
 plotParams.meanWFsAxis_default = [(sum(abs(gwfparams.wfWin))/2 - 10), (sum(abs(gwfparams.wfWin))/2 + 85), -35, +35]; % % changed from -15 +30, first two figures are start and end time to plot, last two are how 																														% many depths below and above soma to show respectively
elseif isequal(gwfparams.analyse_data, 'LFP')
plotParams.meanWFsAxis_default = [(sum(abs(gwfparams.wfWin))/2 - 10), (sum(abs(gwfparams.wfWin))/2 + 75), -15, +35]; %FK: for extended waveform plots, LFP
	% plotParams.meanWFsAxis_default = [(sum(abs(gwfparams.wfWin))/2 - 10*12), (sum(abs(gwfparams.wfWin))/2 + 62*12), -15, +35]; %FK: for extended waveform plots in AP data
	% plotParams.meanWFsAxis_default = [(sum(abs(gwfparams.wfWin))/2 - 40), (sum(abs(gwfparams.wfWin))/2 + 50), -15, +35]; %FK: centering somatic spike in plot of burst 'seconds'
end

plotParams.all_neurons = unique(sp.clu); % to have neuron IDs for plots

% parameters for AIS
if isequal(NP_probe,'1.0')
plotParams.chans_above = 5;
else
plotParams.chans_above = 2;    
end
plotParams.chans_below = 2;         % since it is only a single column - inter-site interval is 40 um; AIS can be about 200 um below soma
plotParams.scaling = 'none';        % what to scale the trace by - 'no' perserves amplitudes, 'min' normalises by minimum voltage, 'amp' normalises by difference between minimum and maximum (normally upstroke after spike)
plotParams.plot_dendrites = 'y';     % 'y' or 'n' for plotting dendritic traces
plotParams.cMap = jet;              % parula, jet, hsv, colormap_blueblackred
plotParams.plot_networkActivity = false; % plot WF averages of events that are associated with most/least network activity
plotParams.plotWhiskerTouches = false;		% to analyse whisker touches and get the PSTHs + waveforms of associated spikes
plotParams.fig_format = 'jpeg';     % format for saving figures ('fig', 'jpeg' or 'svg').
plotParams.fig_format_2 = 'fig';
% plotParams.columns = 'both';      % TODO: can be 'highest' or 'both' - this refers to columns of channels on the neuropixels probe that have the first and second highest amplitude
plotParams.dot_prod = false;		% do dot product analysis? - takes a lot of additional time
plotParams.closeAll = true;		% to close figures for each neuron (at the end of loop).
plotParams.save_bAP_extent = true;	% to save the gathered bAP extent values to struct
plotParams.save_splines = true;	% to save spline fields used to generate space-time plots for later analysis (e.g. subtraction)
plotParams.save_fig = true;			% save figures
plotParams.plot_ISIs = true;		% this is if you want a plot of the ISI histogram of the given units
plotParams.plotting_LED_starts = false;
disp('Plotting parameters: set')
end