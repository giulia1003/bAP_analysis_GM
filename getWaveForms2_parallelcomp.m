function [wf_shift] = getWaveForms2_parallelcomp(gwfparams,drift_stamps,spikeTimes,templateYpos,nrn,hippocampal_units,NP_probe,info, plotParams)
% function wf = getWaveForms(gwfparams)
%
% Extracts individual spike waveforms from the raw datafile, for multiple
% clusters. Returns the waveforms and their means within clusters.
%
% Contributed by C. Schoonover and A. Fink
%
% % EXAMPLE INPUT
% gwfparams.dataDir = '/path/to/data/';    % KiloSort/Phy output folder
% gwfparams.fileName = 'data.dat';         % .dat file containing the raw 
% gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
% gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
% gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
% gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
% gwfparams.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
%
% % OUTPUT
% wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
% wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
% wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
%                                          % nClu: number of different clusters in .spikeClusters
%                                          % nSWf: number of samples per waveform
%
% % USAGE
% wf = getWaveForms(gwfparams);


% Load .dat and KiloSort/Phy output
fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'});


chMap = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab

if isequal(NP_probe,'1.0')
chMap = 1:385;
nChInMap = numel(chMap);
else
chMap = 1:384;
channels_id = [info.chanCoordsAll(:,2) chMap'];
channels_id = sortrows(channels_id);
chMap = channels_id(:,2);
nChInMap = numel(chMap);
end

% Read spike time-centered waveforms
unitIDs = unique(gwfparams.spikeClusters);
numUnits = size(unitIDs,1);

spikeTimes = spikeTimes';
curUnitInd = numUnits;


numTimePoints = numel(spikeTimes);


% introduced by GM on 19.06.2023 
% first we reduce the memory size of the data (mmf) by only taking the
% spiketimes ±200, saving this new array as selectedData
% we use linearIndices to extract the data we want from selectedData for each time point (=spiketime), 
% including the 200 timestamps before and after the time of spike
% The resulting data is stored in a cell array outputCellArray where each element corresponds to a spiketime
% Finally, we convert the cell array to a 3D matrix outputMatrix with dimensions 385x401xnumTimePoints.

% Extend the time points by adding ±200 to each one (mapping extracellular
% space around time of spike)
extendedTimePoints = [];
for i = 1:numel(spikeTimes)
    extendedTimePoints = [extendedTimePoints, (spikeTimes(i)-gwfparams.wfWin(2)):(spikeTimes(i)+gwfparams.wfWin(2))];
end

% Remove any out-of-bounds time points
% extendedTimePoints = unique(extendedTimePoints);
extendedTimePoints = extendedTimePoints(extendedTimePoints > 0 & extendedTimePoints <= size(mmf.Data.x,2));

% Extract the memory map file at all time points in the extendedTimePoints vector
selectedData = [mmf.Data.x(:, extendedTimePoints)];

% Compute linear indices for the memory map data
linearIndices = bsxfun(@plus, (1:gwfparams.nCh)', (extendedTimePoints - 1) * gwfparams.nCh);

% Extract the matrices at each time point using linear indexing
selectedData = mmf.Data.x(linearIndices);

% Reshape selectedData to a cell array with each element corresponding to a time point
 windowExtent = gwfparams.wfWin(2)*2+1;
outputCellArray = mat2cell(selectedData, [gwfparams.nCh], [ones(1, numTimePoints) * windowExtent]);

% in the cell array:
% median subtraction (getting median at every time point across channels)
medianArray = cellfun(@(x) repmat(median(x, 1), size(x,1), 1), outputCellArray, 'UniformOutput', false);
voltageMapArray = cellfun(@(x, y) x - y, outputCellArray, medianArray, 'UniformOutput', false);

% flip matrix if looking at hippocampal units
if hippocampal_units
voltageMapArray = cellfun(@(x) flip(x), voltageMapArray, 'UniformOutput', false);
end

switch NP_probe

    case '1.0'
% apply drift correction
drift_um = drift_stamps;
driftCell = num2cell(drift_um)';

if isequal(NP_probe,'1.0')
x = [20:10:3850]; x = [0 x]'; % x = repmat(x,2)'; x = x(1:size(x,1)/2,:); % it's 192 rows
xCell = repmat({x}, 1, numel(spikeTimes));

else
x = [6:6:1152]; x = [0 x]'; 
xCell = repmat({x}, 1, numel(spikeTimes));
end

extendedTimePoints = [];
for i = 1:numel(spikeTimes)
    extendedTimePoints = [extendedTimePoints; ([spikeTimes(i)-gwfparams.wfWin(2):spikeTimes(i)+gwfparams.wfWin(2)])];
    
end

yCell = mat2cell(extendedTimePoints, ones(1, size(extendedTimePoints, 1)), size(extendedTimePoints, 2))';
vCell = voltageMapArray;
convertToDouble = @(x) double(x);
vCell = cellfun(convertToDouble, vCell, 'UniformOutput', false);

xqCell = cellfun(@(x, y) x + y, xCell, driftCell, 'UniformOutput', false);

yqCell = yCell;


formula = @(yCell, xCell, vCell, yqCell, xqCell) interp2(yCell, xCell, vCell, yqCell, xqCell);

interpolatedVoltageMapArray = cellfun(formula, yCell, xCell, vCell, yqCell, xqCell, 'UniformOutput', false);   
 
        
% Convert the cell array to a 3D matrix with dimensions 385x401xnumTimePoints
outputMatrix = cell2mat(permute(interpolatedVoltageMapArray, [1, 3, 2]));
 
case 'UHD'

% to do: implement drift correction on UHD recordings
   outputMatrix =  cell2mat(permute(voltageMapArray, [1, 3, 2]));

end

waveFormsMean = mean(outputMatrix,3, 'omitnan');

% if you want to make movie of single-spike voltage maps
all_spikes_map_tab = []; all_spikes_shuf_map_tab = []; spike_time_tab = [];

if info.single_spike_maps

if isequal(NP_probe, '1.0')

extractMapFun = @(x) computeSplineLFP_NoHistology(x);
resultMat = arrayfun(extractMapFun, outputMatrix);

elseif isequal(NP_probe, 'UHD')

all_spikes_map_tab = outputMatrix;
spike_time_tab = spikeTimes;

parpool;

p = gcp('nocreate');    addAttachedFiles(gcp,["get_WFs_NO_CSD_spline.m"]); 
addAttachedFiles(gcp,["getWaveForms_shuf.m"]);  addAttachedFiles(gcp,["computeSplineLFP_NoHistology.m"]);   

n = 100;

v_map = cell(1, n);  
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

parfor idx = 1 : numel(shuffle_vec)
    
   	spiketimes_shuf = spikeTimes + shuffle_nums(idx);
    spiketimes_shuf = round(spiketimes_shuf);
     
	[voltage_map_shuf] = get_WFs_NO_CSD_spline(info, plotParams, gwfparams, NP_probe,spiketimes_shuf);
	       
 
     v_map{idx} = voltage_map_shuf;
     
	
end
 delete(gcp('nocreate'));


reshapedArrays = reshape(v_map, [1, 1, size(v_map,2)]);

elseif isequal(NP_probe, 'UHD_4')

extractMapFun = @(x) computeSplineLFP_UHD(x);   
resultMat = arrayfun(extractMapFun, outputMatrix);

vmap_1 = resultMat(1:48,:);
vmap_2 = resultMat(49:end,:);
resultMat = cat(1,vmap_2, vmap_1);

end


end

 %% OLD METHOD, using a loop

% % IMPORTANT - THIS IS FIDDLING AROUND
% if isequal(NP_probe,'1.0')
% chMap = 1:385;
% nChInMap = numel(chMap);
% else
% chMap = 1:384;
% channels_id = [info.chanCoordsAll(:,2) chMap'];
% channels_id = sortrows(channels_id);
% chMap = channels_id(:,2);
% nChInMap = numel(chMap);
% end
% 
% 
% % Read spike time-centered waveforms
% unitIDs = unique(gwfparams.spikeClusters);
% numUnits = size(unitIDs,1);
% 
% spikeTimes = spikeTimes';
% curUnitInd = numUnits;
% 
% parpool;
% 
% p = gcp('nocreate');     addAttachedFiles(gcp,["getWaveForms2.m"]);
%   
%   
%    waveforms_num = [1:1:1:gwfparams.nWf];
%   
%    shifted_map = cell(1, numel(waveforms_num));
%    v_shifted_min_tab = cell(1,numel(waveforms_num));
%    
%   parfor idx = 1:numel(waveforms_num)
% 	
%    curUnitnSpikes = size(drift_stamps,1);
%      
%   curSpikeTime = waveforms_num(idx);
% 	  	% curSpikeTime = 1:min([wfs_num curUnitnSpikes])
% 
% 		
% 	tmpWf = mmf.Data.x(1:gwfparams.nCh,spikeTimes(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimes(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
% 				
% % 		median subtraction (getting median at every time point across channels)
% 		medianwv = median(tmpWf, 1, 'omitnan');
%         median_mat = repmat(medianwv, size(tmpWf, 1), 1);
% 		tmpWf = tmpWf - median_mat;
% 				
% 		
% 		if hippocampal_units
% 		% remember to FLIP MATRIX if hippocampal_units is selected
% 		tmpWf = flip(tmpWf);
% 		end
% 	
% 	
% % 		waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf(chMap,:);
% 		% using interpolation to shift waveform by drift_stamps for each spike
% 		% v_map_shifted = interp2(x,y,v,xq,yq)
% 		% x = position of every channel (1:385) in um
% 		% y = timestamps around spike used to make voltage map (normally 401)
% 		% v = n_Ch x timestamps (385 x 401), tmpWf
% 		% xq = n_Ch around soma +- drift (depending on how probe shifted)
% 		% yq = like y, timestamps in voltage map
% 		
% 		% for the moment no drift correction on the UHD data
% 		
% 		drift_um = drift_stamps(curSpikeTime);
% 		x = []; 
% 		if isequal(NP_probe,'1.0')
% 		x = [20:10:3850]; x = [0 x]'; % x = repmat(x,2)'; x = x(1:size(x,1)/2,:); % it's 192 rows
% 		else
% 		x = [6:6:1152]; x = [0 x]'; 
% 		end
% % 		x = x(:,1);
% 		y = [spikeTimes(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimes(curUnitInd,curSpikeTime)+gwfparams.wfWin(end)]; 
% 		if size(chMap,2) == 1
% 		v = tmpWf(1:size(chMap,1),:); v = double(v);
% 		else
% 		v = tmpWf(1:size(chMap,2),:); v = double(v);
% 		end
% 		unit_pos = templateYpos(nrn); ch_c = x(:,1); % channels' positions are paired (i.e. ch1 and ch2 both at 20um)
% 		n = repmat(unit_pos,size(ch_c,1),1); [minValue, closestIndex] = min(abs(ch_c - n.'));
% 		 % find closest channel to your unit in um
%         closestValue = ch_c(closestIndex); ch_soma = find(ch_c == closestValue); ch_soma_um = ch_c(ch_soma);
% 
% 		if drift_um < 0 
% 		xq = [];
% 		xq = bsxfun(@plus, x,drift_um);
% 		else 
% 		xq = [];
% 		xq = bsxfun(@minus, x,drift_um);
% 		end
% 		yq = y; 
%         
%         v_shifted = interp2(y,x,v,yq,xq);
% 		
% % 		shifted_map(curUnitInd,curSpikeTime,:,:) = v_shifted(chMap,:);
% 		shifted_map{idx} = v_shifted;
% 	 		
% % 		min_shift = abs(min(v_shifted(:))); % normalise every spikemap by negative peak (= soma)
% % 		v_shifted_norm = v_shifted./min_shift;
% % 		shifted_map(curUnitInd,curSpikeTime,:,:) = v_shifted_norm(chMap,:);
% 		
% 		% compare spike amplitude at the soma across maps
% 		% to check for horizontal probe drift 
% 		% put all amplitudes in one structure
% 			
% %  		min_shift = abs(min(tmpWf(:))); 
% %         v_shifted_min_tab{idx} = min_shift;
% % 		
% % 	    v_shifted_min = mean(v_shifted_min_tab, 'omitnan'); 
% % 	[amp_soma_across_time,mean_soma_depth,mean_soma_time] = get_amplitude_and_position_soma(amp_soma_across_time,soma_time_tab,soma_depth_tab,...
% % 	curSpikeTime,shifted_map,curUnitInd);
% 
% 
%  end   % end for loop
% % 
% % wfs_mean_shifted = nan(nChInMap,wfNSamples,size(shifted_map,2));
% 
% reshapedArrays = reshape(shifted_map, [1, 1, size(shifted_map,2)]);
%  
% wfs_mean_shifted = cat(3, reshapedArrays{:});
% 
%  wfs_mean_shifted = mean(wfs_mean_shifted,3, 'omitnan');
% 
% 
% delete(gcp('nocreate'));



%% Package in wf struct
wf_shift.unitIDs = unitIDs;
wf_shift.spikeTimeKeeps = spikeTimes;
wf_shift.waveFormsMean = waveFormsMean;
wf_shift.all_spikes_map = all_spikes_map_tab;
wf_shift.all_spikes_times = spike_time_tab;
wf_shift.shuffled_map = all_spikes_shuf_map_tab;

% try
% wf_shift.min_value_soma = - v_shifted_min_tab;
% catch
% wf_shift.min_value_soma = [];
% end

% wf_shift.soma_amplitude_over_time = amp_soma_across_time;

 

end
