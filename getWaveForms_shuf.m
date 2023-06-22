function [wf_shuf] = getWaveForms_shuf(gwfparams)
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


%% IMPORTANT - THIS IS FIDDLING AROUND
chMap = 1:385;
nChInMap = numel(chMap);


spikeTimeKeeps = nan(1,gwfparams.nWf);
% waveForms = nan(numUnits,gwfparams.nWf,nChInMap,wfNSamples);
% waveFormsMean = nan(numUnits,nChInMap,wfNSamples);


% curSpikeTimes has to be a column vector
if size(gwfparams.spikeTimes,1) == 1
curSpikeTimes = gwfparams.spikeTimes';
else
curSpikeTimes = gwfparams.spikeTimes;

end
curUnitnSpikes = size(curSpikeTimes,1);
spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
spikeTimeKeeps(1,1:min([gwfparams.nWf curUnitnSpikes])) = sort(spikeTimesRP(1:min([gwfparams.nWf curUnitnSpikes])));
    
% 	  v_map = cell(1, curUnitnSpikes);   

numTimePoints = numel(spikeTimeKeeps);

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
for i = 1:numel(spikeTimeKeeps)
    extendedTimePoints = [extendedTimePoints, (spikeTimeKeeps(i)-gwfparams.wfWin(2)):(spikeTimeKeeps(i)+gwfparams.wfWin(2))];
end

% Remove any duplicate or out-of-bounds time points
% extendedTimePoints = unique(extendedTimePoints);
extendedTimePoints = extendedTimePoints(extendedTimePoints > 0 & extendedTimePoints <= size(gwfparams.mmf.Data.x,2));

% Extract the memory map file at all time points in the extendedTimePoints vector
% selectedData = [gwfparams.mmf.Data.x(:, extendedTimePoints)];

% Compute linear indices for the memory map data
linearIndices = bsxfun(@plus, (1:gwfparams.nCh)', (extendedTimePoints - 1) * gwfparams.nCh);

% Extract the matrices at each time point using linear indexing
selectedData = gwfparams.mmf.Data.x(linearIndices);

% Reshape selectedData to a cell array with each element corresponding to a time point
outputCellArray = []; windowExtent = gwfparams.wfWin(2)*2+1;
outputCellArray = mat2cell(selectedData, gwfparams.nCh, ones(1, numTimePoints) * windowExtent);

% Convert the cell array to a 3D matrix with dimensions 385x401xnumTimePoints
outputMatrix = cell2mat(permute(outputCellArray, [1, 3, 2]));

waveFormsMean = mean(outputMatrix,3);

      
%% old method, using parallel for loop     
% to extract voltage map at every spike time

% IMPORTANT - THIS IS FIDDLING AROUND
% chMap = 1:385;
% nChInMap = numel(chMap);
% 
% % Read spike time-centered waveforms
% unitIDs = unique(gwfparams.spikeClusters);
% numUnits = size(unitIDs,1);
% spikeTimeKeeps = nan(numUnits,gwfparams.nWf);
% waveForms = nan(numUnits,gwfparams.nWf,nChInMap,wfNSamples);
% waveFormsMean = nan(numUnits,nChInMap,wfNSamples);

% for curUnitInd=1:numUnits

% % curUnitInd = numUnits;

% parpool;
% 
% p = gcp('nocreate');     addAttachedFiles(gcp,["getWaveForms_shuf.m"]);
% % 	  v_map = cell(1, curUnitnSpikes);         
%    	parfor curSpikeTime = 1:curUnitnSpikes
% 		try

% 			
% 		tmpWf = mmf.Data.x(1:gwfparams.nCh,spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
%         waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf(chMap,:);
% 		min_value = abs(min(tmpWf(:))); 
%  		tmpWf_min = tmpWf./min_value;
% 		waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf_min(chMap,:);
% 		
% 			v_map{curSpikeTime} = tmpWf_min;
% 		catch
% 		continue
% 		end
% 	end
% 	

% reshapedArrays = reshape(v_map, [1, 1, size(v_map,2)]);
%  
% wfs_mean = cat(3, reshapedArrays{:});
% 
%  waveFormsMean = mean(wfs_mean,3, 'omitnan');
% 
%  delete(gcp('nocreate'));

%% Package in wf struct
% wf_shuf.unitIDs = unitIDs;
wf_shuf.spikeTimeKeeps = spikeTimeKeeps;
% wf_shuf.waveForms = waveForms;
wf_shuf.waveFormsMean = waveFormsMean;

end