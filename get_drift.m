function [drift_stamps,spikeTimeKeeps,nrn_idx] = get_drift(nrn,nrn_ID,sp,gwfparams,drift_nrn,batch_id,templateYpos,NP_probe)


drift_tab = []; centre_batch = [];
batch_time = [];


% Load .dat and KiloSort/Phy output
fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'});
chMap = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab


chMap = 1:385;
nChInMap = numel(chMap);

% Read spike time-centered waveforms
unitIDs = unique(gwfparams.spikeClusters);
numUnits = size(unitIDs,1);
spikeTimeKeeps = nan(numUnits,gwfparams.nWf);
waveForms = nan(numUnits,gwfparams.nWf,nChInMap,wfNSamples);
waveFormsMean = nan(numUnits,nChInMap,wfNSamples);
%% extracting waveforms for each spike time, correcting for vertical probe drift 

for curUnitInd=1:numUnits
    curUnitID = unitIDs(curUnitInd);
    curSpikeTimes = gwfparams.spikeTimes(gwfparams.spikeClusters==curUnitID);
    curUnitnSpikes = size(curSpikeTimes,1);
    spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
    spikeTimeKeeps(curUnitInd,1:min([gwfparams.nWf curUnitnSpikes])) = sort(spikeTimesRP(1:min([gwfparams.nWf curUnitnSpikes])));
 
% 	nrn_idx = find(batch_id(:,2)== nrn);  
% 	nrn_spikes = batch_id(nrn_idx,1);

	spikeTimeKeeps = spikeTimeKeeps';


nrn_idx = find(batch_id(:,2) == nrn_ID);
batch_unit = batch_id(nrn_idx,5); % find where unit is in time (batch), batches are saved in column 5 of batch_id
max_batch = max(batch_unit); 
% if max_batch > size(drift_nrn,1) % discard batches that are not computed in drift array
% t = find(batch_unit > size(drift_nrn,1));
% batch_unit(t) =[];
% nrn_idx(t) = [];
% nrn_spikes(t) = [];
% end

% get drift for every segment of the probe, apply correction depending on where unit sits 

if isequal(NP_probe,'1.0')
	
if templateYpos(nrn) <= 390
a = mean(drift_nrn(:,1));
drift_tab = repmat(a,size(drift_nrn,1),1);
 elseif templateYpos(nrn) <= 765
a = [drift_nrn(:,1),drift_nrn(:,2)]; % unit falls between two blocks, get mean drift
drift_tab = mean(a,2);
elseif templateYpos(nrn) <= 1150
a = [drift_nrn(:,2),drift_nrn(:,3)];
drift_tab = mean(a,2);
elseif templateYpos(nrn) <= 1530
a = [drift_nrn(:,3),drift_nrn(:,4)];
drift_tab = mean(a,2);			
elseif round(templateYpos(nrn)) <= 1915
a = [drift_nrn(:,4),drift_nrn(:,5)];
drift_tab = mean(a,2);			
elseif templateYpos(nrn) < 2300
a = [drift_nrn(:,5),drift_nrn(:,6)];
drift_tab = mean(a,2);					
elseif templateYpos(nrn) <= 2680
a = [drift_nrn(:,6),drift_nrn(:,7)];
drift_tab = mean(a,2);	
elseif templateYpos(nrn) <= 3065
a = [drift_nrn(:,7),drift_nrn(:,8)];
drift_tab = mean(a,2);			
elseif templateYpos(nrn) <= 3445
a = [drift_nrn(:,8),drift_nrn(:,9)];
drift_tab = mean(a,2);				
elseif templateYpos(nrn) <= 3830
a = mean(drift_nrn(:,9));
drift_tab = repmat(a,size(drift_nrn,1),1);
 else
end
 
elseif isequal(NP_probe, 'UHD') || isequal(NP_probe, 'UHD_4') 
	
if templateYpos(nrn) <= 390
a = mean(drift_nrn(:,1));
drift_tab = repmat(a,size(drift_nrn,1),1);
 elseif templateYpos(nrn) <= 765
a = [drift_nrn(:,1),drift_nrn(:,2)]; % unit falls between two blocks, get mean drift
drift_tab = mean(a,2);
elseif templateYpos(nrn) <= 1150
a = [drift_nrn(:,2),drift_nrn(:,3)];
drift_tab = mean(a,2);
elseif templateYpos(nrn) <= 1530
a = [drift_nrn(:,3),drift_nrn(:,4)];
drift_tab = mean(a,2);			
elseif round(templateYpos(nrn)) <= 1915
a = [drift_nrn(:,4),drift_nrn(:,5)];
drift_tab = mean(a,2);			
elseif templateYpos(nrn) < 2300
a = [drift_nrn(:,5),drift_nrn(:,6)];
drift_tab = mean(a,2);					
elseif templateYpos(nrn) <= 2680
a = [drift_nrn(:,6),drift_nrn(:,7)];
drift_tab = mean(a,2);	
elseif templateYpos(nrn) <= 3065
a = [drift_nrn(:,7),drift_nrn(:,8)];
drift_tab = mean(a,2);			
elseif templateYpos(nrn) <= 3445
a = [drift_nrn(:,8),drift_nrn(:,9)];
drift_tab = mean(a,2);				
elseif templateYpos(nrn) <= 3830
a = mean(drift_nrn(:,9));
drift_tab = repmat(a,size(drift_nrn,1),1);
 else
end

end

%% interpolating through drift batches to get drift (um) at each spike time

xdrift = drift_tab;
xdrift_rep =[];
xdrift_rep = repmat(xdrift,1,10)';
xdrift = xdrift_rep(:);

% using drift()
% get yourself the average length of each batch in timestamps = max(sp.st)*30000/n of batches
% repeat drift value for each batch * timestamps in one batch
% take average drift value if units sits in between overlapping chunks of the probe
% x = length of recording in timestamps (size(batch_id,1))
% y = vector with drift across entire recording (check that length in timestamps matches)
% p = 0.001 smoothing factor
% drift_rec = csaps(x,y,p);
p = 0.7;

% converting spike times to the corresponding point in the interpolation 
stamps2chunk = size(batch_id,1)/size(xdrift,1);
spiketimes_chunk = spikeTimeKeeps./stamps2chunk/10;
% 
 x = 1:1:size(xdrift,1);
  
drift_rec = csaps(x,xdrift,p); % drift interpolation across the recording, you need to extract values for each time point
drift_stamps = fnval(drift_rec,spiketimes_chunk); % interpolated drift for every spike time
% drift_stamps_x = fnval(drift_rec,x);

% NEW METHOD: linear interpolation strategy to do drift correction (currently has a lot of background noise)

% x = length of recording in timestamps (sp.st)
% batch_length = x/n_of_batches
% drift_timestamps = (1:n_of_bacthes)* batch_length - (batch_length/2)
% drift_spikes = interp1(drift_timestamps,drift_tab,spiketimes);

% x = size(sp.st,1);
% n_of_batches = size(drift_nrn,1);
% batch_length = x/n_of_batches;
% drift_timestamps = (1:n_of_batches)*batch_length - (batch_length/2);
% drift_spikes = interp1(drift_timestamps,drift_tab,spikeTimeKeeps,'linear','extrap'); 
% pchip extrapolates by default. otherwise use 'linear','extrap'

% drift_stamps = drift_spikes;


%% plotting drift
% figure,plot(drift_stamps_x) 
%  hold on, scatter(spiketimes_chunk, ones(length(spiketimes_chunk),1), 'r.')
%  xlabel('batch')
%  ylabel('drift (\mum)')
%  legend('drift (\mum)','spike time')
%  xlim([0 23052])
%  box off
 

% to plot drift traces along the probe:
% drift1 = drift(:,1); 
% drift2 = drift(:,2); drift2 = bsxfun(@plus, drift2, 10);
% drift3 = drift(:,3); drift3 = bsxfun(@plus, drift3, 20);
% drift4 = drift(:,4); drift4 = bsxfun(@plus, drift4,30);
% drift5 = drift(:,5); drift5 = bsxfun(@plus, drift5, 40);
% drift6 = drift(:,6); drift6 = bsxfun(@plus, drift6, 50);
% drift7 = drift(:,7); drift7 = bsxfun(@plus, drift7, 60);
% drift8 = drift(:,8); drift8 = bsxfun(@plus, drift8, 70);
% drift9 = drift(:,9); drift9 = bsxfun(@plus, drift9, 80);
% 
%  figure,plot(drift1), hold on, 
%  plot(drift2)
%  plot(drift3)
%  plot(drift4)
 

 
 
end
end