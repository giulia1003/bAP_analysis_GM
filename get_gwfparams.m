function [gwfparams] = get_gwfparams(myKsDir,gwfparams)

% setting parameters
gwfparams.dataDir = myKsDir;
if isequal(gwfparams.analyse_data, 'field')
	datDir = dir(fullfile(myKsDir, '*ap*.bin')); 
% 	datDir_LFP = dir(fullfile(myKsDir, '*lf*.bin'));
elseif isequal(gwfparams.analyse_data, 'LFP') 
	datDir = dir(fullfile(myKsDir, '*lf*.bin')); 
	
elseif isequal(gwfparams.analyse_data, 'both')  
	datDir = dir(fullfile(myKsDir, '*ap*.bin')); 
 	datDirLFP = dir(fullfile(myKsDir, '*lf*.bin'));
end

gwfparams.dataType = 'int16'; % NP data is naturally in int16 format but preprocessing changes this to single 

% if there are different versions of data, get the right one
if numel(datDir) > 1
	disp('There are more than one results when determining the data file.')
	datFiles = [];
	for i = 1:numel(datDir)
		datFiles{i,1} = datDir(i).name;
	end
	datFiles
	corr_datDir = input('Please enter the index of the correct file (see above): ');
	datDir = datDir(corr_datDir);
	disp(['chosen file is: ' datDir.name])
	if contains(datDir.name, '_single')  % adjust data type according to how .lf file was saved
		gwfparams.dataType = 'single';				% Data type of .dat file
	end
	gwfparams.fileName = datDir.name;
else 
	gwfparams.fileName = datDir.name;
	if isequal(gwfparams.analyse_data, 'both')  
    gwfparams.fileNameLFP = datDirLFP.name;
	end
end
gwfparams.nCh = 385;                        % Number of channels

if isequal(gwfparams.analyse_data, 'field')
	gwfparams.wfWin = [-200 200];			% Number of samples before and after spiketime to include in waveform (recommended [-200 200])
	gwfparams.use_nWfs = 2000; %'all';      % 'all' or number of waveforms per unit to pull out (recommended ~ 2000)
elseif isequal(gwfparams.analyse_data, 'LFP')
	gwfparams.wfWin = [-75 75];				% must be symmetrical
	gwfparams.use_nWfs = 3000; %'all';
 elseif isequal(gwfparams.analyse_data, 'both')
		gwfparams.wfWinLFP = [-75 75];		 
	    gwfparams.wfWin = [-200 200];	
end

fprintf(['\nAnalysing data of type:  ', gwfparams.analyse_data,'  - CHECK if correct\n\n'])
end