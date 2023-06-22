
function [gwfparams] = choose_burst_spikes(ISI_range, prespike_restriction, postspike_restriction, take_burst_spikes, take_burst_type, gwfparams, nrn_ID)


	if ~isequal(ISI_range, 'all') % if we do not want all ISI spikes - expect a two element vector and specify times
% 		
 		ISI_range = ISI_range .* 30; % conversion to samples for extracting spikes (st.sp is in seconds) - each ms is 30 samples
		prespike_restriction = prespike_restriction .*30;
		postspike_restriction = postspike_restriction .*30;
		
		all_spikeTimes = gwfparams.spikeTimes;
		all_spikeClusters = gwfparams.spikeClusters;
		
		if ~isequal(take_burst_spikes, 'singles')
		gwfparams.spikeTimes = []; % clearing so spikes with specific ISIs can be assigned
		gwfparams.spikeClusters = [];
		end
		
		ISIs_before = [];
		ISIs_after = [];
		bursts_list = [];
		this_ISI_after_list = [];

		for i = 4:(numel(all_spikeTimes)-5) % iterating through each spike and checking the interval before and after
			
			this_ISI_before = all_spikeTimes(i) - all_spikeTimes(i - 1); % time (ISI) between i-th spike and the spike before
			this_ISI_2_before = all_spikeTimes(i - 1) - all_spikeTimes(i - 2); % time (ISI) between the spike after the i-th spike and the spike after that
			this_ISI_3_before = all_spikeTimes(i - 2) - all_spikeTimes(i - 3); % time (ISI) between the spike after the i-th spike and the spike after that
			this_ISI_after = all_spikeTimes(i + 1) - all_spikeTimes(i); % time (ISI) between i-th spike and the spike after
			this_ISI_2_after_the_first = all_spikeTimes(i + 2) - all_spikeTimes(i + 1); 
			this_ISI_3_after_the_first = all_spikeTimes(i + 3) - all_spikeTimes(i + 2);
			this_ISI_4_after_the_first = all_spikeTimes(i + 4) - all_spikeTimes(i + 3);
			this_ISI_5_after_the_first = all_spikeTimes(i + 5) - all_spikeTimes(i + 4); % did this to check if there are bursts of 4 or 5 spikes
			
			this_ISI_after_list = [this_ISI_after_list this_ISI_after];

			if ... % all the below conditions are met - the current spike will be used when calculating spike-triggered averages % 'all' doesn't restrict the interval before the spike such taht first spikes are also in this criterion
				isequal(take_burst_spikes, 'all') && ...
				((this_ISI_before > ISI_range(1)   && ... % time before the i-th spike is larger than the lower bound of specified ISI
				this_ISI_after > ISI_range(1) && ... % same but for the spike after
				this_ISI_after < ISI_range(2)) ... % same but for spike after
				|| ...
				(this_ISI_before > ISI_range(1)   && ... % conditions to include last spikes of bursts
				this_ISI_before < ISI_range(2) && ... 
				this_ISI_after > ISI_range(1))) ... 

				ISIs_before = [ISIs_before this_ISI_before];
				ISIs_after = [ISIs_after this_ISI_after];
				
				gwfparams.spikeTimes =  [gwfparams.spikeTimes all_spikeTimes(i)];
				gwfparams.spikeClusters = [gwfparams.spikeClusters all_spikeClusters(i)];
			elseif ...
				isequal(take_burst_spikes, 'firsts') && ... 
				this_ISI_before > prespike_restriction && ... % if no spike are present in a defined amount of ms before the spike (50ms = 1500)
				this_ISI_after > ISI_range(1) && ...
				this_ISI_after < ISI_range(2) %&& ... 
				
				if isequal(take_burst_type, 'all') % get the correct type of bursts
					ISIs_before = [ISIs_before, this_ISI_before];
					if this_ISI_2_after_the_first < ISI_range(2) && this_ISI_3_after_the_first > postspike_restriction   % only get the third spike if it is a triple burst
						ISIs_after = [ISIs_after; [this_ISI_after, this_ISI_2_after_the_first]];
						gwfparams.spikeTimes =  [gwfparams.spikeTimes all_spikeTimes(i)];
						gwfparams.spikeClusters = [gwfparams.spikeClusters all_spikeClusters(i)];
						gwfparams.ISIs_after = ISIs_after;
					elseif this_ISI_2_after_the_first > postspike_restriction
						ISIs_after = [ISIs_after; [this_ISI_after, 0]];
						gwfparams.spikeTimes =  [gwfparams.spikeTimes all_spikeTimes(i)];
						gwfparams.spikeClusters = [gwfparams.spikeClusters all_spikeClusters(i)];
						gwfparams.ISIs_after = ISIs_after;
					end
				elseif isequal(take_burst_type, 'triple') && this_ISI_2_after_the_first < ISI_range(2) && this_ISI_3_after_the_first > postspike_restriction
					ISIs_before = [ISIs_before, this_ISI_before];
					ISIs_after = [ISIs_after; [this_ISI_after, this_ISI_2_after_the_first]];
					
					gwfparams.spikeTimes = [gwfparams.spikeTimes all_spikeTimes(i)];
					gwfparams.spikeClusters = [gwfparams.spikeClusters all_spikeClusters(i)];
					gwfparams.ISIs_after = ISIs_after;	
				
				elseif isequal(take_burst_type, 'double') && this_ISI_2_after_the_first > prespike_restriction
					ISIs_before = [ISIs_before, this_ISI_before];
					ISIs_after = [ISIs_after; [this_ISI_after, 0]]; % 0 index is for compatibility
					
					gwfparams.spikeTimes =  [gwfparams.spikeTimes all_spikeTimes(i)];
					gwfparams.spikeClusters = [gwfparams.spikeClusters all_spikeClusters(i)];
					gwfparams.ISIs_after = ISIs_after;
				end 
				
				bursts_list = [bursts_list; [this_ISI_after, 0, 0, 0, 0]]; % make list with all bursts 
				if this_ISI_2_after_the_first < ISI_range(2) && ~this_ISI_3_after_the_first < ISI_range(2)
					bursts_list = [bursts_list; [this_ISI_after this_ISI_2_after_the_first, 0, 0, 0]];
				elseif this_ISI_2_after_the_first < ISI_range(2) && this_ISI_3_after_the_first < ISI_range(2) && ~this_ISI_4_after_the_first < ISI_range(2)
					bursts_list = [bursts_list; [this_ISI_after this_ISI_2_after_the_first, this_ISI_3_after_the_first, 0, 0]];
				elseif this_ISI_2_after_the_first < ISI_range(2) && this_ISI_3_after_the_first < ISI_range(2) && this_ISI_4_after_the_first < ISI_range(2) && ~this_ISI_5_after_the_first < ISI_range(2)
					bursts_list = [bursts_list; [this_ISI_after this_ISI_2_after_the_first, this_ISI_3_after_the_first, this_ISI_4_after_the_first, 0]];
				elseif this_ISI_2_after_the_first < ISI_range(2) && this_ISI_3_after_the_first < ISI_range(2) && this_ISI_4_after_the_first < ISI_range(2) && this_ISI_5_after_the_first < ISI_range(2) 
					bursts_list = [bursts_list; [this_ISI_after this_ISI_2_after_the_first, this_ISI_3_after_the_first, this_ISI_4_after_the_first, this_ISI_5_after_the_first]];
				end
			elseif ...
				isequal(take_burst_spikes, 'seconds') && ...
				this_ISI_before > ISI_range(1) && ... 
				this_ISI_after > ISI_range(1) && ...
				this_ISI_before < ISI_range(2) && ...
				this_ISI_2_before > prespike_restriction && ... % if the spike before is the first spike of the burst
				this_ISI_after < ISI_range(2) % this can be useful when you want to look at seconds spikes of bursts with at least 3 spikes 

				ISIs_before = [ISIs_before this_ISI_before];
				ISIs_after = [ISIs_after this_ISI_after];
				
				gwfparams.spikeTimes =  [gwfparams.spikeTimes all_spikeTimes(i)];
				gwfparams.spikeClusters = [gwfparams.spikeClusters all_spikeClusters(i)];
			elseif ...
				isequal(take_burst_spikes, 'thirds') && ...
				this_ISI_before > ISI_range(1) && ...
				this_ISI_before < ISI_range(2) && ...
				this_ISI_2_before > ISI_range(1) && ...
				this_ISI_2_before < ISI_range(2) && ...
				this_ISI_3_before > prespike_restriction %&& ... 
% 				this_ISI_after > postspike_restriction  

				ISIs_before = [ISIs_before this_ISI_before];
				ISIs_after = [ISIs_after this_ISI_after];
				
				gwfparams.spikeTimes =  [gwfparams.spikeTimes all_spikeTimes(i)];
				gwfparams.spikeClusters = [gwfparams.spikeClusters all_spikeClusters(i)];
			elseif ...
				isequal(take_burst_spikes, 'inside') && ...
				this_ISI_before > ISI_range(1) && ... % time before the i-th spike is larger than the lower bound of specified ISI
				this_ISI_after > ISI_range(1) && ... % same but for the spike after
				this_ISI_before < ISI_range(2) && ... % time before the i-th spike is smaller than the upper bound of specified ISI
				this_ISI_after < ISI_range(2) % same but for spike after
				
				ISIs_before = [ISIs_before this_ISI_before];
				ISIs_after = [ISIs_after this_ISI_after];
				
				gwfparams.spikeTimes =  [gwfparams.spikeTimes all_spikeTimes(i)];
				gwfparams.spikeClusters = [gwfparams.spikeClusters all_spikeClusters(i)];
			elseif ...
				isequal(take_burst_spikes, 'singles') && ... % criterion is to not be in the burst spike interval 
				this_ISI_before < prespike_restriction(1) || this_ISI_before > prespike_restriction(end)   % && ... 
				%this_ISI_after ~= prespike_restriction
				
				ISIs_before = [ISIs_before this_ISI_before];
% 				ISIs_after = [ISIs_after this_ISI_after]
				
% 				gwfparams.spikeTimes =  [gwfparams.spikeTimes all_spikeTimes(i)];
% 				gwfparams.spikeClusters = [gwfparams.spikeClusters all_spikeClusters(i)];
				gwfparams.spikeTimes(i) =  0; % marking all spikes that are part of bursts to leave singlets
				gwfparams.spikeClusters(i) = 0;
			end
		end		
		
		if isequal(take_burst_spikes, 'singles') % when loop is done, get new spike list with only single spikes, take out all spikes, whcih have been marked as burst spikes
			gwfparams.spikeTimes(gwfparams.spikeTimes==0) = [];
			gwfparams.spikeClusters(gwfparams.spikeClusters==0) = [];
		end
		
		if size(gwfparams.spikeTimes,1) == 1 % clean up formatting which can vary between conditions
			gwfparams.spikeTimes = gwfparams.spikeTimes';
			gwfparams.spikeClusters = gwfparams.spikeClusters';
		end
		
% 		if ~isequal(take_burst_spikes, 'singles')
% 			fprintf(['Number of spikes for ' num2str(nrn_ID) ' in ISI interval: ' num2str(ISI_range(1)/30) 'ms to ' num2str(ISI_range(2)/30) 'ms is: ' num2str(numel(gwfparams.spikeTimes)) '\n'])
% 		else
% 			fprintf(['Number of single spikes for nrn ' num2str(nrn_ID) ' is: ' num2str(numel(gwfparams.spikeTimes)) '\n'])
% 		end
% 		gwfparams.ISIs_after = ISIs_after;	
	end
	
% 	if numel(gwfparams.spikeTimes) < gwfparams.nWf          % to get all spikes of a particular neuron (mostly for AIS initiation)
% 		gwfparams.nWf = numel(gwfparams.spikeTimes);        % just take all spikes if there are fewer than the number requested
% % 		fprintf(['Number of requested spikes is too high, taking only ' num2str(numel(gwfparams.spikeTimes)) '\n'])
% 	end
	
% truncate the start of spike times so you don't index negatively when the unit has very early spike times 
% this happens mostly when you take 'all' spikes of 'firsts' or especially 'singles'
% gwfparams.spikeTimes = gwfparams.spikeTimes(gwfparams.spikeTimes>(gwfparams.wfWin(2)*12)); %need to multiply window since it is in AP stime stamps here
% gwfparams.spikeClusters(size(gwfparams.spikeTimes,1)+1:end) = []; 
