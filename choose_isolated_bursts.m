function  [gwfparams,this_ISI_before_list] = choose_isolated_bursts(info,ISI_range, prespike_restriction, postspike_restriction, take_burst_type, gwfparams, nrn_ID)

ISI_range = ISI_range .* 30; % conversion to samples for extracting spikes - each ms is 30 samples
% prespike_restriction_1: lower limit of ISI permitted (e.g. 1 ms)
% prespike_restriction_2: higher limit of ISI permitted (e.g. 10 ms)

prespike_restriction_1 = prespike_restriction(1)*30; % convert extremities of ISI range from ms to timestamps
prespike_restriction_2 = prespike_restriction(end)*30;

prespike_restriction_first_spikelets = prespike_restriction_2*10;   % first spikelet should be 100 ms from any previous spike
prespike_restriction_third_spikelets = prespike_restriction(5)*30;

all_spikeTimes = gwfparams.spikeTimes;
all_spikeClusters = gwfparams.spikeClusters;
		

this_ISI_list = diff(all_spikeTimes);
this_ISI_before_list = []; 

if ~(info.single_spike_maps) && ~(isequal(take_burst_type, 'all'))
	
% take only spikes fired at least 30 ms after the previous spikes that also have
% the next spike within 1-10 ms
if isequal(take_burst_type,'firsts')

all_spikeTimes = cat(1,0,all_spikeTimes,0);
all_spikeClusters = cat(2,0,all_spikeClusters,0);

gwfparams.spikeTimes = cat(1,0,gwfparams.spikeTimes,0);
gwfparams.spikeClusters = cat(2,0,gwfparams.spikeClusters,0);

for i = 2:size(all_spikeTimes,1)-1 % iterating through each spike and checking the interval before and after
			
			this_ISI_before = all_spikeTimes(i) - all_spikeTimes(i - 1); % time (ISI) between i-th spike and the spike before
			this_ISI_after = all_spikeTimes(i+1) - all_spikeTimes(i);
			this_ISI_before_list = [this_ISI_before_list this_ISI_before];
			
		if this_ISI_before < prespike_restriction_first_spikelets || this_ISI_after > prespike_restriction_2
			% if the interval between two spikes does not fall into the specified ISI
				
				gwfparams.spikeTimes(i) =  0; % marking all spikes that are part of bursts to leave singlets
				gwfparams.spikeClusters(i) = 0;
	 end
end


gwfparams.spikeTimes(gwfparams.spikeTimes==0) = [];

gwfparams.spikeClusters(gwfparams.spikeClusters==0) = [];

%%
elseif isequal(take_burst_type,'seconds')  

% select only spikes with ISI 1-10ms before and after
all_spikeTimes = cat(1,0,0,all_spikeTimes);
all_spikeClusters = cat(2,0,0,all_spikeClusters);

gwfparams.spikeTimes = cat(1,0,0,gwfparams.spikeTimes);
gwfparams.spikeClusters = cat(2,0,0,gwfparams.spikeClusters);

for i = 3:size(all_spikeTimes,1) % iterating through each spike and checking the interval before and after
			
	% this_ISI_before_1: time between second spikelet and first spikelet
	% this_ISI_before_2: time between first spieklet and whatever came before 
	
			this_ISI_before_1 = all_spikeTimes(i) - all_spikeTimes(i - 1); % time (ISI) between i-th spike and the spike before
			this_ISI_before_2 = all_spikeTimes(i-1) - all_spikeTimes(i - 2);

			this_ISI_before_list = [this_ISI_before_list this_ISI_before_1];
			
		if this_ISI_before_1 > prespike_restriction_2 || this_ISI_before_2 < prespike_restriction_first_spikelets
			% if the previous spike was more than 10 ms away
			% and the spike before the previous one is NOT at least 60 ms away
			
			% discard spike(i)
				
				gwfparams.spikeTimes (i) =  0; % marking all spikes that are part of bursts to leave singlets
				gwfparams.spikeClusters(i) = 0;
			end
end

% when loop is done, get new spike list with only single spikes,
% take out all spikes whcih have been marked as non-burst spikes

gwfparams.spikeTimes(gwfparams.spikeTimes==0) = [];
gwfparams.spikeClusters(gwfparams.spikeClusters==0) = [];

	
	 	
elseif isequal(take_burst_type,'thirds')
	% take spikes of a burst (1-5 ms away from previous spike) except for the first one	first_sp_indices = find(this_ISI_list > 300);

all_spikeTimes = cat(1,0,0,0,all_spikeTimes);
all_spikeClusters = cat(2,0,0,0,all_spikeClusters);

gwfparams.spikeTimes = cat(1,0,0,0,gwfparams.spikeTimes);
gwfparams.spikeClusters = cat(2,0,0,0,gwfparams.spikeClusters);

for i = 4:size(all_spikeTimes,1) % iterating through each spike and checking the interval before and after
			
	% this_ISI_before_1: time between third spikelet and second spikelet
	% this_ISI_before_2: time between second spikelet and first spikelet
	% this_ISI_before_3: time between first spikelet and whatever came before 
	
			this_ISI_before_1 = all_spikeTimes(i) - all_spikeTimes(i - 1); % time (ISI) between i-th spike and the spike before
			this_ISI_before_2 = all_spikeTimes(i-1) - all_spikeTimes(i - 2);
			this_ISI_before_3 = all_spikeTimes(i-2) - all_spikeTimes(i - 3);

			this_ISI_before_list = [this_ISI_before_list this_ISI_before_1];
			
		if this_ISI_before_1 > prespike_restriction_third_spikelets || this_ISI_before_2 > prespike_restriction_2 || this_ISI_before_3 < prespike_restriction_first_spikelets
			% if the previous spike was more than 10 ms away
			% and the one before (first spikelet) was more than 10 ms away
			% and the spike before the previous one is NOT at least 60 ms away
			
			% discard spike(i)
				
				gwfparams.spikeTimes (i) =  0; % marking all spikes that are part of bursts to leave singlets
				gwfparams.spikeClusters(i) = 0;
			end
end

% when loop is done, get new spike list with only single spikes,
% take out all spikes whcih have been marked as non-burst spikes

gwfparams.spikeTimes(gwfparams.spikeTimes==0) = [];
gwfparams.spikeClusters(gwfparams.spikeClusters==0) = [];

else
	
end
	%% TO GET FIRST, SECOND AND THIRD SPIKELETS IN BURSTS 
	% AND MAKE SINGLE-SPIKE VOLTAGE MAPS FOR EVERY NEURON	
	
elseif info.single_spike_maps && isequal(take_burst_type, 'all')
	
	% to get a snapshot of entire bursts, get only first spikelet	
% 1ST SPIKELET

gwfparams.spikeTimes_1 = gwfparams.spikeTimes;
gwfparams.spikeClusters_1 = gwfparams.spikeClusters;

all_spikeTimes = cat(1,0,all_spikeTimes,0);
all_spikeClusters = cat(2,0,all_spikeClusters,0);

gwfparams.spikeTimes_1 = cat(1,0,gwfparams.spikeTimes_1,0);
gwfparams.spikeClusters_1 = cat(2,0,gwfparams.spikeClusters_1,0);


for i = 2:size(all_spikeTimes,1)-1 % iterating through each spike and checking the interval before and after
			
			this_ISI_before = all_spikeTimes(i) - all_spikeTimes(i - 1); % time (ISI) between i-th spike and the spike before
			this_ISI_after = all_spikeTimes(i+1) - all_spikeTimes(i);
			this_ISI_before_list = [this_ISI_before_list this_ISI_before];
			
		if this_ISI_before < prespike_restriction_first_spikelets || this_ISI_after > prespike_restriction_2
			% if the interval between two spikes does not fall into the specified ISI
				
				gwfparams.spikeTimes_1(i) =  0; % marking all spikes that are part of bursts to leave singlets
				gwfparams.spikeClusters_1(i) = 0;
	 end
end


gwfparams.spikeTimes_1(gwfparams.spikeTimes_1==0) = [];

gwfparams.spikeClusters_1(gwfparams.spikeClusters_1 ==0) = [];

% % 2ND SPIKELET 
% % select only spikes with ISI 1-10ms before and after
% all_spikeTimes = cat(1,0,0,all_spikeTimes);
% all_spikeClusters = cat(2,0,0,all_spikeClusters);
% 
% gwfparams.spikeTimes_2 = gwfparams.spikeTimes;
% gwfparams.spikeClusters_2 = gwfparams.spikeClusters;
% 
% gwfparams.spikeTimes_2 = cat(1,0,0,gwfparams.spikeTimes_2);
% gwfparams.spikeClusters_2 = cat(2,0,0,gwfparams.spikeClusters_2);
% 
% for i = 3:size(all_spikeTimes,1) % iterating through each spike and checking the interval before and after
% 			
% 	% this_ISI_before_1: time between second spikelet and first spikelet
% 	% this_ISI_before_2: time between first spieklet and whatever came before 
% 	
% 			this_ISI_before_1 = all_spikeTimes(i) - all_spikeTimes(i - 1); % time (ISI) between i-th spike and the spike before
% 			this_ISI_before_2 = all_spikeTimes(i-1) - all_spikeTimes(i - 2);
% 
% 			this_ISI_before_list = [this_ISI_before_list this_ISI_before_1];
% 			
% 		if this_ISI_before_1 > prespike_restriction_2 || this_ISI_before_2 < prespike_restriction_first_spikelets
% 			% if the previous spike was more than 10 ms away
% 			% and the spike before the previous one is NOT at least 60 ms away
% 			
% 			% discard spike(i)
% 				
% 				gwfparams.spikeTimes_2 (i) =  0; % marking all spikes that are part of bursts to leave singlets
% 				gwfparams.spikeClusters_2(i) = 0;
% 			end
% end
% 
% % when loop is done, get new spike list with only single spikes,
% % take out all spikes whcih have been marked as non-burst spikes
% 
% gwfparams.spikeTimes_2(gwfparams.spikeTimes_2 == 0) = [];
% gwfparams.spikeClusters_2(gwfparams.spikeClusters_2==0) = [];
% 
% 	
% % 3RD SPIKELET
% 
% gwfparams.spikeTimes_3 = gwfparams.spikeTimes;
% gwfparams.spikeClusters_3 = gwfparams.spikeClusters;	
% 
% all_spikeTimes = cat(1,0,0,0,all_spikeTimes);
% all_spikeClusters = cat(2,0,0,0,all_spikeClusters);
% 
% gwfparams.spikeTimes_3 = cat(1,0,0,0,gwfparams.spikeTimes_3);
% gwfparams.spikeClusters_3 = cat(2,0,0,0,gwfparams.spikeClusters_3);
% 
% for i = 4:size(all_spikeTimes,1) % iterating through each spike and checking the interval before and after
% 			
% 	% this_ISI_before_1: time between third spikelet and second spikelet
% 	% this_ISI_before_2: time between second spikelet and first spikelet
% 	% this_ISI_before_3: time between first spikelet and whatever came before 
% 	
% 			this_ISI_before_1 = all_spikeTimes(i) - all_spikeTimes(i - 1); % time (ISI) between i-th spike and the spike before
% 			this_ISI_before_2 = all_spikeTimes(i-1) - all_spikeTimes(i - 2);
% 			this_ISI_before_3 = all_spikeTimes(i-2) - all_spikeTimes(i - 3);
% 
% 			this_ISI_before_list = [this_ISI_before_list this_ISI_before_1];
% 			
% 		if this_ISI_before_1 > prespike_restriction_third_spikelets || this_ISI_before_2 > prespike_restriction_2 || this_ISI_before_3 < prespike_restriction_first_spikelets
% 			% if the previous spike was more than 10 ms away
% 			% and the one before (first spikelet) was more than 10 ms away
% 			% and the spike before the previous one is NOT at least 60 ms away
% 			
% 			% discard spike(i)
% 				
% 				gwfparams.spikeTimes_3 (i) =  0; % marking all spikes that are part of bursts to leave singlets
% 				gwfparams.spikeClusters_3(i) = 0;
% 			end
% end
% 
% % when loop is done, get new spike list with only single spikes,
% % take out all spikes whcih have been marked as non-burst spikes
% 
% gwfparams.spikeTimes_3(gwfparams.spikeTimes_3==0) = [];
% gwfparams.spikeClusters_3(gwfparams.spikeClusters_3==0) = [];
% 
% gwfparams.spikeTimes = cat(1,gwfparams.spikeTimes_1, gwfparams.spikeTimes_2,gwfparams.spikeTimes_3);
% gwfparams.spikeClusters = cat(2, gwfparams.spikeClusters_1, gwfparams.spikeClusters_2, gwfparams.spikeClusters_3);

% gwfparams.spikeTimes_1 = [];
% gwfparams.spikeTimes_2 = [];
% gwfparams.spikeTimes_3 = [];

% gwfparams.spikeClusters_1 = [];
% gwfparams.spikeClusters_2 = [];
% gwfparams.spikeClusters_3 = [];

gwfparams.spikeTimes = gwfparams.spikeTimes_1;
gwfparams.spikeClusters_1 = gwfparams.spikeClusters_1;

end
	

end % end of function




	
		
		
		
		
		
		
		
		
		
		
		
		
