function [soma_trace, above_soma_trace,high_trace, correlation_matrix] = get_amplitude_above_soma(info,sp,select_height,wfMean_mean,isi,w_drive,rec_date,nrn,NP_probe)

% get the channel that corresponds to the wanted height
% [~, chan] = min(abs(sp.ycoords - select_height));

chan = select_height;

chan_min = min(wfMean_mean, [], 2);
min_chan = find(chan_min == min(chan_min));

if isequal(NP_probe, 'UHD') || isequal(NP_probe, 'UHD_4')
soma_trace = wfMean_mean(chan, :);	
above_soma_trace = wfMean_mean(chan+20,:);
high_trace = wfMean_mean(chan+30,:);
else

soma_trace = wfMean_mean(min_chan, :);
try
above_soma_trace = wfMean_mean(min_chan+8,:);

high_trace = wfMean_mean(min_chan+15,:);
catch
	
above_soma_trace = [];
high_trace = [];
end
end


%% plot to get the raw signal from the somatic channel upwards
% after spike time
% % 
figure_dir = strcat(w_drive,'\',rec_date,'\','nrn_',num2str(nrn),'\','drift_correction','\',date);
if ~exist(figure_dir,'dir')
	mkdir(figure_dir)
end

chan_tab = [];
fig = figure; hold on

switch NP_probe
	
	case '1.0'
	before_soma = min_chan-30;
    
% plot raw traces from soma upwards
for i = 1:80

	chan = before_soma+i;
	x = wfMean_mean(chan,:);
	chan_tab = [chan_tab; x];
	x1 = bsxfun(@plus,x,10*i);
	x1 = medfilt1(x1,3);
	plot(x1,'Color','k','LineWidth',1.2)
   
	
end

	xlim([140 300]); ylim([72 580]);
	xticks([140:30:400]);
	yticks([0:50:850]); yticklabels([-600:100:2000]);
	xticklabels([-2:4]);
	ylabel('distance from soma (\mum)')
	xlabel('time (ms)')

	case 'UHD' 
	before_soma = min_chan-20;	
	for i = 1:40

	chan = before_soma+i;
	x = wfMean_mean(chan,:);
	
	x1 = bsxfun(@plus,x,10*i);
	x1 = medfilt1(x1,3);
	plot(x1,'Color','k','LineWidth',1.2)
	end
	
	case 'UHD_4'
		before_soma = min_chan-20;	
	for i = 1:40

	chan = before_soma+i;
	x = wfMean_mean(chan,:);
	
	x1 = bsxfun(@plus,x,10*i);
	x1 = medfilt1(x1,3);
	plot(x1,'Color','k','LineWidth',1.2)
	end

end
		
	if ~(info.analyse_LED_spikes) && ~isempty(info.analyse_running)
	if isequal(info.analyse_running,'stationary') 
exportgraphics(fig,[figure_dir,'\','Raw_trace_stationary_',num2str(isi),'_ISI.jpg'], 'Resolution',600)
	elseif isequal(info.analyse_running, 'running')
exportgraphics(fig,[figure_dir,'\','Raw_trace_running_',num2str(isi),'_ISI.jpg'], 'Resolution',600)
	
	end
	
	elseif info.analyse_LED_spikes && ~isempty(info.analyse_running)
		
		if isequal(info.analyse_running,'stationary') 
exportgraphics(fig,[figure_dir,'\','Raw_trace_stationary_LED_',num2str(isi),'_ISI.jpg'], 'Resolution',600)
	elseif isequal(info.analyse_running, 'running')
exportgraphics(fig,[figure_dir,'\','Raw_trace_running_LED',num2str(isi),'_ISI.jpg'], 'Resolution',600)
		
		end	
	
exportgraphics(fig,[figure_dir,'\','Raw_trace_',num2str(isi),'_ISI.jpg'], 'Resolution',600)
		
			
    end

% compute correlation coefficient between channels 

% restricting signal to -1ms before spike time and + 2ms after spike time, only taking timepoints around spike time (which is at
% time 200)
chan_tab(:,1:170) = [];
chan_tab(:,90:end) = [];

data_transposed = chan_tab';

correlation_matrix = corrcoef(data_transposed);

% figure;
% heatmap(correlation_matrix, 'Colormap', parula, 'ColorbarVisible', 'on', 'XLabel', 'Channel', 'YLabel', 'Channel');

end