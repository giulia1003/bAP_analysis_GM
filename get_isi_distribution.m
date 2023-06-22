function [isi_properties,bursty_neurons,non_bursty_neurons] = get_isi_distribution(analyse_nrn,sp,infodir,plot_figs)

% plot ISI intervals before and after spikes 
% for every neuron
% on log scale (useful for outliers with long isi)
% set threshold for bursts (e.g. 1-5 ms)

isi_properties = []; % structure to collect data
perc_bursts_tab = []; 

%%

prespike_isi = [];
postspike_isi = [];

% units_bursts = [analyse_nrn(1):1:analyse_nrn(end)];

for j = 1:length(analyse_nrn)
	
nrn = analyse_nrn(j);

nrn_ID = double(sp.cids(nrn));
spiketimes = ceil(sp.st(sp.clu==nrn_ID)*30000);  % timestamps

isi_list = diff(spiketimes); isi_list = isi_list./30; % isi in ms

for i = 2:(numel(spiketimes)-2)
	
	prespike_ms = isi_list(i-1); 
	postspike_ms = isi_list(i+1);
	
	prespike_isi = [prespike_isi; prespike_ms];
	postspike_isi = [postspike_isi; postspike_ms];
	
end

if plot_figs
	
fig1 = figure; hold on

loglog(postspike_isi,prespike_isi,'o'), hold on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([10^0 10^6])
ylim([10^0 10^6])
% yline(10^1/2, '--r')
yline(10^1/2, '--r','LineWidth',1.5)
box off
ylabel('interval before (ms)')
xlabel('interval after (ms)')
title(strcat('Interspike intervals nrn: ',num2str(nrn)))

cd(infodir)

saveas(fig1,strcat('log scale isi nrn_',num2str(nrn),'.fig'));
close(fig1)
end

% plot neuron AUTOCORRELOGRAM 

binsize =  0.001;  % in s
maxlag = 0.05; % in s

[ac,xbin] = acf(spiketimes,binsize,maxlag,nrn,plot_figs);


bursts_sp = find(isi_list <= 10);
perc_bursts = (numel(bursts_sp)*100)/numel(spiketimes);

perc_bursts_tab = [perc_bursts_tab; perc_bursts];

isi_properties(j).nrn = nrn;
isi_properties(j).isi_distribution = isi_list;
isi_properties(j).percentage_of_bursts = perc_bursts;
isi_properties(j).ac_estimate = ac;

end

% figure,histogram(perc_bursts_tab,10)

thresh = multithresh(perc_bursts_tab,2);
thresh = max(thresh); % get threshold for bursty neurons

id_bursty_neurons = find(perc_bursts_tab >= thresh);

bursty_neurons = analyse_nrn(id_bursty_neurons);
non_bursty_neurons = analyse_nrn; non_bursty_neurons(id_bursty_neurons) = [];



end