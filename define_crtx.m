function [bAP_nrn_ind,templateYpos, crtx_bottCh, crtx_topCh,driftEvents] = define_crtx(myKsDir, sp, linux, procDataDir_General, previously_run)

    if ~previously_run
% 		try		%FK, 27.06.21, KS 3 data loading issues
			[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
% 		catch
% 			warning('not plotting drift map')
% 			warning('pc_features.npy and pc_feature_ind.npy were NOT found in the ksDir, this may be due to spike sorting with KS 3')
% 		end


        % adding a fake spike to 'spikeDepths' as a simple way to plot along
        % the whole probe - DOESNT WORK

        spikeTimes([end+1 end+2]) = [0 0]; % at the start of the recording
        spikeAmps([end+1 end+2]) = [150 150]; % approximate amplitude for high amp neurons
        spikeDepths([end+1 end+2]) = [0 3840]; % lowest and highest channel on probe

        if linux
            drift_fig = figure('Units','normalized', 'position', [0 1 1/3 1/2]);
        else
            drift_fig = figure;
        end

try
       [driftEvents] = plotDriftmap2(spikeTimes, spikeAmps, spikeDepths, 'show');
         plotDriftmap(spikeTimes, spikeAmps, spikeDepths, 'show');
		try
        saveas(drift_fig, [procDataDir_General, '/drift']);
		catch
		end
        % basic quantification of spiking plot
catch
	driftEvents = [];
end
        depthBins = 0:40:3840;
        ampBins = 0:30:min(max(spikeAmps),800);
        recordingDur = sp.st(end);


        [pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
        plotWFampCDFs(pdfs, cdfs, ampBins, depthBins, linux, procDataDir_General);


        % Plotting some basics about LFPs
       
 
        plotLFP(myKsDir, linux, procDataDir_General) %will plot lfp for *.lf.bin file in myKsDir with parameters set in function
        % Computing some useful details about spikes/neurons (like depths)

        [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
            templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);



        %% IMPORTANT!!! this section of the code will ask you to define cortex in console window
        % it also calculates depths and assigns cortical units

        fprintf('\nBased on the plots above, please provide location of cortex (approximately 1000 um thick):\n\n')

        crtx_bottCh = input('Please enter the depth for bottom bit of cortex: ')/10;

        fprintf(['Channel ' num2str(crtx_bottCh) ' assigned as lower boundary of cortex.\n\n']);

        crtx_topCh = input('Please enter the depth for pia surface: ')/10;

        fprintf(['Channel ' num2str(crtx_topCh) ' assigned as upper boundary of cortex.\n\n']);


        if not (size(templateYpos,1) == max(size(sp.cids)))

            templateYpos = templateYpos(1:max(size(sp.cids)));
            fprintf(2, '\nWARNING: Truncating templateYpos variable\n')

        end

    elseif previously_run
		
        crtx_topCh = crtx_topCh_prev;
        crtx_bottCh = crtx_bottCh_prev;
		
% 		[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
%          [driftEvents] = plotDriftmap2(spikeTimes, spikeAmps, spikeDepths, 'show');
        [~, ~, templateYpos, ~, ~, ~, ~] = templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
		
		if not (size(templateYpos,1) == size(sp.cids,1))

            templateYpos = templateYpos(1:size(sp.cids,1));
            fprintf(2, '\nWARNING: Truncating templateYpos variable\n')

        end
    end
    
    
    bAP_nrn_ID = sp.cids(templateYpos < 10*crtx_topCh & templateYpos > 10*crtx_bottCh);

    bAP_nrn_ind = zeros(size(bAP_nrn_ID));

    all_neurons = unique(sp.clu);

    for i = 1: length(bAP_nrn_ind)
% 		try  % FK, 27.06.21, implemented this because of never before seen error that some neurons wouldn't show up in all_neurons
			bAP_nrn_ind(i) = find(all_neurons == bAP_nrn_ID(i));
% 		catch 
% 			bAP_nrn_ind(i) = 0;
% 			warning(['couldn''t find wanted neuron ',  num2str(bAP_nrn_ID(i)), ' in neuron IDs - skipping!'])
% 			continue
% 		end
	end
% 	bAP_nrn_ind = bAP_nrn_ind(find(bAP_nrn_ind)); % also 27.06.21
end

