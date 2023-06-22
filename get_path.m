function [mouse,session,linux,rec_date,insertion_angle_name,myKsDir,procDataDir_General,infodir] = get_path(data_path,NP_probe)

% 
% rep_path = 'C:\Users\tex_analysis\Desktop\Dendrites_project\NP2_bAP_analysis_GM_v20220104';
% addpath(genpath(rep_path)) % repository includes some functions from Nick's Spikes and npy-matlab,
% some are modified from the orginal

linux = 0;            %if on the linux machine; otherwise 0
fprintf(string(strcat('\nSetting configurations for probe:', {' '}, NP_probe,{' '},'- CHECK if correct','.\n\n')))


split_data_path = split(data_path, '\');
 % the naming convention is not always the same because different people + different probes
 % so wI added these lines
if size(split_data_path,1) == 5
session = split_data_path{5};
mouse = split_data_path{3};
rec_date = session(7:14);
elseif size(split_data_path,1) == 4
session = split_data_path{4};
mouse = split_data_path{2}
rec_date = session(7:14);
elseif size(split_data_path,1) == 6 
	mouse = split_data_path{4};
    rec_date = split_data_path{5};
	session = split_data_path{6};
elseif size(split_data_path,1) == 7
mouse = split_data_path{4};
rec_date = split_data_path{5};
	session = split_data_path{7};
	session = session(1:end-3);
end
insertion_angle_num = data_path(end-1:end);
insertion_angle_num = insertion_angle_num(regexp(insertion_angle_num, '\d'));
insertion_angle_name = strcat('angle_', insertion_angle_num);

fprintf(string(strcat('\nNow analysing data of mouse', {' '}, mouse, ' with angle', {' '}, insertion_angle_num, '.\n\n')))
procDataDir_General = strcat(data_path, '/general');                     % this is where variables such as cortical channels will be stored

myKsDir = data_path; %% set Kilosort path and figure output directory

infodir = strcat(data_path, '\info_nrn\');
% if ~exist('infodir','dir')
% mkdir(infodir);
% end

disp(['Kilosort output directory: ', myKsDir]) % print out
fprintf('\n')

end