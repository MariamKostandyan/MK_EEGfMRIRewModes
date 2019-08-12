% batch for running marsbar extractin betas from roi
% modified from http://marsbar.sourceforge.net/faq.html#marsbar_batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% requires function "keep" (put it into spm folder) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
marsbar('on');

% rois
rois={'inc_con_-27_-64_44_roi.mat'}


rootdir = 'D:\MK_MREEGRewModes_2016\fMRI\FirstLevel_Mixed_NoDeriv\Concatenated_NoBlocks_SepCond\';
cd(rootdir);

    
for r=1:length(rois) 
        roi=rois{r}
        roi_file = ['D:\MK_MREEGRewModes_2016\fMRI\ROIAnanlysis\ROIs_Inc_vs_Con\v6_SPMcluster\roi_clusterINCvsCON\', roi];

  
    
        
sub={'sub09', 'sub11', 'sub12', 'sub13', 'sub14', 'sub15','sub16','sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24','sub25','sub26','sub27','sub28','sub29','sub32','sub33','sub34','sub35','sub36'};

        
% subject loop        
for z=1:length(sub)
    subject=sub{z}
    clear ('SPM', 'E', 'D', 'xCon', 'xSPM', 'Bcov', 'Y');
    spm_name = ['D:\MK_MREEGRewModes_2016\fMRI\FirstLevel_Mixed_NoDeriv\Concatenated_NoBlocks_SepCond\', subject, '\model01\spm.mat\'];
    
 
% Make marsbar design object
D  = mardo(spm_name);

% % Set fmristat AR modelling
 D = autocorr(D, 'fmristat', 2);

% Make marsbar ROI object
R  = maroi(roi_file);

% Fetch data into marsbar data object
Y  = get_marsy(R, D, 'mean');

% Get contrasts from original design
xCon = get_contrasts(D);

% Estimate design on ROI data
E = estimate(D, Y);

% Put contrasts from original design back into design object
E = set_contrasts(E, xCon);

% Get definitions of all events in model
[e_specs, e_names] = event_specs(E);
n_events = size(e_specs, 2);
dur = 0;

% get design betas
% z is the variable for subjects
b{:,z} = betas(E); 
tmp=['beta_',roi];
cd 'D:\MK_MREEGRewModes_2016\fMRI\ROIAnanlysis\ROIs_Inc_vs_Con\v6_SPMcluster\Batch_ROI_results';


% pct_ev: extract percent signal change for all events in design:
% z is the variable for subjects
% for e_s = 1:n_events
%   pct_ev(e_s,z) = event_signal(E, e_specs(:,e_s), dur);
% end

% tmp2=['pct_ev',roi];
% 
% % fir_tc: extract FIR time courses of all bins averaged across runs:
% % Get compound event types structure
% ets = event_types_named(E);
% n_event_types = length(ets);
% 
% % Bin size in seconds for FIRclear
% 
% bin_size = tr(E);
% 
% % Length of FIR in seconds 
% fir_length = 18;
% 
% % Number of FIR time bins to cover length of FIR
% bin_no = fir_length / bin_size;
% 
% % Options - here 'single' FIR model, return estimated % signal change
% % z is used to calculate column number
% opts = struct('single', 1, 'percent', 1);
% for e_t = 1:n_event_types
%     fir_tc(:, (((z-1)*n_event_types)+e_t)) = event_fitted_fir(E, ets(e_t).e_spec, bin_size, ...
%     bin_no, opts);
% end
% 
% tmp3=['fir_',roi];
% 
% % get stats and stuff for all contrasts into statistics structure
% %marsS = compute_contrasts(E, 1:length(xCon)); 
%   
% 
 end;



save (tmp, 'b');
% save (tmp2, 'pct_ev');
% save (tmp3, 'fir_tc');

% keep ('rois', 'r', 'roi', 'roi_file');

end;

