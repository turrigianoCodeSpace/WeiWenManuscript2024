% This script takes the analyzed f-I data struct and calculates kinetic
% features for the first action potential at the rheobase

%%%% LOAD THE ANALYZED FI DATA IN THE WORKSPACE BEFORE RUNNING %%%%
%%%% NEED CELL CAPACITANCE VALUES AS WELL %%%%
%% Experimental parameters

%location of the analyzed fI files
save_file_path = '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/analyzed_AP_kinetics/';

%name of the saved file
save_file_name = 'AP_kinetics_220805';

%sampling rate (for fI protocol should be 10k, but should double check)
sp_rate = 10000;

%save results
save_result = 1;

%% calculations


%pre-allocations

%for all action potentials detected in all current steps
ap_rise = cell(numel(cell_id),1);
ap_fall = cell(numel(cell_id),1);
ap_max_rise_rate = cell(numel(cell_id),1);
ap_max_fall_rate = cell(numel(cell_id),1);
max_Na_current = cell(numel(cell_id),1);

%for the first action potential at the rheobase
ap_rise_1st = NaN(numel(cell_id),1);
ap_fall_1st = NaN(numel(cell_id),1);
ap_max_rise_rate_1st = NaN(numel(cell_id),1);
ap_max_fall_rate_1st = NaN(numel(cell_id),1);
max_Na_current_1st = NaN(numel(cell_id),1);


%start calculations
for ci = 1:numel(cell_id)
    if isempty(cell_id{1,ci})
        ap_rise{1,ci} = [];
        ap_fall{1,ci} = [];
        ap_max_rise_rate{1,ci} = [];
        ap_max_fall_rate{1,ci} = [];
        max_Na_current{1,ci} = [];
        
    else
        for ti = 1:size(cell_id{1,ci},1)
            ti_ind = cell_id{1,ci}(ti,1);

            for si = 1:size(AP_peak{1,ci}{1,ti_ind},1)

                if isnan(AP_peak{1,ci}{1,ti_ind}(1,1))
                    ap_rise{1,ci}{1,ti_ind} = [];
                    ap_fall{1,ci}{1,ti_ind} = [];
                    ap_max_rise_rate{1,ci}{1,ti_ind} = [];
                    ap_max_fall_rate{1,ci}{1,ti_ind} = [];
                    max_Na_current{1,ci} = [];
                else
                    
                    % ap upstroke rise time (from threshold to peak, in ms)
                    ap_rise{1,ci}{1,ti_ind}(si,1) = ...,
                        (AP_peak{1,ci}{1,ti_ind}(si,1) - V_th{1,ci}{1,ti_ind}(si,1))/sp_rate*1000;

                    % ap downstroke fall time (from peak to fast trough, in
                    % ms)
                    ap_fall{1,ci}{1,ti_ind}(si,1) = ...,
                        (f_trough{1,ci}{1,ti_ind}(si,1) - AP_peak{1,ci}{1,ti_ind}(si,1))/sp_rate*1000;

                    %maximum dV/dt during upstroke (V/s)
                    ap_max_rise_rate{1,ci}{1,ti_ind}(si,1) = ...,
                        max(dV_sec{1,ci}{1,ti_ind}(V_th{1,ci}{1,ti_ind}(si,1):AP_peak{1,ci}{1,ti_ind}(si,1)));

                    %maximum dV/dt during downstroke (V/s)
                    ap_max_fall_rate{1,ci}{1,ti_ind}(si,1) = ...,
                        -min(dV_sec{1,ci}{1,ti_ind}(AP_peak{1,ci}{1,ti_ind}(si,1):f_trough{1,ci}{1,ti_ind}(si,1)));

                    %maximum sodium current estimation (max(dV/dt)*Cm, in
                    %nA) (might be an over-estimation because there's
                    %current injection)
                    max_Na_current{1,ci}{1,ti_ind}(si,1) = ap_max_rise_rate{1,ci}{1,ti_ind}(si,1)*mean(Cm{1,ci}(:,1))/1000;

                end
            end % per AP
        end % per trace
    end

    for tii = 1:numel(ap_rise{1,ci})
        if isempty(ap_rise{1,ci}{1,tii})
            continue
        else
            ap_rise_1st(ci,1) = ap_rise{1,ci}{1,tii}(1,1);
            ap_fall_1st(ci,1) = ap_fall{1,ci}{1,tii}(1,1);
            ap_max_rise_rate_1st(ci,1) = ap_max_rise_rate{1,ci}{1,tii}(1,1);
            ap_max_fall_rate_1st(ci,1) = ap_max_fall_rate{1,ci}{1,tii}(1,1);
            max_Na_current_1st(ci,1) = max_Na_current{1,ci}{1,tii}(1,1);
            break
        end
    end

end % per cell

%% save results
if save_result == 1

    cd(save_file_path)
    save(strcat(save_file_name,'.mat'), 'cell_id', 'dV_sec', 'AP_peak', 'AP_peak_1st',...,
        'f_trough','V_th','ap_rise','ap_fall','ap_max_fall_rate','ap_max_rise_rate',...,
        'max_Na_current','ap_rise_1st','ap_fall_1st','ap_max_fall_rate_1st',...,
        'ap_max_rise_rate_1st','max_Na_current_1st')
end
