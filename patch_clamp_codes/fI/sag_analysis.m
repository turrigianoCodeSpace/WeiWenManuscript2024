% This script takes hyperpolarization traces obtained under the current clamp
% and calculate the sag amplitude.
% The sag amplitude is defined as the difference between the peak
% hyperpolarized membrane potential and the following steady-state potential 
% during the pulse.
% Steady-state is defined as the region where membrance potential
% is less than 10% of the peak potential.

% For cell-to-cell comparison, normalzied sag amp ([0,1])is also reported: 
% the sag amplitude of each cell is divided by the peak deflection
% (difference between resting and peak hyperpolarized potential)

%% General analysis settings

%loop through all the folders
experiment = {'221202'};

%subfolders under a certain experiment if any (e.g. before/after)
sub = '';

%filepath where you keep data folders
fp_data = '/Users/wwneuro/My_Drive/Lab/Data/culture_experiments/hyperpolarization/';

%location of the analyzed sag data
fp_analyzed_data = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/analyzed_sag';

%save results
save_results = 1;

%whether to plot sag calculation
figure_sag_on = 1;

%whether to plot pp calculation
figure_pp_on = 0;

%%%% recording settings
%hyperpolarizing current increment (in A)
hp_curr = -40 * 10^-12;

%latency of pulse onset (in s)
step_start = 1;

%duration of the current injection pulse (in s)
pulse = 1;

%sampling rate (in Hz)
samprate = 10000;

%% data readout and analysis

for jj = 1:1%numel(experiment)
    current_experiment = experiment{jj};
    
    %file names of the analyzed fI and pp data were they to be saved
    save_file_name = strcat('sag_', current_experiment, sub, '.mat');

    current_data_folder = strcat(fp_data, current_experiment, '/', sub);
    [raw_h5_files, aDAT_sag, cell_id, cell_num, filename] = h5_file_readout(current_data_folder);


    %%% pre-allocate cell arrays

    Rin = cell(1,max(cell_num));
    Vm = cell(1,max(cell_num));
    Cm = cell(1,max(cell_num));
    peakAmp = cell(1,max(cell_num));
    sagAmp = cell(1,max(cell_num));
    sagAmp_norm = cell(1,max(cell_num));
    sagTau = cell(1,max(cell_num));
    curr_inj = cell(1,max(cell_num));
    rebound_amp = cell(1,max(cell_num));
    
    
    %% passive property and sag parameter calculation
    for ci = 1:max(cell_num)
        if isempty(cell_id{1,ci})
            continue
        else
            trace_start = cell_id{1,ci}(1,1);
            trace_end = cell_id{1,ci}(end,1);

            for ti = trace_start : trace_end
                if ~ismember(ti,cell_id{1,ci})
                    Rin{1,ci}(ti,1) = NaN;
                    Vm{1,ci}(ti,1) = NaN;
                    Cm{1,ci}(ti,1) = NaN;
                    peakAmp{1,ci}(ti,1) = NaN;
                    sagAmp{1,ci}(ti,1) = NaN;
                    sagAmp_norm{1,ci}(ti,1) = NaN;
                    sagTau{1,ci}(ti,1) = NaN;
                    curr_inj{1,ci}(ti,1) = NaN;

                    continue
                else
                    data = aDAT_sag{1,ci}(:,ti);
                    curr_inj_curr = hp_curr*(ti-trace_start+1);

                    [Rin_curr, Cm_curr, mTau_curr, fit_vals_curr, trace_time_curr, exp_fit_curr] = ...,
                        get_PP_I_clamp(data, step_start, pulse, curr_inj_curr, samprate, figure_pp_on);

                    Rin{1,ci}(ti,1) = Rin_curr;
                    Cm{1,ci}(ti,1) = Cm_curr;
                    Vm{1,ci}(ti,1) = mean(data(0.2*samprate:0.5*samprate),'omitnan');
                    sagTau{1,ci}(ti,1) = mTau_curr;
                    curr_inj{1,ci}(ti,1) = curr_inj_curr*10^12; %in pA
 
                    %region for sag analysis
                    sag_val = data(step_start*0.8*samprate : (step_start+pulse*1.2)*samprate, 1);

                    %region for rebound depolarization analysis
                    reb_val = data((step_start+pulse-0.1)*samprate:(step_start+pulse+0.5)*samprate,1);
%                     
%                     figure
%                     plot(sag_val)
%                     hold on

                    % rebound deploarization analysis
                    [reb_peak,reb_peak_ind] = max(reb_val(1:0.3*samprate));
                    reb_peak_amp = reb_peak-Vm{1,ci}(ti,1);

                    rebound_amp{1,ci}(ti,1) = reb_peak_amp;

                    if figure_sag_on == 1
                        figure('position',[56 200 1200 490]); 
                         subplot(1,3,3)
                         plot(reb_val)
                         hold on
                         

                         reb_val_fit = reb_val(0.3*samprate:end);
                         if reb_peak < mean(reb_val_fit,'omitnan') + 4*std(reb_val_fit,'omitnan')
                             scatter(reb_peak_ind,reb_peak,14,'red','filled')
                         else
                             scatter(reb_peak_ind,reb_peak,14,'black','filled')
                         end

                         title('Rebound')
                         legend({'data','rebound peak'},'Location','southeast')
                         hold off
                    end


                    [peak, ind] = min(sag_val);
                    
                    if ind >= 4000
                        peakAmp{1,ci}(ti,1) = 0;
                        sagAmp{1,ci}(ti,1) = 0;
                        sagAmp_norm{1,ci}(ti,1) = 0;
                        
                        disp(strcat('No peak found in Cell',num2str(ci),' Trace',num2str(ti),'!'))
                        
                    else
                        %define steady state
                        ss_est_range = ind+0.5*pulse*samprate : ind+0.8*pulse*samprate;
                        ss_amp = mean(sag_val(ss_est_range),'omitnan');
                        %plot(ss_range,data(ss_range))
                        
                        %peak deflection
                        peak_def = peak - Vm{1,ci}(ti,1);


                        %sag amplitude
                        sag_amp = ss_amp - peak;
                        sag_amp_norm = abs(sag_amp/peak_def); %normalized sag amplitude (to peak delfection)

                        peakAmp{1,ci}(ti,1) = peak_def;
                        sagAmp{1,ci}(ti,1) = sag_amp;
                        sagAmp_norm{1,ci}(ti,1) = sag_amp_norm;
    

    
                        if figure_sag_on == 1
                             subplot(1,3,1)
                             plot(sag_val)
                             hold on
                             plot(ss_est_range, sag_val(ss_est_range),'r','LineWidth',2)
              
                             plot(sag_val(1:ind),'c','LineWidth',2)
                             scatter(ind,peak,14,'red','filled')
                             title(strcat('Cell',num2str(ci),', Trace',num2str(ti)))
                             legend({'data','steady-state','hyperpolarization','peak'},'Location','northeast')
                             hold off
    
                             subplot(1,3,2)
                             plot(trace_time_curr,fit_vals_curr,'LineWidth',2)
                             hold on
                             plot(exp_fit_curr)
                             title('Fitting')
                             legend({'data','first-exp fitting'},'Location','southeast')
                             hold off

                        end
                    
                    end
                end
            end %trace
        end
    end %cell
            
            
end %experiment

%% save to file
if save_results == 1
    cd(fp_analyzed_data)
    save(save_file_name, 'Rin','Cm','Vm','peakAmp','sagAmp','sagAmp_norm','sagTau','curr_inj','rebound_amp',...
        'cell_id','aDAT_sag')
end