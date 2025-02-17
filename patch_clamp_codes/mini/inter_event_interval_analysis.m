%%% This script takes the analyzed mini data structure and returns the
%%% inter-event-interval vectors for each trace of each cell

%%
fp_analyzed_mini = '/Users/wwneuro/My_Drive/Lab/Data_analysis/slice_NT/analyzed_mini_results/frq_only/rise_1/';

experiment = '231117';

load(strcat(fp_analyzed_mini,'MINIANALYSIS_',experiment,'.mat'))

fp_result = '/Users/wwneuro/My_Drive/Lab/Data_analysis/slice_NT/analyzed_mini_results/frq_only/IEI_analysis/';

samp_rate = 5000; %Hz

plot_on = 0;

figure_window = [500 100 800 420];

%% calculate IEI for all events in each cell
IEI = cell(1,numel(AMP_ALL{1})); %results saved in mseconds
IEI_clean = cell(1,numel(AMP_ALL{1})); % NaNs removed, IEI values for each cell saved in a single column

for ci = 1:numel(AMP_ALL{1})
    curr_cell = AMP_ALL{1}{ci};
    timeindx = TIME_INDICES{1}{ci};
    event_start_ind = EVENT_START_INDICES{1}{ci};
    decay_ind = DECAY_INDICES{1}{ci};
    sm_pk_ind = SMOOTH_PEAK_INDICES{1}{ci};

    for ti = 1:size(curr_cell,2)
         if plot_on == 1
             curr_Dat = nonzeros(aDAT_all{1}{ci}(:,ti));

             figure(figure('Position',figure_window))
             hold on
             plot(curr_Dat,'k')
         end


        for ei = 1:size(curr_cell,1)
            if ei+1 >= size(curr_cell,1)
                continue
            else
                if ~isnan(curr_cell(ei,ti)) && curr_cell(ei+1,ti) ~= 0                  
                   if isnan(decay_ind(ei,ti))
                       pre_end = imeindx(ei,ti)+sm_pk_ind(ei,ti)+10; %if no decay ind, then default decay set to 2 ms;
                   else
                       pre_end = timeindx(ei,ti)+sm_pk_ind(ei,ti)+decay_ind(ei,ti);
                   end

                   if isempty(event_start_ind{ei+1,ti})
                       IEI{1}{ci}(ei,ti) = NaN;
                       continue
                   else
                       post_start = timeindx(ei+1,ti)+event_start_ind{ei+1,ti}-1;
                   end

                   c_iei = (post_start - pre_end)/samp_rate*1000; %in ms

                   IEI{1}{ci}(ei,ti) = c_iei;

                   if plot_on == 1
                       plot(pre_end,curr_Dat(pre_end),'ro') % end of pre-event
                       plot(post_start,curr_Dat(post_start),'go') % start of post-event
                   end

                else
                    IEI{1}{ci}(ei,ti) = NaN;
                end
            end
        end %per event
    end % per trace
end % per cell

%clean up IEI
for cj = 1:size(IEI{1},2)
    nan_ind = ~isnan(IEI{1}{cj});
    IEI_clean{1}{cj} = nonzeros(IEI{1}{cj}(nan_ind));
end


%% save analyzed files
cd(fp_result)
save_file_name = strcat('IEI_',experiment,'.mat');
save(save_file_name, 'IEI','IEI_clean')
