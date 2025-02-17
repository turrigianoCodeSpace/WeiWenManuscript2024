%% Experimental settings

%save results
save_results = 1;

%figure on
figure_on_unscaled = 1;
figure_on_scaled = 1;

%experiment name 
experiment = {'230331'};

%appendix name
appendix = 'minus40';

%subfolder
sub = 'rise';

%where to save grouped files
fp_wa_group = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/waveform_average_by_groups/';

%experimental conditions
exp_con = {'ctrl','apv'};

%location of waveform average files (per cell)
fp_wak = '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/waveform_average_kinetics/';

%%
%pre-allocation
all_wavg_temp = cell(1,numel(experiment));

%loop through folder
cd(strcat(fp_wak, sub))
all_files_temp = dir;
for fi = 1:size(all_files_temp,1)
    if strcmp('.DS_Store',all_files_temp(fi).name)
        delete '.DS_Store'
        continue
    end
end

all_files = dir;
file_num = numel(all_files)-2;

for fi = 1:file_num
    curr_name = all_files(fi+2).name;
    curr_prefix = curr_name(1:6);


    for ei = 1:numel(experiment)
        if ~strcmp(curr_prefix,experiment{ei})
            continue
        else

            for cdi = 1:numel(exp_con)

                curr_con = curr_name(15:15+numel(exp_con{cdi})-1);

                if ~strcmp(curr_con,exp_con{cdi})
                    continue
                else
                    load(all_files(fi+2).name)

                    for ci = 1:numel(cell_id)

                        all_wavg_temp{1,ei}{1,cdi}(:,ci) = wavg_per_cell{1,1}(:,ci);
                    end
                end
            end
        end
    end
end

%% scale peak and plot
%cultured neuron + APV color scheme
color1 = '#D095DB'; 
color2 = '#573E5C';

for ci = 4:4%numel(cell_id)
    condition1 = all_wavg_temp{1,1}{1,1}(:,ci);
    condition2 = all_wavg_temp{1,1}{1,2}(:,ci);

    [peak_cond1, peak_ind_cond1] = min(condition1);
    [peak_cond2, peak_ind_cond2] = min(condition2);
    
    %the entire waveform scaled by peak amplitude (will get similar results if 
    % using the scaling factor calcuated from the cumulative population, faster this way)
    
    scale_factor = mean(condition1((peak_ind_cond1-1):(peak_ind_cond1+1),1))/...,
        mean(condition2((peak_ind_cond2-1):(peak_ind_cond2+1),1));
    
    condition2_scaled = NaN(size(condition1,1),1);
    diff = size(condition1,1)-size(condition2,1);
    
    if diff >=0
        condition2_scaled(diff+1:end) = condition2 .* scale_factor;
    else
        condition2_scaled = condition2(abs(diff)+1 : end) .* scale_factor;
    end
        
    
    %plotting unscaled
    
    plot_range = (10:150);
    y_limit = [-Inf Inf];
    
    if figure_on_unscaled == 1
        figure();
        plot(condition1(plot_range),'Color',color1,'LineWidth',4)
        hold on
        %plot(condition3(plot_range),'Color',color3,'LineWidth',4)   
        plot(condition2(plot_range),'Color',color2,'LineWidth',4)
        %hold on
        %plot(meta_ave.DR_CNO_24(20:110),'Color',[0.27 0.51 0.71],'LineWidth',4)
        %draw scale
        scale_x_start = plot_range(end)-30;
        plot([scale_x_start; scale_x_start+10],[-10; -10], '-k',[scale_x_start;scale_x_start],[-10; -5], '-k', 'LineWidth',2)
        text(scale_x_start, -8, '5 pA', 'HorizontalAlignment', 'right')
        text(scale_x_start+5, -11, '2 ms', 'HorizontalAlignment', 'center')
        
        title(strcat('cell ',num2str(ci), ' unscaled'))
        ylim(y_limit)
        legend(exp_con{1,:},'Location','southeast')
        fontsize(gcf,scale = 1.5)
        
        box off
        hold off
    end
    
    
    %plotting scaled
    if figure_on_scaled == 1    
        figure();
        plot(condition1(plot_range),'Color',color1,'LineWidth',4)
        hold on
        plot(condition2_scaled(plot_range),'Color',color2,'LineWidth',4)
    %     hold on
    %     plot(shk3_scaled(20:110),'Color',[0.27 0.51 0.71],'LineWidth',4)
        title(strcat('cell ',num2str(ci), ' scaled'))
        ylim(y_limit)
        legend(exp_con{1,:},'Location','southeast')
        fontsize(gcf,scale = 1.5)


        box off
        hold off
    end

    
end

