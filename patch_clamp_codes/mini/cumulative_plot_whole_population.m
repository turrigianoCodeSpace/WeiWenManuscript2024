% This script takes grouped mini events, randomly picks 100 events from
% each cell under each experimental condition, and then calculates the
% cumulative distribution

%% load mini populations (grouped by experimental condition)

%whether to save results
save_results = 1;

%location of grouped mini events
fp_all_mini_group = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/slice_NT/all_mini_by_group/';

%sub folder
sub = 'rise_1';

%experiment name (correpsonding to the sheet name in the cell_id_index
% excel file)
exp_name = 'WT_Light_Dark_mini_rise_1';

%load grouped mini events
cd(strcat(fp_all_mini_group, sub))
load(strcat(exp_name,'_2024','.mat'))

%where to save the analyzed cumulative data
fp_cumu = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/slice_NT/cumulative_data/cumulative_stats/';

%name of the saved cumulative data
save_file_name = strcat(exp_name, '_2024','_cumu_stats.mat');


%% prepare data
% choose experimental conditions you want to scale and plot for comparison
% refer to the exp_con for the order of cell arrays in cumulative
exp_ind = 1; %DR+CNO %Dark
ctrl_ind = 2; %CNO %Light
add_ind = 3; %TTX+GLYX %Dark 6h

%color code (darker analogous purple/blue)
% color{2} = '#F8C2D2';
% %color{2} = '#573E5C'; 
% color{1} = '#63B6FF';
% color{3} = '#D095DB';

%saline
% color{2} = '#A3A194'; 
% color{1} = '#B679F2'; 
% %xpro
% color{2} = '#669999'; 
% color{1} = '#732666';

%Light/Dark color scheme
color{1} = '#309DD9';
color{2} = '#8FA4BF';
color{3} = '#295ABC';
%%

%amp cutoff (high) (if no cutoff, then h_cut = 0)
h_cut = 0;

exp_amps = selected_mini_events{1,exp_ind}(:,1); 
ctrl_amps = selected_mini_events{1,ctrl_ind}(:,1); 
add_amps = selected_mini_events{1,add_ind}(:,1); 

% make sure both groups have same number of selected mini events
% the group with more events is normalized to the one with less events
% (random selection)

if size(exp_amps,1) > size(ctrl_amps,1)
    rand_index = randi([1,size(exp_amps,1)],size(ctrl_amps,1),1);
    exp_amps_final = exp_amps(rand_index);
    ctrl_amps_final = ctrl_amps;
    add_amps_final = add_amps(rand_index);
    
elseif size(exp_amps,1) < size(ctrl_amps,1)
    rand_index = randi([1,size(ctrl_amps,1)],size(exp_amps,1),1);
    exp_amps_final = exp_amps;
    ctrl_amps_final = ctrl_amps(rand_index);
    add_amps_final = add_amps(rand_index);
    
else
    exp_amps_final = exp_amps;
    ctrl_amps_final = ctrl_amps;
    rand_index = 'NaN';
    add_amps_final = add_amp;

end    

%sort event amplitudes from low to high
exp_amps_sort_t = sort(exp_amps_final);
ctrl_amps_sort_t = sort(ctrl_amps_final);
add_amps_sort_t = sort(add_amps_final);

if h_cut == 0
    exp_amps_sort = exp_amps_sort_t;
    ctrl_amps_sort = ctrl_amps_sort_t;
    add_amps_sort = add_amps_sort_t;
else

    if isempty(find(exp_amps_sort_t>=h_cut,1,'first'))
        e_cut_ind = numel(exp_amps_sort_t);
    else
        e_cut_ind = find(exp_amps_sort_t>=h_cut,1,'first');
    end

    if isempty(find(ctrl_amps_sort_t>=h_cut,1,'first'))
        c_cut_ind = numel(ctrl_amps_sort_t);
    else
        c_cut_ind = find(ctrl_amps_sort_t>=h_cut,1,'first');
    end

    if isempty(find(add_amps_sort_t>=h_cut,1,'first'))
        a_cut_ind = numel(add_amps_sort_t);
    else
        a_cut_ind = find(add_amps_sort_t>=h_cut,1,'first');
    end

    min_pop = min([c_cut_ind e_cut_ind]);
    exp_amps_sort = exp_amps_sort_t(1:min_pop);
    ctrl_amps_sort = ctrl_amps_sort_t(1:min_pop);
    add_amps_sort = add_amps_sort_t(1:min_pop);


%     if c_cut_ind > e_cut_ind
%         exp_amps_sort = exp_amps_sort_t(1:e_cut_ind);
%         ctrl_amps_sort = ctrl_amps_sort_t(1:e_cut_ind);
%     else
%         exp_amps_sort = exp_amps_sort_t(1:c_cut_ind);
%         ctrl_amps_sort = ctrl_amps_sort_t(1:c_cut_ind);
%     end
% 
%     add_amps_sort = add_amps_sort_t(1:a_cut_ind);

end

%% Scaling

%use linear fitting to get scale factor (scale exp to ctrl)
f = fittype('a*x+b','coefficients',{'a','b'});
slope_start = (max(exp_amps_sort)-min(exp_amps_sort))/(max(ctrl_amps_sort)-min(ctrl_amps_sort));
[lnfit,gof] = fit(ctrl_amps_sort,exp_amps_sort,f,'StartPoint',[slope_start 0]);

exp_scaled = (exp_amps_sort-lnfit.b) ./ lnfit.a;

% if isempty(find(exp_scaled_t < 5, 1, 'last'))
%     exp_scaled = exp_scaled_t;
% else
%     scaled_cut_ind = find(exp_scaled_t < 5, 1, 'last');
% 
%     exp_scaled = exp_scaled_t(scaled_cut_ind+1:end);
% end

%plot scaled data
figure()
plot(ctrl_amps_sort,ctrl_amps_sort,'k', 'Linewidth',2)
hold on
plot(ctrl_amps_sort,exp_amps_sort,'o','MarkerSize',12)
plot(lnfit)

ax1 = gca;
ax1.FontSize = 14;
ax1.LineWidth = 2;
ax1.YLabel.String = strcat(exp_con{exp_ind}, ' (pA)');
ax1.YLabel.FontSize = 14;
%ax1.XTick = [0 20 40 60 80];
ax1.XLabel.String = strcat(exp_con{ctrl_ind}, ' (pA)');
ax1.XLabel.FontSize = 14;
legend(strcat(exp_con{ctrl_ind}, ' vs. ', exp_con{ctrl_ind}),...
    strcat(exp_con{exp_ind}, ' vs. ', exp_con{ctrl_ind}),'FontSize',12,'Location','southeast','Interpreter','none')

legend('boxoff')
text(10,max(ctrl_amps_sort),strcat('y=',num2str(lnfit.a),'x',num2str(lnfit.b)),'FontSize',12)

box off

hold off

%% generate cumulative distribution plots

%cumulative distribution for scaled exp data
[cumu_exp_X_scaled, cumu_exp_Y_scaled] = cumhist(exp_scaled,[min(exp_scaled) max(exp_scaled)],0.1);

%cumulative distribution for sorted exp and control data
[cumu_ctrl_X,cumu_ctrl_Y] = cumhist(ctrl_amps_sort, [min(ctrl_amps_sort) max(ctrl_amps_sort)],0.1);
[cumu_exp_X,cumu_exp_Y] = cumhist(exp_amps_sort, [min(exp_amps_sort) max(exp_amps_sort)],0.1);
[cumu_add_X,cumu_add_Y] = cumhist(add_amps_sort,[min(add_amps_sort) max(add_amps_sort)],0.1);

%two-sample Kolmogorov-Smirnov test
%after scaling
[h1,p1,ks2stat1] = kstest2(exp_scaled,ctrl_amps_sort);
%before scaling
[h2,p2,ks2stat2] = kstest2(exp_amps_sort,ctrl_amps_sort);
%additional
[h3,p3,ks2stat3] = kstest2(exp_amps_sort,add_amps_sort);
[h4,p4,ks2stat4] = kstest2(ctrl_amps_sort,add_amps_sort);

figure();
hold on

plot(cumu_ctrl_X, cumu_ctrl_Y,'Color',color{ctrl_ind},'LineWidth',2)

%plot(cumu_exp_X_scaled, cumu_exp_Y_scaled,'k','LineWidth',2)

plot(cumu_exp_X, cumu_exp_Y,'Color',color{exp_ind},'LineWidth',2)

plot(cumu_add_X, cumu_add_Y,'Color',color{add_ind},'LineWidth',2)

ax = gca;
ax.FontSize = 12;
ax.LineWidth = 2;
ax.YLabel.String = 'Cumulative %';
ax.YLabel.FontSize = 14;
ax.XLim = [0 40];
ax.YLim = [0 110];
ax.YTick = [0 50 100];
ax.XLabel.String = 'mEPSC amplitude (pA)';
ax.XLabel.FontSize = 14;
legend(exp_con{ctrl_ind},'Scaled',exp_con{exp_ind},'FontSize',14,'Location','southeast','Interpreter','none')
legend('boxoff')
text(ax.XLim(2)-25,70,strcat(exp_con{exp_ind}, ' vs. ',exp_con{ctrl_ind},': p=',num2str(p2)),'FontSize',12,'Interpreter','none')
text(ax.XLim(2)-25,60,strcat('Scaled vs. ',exp_con{ctrl_ind},': p=',num2str(p1)),'FontSize',12,'Interpreter','none')
%title(strcat(num2str(ax.XLim(2)),'pA cutoff'))
box off
hold off

%% save files
if save_results == 1
    cd(strcat(fp_cumu, sub))
    save(save_file_name,'exp_amps_sort','ctrl_amps_sort','exp_scaled','cumu_exp_X_scaled','cumu_exp_Y_scaled',...
        'h1','p1','ks2stat1','h2','p2','ks2stat2','rand_index')

%     save(save_file_name,'exp_amps_sort','ctrl_amps_sort','add_amps_sort','exp_scaled','cumu_exp_X_scaled','cumu_exp_Y_scaled',...
%     'h1','p1','ks2stat1','h2','p2','ks2stat2','p3','p4','rand_index')
end
