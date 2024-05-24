%% Select which cell and which trace to plot

% where to save the selected trace data
save_data_path = '/Users/wwneuro/My_Drive/Lab/Data_analysis/chronic_DREADDs/chronic_hm4di/representative_traces/fI';

% name of the saved data
save_data_name = 'TTX_2h_6h_ActD_IHP_175p.mat';

%save results?
save_results = 0;

% analyzed fI results folder
fp_data = '/Users/wwneuro/My_Drive/Lab/Data_analysis/chronic_DREADDs/chronic_hm4di/analyzed_fI_results';

% experimental conditions
cond = {'CNO_saline','DR_CNO_saline','CNO_xpro','DR_CNO_xpro'};

% Choose f_I analysis mat file
experiment = {'fI_211112','fI_211111','fI_211014','fI_211014'};

% Choose a cell in each experiment for plotting
cell_num = [6 1 12 3];

% Plotting mode (0 for non-normarlized traces, 1 for normalized traces) 
plot_mode = 0;

%data strcuture mode (1 for fI data only protocol [old], 2 for alternating
%protocol [new])
dt_mode = 1;

% current step
stim = 200;

% current increment
cur_inc = 20;

%tint factors (btw 0 and 1, closer to 1 = more tint)
tint_factor = [0 0 0 0];

%tint and convert function
    %x-RGB value, y-tint factor
fun = @(x,y) (x+(255-x)*y)./255;

% color1 = '#E88BE2';
% color2 = '#78BEEB';
%darker color transition (still purple)
colorcode = cell(1,3);
% colorcode{1} = '#F9C3D2'; 
% colorcode{2} = '#913DA2';
% colorcode{3} = '#573E5C';
% colorcode{4} = '#7F7CD6';

%saline
colorcode{1} = '#A3A194'; 
colorcode{2} = '#B679F2'; 
%xpro
colorcode{3} = '#669999'; 
colorcode{4} = '#732666';
%cpp
% colorcode{1} = '#000000';
% colorcode{2} = '#C95142';

% %CT color scheme
% colorcode{1} = '#000000'; %black
% colorcode{2} = '#DB7093'; %pink
% colorcode{3} = '#4680B2'; %steel blue

%Light/Dark color scheme
% colorcode{1} = '#309DD9';
% colorcode{2} = '#8FA4BF';
% colorcode{3} = '#295ABC';

%Light/Dark, CPP/Saline color scheme
% colorcode{2} = '#BF8C60';
% colorcode{1} = '#DD6F16';

%culture- NT,APV,TTX,TTX_GLYX
% colorcode{1} = '#F8C2D2'; 
% % col_bar{2} = '#913DA2'; 
% colorcode{2} = '#573E5C';
% colorcode{3} = '#AFCEE6';

% %fig3 2/6h and ActD data- green/orange
% colorcode{1} = '#95E29E'; %green
% colorcode{2} = '#2F8B3A'; %green
% colorcode{3} = '#F8C2D2'; %orange
% colorcode{4} = '#D973BE'; %orange

%% data readout
% Variable initializations
cell_data = cell(1,numel(experiment));
trace_data = cell(1,numel(experiment));
trace_id = NaN(1,numel(experiment));
rheo_id = NaN(1,numel(experiment));

for fi = 1:numel(experiment)
    curr_file = matfile(strcat(fp_data,'/',experiment{fi}, '.mat'));

    if dt_mode == 1
        curr_file_data = curr_file.aDAT;
    elseif dt_mode == 2
        curr_file_data = curr_file.aDAT_fI;
    end

    curr_cell_id = curr_file.cell_id;
    curr_rheobase_ind = curr_file.rheobase_ind;
    
    cell_data{1,fi} = curr_file_data{1,cell_num(fi)};
    rheo_id(1,fi) = curr_rheobase_ind(fi,1);
    
    if plot_mode == 0
        trace_id(1,fi) = curr_cell_id{1,cell_num(fi)}(1,1)+stim/cur_inc-1;
    elseif plot_mode == 1
        trace_id(1,fi) = curr_cell_id{1,cell_num(fi)}(1,1)+rheo_id(1,fi)-1+stim/cur_inc;
    end
    
    trace_data{1,fi} = cell_data{1,fi}(:,trace_id(1,fi));
end

%% Plotting 

plot_range =8000:23000;
figure('position',[56 200 800 400])
hold on
for ci = 1:numel(cell_num)
    plot(trace_data{1,ci}(plot_range,1),'Color',colorcode{ci},'LineWidth',3)
end

scale_coor = plot_range(end)-plot_range(1);
%draw scale
plot([scale_coor-500; scale_coor+1500],[0; 0], '-k',[scale_coor-500;scale_coor-500],[0; 20], '-k', 'LineWidth',2)
text(scale_coor-600, 5, '20 mV', 'HorizontalAlignment', 'right')
text(scale_coor+1400, -5, '0.2 s', 'HorizontalAlignment', 'center')

legend(cond,'Interpreter','none')

if plot_mode == 0
    title(strcat(num2str(stim), 'pA injected'))
elseif plot_mode == 1
    title(strcat(num2str(stim),'pA over rheobase'))
end

hold off

%% Drawing current step

% figure('position',[56 200 1000 490])
% bl_before = zeros(1,1999);
% current_step = repmat(stim/1000, 1, 5000);
% bl_after = zeros(1,1999);
% 
% current_all = horzcat(bl_before,current_step,bl_after);
% 
% plot(current_all,'k', 'LineWidth',3)
% ylim([-0.2 1])
    
%% save data

if save_results == 1
    cd(save_data_path)
    
    save(save_data_name, 'cond','experiment','cell_num','trace_data')
end
