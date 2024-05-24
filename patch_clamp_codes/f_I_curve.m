%%%% This script takes the imported f-I mat file and calculate basic stats
%%%% for each property

%%%% Under each field, data from each experimental condition (i.e. WT,
%%%% WT_CNO, DR_CNO, and, WT_saline), are saved in corresponding columns 1, 2, 3,
%%%% and 4.

%% Initialization

%location where the mat file will be saved
fp_analyzed_data = '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/fI_data_by_groups';

%name of the mat file
filename = {'fI_WT_Light_Dark_IE_stats.mat'};

%experimental conditions (stats in each column follow this order)
cond = {'CNO_saline','DR_CNO_saline'};

%functions for calculating stats 
func = {@mean, @std, @nansem};

%field names (corresponding to properties)
fields = {'MFR','IFR','mean_IFR'};

%current steps mode
curr_mo = 2;
%1- all current steps saved in the "current_inj" vector
%2- current steps saved under each condition in the "curr_c.curr_inj" field

%plot
plot_on = 1;

%save results
save_results = 0;

%% calculate stats 

%calculate stats for each current step
struct_all = cell(1,3); %results from each function are saved in corresponding cells

for sti = 1:numel(func) % per stat
    for cdi = 1:numel(cond) % per experimental condition
        curr_c = eval(cond{cdi});
        
        for fi = 1:numel(fields) % per field
            if size(curr_c.(fields{fi}),2)<2 % ignore fields that are not current-step related
                continue
            else
                for cui = 1:size(curr_c.(fields{fi}),1) % per current step
                    struct_all{sti}.(fields{fi})(cui,cdi) = func{sti}(curr_c.(fields{fi})(cui,:));
                end
            end
        end
    end
end

%current injection vector (just pick one vector, should all be the same in
%terms of plotting

if curr_mo == 2
    current_inj = curr_c.curr_inj(:,1);
end

%transfer results to corresponding data structs for storage
ave = struct_all{1};
std = struct_all{2};
sem = struct_all{3};

%calculate median for each cell (med)
% and
%get values from the rheobase step for each cell (rheo)
%get values from 100 pA over rheobase for each cell (above100)

for cdii = 1:numel(cond) % per experimental condition
        curr_c = eval(cond{cdii});
        
        for fii = 1:numel(fields) % per field
            if size(curr_c.(fields{fii}),2)<2 % ignore fields that are not current-step related
                continue
            else
                for cuii = 1:size(curr_c.(fields{fii}),2) % per cell
                    med.(fields{fii})(cuii,cdii) = median(curr_c.(fields{fii})(:,cuii),'omitnan');
                end
                
%                 vals = get_step_values(curr_c.(fields{fii}), 2, 0);
%                 rheo.(fields{fii})(1:size(vals,1),cdii) = vals;
%                 
%                 vals_1 = get_step_values(curr_c.(fields{fii}), 2, 100);
%                 above100.(fields{fii})(1:size(vals_1,1),cdii) = vals_1;
                
            end
        end
end
 
%current injection vector (just pick one vector, should all be the same in
%terms of plotting

%% plotting

% data range for plotting
plot_range = 1:20;
plot_variable = 'mean_IFR';
plot_stat = 'ave';
err = 'sem';
y_label = 'mean IFR (Hz)';
mksz = 12; %marker size
figure_window = [500 100 500 400];

%marker/line type
marker{1} = '-o';
%marker{2} = '-o';
marker{2} = '-o';
marker{3} = ':^';
marker{4} = ':^';

%color code

%Light/Dark color scheme
% color_edge{1} = '#309DD9';
% color_edge{2} = '#8FA4BF';
% color_edge{3} = '#295ABC';

%Light/Dark, CPP/Saline color scheme
% color_edge{1} = '#BF8C60';
% color_edge{2} = '#DD6F16';

% %CT color scheme
% color_edge{1} = '#000000'; %black
% color_edge{2} = '#DB7093'; %pink
% color_edge{3} = '#4680B2'; %steel blue

%CPP color scheme
% color1 = '#F17F73'; %salmon
% color2 = '#4680B2'; %steel blue
% color3 = '#F9A36B'; %light orange
% color4 = '#63B6FF'; %light blue

%culture TTX/APV data color scheme
%darker color transition (still purple)
% color_edge{1} = '#F8C2D2'; 
% color_edge{2} = '#573E5C';
% color_edge{3} = '#913DA2';
% color_edge{4} = '#7F7CD6'; 

% %PhTx/TTX color scheme
% color_edge{1} = '#000000'; %NT
% color_edge{2} = '#D973BE'; %TTX
% color_edge{3} = '#F2A413'; %PhTX
% color_edge{4} = '#5D71B3'; %TTX+PhTX

% %fig3 2/6h and ActD data- green/orange
% color_edge{1} = '#95E29E'; %green
% color_edge{2} = '#2F8B3A'; %green
% color_edge{3} = '#F8C2D2'; %orange
% color_edge{4} = '#D973BE'; %orange
 
 %GLYX data
% color_edge{1} = '#F8C2D2'; 
% color_edge{2} = '#573E5C';
% color_edge{3} = '#AFCFE8';

% %sTNFR data- monochromatic orange 20h
% color_edge{1} = '#97B09A'; 
% color_edge{2} = '#B82A57';
% color_edge{3} = '#B3812D'; 
% color_edge{4} = '#665538';

%sTNFR data- monochromatic green 6h
% color_edge{1} = '#97B09A'; 
% color_edge{2} = '#BF5841';
% color_edge{3} = '#B3812D'; 
% color_edge{4} = '#334741';

%xpro, saline, and cpp color scheme
%saline
color_edge{1} = '#A3A194'; 
color_edge{2} = '#B679F2'; 
%xpro
% color_edge{1} = '#669999'; 
% color_edge{2} = '#732666'; 
%cpp
% color_edge{3} = '#000000';
% color_edge{4} = '#C95142';

color_face{1} = 'none'; 
color_face{2} = 'none';
% color_face{1} = color_edge{1};
% color_face{2} = color_edge{2};

color_face{3} = 'none';
color_face{4} = 'none';


curr_y = eval(strcat(plot_stat,'.',plot_variable));
curr_y_err = eval(strcat(err,'.',plot_variable));

if plot_on == 1

    figure('Position',figure_window)
    hold on
    
    for gi = 1:numel(cond)
        %loop all conditions
        errorbar([0,current_inj(plot_range)']',[0,curr_y(plot_range,gi)']',[0,curr_y_err(plot_range,gi)']',marker{gi},'MarkerSize',mksz,...
                'MarkerEdgeColor',color_edge{gi},'MarkerFaceColor',color_face{gi},'LineWidth',2,'Color',color_edge{gi})
    end
%         %cond2
%     errorbar(current_inj(plot_range),curr_y(plot_range,2),curr_y_err(plot_range,2),'-o','MarkerSize',mksz,...
%             'MarkerEdgeColor',color2, 'MarkerFaceColor',color2,'LineWidth',2,'Color',color2)    
%         
%     %cond3
%     errorbar(current_inj(plot_range),curr_y(plot_range,3),curr_y_err(plot_range,3),':^','MarkerSize',mksz,...
%             'MarkerEdgeColor',color3, 'MarkerFaceColor','none','LineWidth',2,'Color',color3)
%     %cond4
%     errorbar(current_inj(plot_range),curr_y(plot_range,4),curr_y_err(plot_range,4),':^','MarkerSize',mksz,...
%             'MarkerEdgeColor',color4, 'MarkerFaceColor','none','LineWidth',2,'Color',color4)
        
%     %cond5
%     errorbar(current_inj(plot_range),curr_y(plot_range,5),curr_y_err(plot_range,5),':^','MarkerSize',mksz,...
%             'MarkerEdgeColor',color5, 'MarkerFaceColor','none','LineWidth',2,'Color',color5)
    %WT+saline
%     errorbar(current_inj(plot_range),curr_y(plot_range,2),curr_y_err(plot_range,2),'-o','MarkerSize',8,...
%             'MarkerEdgeColor','k', 'MarkerFaceColor','k','LineWidth',2,...
%             'Color','k')

    ax = gca;
    ax.FontSize = 14;
    ax.LineWidth = 2;
    ax.YLabel.String = y_label;
    ax.YLabel.FontSize = 14;
    ax.XLabel.String = 'Injected Current (pA)';
    ax.YLabel.FontSize = 14;
    ax.YLabel.Interpreter = 'none';
    ax.XLim = [0 420];
    ax.YLim = [0 140];

    legend(cond{1,:},'FontSize',14,'Location','northwest','Interpreter','none','Box','off')
    %title('Bursting')

end

%% save results
if save_results == 1
    cd (fp_analyzed_data)
    save(filename{1},'ave','std','sem','med','cond')
end
