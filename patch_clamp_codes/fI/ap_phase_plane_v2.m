%%% This script take the first APs at rheobase for all cells in a given
%%% conditon and generate phase plane plot for each of them
%%% Also generate the meta average phase plane plot of all cells

%% load first_ap data struct
fp_fads = '/Users/wwneuro/My_Drive/Lab/Data_analysis/chronic_DREADDs/chronic_hm4di/fI_data_by_groups/firstAP_by_group/';

% should match the sheet name in the cell_id_index excel sheet
exp_name = 'WT_CP_48saline_24CNO_pooled';

load(strcat(fp_fads,exp_name,'.mat'))

figure_on = 1;

save_results = 1;

% save results to path
fp_results = '/Users/wwneuro/My_Drive/Lab/Data_analysis/chronic_DREADDs/chronic_hm4di/fI_data_by_groups/phase_plane_data';

% color codes
%saline
color{2} = '#A3A194'; 
color{1} = '#B679F2'; 

%xpro
% color{1} = '#669999'; 
% color{2} = '#732666'; 

%% generate phase plane plot for each cell
dVdt = cell(1,numel(exp_con)); %in V/s
dVdt_meta = cell(1,numel(exp_con)); % in V/s
V_meta = cell(1,numel(exp_con)); %in mV

for cdi = 1:numel(exp_con)
    c_V = all_ap{cdi}; % in mV
    c_dVdt = NaN(size(c_V,1)-1,1);
    c_V_meta =  NaN(size(c_V,1),1);
    c_dVdt_meta = NaN(size(c_V,1)-1,1);

    if figure_on == 1
        figure(cdi)
        hold on
    
        for ci = 1:size(c_V,2)
            c_dVdt(1:size(c_V,1)-1,ci) = diff(c_V(1:size(c_V,1),ci)).*10; %in V/s
            plot(c_V(1:size(c_V,1)-1,ci),c_dVdt(1:size(c_V,1)-1,ci),...,
                'Color',[0.8 0.8 0.8],'LineWidth',0.5)
        end
        
    dVdt{cdi} = c_dVdt;


        for ddi = 1:size(c_V,1)
            c_V_meta(ddi,1) = mean(c_V(ddi,1:size(c_V,2)),'omitnan');
        end

        c_dVdt_meta(1:size(c_V,1)-1,1) = diff(c_V_meta(1:size(c_V_meta,1),1)).*10;
        plot(c_V_meta(1:size(c_V_meta,1)-1),c_dVdt_meta(1:size(c_dVdt_meta,1),1),...,
            'Color',color{cdi},'LineWidth',3)

        ax = gca;
        ax.FontSize = 14;
        ax.LineWidth = 2;
        ax.YLabel.String = 'dV/dt (V/s)';
        ax.YLabel.FontSize = 14;
        ax.XLabel.String = 'V (mV)';
        ax.YLabel.FontSize = 14;
        ax.YLabel.Interpreter = 'none';
        
        hold off

        dVdt_meta{cdi} = c_dVdt_meta;
        V_meta{cdi} = c_V_meta;
    end

end

%% plot conditions in one figure
figure()
hold on
for cdii = 1:numel(exp_con)
    plot(V_meta{cdii}(1:size(V_meta{cdii},1)-1,1),dVdt_meta{cdii}(1:size(dVdt_meta{cdii},1),1),...,
    'Color',color{cdii},'LineWidth',3)
end

ax = gca;
ax.FontSize = 14;
ax.LineWidth = 2;
ax.YLabel.String = 'dV/dt (V/s)';
ax.YLabel.FontSize = 14;
ax.XLabel.String = 'V (mV)';
ax.YLabel.FontSize = 14;
ax.YLabel.Interpreter = 'none';
legend(exp_con{1,:},'FontSize',14,'Location','northwest','Interpreter','none','Box','off')

hold off

%% save results
if save_results == 1
    cd(fp_results)
    save(strcat(exp_name,'.mat'),'exp_con','all_ap','dVdt','V_meta','dVdt_meta','cell_id_index')
end
