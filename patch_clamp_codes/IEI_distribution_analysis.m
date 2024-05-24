%%% This script takes analyzed IEI data structs, groups IEI values by
%%% experimental conditions, and perform distribution analyses

%% Experimental settings
save_results = 1;

%experimental conditions (should match the numbers assigned to each group
%in the cell_id excel sheet)
exp_con = {'Dark','Light_6h','Dark_6h'};

%experiment name (correpsonding to the sheet name in the cell_id_index
% excel file)
exp_name = 'WT_Light_Dark_mini_rise_1';

%file name of results
save_file_name = 'WT_Light_Dark_mini_IEI_distribution.mat';

%where to save grouped files
fp_mini_IEI_group = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/slice_NT/mini_IEI_by_groups/';

%location of IEI files
fp_IEI = '/Users/wwneuro/My_Drive/Lab/Data_analysis/slice_NT/analyzed_mini_results/frq_only/IEI_analysis/';

%cell_id index readout
fp_mini = '/Users/wwneuro/My_Drive/Lab/Data_analysis/slice_NT/mini_data_by_groups/';

cell_id_index = readtable(strcat(fp_mini,'cell_id_index.xlsx'),'Sheet',exp_name);

%histogram?
histo_on = 1;

h_bin = 10; % in ms

%cumulative distribution?
cumu_on = 1;

c_bin = 0.01;

%color schemes

%Light/Dark color scheme
col{1} = '#309DD9';
col{2} = '#8FA4BF';
col{3} = '#295ABC';
%% group IEIs by condition

%pre-allocation
all_IEI = cell(1,numel(exp_con));

%loop through folder
cd(fp_IEI)
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
    curr_date = curr_name(end-9:end-4);
    
    cy = strcat('20',curr_date(1:2));
    cM = curr_date(3:4);
    cday = curr_date(5:6);
    
    date = datetime(strcat(cy,'-',cM,'-',cday),'InputFormat','y-M-d',...
        'Format', 'M/d/y');
    if ~(ismember(date,cell_id_index.Date))
        continue
    else
        load(all_files(fi+2).name)
        
        for ci = 1:numel(IEI_clean{1})
            if isempty(IEI_clean{1}{ci})
                continue
            else
                if isempty(find(cell_id_index.Date == date & cell_id_index.Cell_num == ci,1))
                    continue
                else
                
                    row = find(cell_id_index.Date == date & cell_id_index.Cell_num == ci);
                    cond_i = cell_id_index{row,'Cat'};
                    cell_ID = cell_id_index{row,'Cell_ID'};
                    all_IEI{1,cond_i}{1,cell_ID} = IEI_clean{1}{ci};
                end
            end
        end
    end
        
end

%% generate histograms for all events (no sampling)

% gather IEIs for all cells into one column

meta_IEI = cell(1,numel(exp_con));

for gi = 1:numel(exp_con)
    lgthct = 1;

    for ci = 1:numel(all_IEI{gi})
        lgth = size(all_IEI{gi}{ci},1);
        meta_IEI{gi}(lgthct:lgthct+lgth-1,1) = all_IEI{gi}{ci}(1:lgth,1);
        lgthct = lgthct+lgth;
    end
end

if histo_on == 1

    figure('Position',[500 100 1000 500])
    hold on
    
    for hi = 1:numel(meta_IEI)
        histogram(meta_IEI{hi}','BinWidth',h_bin,...
            'FaceColor',col{hi},'EdgeColor',col{hi},...
            'FaceAlpha',0.5)
    end
end

%% sample IEI values to have equal number of IEI values for each cell within each condition

%pre-allocation
selected_IEI = cell(1,numel(all_IEI)); % 100 values from each cell
cell_num = NaN(numel(all_IEI),1);

for ei = 1:numel(all_IEI)
    
    i_ct = 1;
    cellct = 0;
    
    %if less than 100 events in this cell, select all events
    %if more than 100 events, then randomly pick 100 events
    for cj = 1:numel(all_IEI{ei})
        curr_cell = all_IEI{ei}{cj};
        cellct = cellct+1;
        
        if size(curr_cell,1) < 100
            sel_IEI = curr_cell;
        else
            rand_index_exp = randi([1 size(curr_cell,1)],100,1);
            sel_IEI = curr_cell(rand_index_exp);
        end
        
        selected_IEI{ei}(i_ct:i_ct+numel(sel_IEI)-1,1) = sel_IEI;
        i_ct = i_ct + numel(sel_IEI);
    end
    cell_num(ei,1) = cellct;
end

IEI_size = min(cell_num)*100;
eql_sel_IEI = NaN(IEI_size,numel(cell_num)); %equal number of IEI values for each condition

for cei = 1:numel(cell_num)
    curr_cond = selected_IEI{cei};

    if size(curr_cond,1)>IEI_size
        rand_ind = randi([1 size(curr_cond,1)],IEI_size,1);
        hit_IEI = curr_cond(rand_ind);
    else
        hit_IEI = curr_cond;
    end

    eql_sel_IEI(:,cei) = hit_IEI;
end

%% K-S test

%cond 1-2
[h1,p1,ks2stat1] = kstest2(eql_sel_IEI(:,1),eql_sel_IEI(:,2));

%cond 1-3
[h2,p2,ks2stat2] = kstest2(eql_sel_IEI(:,1),eql_sel_IEI(:,3));

%cond 2-3
[h3,p3,ks2stat3] = kstest2(eql_sel_IEI(:,2),eql_sel_IEI(:,3));


%% histogram for selected events

if histo_on == 1

    figure('Position',[500 100 1000 500])
    hold on
    
    for hi = 1:size(eql_sel_IEI,2)
        histogram(eql_sel_IEI(:,hi)','BinWidth',h_bin,...
            'FaceColor',col{hi},'EdgeColor',col{hi},...
            'FaceAlpha',0.5)
    end
end

legend(exp_con,'Interpreter','none')

%% cumulative density function
 
cumu_IEI = cell(1,size(eql_sel_IEI,2));

if cumu_on == 1

    for cmi = 1:size(eql_sel_IEI,2)
        curr_cond = eql_sel_IEI(:,cmi);
        min_IEI = min(curr_cond);
        max_IEI = max(curr_cond);

        [cu_X,cu_Y] = cumhist(curr_cond,[min_IEI max_IEI],c_bin);
        cumu_IEI{cmi}(:,1) = cu_X;
        cumu_IEI{cmi}(:,2) = cu_Y;
    end

    figure();
    hold on

    for cui = 1:numel(cumu_IEI)
        plot(cumu_IEI{cui}(:,1),cumu_IEI{cui}(:,2),'Color',col{cui},'LineWidth',2)
    end

    ax = gca;
    ax.FontSize = 12;
    ax.LineWidth = 2;
    ax.YLabel.String = 'Cumulative %';
    ax.YLabel.FontSize = 14;
    ax.XLim = [0 1200];
    ax.YLim = [0 110];
    ax.YTick = [0 50 100];
    ax.XLabel.String = 'mEPSC IEI (ms)';
    ax.XLabel.FontSize = 14;
    legend(exp_con{1,:},'FontSize',14,'Location','southeast','Interpreter','none')
    legend('boxoff')
    %text(ax.XLim(2)-1500,70,strcat(exp_con{1}, ' vs. ',exp_con{2},': p=',num2str(p1)),'FontSize',12,'Interpreter','none')
    %text(ax.XLim(2)-1500,60,strcat(exp_con{1}, ' vs. ',exp_con{3},': p=',num2str(p2)),'FontSize',12,'Interpreter','none')
    %text(ax.XLim(2)-1500,50,strcat(exp_con{2}, ' vs. ',exp_con{3},': p=',num2str(p3)),'FontSize',12,'Interpreter','none')
    %title(strcat(num2str(ax.XLim(2)),'pA cutoff'))
    box off
    hold off

end

%% save data

if save_results == 1
    cd(fp_mini_IEI_group)
    save(save_file_name,'exp_con','all_IEI','cumu_IEI','eql_sel_IEI','meta_IEI',...
        'p1','p2','p3','selected_IEI')
end


