%%% analysis of width at half height for a given AP in a spike train %%%


%% experimental settings

%save results
save_results = 1;

%experiment name (correpsonding to the sheet name in the cell_id_index
% excel file)
exp_name = 'TTX_6h_24h_ActD';

%where to save grouped width mat files
fp_width_group = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/fI_data_by_groups/width_by_group';

%location of analyzed fI data folder
fp_fI = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/analyzed_fI_results/';

%location of grouped fI data folder
fp_fI_g = ...,
    '/Users/wwneuro/My_Drive/Lab/Data_analysis/culture_experiments/fI_data_by_groups/';

%experimental conditions
exp_con = {'Ctrl_TTX_6h','ActD_TTX_6h','TTX_6h','TTX_24h'};

%number of current injections
cj = 7;

%import cell_id_index table 
cd(fp_fI_g)
cell_id_index = readtable('cell_id_index.xlsx','Sheet',exp_name);

%% group all width values for all cells by experimental conditions

all_width = cell(1,numel(exp_con));%same order as exp_con
dt_t = cell(1,numel(exp_con));

cd(fp_fI)
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
    if contains(curr_name,'test')
        continue
    else
        curr_date = curr_name(4:9);
    end
    
    cy = strcat('20',curr_date(1:2));
    cM = curr_date(3:4);
    cday = curr_date(5:6);
    
    date = datetime(strcat(cy,'-',cM,'-',cday),'InputFormat','y-M-d',...
        'Format', 'M/d/y');
    if ~(ismember(date,cell_id_index.Date))
        continue
    else
        load(all_files(fi+2).name)

        for ci = 1:size(width,2)
            if isempty(width{ci})
                continue
            else
                if isempty(find(cell_id_index.Date == date & cell_id_index.Cell_num == ci,1))
                    continue
                else
                    row = find(cell_id_index.Date == date & cell_id_index.Cell_num == ci);
                    cond_i = cell_id_index{row,'Cat'};
                    cell_ID = cell_id_index{row,'Cell_ID'};

                    curr_width = width{ci};
                    dt_t{cond_i}.Rheobase(cell_ID,1) = rheobase(ci,1);
                    dt_t{cond_i}.Rheobase_ind(cell_ID,1) = rheobase_ind(ci,1);

                    %convert zeros to NaNs
                    for wi = 1:size(curr_width,1)
                        for wj = 1:size(curr_width,2)
                            if curr_width(wi,wj) == 0
                                curr_width(wi,wj) = NaN;
                            end
                        end
                    end
                           
                    all_width{cond_i}{cell_ID} = curr_width;
                end
            end
        end
    end
end

%% Width analyses

firstAP_width = cell(1,numel(exp_con)); %first AP at rheobase
firstAP_lastInj_width = cell(1,numel(exp_con)); %first AP at the last current injection (max)
firstAPs_med_width = cell(1,numel(exp_con)); %median of first APs across all current injections
firstAPs_delta_width = cell(1,numel(exp_con)); %rate of change in width across first APs for a given cell

allAPs_med_width = cell(1,numel(exp_con));%median for each spike train
allAPs_delta_width = cell(1,numel(exp_con)); %rate of change in width across a spike train

lastAPs_width = cell(1,numel(exp_con)); %last APs of each spike train
lastAPs_med_width = cell(1,numel(exp_con)); %median of last APs across all current injections
lastAPs_delta_width = cell(1,numel(exp_con)); %rate of change in width across last APs for a given cell

for cdi = 1:numel(all_width)
    for ci = 1:numel(all_width{cdi})
        firstAP_width{cdi}(ci,1) = all_width{cdi}{ci}(find(~isnan(all_width{cdi}{ci}),1,'first'));
        firstAP_lastInj_width{cdi}(ci,1) = all_width{cdi}{ci}(size(all_width{cdi}{ci},1),1);
        firstAPs_med_width{cdi}(ci,1) = median(all_width{cdi}{ci}(:,1),'omitnan');

        [fap_num,f_ind] = count_non_nan(all_width{cdi}{ci}(:,1));
        firstAPs_delta_width{cdi}(ci,1) = sum(diff(all_width{cdi}{ci}(:,1)),'omitnan')/fap_num;


        
        for ti = 1:size(all_width{cdi}{ci},1)
            allAPs_med_width{cdi}{ci}(ti,1) = median(all_width{cdi}{ci}(ti,:),'omitnan');

            [afap_num,af_ind] = count_non_nan(all_width{cdi}{ci}(ti,:));
            allAPs_delta_width{cdi}{ci}(ti,1) = sum(diff(all_width{cdi}{ci}(ti,:)),'omitnan')/afap_num;

            if find(isnan(all_width{cdi}{ci}(ti,:)),1,'first') == 1
                lastAPs_width{cdi}{ci}(ti,1) = NaN;
            else
                nan_ind = find(isnan(all_width{cdi}{ci}(ti,:)),1,'first');
                if isempty(nan_ind)
                    lastAPs_width{cdi}{ci}(ti,1) = all_width{cdi}{ci}(ti,end);
                else
                    lastAPs_width{cdi}{ci}(ti,1) = all_width{cdi}{ci}(ti,nan_ind-1);
                end
            end
        end

        lastAPs_med_width{cdi}(ci,1) = median(lastAPs_width{cdi}{ci},'omitnan');

        [lfap_num,lf_ind] = count_non_nan(lastAPs_width{cdi}{ci});
        lastAPs_delta_width{cdi}(ci,1) = sum(diff(lastAPs_width{cdi}{ci}),'omitnan')/lfap_num;
    end
end


%% save in struct

% vars = {'firstAP','firstAP_lastInj','firstAPs_med','firstAPs_delta',...
%     'allAPs_med','allAPs_delta','lastAPs','lastAPs_med','lastAPs_delta'};

for dti = 1:numel(exp_con)
    for cti = 1:numel(all_width{dti})
        dt_t{dti}.firstAP = firstAP_width{dti};
        dt_t{dti}.firstAP_lastInj = firstAP_lastInj_width{dti};
        dt_t{dti}.firstAPs_med = firstAPs_med_width{dti};
        dt_t{dti}.firstAPs_delta = firstAPs_delta_width{dti};

        dt_t{dti}.lastAPs_med = lastAPs_med_width{dti};
        dt_t{dti}.lastAPs_delta = lastAPs_delta_width{dti};
        dt_t{dti}.all_width = all_width{dti};

        for tti = 1:numel(allAPs_med_width{dti})

            dt_t{dti}.allAPs_med(1:cj,tti) = allAPs_med_width{dti}{tti}(end-(10-cj)-cj+1:end-(10-cj),1);
            dt_t{dti}.allAPs_delta(1:cj,tti) = allAPs_delta_width{dti}{tti}(end-(10-cj)-cj+1:end-(10-cj),1);
            dt_t{dti}.lastAPs(1:cj,tti) = lastAPs_width{dti}{tti}(end-(10-cj)-cj+1:end-(10-cj),1);
        end
    end
end


Ctrl_TTX_6h_w = dt_t{1};
ActD_TTX_6h_w = dt_t{2};
TTX_6h_w = dt_t{3};
TTX_24h_w = dt_t{4};

%% save to file
cd(fp_width_group)

exp_con_w = cell(1,numel(exp_con));
for exi = 1:numel(exp_con)
    exp_con_w{exi} = strcat(exp_con{exi},'_w');
end
save(strcat(exp_name,'_width.mat'),exp_con_w{1,:},'exp_con','cell_id_index')
