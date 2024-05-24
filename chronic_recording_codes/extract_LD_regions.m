%%% Sorted in vivo data from ATP 2021 Neuron paper %%%
% all can be found on Turrigiano lab server- ATP- Main

%%
save_sorted_data_path = '/Users/wwneuro/My_Drive/Lab/MATLAB/invivo';

save_sorted_data_name = 'ER_ctrl_withFS_v2';

OG_data_name = 'CONTCELL_ER_MLS_v10_under3Quality';

cd(save_sorted_data_path)
load(strcat(OG_data_name,'.mat'))


master_data = CellStrct;
ext_data = struct();

%whether pick a specific time window
time_int = []; %in hours

%%
%define light/dark in hours (starting from light and alternate between light/dark)
% 1 = light; 0 = dark. (e.g. 12h 0 indicates ZT12 to ZT24 is dark) 

LD_index(1:26,1) = (0:12:300)';
for ii = 1:size(LD_index,1)
    if mod(LD_index(ii,1)/12,2) == 0
        LD_index(ii,2) = 1;
    else
        LD_index(ii,2) = 0;
    end
end


%%
cell_ct = 0;

for ci = 1:size(master_data,2)
    if master_data(ci).deprived == 0
        cell_ct = cell_ct + 1;

        if master_data(ci).tailSlope >= 0.004 && master_data(ci).neg_pos_time >= 0.39
            ext_data(cell_ct).ctype = 'RS';
        else
            ext_data(cell_ct).ctype = 'FS';
        end

        ext_data(cell_ct).time = master_data(ci).time;
        ext_data(cell_ct).OGcellID = ci;
        ext_data(cell_ct).animal = master_data(ci).animal;
        ext_data(cell_ct).channel = master_data(ci).channel;

        ext_data(cell_ct).onTime = master_data(ci).onTime;
        ext_data(cell_ct).offTime = master_data(ci).offTime;
        ext_data(cell_ct).EXPTSTART = master_data(ci).EXPTSTART;
        ext_data(cell_ct).trem = master_data(ci).trem;
        

        %hours from when lights come on first day (7:30a)
        ext_data(cell_ct).onHour = ext_data(cell_ct).onTime/3600; 
        ext_data(cell_ct).offHour = ext_data(cell_ct).offTime/3600;
        ext_data(cell_ct).duration = ext_data(cell_ct).offHour - ext_data(cell_ct).onHour;

    end
end


%% categorize the duration of each cell into light/dark
cell_num = size(ext_data,2);

for ci = 1:cell_num
    on_ind = NaN(numel(ext_data(ci).onHour),1);
    off_ind = NaN(numel(ext_data(ci).offHour),1);
    ri_ct = 1;

    clear LD_profile

    for ri = 1:numel(ext_data(ci).onHour)
        
        on_ind(ri) = find(ext_data(ci).onHour(ri) < LD_index(:,1),1,'first')-1;
        off_ind(ri) = find(ext_data(ci).offHour(ri) <= LD_index(:,1),1,'first')-1;
    
        if ri > 1 && on_ind(ri) == off_ind(ri-1)
            c_region = off_ind(ri) - off_ind(ri-1);
        else
            c_region = off_ind(ri) - on_ind(ri);
        end

        LD_profile(ri_ct:ri_ct+c_region,1:2) = LD_index(on_ind(ri):off_ind(ri),1:2);
        ri_ct = ri_ct+c_region;
    end

    ext_data(ci).LD_profile = LD_profile;

end

%% create animal/cell index (original cell IDS saved by animal and categoried by RS/FS)

%save by animal
animal_id_t = cell(size(ext_data,2),1);

for cci = 1:size(ext_data,2)
    animal_id_t{cci} = ext_data(cci).animal;
end

animal_id = unique(animal_id_t);
animal_t = cell(1,size(animal_id,1));

%first col saves original cell ID, second col saves cell type
% (1- RS; 0-FS)
cell_ogID = cell(1,size(animal_id,1));

cell_ct_2 = 0;

for ai = 1:size(animal_id,1)
    for aci = 1:size(ext_data,2)
        if strcmp(ext_data(aci).animal, animal_id{ai})
            cell_ct_2 = cell_ct_2 + 1;

            animal_t{ai}{cell_ct_2} = ext_data(aci).LD_profile;
            cell_ogID{ai}(cell_ct_2,1) = ext_data(aci).OGcellID;
            if strcmp(ext_data(aci).ctype,'RS')
                cell_ogID{ai}(cell_ct_2,2) = 1;
            elseif strcmp(ext_data(aci).ctype, 'FS')
                cell_ogID{ai}(cell_ct_2,2) = 0;
            end
            
        end
    end

    cell_ct_2 = 0;
end

%% pick cells that meet the criteria (if needed)
ext_data_final = struct();
animal_status = cell(1,size(animal_t,2));

if ~isempty(time_int)

    ct1 = 0;
    for cdi = 1:size(ext_data,2)
        for cai = 1:size(animal_id,1)
            if strcmp(ext_data(cdi).animal, animal_id{cai}) == 0
                continue
            else
                if ext_data(cdi).onHour <= time_int(1) && ext_data(cdi).offHour >= time_int(2)
                    ct1 = ct1+1;
            
                    ext_data_final(ct1).ctype = ext_data(cdi).ctype;
                    ext_data_final(ct1).animal = ext_data(cdi).animal;
                    ext_data_final(ct1).OGcellID = ext_data(cdi).OGcellID;
                    ext_data_final(ct1).time = ext_data(cdi).time;
                    ext_data_final(ct1).channel = ext_data(cdi).channel;
                    ext_data_final(ct1).onTime = ext_data(cdi).onTime;
                    ext_data_final(ct1).offTime = ext_data(cdi).offTime;
                    ext_data_final(ct1).onHour = ext_data(cdi).onHour;
                    ext_data_final(ct1).offHour = ext_data(cdi).offHour;
                    ext_data_final(ct1).duration = ext_data(cdi).duration;
                    ext_data_final(ct1).LD_profile = ext_data(cdi).LD_profile;
    
                    for cci = 1:size(cell_ogID{cai},1)
                        if ext_data(cdi).OGcellID == cell_ogID{cai}(cci,1)
                            animal_status{cai}(cci,1) = 1;
                        end
                    end
                else
                        for cci = 1:size(cell_ogID{cai},1)
                            if ext_data(cdi).OGcellID == cell_ogID{cai}(cci,1)
                                animal_status{cai}(cci,1) = 0;
                            end
                        end
                end
            end
        end
    end

else

    ext_data_final = ext_data;
end

%% 

%% save data 
cd(save_sorted_data_path)
if ~isempty(time_int)
    final_save_name = strcat(save_sorted_data_name,'_',num2str(time_int(1)),'_',num2str(time_int(2)),'.mat');
else
    final_save_name = strcat(save_sorted_data_name,'.mat');
end

save(final_save_name,...,
    'LD_index','ext_data_final','animal_id',...,
    'cell_ogID','time_int','animal_status','-v7.3');