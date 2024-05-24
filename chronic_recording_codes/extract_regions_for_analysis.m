%%% Sorted in vivo data from ATP 2021 Neuron paper %%%
% all can be found on Turrigiano lab server- ATP- Main

%%
save_sorted_data_path = '/Users/wwneuro/My_Drive/Lab/MATLAB/invivo';

save_sorted_data_name = 'ER_ctrl_withFS.mat';

OG_data_name = 'CONTCELL_ER_MLS_v10_under3Quality';

cpp_inj = 0;

cd(save_sorted_data_path)
load(strcat(OG_data_name,'.mat'))

%%
master_data = CONTCELL_CPP2recov.MASTER;
ext_data = struct();

%%
cell_ct = 0;

for ci = 1:size(master_data,2)
    if master_data(ci).deprived == 0 && numel(master_data(ci).onTime) == 1
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



%%
%define light/dark in hours (starting from light and alternate between light/dark)
% 1 = light; 0 = dark.

LD_index(1:24,1) = (0:12:276)';
for ii = 1:size(LD_index,1)
    if mod(LD_index(ii,1)/12,2) == 0
        LD_index(ii,2) = 1;
    else
        LD_index(ii,2) = 0;
    end
end

%%
%categorize the duration of each cell into light/dark, and also label +/-
%CPP injection 
%(LD_profile: 1st col- hours; 2nd col- light/dark; 3rd col- +/- CPP)

%for animals that had one CPP injection at the time of ER, CPP session is
%from 168h to 180h (LD_index = 15)

%for animals that had two CPP injections, CPP sessions are: 192h to 204h,
%216h to 228h (LD_index = 17 and 19)

if cpp_inj == 1
    cpp_ind = 168;
elseif cpp_inj == 2
    cpp_ind = [192 216];
elseif cpp_inj == 0
    cpp_ind = 0;
end

cell_num = size(ext_data,2);
%LD_profile = cell(1,cell_num);

for ci = 1:cell_num
    on_ind = find(ext_data(ci).onHour < LD_index(:,1),1,'first')-1;
    off_ind = find(ext_data(ci).offHour <= LD_index(:,1),1,'first')-1;
    LD_profile = LD_index(on_ind:off_ind,1:2);

    if cpp_ind ~= 0

        for ldi = 1:size(LD_profile,1)
            if ismember(LD_profile(ldi,1),cpp_ind)
                LD_profile(ldi,3) = 1;
            else
                LD_profile(ldi,3) = 0;
            end
        end
    end

    ext_data(ci).LD_profile = LD_profile;

end

%% Further exclusion 
% 
% For an animal to be included for correlation analysis, it would need to
% have at least two cells that have at least 1 overlapping session with CPP


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


animal_status = cell(1,size(animal_t,2));
animal_hit = NaN(size(animal_t,2),1);

for asi = 1:size(animal_t,2)
    for asci = 1:size(animal_t{asi},2)
       
        animal_status{asi}(asci,1) = sum(animal_t{asi}{asci}(:,3));
        
    end

    if sum(animal_status{asi}) > 1
        animal_hit(asi,1) = 1;
    else
        animal_hit(asi,1) = 0;
    end
    
end


cell_ct_3 = 0;
ext_data_final = struct();

for cdi = 1: size(ext_data,2)
    for cai = 1:size(animal_id,1)

        if strcmp(ext_data(cdi).animal, animal_id{cai}) == 0
            continue
        else
            if animal_hit(cai,1) == 1
                cell_ct_3 = cell_ct_3 + 1;

                ext_data_final(cell_ct_3).ctype = ext_data(cdi).ctype;
                ext_data_final(cell_ct_3).animal = ext_data(cdi).animal;
                ext_data_final(cell_ct_3).OGcellID = ext_data(cdi).OGcellID;
                ext_data_final(cell_ct_3).time = ext_data(cdi).time;
                ext_data_final(cell_ct_3).channel = ext_data(cdi).channel;
                ext_data_final(cell_ct_3).onTime = ext_data(cdi).onTime;
                ext_data_final(cell_ct_3).offTime = ext_data(cdi).offTime;
                ext_data_final(cell_ct_3).onHour = ext_data(cdi).onHour;
                ext_data_final(cell_ct_3).offHour = ext_data(cdi).offHour;
                ext_data_final(cell_ct_3).duration = ext_data(cdi).duration;
                ext_data_final(cell_ct_3).LD_profile = ext_data(cdi).LD_profile;

            end
        end
    end
end


%% save extracted data

save(save_sorted_data_name,'ext_data','LD_index','cpp_ind','animal_hit',...,
    'animal_status','ext_data_final','animal_id','cell_ogID')
