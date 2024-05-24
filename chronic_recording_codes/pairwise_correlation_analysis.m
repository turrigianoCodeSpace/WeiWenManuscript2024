%%%%% Pairwise correlation analysis %%%%%

% Algorithm from ATP2019 (PNAS)

% Each spike train is binned into spike counts of bin size 100 ms,
% generating a vector of spike counts for each cell. Then correlation
% coefficient p for a pair of neurons is calculated in 30 min episodes
% using a sliding window of 5 min. For a 12h light or dark session, this
% will produce 139 p values, which are then averaged to produce a single
% value as an indicator of the correlation in this 12h.

%%
save_data_path = '/Users/wwneuro/My_Drive/Lab/MATLAB/invivo';

save_data_name = 'CPP2recov_1st_inj_ctrl_pr_correlation_v2.mat';

%define dark, light, light+CPP (12h session)
dark = [156 168]; %hours
light = [144 156];
light_cpp = [192 204];

%bin size (in seconds)
bin = 0.1; 

%episode length (in minutes)
epi = 30;

%sliding window size (in minutes)
sw = 5;

%whether to include FS cells 
FS_on = 0;

%%

%data structure for spike vectors
%animal -> pairs of cell -> two fields:trains of bin vectors; trains of
%spike vectors generated from the corresponding bin vectors

spike_bnd_dark = cell(1,numel(cell_ogID)); % save spike vectors
spike_bnd_light = cell(1,numel(cell_ogID));
spike_bnd_light_cpp = cell(1,numel(cell_ogID));

%data structure for pairwise correlation values
spike_corr_dark = cell(1,numel(cell_ogID));
spike_corr_light = cell(1,numel(cell_ogID));
spike_corr_light_cpp = cell(1,numel(cell_ogID));

%data structure for normalzied pairwise correlation values
spike_corr_dark_norm = cell(1,numel(cell_ogID));
spike_corr_light_norm = cell(1,numel(cell_ogID));
spike_corr_light_cpp_norm = cell(1,numel(cell_ogID));

%mean correlation for each pair (saved by animals)
corr_ave = cell(1,3);

%save pairs of cell within each animal
cell_pairs = cell(1,numel(cell_ogID));

%identity of cell pairs
%1- RS to RS; 2- FS to FS; 3- RS to FS or FS to RS
cell_pair_type = cell(1,numel(cell_ogID));

%% 
for ai = 1:numel(cell_ogID)
    if animal_hit(ai) == 0
        continue
    else

        if FS_on == 0
            fs_hit = cell_ogID{ai}(:,2) == 1 ;
            cell_fs_hit = NaN(numel(fs_hit),1);

            for chi = 1: size(cell_ogID{ai},1)
                if animal_status{ai}(chi,1) ~= 0 && fs_hit(chi,1) ~= 0
                    cell_fs_hit(chi,1) = 1;
                else
                    cell_fs_hit(chi,1) =0;
                end
            end
            
            cell_hit = cell_ogID{ai}(cell_fs_hit ~= 0);
            
        else
            cell_hit = cell_ogID{ai}(animal_status{ai} ~= 0);
        end
        
        all_pairs = nchoosek(cell_hit,2);
        for api = 1:size(all_pairs,1)
            for aci = 1:size(cell_ogID{ai})
                c1_ind = find(cell_ogID{ai}(:,1) == all_pairs(api,1));
                c2_ind = find(cell_ogID{ai}(:,1) == all_pairs(api,2));

                if cell_ogID{ai}(c1_ind,2) == 1
                    if cell_ogID{ai}(c2_ind,2) == 1
                        cell_pair_type{ai}(api,1) = 1; %RS to RS
                    elseif cell_ogID{ai}(c2_ind,2) == 0
                        cell_pair_type{ai}(api,1) = 3; %RS to FS
                    end
                elseif cell_ogID{ai}(c1_ind,2) == 0
                     if cell_ogID{ai}(c2_ind,2) == 1
                        cell_pair_type{ai}(api,1) = 3; %FS to RS
                    elseif cell_ogID{ai}(c2_ind,2) == 0
                        cell_pair_type{ai}(api,1) = 2; %FS to FS
                     end
                end

            end
        end
                 
        cell_pairs{ai} = all_pairs;

        for pi = 1:size(all_pairs,1)
            cell1 = all_pairs(pi,1);
            cell2 = all_pairs(pi,2);

            for ci = 1:size(ext_data_final,2)
                if ext_data_final(ci).OGcellID == cell1
                    spike_train_1 = ext_data_final(ci).time;
                end

                if ext_data_final(ci).OGcellID == cell2
                    spike_train_2 = ext_data_final(ci).time;
                end
            end

            %session indices
            dark_on_1 = find(spike_train_1>=dark(1)*3600,1,'first');
            dark_off_1 = find(spike_train_1>=dark(2)*3600,1,'first');

            light_on_1 = find(spike_train_1>=light(1)*3600,1,'first');
            light_off_1 = find(spike_train_1>=light(2)*3600,1,'first');           

            light_cpp_on_1 = find(spike_train_1>=light_cpp(1)*3600,1,'first');
            if isempty(find(spike_train_1>=light_cpp(2)*3600,1,'first'))
                light_cpp_off_1 = numel(spike_train_1);
                session1 = spike_train_1(end)/3600 - light_cpp(1);
            else
                light_cpp_off_1 = find(spike_train_1>=light_cpp(2)*3600,1,'first');  
                session1 = 12;
            end
            
            dark_on_2 = find(spike_train_2>=dark(1)*3600,1,'first');
            dark_off_2 = find(spike_train_2>=dark(2)*3600,1,'first');

            light_on_2 = find(spike_train_2>=light(1)*3600,1,'first');
            light_off_2 = find(spike_train_2>=light(2)*3600,1,'first');           

            light_cpp_on_2 = find(spike_train_2>=light_cpp(1)*3600,1,'first');
            if isempty(find(spike_train_2>=light_cpp(2)*3600,1,'first'))
                light_cpp_off_2 = numel(spike_train_2);
                session2 = spike_train_2(end)/3600 - light_cpp(1);
            else
                light_cpp_off_2 = find(spike_train_2>=light_cpp(2)*3600,1,'first');  
                session2 = 12;
            end            

            %3rd level- first cell saves bin vectors
            spike_bnd_dark{ai}{pi}{1} = generate_sliding_bins(dark(1),12,bin,epi,sw);
            spike_bnd_light{ai}{pi}{1} = generate_sliding_bins(light(1),12,bin,epi,sw);
            if session1<session2
                spike_bnd_light_cpp{ai}{pi}{1} = generate_sliding_bins(light_cpp(1),session1,bin,epi,sw);
            else
                spike_bnd_light_cpp{ai}{pi}{1} = generate_sliding_bins(light_cpp(1),session2,bin,epi,sw);
            end
                


            %3rd level- second cell saves spike vectors (saved by cell)
            spike_bnd_dark{ai}{pi}{2}{1} = spike_train_1(dark_on_1:dark_off_1);
            spike_bnd_dark{ai}{pi}{2}{2} = spike_train_2(dark_on_2:dark_off_2);

            spike_bnd_light{ai}{pi}{2}{1} = spike_train_1(light_on_1:light_off_1);
            spike_bnd_light{ai}{pi}{2}{2} = spike_train_2(light_on_2:light_off_2);

            spike_bnd_light_cpp{ai}{pi}{2}{1} = spike_train_1(light_cpp_on_1:light_cpp_off_1);
            spike_bnd_light_cpp{ai}{pi}{2}{2} = spike_train_2(light_cpp_on_2:light_cpp_off_2);

            %3rd level- third cell saves binned spike vectors (saved by
            %cell)
            for bvi = 1:numel(spike_bnd_dark{ai}{pi}{1})

                %cell1
                spike_bnd_dark{ai}{pi}{3}{1}{bvi} = ...,
                    histcounts(spike_bnd_dark{ai}{pi}{2}{1},'BinEdges',spike_bnd_dark{ai}{pi}{1}{bvi});

                spike_bnd_light{ai}{pi}{3}{1}{bvi} = ...,
                    histcounts(spike_bnd_light{ai}{pi}{2}{1},'BinEdges',spike_bnd_light{ai}{pi}{1}{bvi});

                 if bvi > numel(spike_bnd_light_cpp{ai}{pi}{1})
                     spike_bnd_light_cpp{ai}{pi}{3}{1}{bvi} = NaN;
                 else
                    spike_bnd_light_cpp{ai}{pi}{3}{1}{bvi} = ...,
                    histcounts(spike_bnd_light_cpp{ai}{pi}{2}{1},'BinEdges',spike_bnd_light_cpp{ai}{pi}{1}{bvi});
                 end
 

                %cell2

                spike_bnd_dark{ai}{pi}{3}{2}{bvi} = ...,
                    histcounts(spike_bnd_dark{ai}{pi}{2}{2},'BinEdges',spike_bnd_dark{ai}{pi}{1}{bvi});

                spike_bnd_light{ai}{pi}{3}{2}{bvi} = ...,
                    histcounts(spike_bnd_light{ai}{pi}{2}{2},'BinEdges',spike_bnd_light{ai}{pi}{1}{bvi});
                 
                if bvi > numel(spike_bnd_light_cpp{ai}{pi}{1})
                     spike_bnd_light_cpp{ai}{pi}{3}{2}{bvi} = NaN;
                 else
                    spike_bnd_light_cpp{ai}{pi}{3}{2}{bvi} = ...,
                    histcounts(spike_bnd_light_cpp{ai}{pi}{2}{2},'BinEdges',spike_bnd_light_cpp{ai}{pi}{1}{bvi});
                end


                %calculate pairwise correlation for this episode (bvi)

                spike_corr_dark{ai}{pi}(bvi,1) = ...
                    pr_correlation_spike(spike_bnd_dark{ai}{pi}{3}{1}{bvi},spike_bnd_dark{ai}{pi}{3}{2}{bvi});

                spike_corr_light{ai}{pi}(bvi,1) = ...
                    pr_correlation_spike(spike_bnd_light{ai}{pi}{3}{1}{bvi},spike_bnd_light{ai}{pi}{3}{2}{bvi});
 
                 spike_corr_light_cpp{ai}{pi}(bvi,1) = ...
                    pr_correlation_spike(spike_bnd_light_cpp{ai}{pi}{3}{1}{bvi},spike_bnd_light_cpp{ai}{pi}{3}{2}{bvi});
                             
            end

            
            %normalize all correlations values of a pair to the selected region
            %during dark (900-1000min from ZT0, stable region)
            curr_mean_corr_dark = mean(spike_corr_dark{ai}{pi}(35:55,1),'omitnan');


            %normalize all correlations values of a pair to the selected region
            %during dark (900-1000min from ZT0, stable region)
            for bvii = 1:numel(spike_bnd_dark{ai}{pi}{1})
                spike_corr_dark_norm{ai}{pi}(bvii,1) = spike_corr_dark{ai}{pi}(bvii,1)/curr_mean_corr_dark;
                spike_corr_light_norm{ai}{pi}(bvii,1) = spike_corr_light{ai}{pi}(bvii,1)/curr_mean_corr_dark;
                spike_corr_light_cpp_norm{ai}{pi}(bvii,1) = spike_corr_light_cpp{ai}{pi}(bvii,1)/curr_mean_corr_dark;
            end

             %calculate mean correlation for a cell pair (each row: one pair)
             %1st cell- dark
             %2nd cell- light
             %3rd cell- light+cpp
             corr_ave{1}{ai}(pi,1) = mean(spike_corr_dark{ai}{pi},'omitnan');
             corr_ave{2}{ai}(pi,1) = mean(spike_corr_light{ai}{pi},'omitnan');
             corr_ave{3}{ai}(pi,1) = mean(spike_corr_light_cpp{ai}{pi},'omitnan');
             
             corr_ave_norm{1}{ai}(pi,1) = mean(spike_corr_dark_norm{ai}{pi},'omitnan');
             corr_ave_norm{2}{ai}(pi,1) = mean(spike_corr_light_norm{ai}{pi},'omitnan');
             corr_ave_norm{3}{ai}(pi,1) = mean(spike_corr_light_cpp_norm{ai}{pi},'omitnan');


        end
    end
end

%% save data
cd(save_data_path)

save(save_data_name,"dark","light","light_cpp","spike_bnd_dark","spike_bnd_light",...
"spike_bnd_light_cpp","spike_corr_dark","spike_corr_light","spike_corr_light_cpp",...
"spike_corr_dark_norm","spike_corr_light_norm","spike_corr_light_cpp_norm",...
"animal_id","animal_hit","animal_status","LD_index","cell_ogID","corr_ave", ...,
"corr_ave_norm","cell_pair_type","cell_pairs","-v7.3")

