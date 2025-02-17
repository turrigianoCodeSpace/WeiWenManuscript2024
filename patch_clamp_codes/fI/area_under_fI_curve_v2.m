%% import mat files, areas of each variable stored under each experimental condition

%location where the mat file will be saved
fp_analyzed_data = '/Users/wwneuro/My_Drive/Lab/Data_analysis/slice_NT/area_under_curve';

%name of the saved file
filename = 'Area_Under_Curve_fI_WT_Light_Dark_IE.mat';

%save results
save_file = 1;

% names of experimental conditions within the analyzed fI mat file
condn = {'Dark','Light_6h','Dark_6h'};
%field names
fn = {'MFR', 'IFR', 'mean_IFR'}; 

%number of traces
cti = 15;

%current steps mode
curr_mo = 2;
%1- all current steps saved in the "current_inj" vector
%2- current steps saved under each condition in the "curr_c.curr_inj" field

    
%% 
data_temp = cell(1,numel(condn));
dst = cell(1,numel(condn));
cell_num = NaN(numel(condn),1);

for gi = 1:numel(condn) %per condition
    %area_names{1,gi} = strcat(condn{1,gi},'_area');
    
   for fii = 1:numel(fn) %per field
       cell_num(gi,1) = size(eval(strcat(condn{gi},'.',fn{fii})),2);
       current_data = eval(strcat(condn{gi},'.',fn{fii}));
       
       if curr_mo == 2
           current_inj = eval(strcat(condn{gi},'.','curr_inj'));
       end
       
       for cii = 1:max(cell_num(gi)) %per cell
      
           if isnan(current_data(:,cii))
               data_temp{1,gi}{1,fii}(cii,1) = NaN;
           else
               %trace_ct = size(current_data,1);
               trace_ct = cti;
               if curr_mo == 2
                  data_temp{1,gi}{1,fii}(cii,1) = areaundercurve(current_inj(1:trace_ct,cii),current_data(1:trace_ct,cii));
               else
                  data_temp{1,gi}{1,fii}(cii,1) = areaundercurve(current_inj(1:trace_ct,1),current_data(1:trace_ct,cii));
               end
           end
          
       end
       
       dst{gi}.(fn{fii}) = data_temp{1,gi}{1,fii}(1:cell_num(gi,1));
       
   end
end

%save data to corresponding data structure (by experimental condition)
%'P_T_NT_Ctrl','PhTX_24h','TTX_24h','PhTX_TTX_24h'
% P_T_NT_Ctrl_area = dst{1};
% PhTX_24h_area = dst{2};
% TTX_24h_area = dst{3};
% PhTX_TTX_24h_area = dst{4};
%sTNFR_TTX_6h_area = dst{4};
Dark_area = dst{1};
Light_6h_area = dst{2};
Dark_6h_area = dst{3};
%% save file
if save_file == 1
    cd(fp_analyzed_data)
    
    condn_area = cell(1,numel(condn));
    for cdi = 1:numel(condn)
        condn_area{cdi} = strcat(condn{cdi},'_area');
    end
        
    save(filename,condn_area{1,:},'condn')
end
