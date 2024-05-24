
%animal 6 and 9 both have 5 RS cells

c_table = pairwise_correlation_LD_CPP;
row_id_6 = zeros(size(c_table,1),1);
row_id_9 = zeros(size(c_table,1),1);

for ri = 1:size(c_table,1)
    if c_table.animalID(ri) == 6 && c_table.pair_type(ri) == '1'
        row_id_6(ri,1) = true;
    else
        row_id_6(ri,1) = false;
    end
end

row_id_6 = logical(row_id_6);

f_table_6 =  pairwise_correlation_LD_CPP(row_id_6,:);

for ri = 1:size(c_table,1)
    if c_table.animalID(ri) == 9 && c_table.pair_type(ri) == '1'
        row_id_9(ri,1) = true;
    else
        row_id_9(ri,1) = false;
    end
end

row_id_9 = logical(row_id_9);

f_table_9 =  pairwise_correlation_LD_CPP(row_id_9,:);

%%  generate matrix for heatmap
pc_mtrx = cell(1,3);
cell_ID = [218 220 223 224 226];
curr_table = f_table_9;
%%

cond = 1;
for pi = 1:5
    c_row = cell_ID(pi);
    for qi = 1:5
        c_col = cell_ID(qi);
        

        if c_row == c_col
            pc_mtrx{cond}(pi,qi) = 0.25;
        else

            f_row_ind = zeros(size(curr_table,1),1);
            for pqi = 1:size(f_table_9,1)
            
             if (curr_table.cell1_ogID(pqi) == c_row && curr_table.cell2_ogID(pqi)== c_col) || ...
                (curr_table.cell2_ogID(pqi) == c_row && curr_table.cell2_ogID(pqi) == c_col)
                 f_row_ind(pqi) = 1;
             else
                 f_row_ind(pqi) = 0;
             end
            end

             
            f_row_ind = logical(f_row_ind);
            if sum(f_row_ind) == 0
                continue
            else

                pc_mtrx{cond}(pi,qi) = curr_table.Dark_raw(f_row_ind,:);
                pc_mtrx{cond}(qi,pi) = curr_table.Dark_raw(f_row_ind,:);
            end
        end
    end
end

Animal9.Dark = pc_mtrx{1};
Animal9.Light = pc_mtrx{2};
Animal9.Light_CPP = pc_mtrx{3};
%% generate heat map from the matrix

figure('Position',[500 100 220 200])  
h = heatmap(Animal9.Light_CPP);
h.Colormap = hot;
h.ColorLimits = [0.01 0.14];
h.ColorMethod = 'median';
h.ColorbarVisible = 'on';



