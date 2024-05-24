function [CI_coor] = generate_sem_CI(data,sem)

% 1st col- lower bound;
% 2nd col- upper bound;
CI_coor = NaN(size(data,1),2);

%  define coordinates for the polygon (one polygon every two points)
for i = 1:size(data,1)
    CI_coor(i,1) = data(i,1)-sem(i,1);
    CI_coor(i,2) = data(i,1)+sem(i,1);
end
end