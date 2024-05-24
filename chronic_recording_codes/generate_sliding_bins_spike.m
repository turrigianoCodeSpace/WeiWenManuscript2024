%%% generate a series of bins with a sliding window

% bin- bin size, in seconds
% sw- sliding window size, in seconds
% s_start- starting time of the spike train, in seconds
% s_end- end time of the spike train, in seconds

function [eg_le_v,eg_tr_v] = generate_sliding_bins_spike(s_start,s_end,bin,sw)

bin_ct = 1;

eg_le = s_start;
eg_tr = s_start+bin;

while eg_tr < s_end
    eg_le = eg_le+sw;
    eg_tr = eg_tr+sw;

    bin_ct = bin_ct+1;
end
 
% leading and trailing edges
eg_le_v = NaN(1,bin_ct-1);
eg_tr_v = NaN(1,bin_ct-1);

for bi = 1:bin_ct-1
    eg_le_v(bi) = s_start+sw*(bi-1);
    eg_tr_v(bi) = eg_le_v(bi)+bin;
end

end