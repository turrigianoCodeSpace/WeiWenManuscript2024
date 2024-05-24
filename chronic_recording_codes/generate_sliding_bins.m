function [all_bins] = generate_sliding_bins(bin_start, total_length, bin_size, episode, window)

%This function generates bins within an episode and then slides across a
%period of time (size = total_length) by a window to produce a train of
%bins with equal length.

%bin_start- starting time, in hours
%total_length- length of the entire time that needs to be binned, in hours
%bin_size- width of the bin, in seconds
%episode- length of an episode, in minutes
%window- size of the sliding window, in minutes


bin_num = 1; 
bin_e = 0:episode:episode;
while bin_e(end) < total_length*60
    bin_e = bin_e+window;
    bin_num = bin_num+1;
end


all_bins = cell(1,bin_num);

bin_ct = 0;
bins = 0:bin_size:episode*60; %for a single episode

for bi = 1:bin_num
    all_bins{bi} = bins+bin_start*3600;
    bins = bins + window*60;
    bin_ct = bin_ct + 1;
end

end