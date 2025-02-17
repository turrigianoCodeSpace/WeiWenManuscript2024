%% Select which cell and which trace to plot

% where to save the selected trace data
save_data_path = '/Users/wwneuro/My_Drive/Lab/Data_analysis/slice_NT/representative_traces/mini';

% name of the saved data
save_data_name = 'Dark_6h.mat';

% data folder
fp_data = '/Users/wwneuro/My_Drive/Lab/Data/slice_NT/mini/';

%whether to save file
save_file = 1;

%Picking traces
experiment = '231213';

trace_file = 'cell1_0009.h5';

sweep = 9;

%Picking regions (in seconds)
start_time = 27.5;

end_time = 29.5;

%% data readout and plotting
extracted_data_file = ws.loadDataFile(strcat(fp_data,experiment,'/',trace_file));
fields = fieldnames(extracted_data_file);  %get the fieldnames of the data struct
sprate = extracted_data_file.header.AcquisitionSampleRate;

sweep_start = str2double(trace_file(7:10));

selected_trace_data = extracted_data_file.(fields{sweep-sweep_start+2}).analogScans(:,1);

time_increment = 1/sprate;
number_of_datapoints = size(selected_trace_data,1); %get the number of data points of the sweep

%pick the interval you want to show (in seconds)
startpoint = start_time * sprate;
endpoint = end_time * sprate;

timepoints = (start_time:time_increment:end_time);

plot_data = detrend(selected_trace_data(startpoint:endpoint,1));

figure('position',[56 200 700 490]);
plot(timepoints,plot_data,'k-','Linewidth',1);
hold on
xlabel('Time (s)');
xlim([start_time end_time]);
ylim([-100 50])
    
%draw scale
plot([start_time+0.5; start_time+0.55],[-80; -80], '-k',[start_time+0.5;start_time+0.5],[-80; -70], '-k', 'LineWidth',2)
text(start_time+0.49, -76, '10 pA', 'HorizontalAlignment', 'right')
text(start_time+0.52, -84, '50 ms', 'HorizontalAlignment', 'center')
%set(gca, 'Visible', 'off')


%% save data
if save_file == 1
    cd(save_data_path)
    save(save_data_name, 'experiment','trace_file','start_time','end_time')
end