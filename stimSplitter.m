%% Description
% This script is for automatically separating stim into different .txt files
% Data should be saved as .txt files with each column being adifferent
% presentation of the stimulus and each row being a sample.
% .Txt file can contain header data before the traces begin, but the number
% of numeric headers must be specified.
% Created by EHazlett, last edited 2018-06-08

close all
clear all

folderold = cd;
%% User editted info
cd('C:\Data Processing\Processing\'); % Look for files in this folder
Files = dir('1212*_All_trace_cleaned.txt'); % Find txt files containing this phrase to batch through
headers = 3; % number of rows containing numeric data in ascii file before the traces start

%% Batch through all files in the folder
for ii = 1:length(Files)
    %% Import data
    filename = Files(ii).name;
    traces = importdata(filename);
    Reps.stim = strrep(traces.textdata(1,2:end),' ','');
    Reps.pulsepolarity = traces.data(1,:);
    Reps.pulsevoltage = traces.data(2,:);
    Reps.timestamp = traces.data(3,:);
    Reps.trace = (traces.data(headers+1:end,:))/10;
    stimList = unique(Reps.stim);
    clear traces
    
    %% Save txt files for each stim
    for pp = 1:length(stimList)
        grouper = strcmp(stimList{pp}, Reps.stim); % drop reps for this stim
        if isempty(Reps.stim(grouper)) %Don't do anything if all reps are dropped
            continue
        end
        
        stimreps = Reps.pulsepolarity(:,grouper);
        stimreps = [stimreps; Reps.pulsevoltage(:,grouper)];
        stimreps = [stimreps; Reps.timestamp(:,grouper)];
        stimreps = [stimreps; Reps.trace(:,grouper)];
        texter = Reps.stim(grouper);
        
        txtname = strrep(filename, 'All_trace', [stimList{pp},'_trace']);
        fid = fopen(txtname, 'w');
        fprintf(fid,[repmat('%s \t',1,length(texter)-1), '%s'],texter{:});
        fprintf(fid,'\n');
        fclose(fid);
        dlmwrite(txtname, stimreps, 'precision','%f', '-append', 'delimiter', '\t','roffset', 0);      
        clear grouper stimreps* texter*
    end
    disp([num2str(ii), '/', num2str(length(Files))]);
    clear Reps reps pp ii samples
end
clear headers *stim background
cd(folderold);