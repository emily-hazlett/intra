%% Description
% This script is for automatically calculating the Input resistance before
% and after a stimulus presentation, and then saving a .txt files with
% IR data. Data should be saved as .txt files with each column being a
% different presentation of the stimulus and each row being a sample.
% .Txt file can contain header data before the traces begin, but the number
% of numeric headers must be specified.  This script will batch through all
% files matching a specified naming structure in a specifies folder.

% Created by EHazlett, last edited 2017-11-06
clear all
folderold = cd;
%% User editted info
% cd('D:\Intracellular\data analysis\processing\sub\'); % Look for files in this folder
cd('C:\Data Processing\Processing\sub\'); % Look for files in this folder
Files = dir('*000005_IRTest.txt'); % Find txt files containing this phrase to batch through

prestim = 100; %ms in sweep before stim onset
poststim = 900; %ms in sweep after stim onset
background = 200; %ms starting from end of poststim window to calculate rmp
headers = 4; % number of rows containing numeric data in ascii file before the traces start

%% Batch through all files in the folder
for ii = 1:length(Files)
    %% Import data - add matching with stim
    % IR txt file
    filename = Files(ii).name;
    trace = dlmread(filename, '\t', 0, 1);
    Pulses = struct('trace', trace(headers+1:end,:), ...
        'voltageInjected', trace(3,:), ...
        'timestamp', trace(1,:));
    clear trace
    
    % Stim trace txt file
    stimfilename = strrep(filename, 'IRTest', 'All_trace');
    if ~isempty(dir(stimfilename))
        traces = importdata(stimfilename);
        Stims.StimName = traces.textdata(1,2:end);
        Stims.PulsePolarity = traces.data(1,:);
        Stims.PulseAmplitude = traces.data(2,:);
        Stims.OutputTimestamp = traces.data(3,:)/1000;
        Stims.trace = traces.data(4:end,:);
        
        testname = ['x', strrep(stimfilename, '.txt','')];
        clear traces
    end
    
    %% Find how pulse changes over course of pulse
    [samples, reps] = size(Pulses.trace);%(headers+1:end, :));
    smoothBin = 3; % How many reps do you want to smooth over
    
    Pulses.traceChange = zeros(samples, reps);
    for kk = 1:reps
        for pp = smoothBin+1:samples
            Pulses.traceSmooth(pp, kk) = mean(Pulses.trace(pp-smoothBin:pp, kk));
            if pp+smoothBin > samples
                continue
            end
            Pulses.traceChange(pp, kk) = mean(Pulses.trace(pp:pp+smoothBin, kk)) - mean(Pulses.trace(pp-smoothBin:pp, kk)); % pos means IR is increasing
        end
    end
    
    %% Cut header and tail from pulse traces
    Pulses.pulser1 = Pulses.traceChange < -0.5 | Pulses.traceChange > 0.5;
    Pulses.pulser2 = Pulses.pulser1;
    Pulses.pulser1(floor(samples/3):end,:) = false; %Only include samples from first third of pulse
    Pulses.pulser2(1:2*floor(samples/3),:) = false; %Only include samples from last third of pulse
    
    %% Find stable state of IR test
    %Need real data to logic
%     Pulses.stable = Pulses.traceChange > .1 | Pulses.traceChange > .1;
%     Pulses.stable1(1:floor(samples/2),:) = false; %Only include increasing samples from second half of pulse
%     Pulses.stable2(floor(samples/2):end,:) = false; %Only include decreasing samples from first half of pulse
    
    % Label samples to drop by rep to mantain matrix shape
    Pulses.todrop = false(samples,reps);
    for kk = 1:reps
        if ~isempty(find(Pulses.pulser1(:,kk), 1, 'last'))
            Pulses.todrop(1:find(Pulses.pulser1(:,kk), 1, 'last'), kk) = true;
        end
        if ~isempty(find(Pulses.pulser2(:,kk), 1, 'first'))
            Pulses.todrop(find(Pulses.pulser2(:,kk), 1, 'first'):end, kk) = true;
        end
%         if ~isempty(find(Pulses.stable1(:,kk), 1, 'first'))
%             Pulses.todrop(find(Pulses.stable1(:,kk), 1, 'first'):end, kk) = true;
%         end
%         if ~isempty(find(Pulses.stable2(:,kk), 1, 'last'))
%             Pulses.todrop(1:find(Pulses.stable2(:,kk), 1, 'last'), kk) = true;
%         end
    end
    
    %% Find mean value of trace
    for kk = 1:reps
        repper = Pulses.trace(:,kk);
        if any(Pulses.todrop(:,kk))
            Pulses.voltage(kk) = mean(repper(Pulses.todrop(:,kk)));
        else
            Pulses.voltage(kk) = mean(repper);
        end
        clear repper
    end
    
    %% Find change in V during test
    count = 1; %%%%%CHECK UNITS%%%%%%%
    for kk = 1:2:reps-1
        IRTest.time1(count) = Pulses.timestamp(1,kk);
        IRTest.time2(count) = Pulses.timestamp(1,kk+1);
        IRTest.voltage(count) = Pulses.voltage(kk+1) - Pulses.voltage(kk); %mV
        IRTest.current(count) = 10/(IRTest.voltage(count)/1000);
        IRTest.resistance(count) = round(10000000*(IRTest.voltage(count)/1000)/(IRTest.current(count)),1);
        count = count+1;
    end
    
    %% Group with stim
    nPulse = find(diff(IRTest.time1)> 800, 1); % number of pulses in each batch of IR test
    for kk = 1:length(Stims.OutputTimestamp)
        if strcmp(Stims.StimName(kk),'None')
            Stims.PreIR(kk) = NaN;
            Stims.PostIR(kk) = NaN;
            continue
        end
        chunker = find(IRTest.time1 > Stims.OutputTimestamp(kk),1);
        Stims.PreIR(kk) = mean(IRTest.resistance(chunker-nPulse:chunker-1));
        Stims.PostIR(kk) = mean(IRTest.resistance(chunker:chunker+nPulse-1));
    end
    
    %% Save txt file with IR data
    if ~isempty(dir(stimfilename))
        txtname = strrep(stimfilename, 'All_trace', 'withIR_All_trace');
        fid = fopen(txtname, 'w');
        fprintf(fid,[repmat('%s \t',1,length(Stims.StimName)-1), '%s'],Stims.StimName{:});
        fprintf(fid,'\n');
        fclose(fid);
        dlmwrite(txtname, Stims.PulsePolarity, '-append', 'delimiter', '\t','roffset', 0);
        dlmwrite(txtname, Stims.PulseAmplitude, '-append', 'delimiter', '\t','roffset', 0);
        dlmwrite(txtname, Stims.OutputTimestamp, '-append', 'delimiter', '\t','roffset', 0);
        dlmwrite(txtname, Stims.PreIR, '-append', 'delimiter', '\t','roffset', 0);
        dlmwrite(txtname, Stims.PostIR, '-append', 'delimiter', '\t','roffset', 0);
        dlmwrite(txtname, Stims.trace, '-append', 'delimiter', '\t','roffset', 0);
    end
end