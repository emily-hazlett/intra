%% Description
% This script is for automatically dropping reps that occur before
% penetration of a cell or after seal degrades, and then saving .txt files
% of each stimulus. Data should be saved as .txt files with each column
% being adifferent presentation of the stimulus and each row being a sample.
% .Txt file can contain header data before the traces begin, but the number
% of numeric headers must be specified.  This script will batch through all
% files matching a specified naming structure in a specifies folder and
% create and save a .tiff image of the summary figure.

% Created by EHazlett, last edited 2017-10-13

clear all
folderold = cd;
%% User editted info
% cd('D:\Intracellular\data analysis\processing\sub\'); % Look for files in this folder
cd('C:\Data Processing\Processing\sub\'); % Look for files in this folder
Files = dir('0000*.txt'); % Find txt files containing this phrase to batch through

prestim = 300; %ms in sweep before stim onset
poststim = 700; %ms in sweep after stim onset
background = 200; %ms starting from end of poststim window to calculate rmp
headers = 4; % number of rows containing numeric data in ascii file before the traces start

%% Batch through all files in the folder
for ii = 1:length(Files)
    %% Import data - add matching with stim
    %     filename = Files(ii).name;
    %     testname = ['x', strrep(filename, '.txt','')];
    %     traces = importdata(filename);
    %     traces.textdata = traces.textdata(1,2:end);
    %     stimList = unique(traces.textdata);
    %     celldata = (traces.data(headers+1:end,:))/10;
    filename = Files(ii).name;
    trace = dlmread(filename, '\t', 0, 1);
    Pulses = struct('trace', trace(headers+1:end,:), ...
        'voltageInjected', trace(4,:), ...
        'timestamp', trace(3,:));
    clear trace
    
    %% Find stable state of IR test
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
    
    Pulses.increasing = Pulses.traceChange > .1;
    Pulses.increasing(1:floor(reps/2),:) = false; %Only include increasing reps from second half of recording
    Pulses.decreasing = Pulses.traceChange < -.1;
    Pulses.decreasing(floor(reps/2):end,:) = false; %Only include decreasing reps from first half
    
    Pulses.todrop = false(samples,reps);
    for kk = 1:reps
        if ~isempty(find(Pulses.increasing(:,kk), 1, 'first'))
            Pulses.todrop(find(Pulses.increasing(:,kk), 1, 'first'):end, kk) = true;
        end
        if ~isempty(find(Pulses.decreasing(:,kk), 1, 'last'))
            Pulses.todrop(1:find(Pulses.decreasing(:,kk), 1, 'last'), kk) = true;
        end
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
    count = 1;
    for kk = 1:2:reps-1
        IRTest.time1(count) = Pulses.timestamp(1,kk);
        IRTest.time2(count) = Pulses.timestamp(1,kk+1);
        IRTest.voltage(count) = Pulses.voltage(kk+1) - Pulses.voltage(kk);
        IRTest.current(count) = 10/(IRTest.voltage(count)/1000);
        IRTest.resistance(count) = 10000000*(IRTest.voltage(count)/1000)/(IRTest.current(count));
        count = count+1;
    end
    
    %% Group with stim
    %     StimIRTest.stimname{1:48} = 'peaches';
    %     StimIRTest.stimtime(1:48) = 1;
    
    count = 1;
    pp = 1;
    for kk = 2:length(IRTest.time2)
        if IRTest.time2(kk)-IRTest.time2(kk-1) > 50000
            StimIRTest.meanDeltaV(count) = mean(IRTest.voltage(pp:kk));
            StimIRTest.meanR(count) = mean(IRTest.resistance(pp:kk));
            %         StimIRTest.IR(count) = IRTest(kk:pp);
            pp = kk;
            count = count+1;
        end
    end
end