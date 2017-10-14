%% Description
% This script is for automatically dropping reps that occur before
% penetration of a cell or after seal degrades, and then saving .txt files
% of each stimulus. Data should be saved as .txt files with each column
% being adifferent presentation of the stimulus and each row being a sample.
% .Txt file can contain header data before the traces begin, but the number
% of numeric headers must be specified.  This script will batch through all
% files matching a specified naming structure in a specifies folder and
% create and save a .tiff image of the summary figure.

% Created by EHazlett, last edited 2017-10-14
clear all
close all

folderold = cd;
%% User editted info
cd('C:\Data Processing\Processing\sub\'); % Look for files in this folder
Files = dir('*All_trace.txt'); % Find txt files containing this phrase to batch through

sweeplength = 1000; %ms in sweep
background = 200; %ms from end of sweep to calculate rmp
headers = 3; % number of rows containing numeric data in ascii file before the traces start

%% Batch through all files in the folder
for ii = 5:5%:length(Files)
    %% Import data
    filename = Files(ii).name;
    traces = importdata(filename);
    Reps.stim = traces.textdata(1,2:end);
    Reps.pulsepolarity = traces.data(1,:);
    Reps.pulsevoltage = traces.data(2,:);
    Reps.trace = (traces.data(headers+1:end,:))/10;
    clear traces
    
    %% RMP distribution
    [samples, reps] = size(Reps.trace);
    Reps.background = Reps.trace(samples - floor(background*samples/1000):end,:);
    
    for pp = 1:reps
        repper = Reps.background(:,pp);
        Reps.rmp(pp) = round(mean2(repper(repper>mean2(repper) - std2(repper) ...
            & repper<mean2(repper) + std2(repper))),1); % round to 1 decimal place
        clear repper
    end
    
    figure
    histogram(Reps.trace)
    title(strrep(filename,'.txt',''),'Interpreter','none','Fontsize',12)
    
    %% Select reps to be dropped
    smoothBin = 3; % How many reps do you want to smooth over
    rmpSmooth = zeros(1,reps-smoothBin);
    for pp = smoothBin+1:reps
        rmpSmooth(pp-smoothBin) = mean(Reps.rmp(pp-smoothBin:pp));
        if pp+smoothBin >reps
            continue
        end
        Reps.rmpChange(pp) = mean(Reps.rmp(pp:pp+smoothBin)) - mean(Reps.rmp(pp-smoothBin:pp)); % pos means rmp is increasing
    end
    
    Reps.increasing = find(Reps.rmpChange > 4);
    Reps.increasing(Reps.increasing<reps/2) = []; %Only include increasing reps from second half of recording
    Reps.decreasing = find(Reps.rmpChange < -4);
    Reps.decreasing(Reps.decreasing>reps/2) = []; %Only include decreasing reps from first half
    Reps.badrmp = Reps.rmp > -45;
    
    Reps.todrop = false(1,reps);
    if ~isempty(Reps.increasing)
        Reps.todrop(min(Reps.increasing):end) = true;
    end
    if ~isempty(Reps.decreasing)
        Reps.todrop(1:max(Reps.decreasing)) = true;
    end
    if ~isempty(Reps.badrmp)
        Reps.todrop(Reps.badrmp) = true;
    end
    
    %% Plot and save figure
    figure
    ax(1) = subplot(3,1,1);
    histogram(Reps.rmp, 'BinWidth', 1)
    title(filename,'Interpreter','none','Fontsize',12)
    
    ax(2) = subplot(3,1,2);
    hold on
    plot(Reps.todrop, 'g', 'LineWidth', 2)
    plot(Reps.rmpChange, 'r')
    ylim([-5,5])
    
    ax(3) = subplot(3,1,3);
    hold on
    plot(Reps.rmp)
    plot(ceil(smoothBin/2):length(rmpSmooth)+ceil(smoothBin/2)-1, rmpSmooth, 'LineWidth', 2)
    
    testname = ['x', strrep(filename, '.txt','')];
    print('-dtiff','-r500',[testname,'_dropreps.tif'])
    
    %% Drop Reps
    stimList = unique(Reps.stim);
    Reps.stim(:,Reps.todrop) = [];
    Reps.trace(:,Reps.todrop) = [];
    
    for pp = 1:length(stimList)
        txtname = strrep(filename, 'All_trace', [stimList{pp},'_trace_1']);
        grouper = strcmp(stimList{pp}, Reps.stim);
        if ~isempty(Reps.stim(grouper))
            dater = Reps.pulsepolarity(:,grouper);
            dater = [dater; Reps.pulsevoltage(:,grouper)];
            dater = [dater; Reps.trace(:,grouper)];
            fid = fopen(txtname, 'w');
            fprintf(fid,[repmat('%s \t',1,sum(grouper)-1), '%s'],string(Reps.stim(grouper)));
            fprintf(fid,'\n');
            fclose(fid);
            dlmwrite(txtname, dater, '-append', 'delimiter', '\t','roffset', 0);
            clear dater grouper txtname fid
        end
    end
    clear Reps reps pp ii samples
end
clear headers *stim background
cd(folderold);