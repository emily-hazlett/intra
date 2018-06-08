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
close all
clear all

folderold = cd;
%% User editted info
cd('C:\Data Processing\Processing\'); % Look for files in this folder
Files = dir('1212*_All*.txt'); % Find txt files containing this phrase to batch through

sweeplength = 1000; %ms in sweep
background = 500; %ms from end of sweep to calculate rmp
smoothBin = 3; % How many reps do you want to smooth over to find dropreps?
badrmp = -40; % Drop reps above this threshold
deltaV = 4; % change in V for figuring out which reps to drop
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
    
%     figure
%     histogram(Reps.trace)
%     title(strrep(filename,'.txt',''),'Interpreter','none','Fontsize',12)
    
    %% Select reps to be dropped
    rmpSmooth = zeros(1,reps-smoothBin);
    for pp = smoothBin+1:reps
        rmpSmooth(pp-smoothBin) = mean(Reps.rmp(pp-smoothBin:pp));
        if pp+smoothBin >reps
            continue
        end
        Reps.rmpChange(pp) = mean(Reps.rmp(pp:pp+smoothBin)) - mean(Reps.rmp(pp-smoothBin:pp)); % pos means rmp is increasing
    end
    
    Reps.decreasing = find(Reps.rmpChange < -deltaV);
    Reps.decreasing = Reps.decreasing(lt(Reps.decreasing,reps/2)); %Only include decreasing reps from first half
    Reps.increasing = find(Reps.rmpChange > deltaV);
    Reps.increasing = Reps.increasing(gt(Reps.increasing,reps/2)); %Only include increasing reps from second half of recording
    
    repper = gt([repmat(mean(Reps.rmp(1:floor(smoothBin/2))),1,floor(smoothBin/2)), ...
        rmpSmooth, repmat(mean(Reps.rmp(reps-ceil(smoothBin/2)+1:end)),1,ceil(smoothBin/2))], badrmp);
    Reps.badrmpFront = find(repper);
    Reps.badrmpFront = Reps.badrmpFront(lt(Reps.badrmpFront,reps/2));
    Reps.badrmpBack = find(repper);
    Reps.badrmpBack = Reps.badrmpBack(gt(Reps.badrmpBack,reps/2));
    clear repper
    
    Reps.todrop = false(1,reps);
    if ~isempty(Reps.increasing)
        Reps.todrop(min(Reps.increasing):end) = true;
    end
    if ~isempty(Reps.decreasing)
        Reps.todrop(1:max(Reps.decreasing)) = true;
    end
    if ~isempty(Reps.badrmpBack)
        Reps.todrop(min(Reps.badrmpBack):end) = true;
    end
    if ~isempty(Reps.badrmpFront)
        Reps.todrop(1:max(Reps.badrmpFront)) = true;
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
        
        % Where do different pulse tests occur
        indPol = find(diff(stimreps(1,:))~=0);
        indVolt = find(diff(stimreps(2,:))~=0);
        indTests = unique(sort([reshape(indPol,[],1);reshape(indVolt,[],1)])); % reshape forces col vector
        
        % Find and save separate tests
        count = 1;
        ff = 1;
        if ~isempty(indTests) % Save a different txt file for each test
            indTests = sort([0;reshape(indTests,[],1); length(texter)]);
            for kk = 1:length(indTests)-1
                stimrepsX = stimreps(:,indTests(kk)+1:indTests(kk+1));
                texterX = texter(:,indTests(kk)+1:indTests(kk+1));
                txtname = strrep(filename, 'All_trace', [stimList{pp},'_trace', '_', num2str(count)]);
                fid = fopen(txtname, 'w');
                fprintf(fid,[repmat('%s \t',1,length(texterX)-1), '%s'],texterX{:});
                fprintf(fid,'\n');
                fclose(fid);
                dlmwrite(txtname, stimrepsX, 'precision','%f', '-append', 'delimiter', '\t','roffset', 0);
                clear stimrepsX* texterX*
                count = count + 1;
            end
        else % save one text file if there's only one test
            txtname = strrep(filename, 'All_trace', [stimList{pp},'_trace', '_', num2str(count)]);
            fid = fopen(txtname, 'w');
            fprintf(fid,[repmat('%s \t',1,length(texter)-1), '%s'],texter{:});
            fprintf(fid,'\n');
            fclose(fid);
            dlmwrite(txtname, stimreps, 'precision','%f', '-append', 'delimiter', '\t','roffset', 0);
        end
        clear grouper stimreps* texter* ind*
    end
    disp([num2str(ii), '/', num2str(length(Files))]);
    clear Reps reps pp ii samples
end
clear headers *stim background
cd(folderold);



%% Old junk

% % Split based on differences in timestamp
% count = 1; % Save txt file for each test run
% 
% if count == 1 %Still save single txt file if only 1 test was run
%     txtname = strrep(filename, 'All_trace', [stimList{pp},'_trace', '_', num2str(count)]);
%     fid = fopen(txtname, 'w');
%     fprintf(fid,[repmat('%s \t',1,length(texter)-1), '%s'],string(texter));
%     fprintf(fid,'\n');
%     fclose(fid);
%     dlmwrite(txtname, stimreps, '-append', 'delimiter', '\t','roffset', 0);
% end
% clear dater grouper txtname fid