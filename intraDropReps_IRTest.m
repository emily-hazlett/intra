%% Description
% This script is for automatically dropping reps that occur before
% penetration of a cell or after seal degrades, and then saving .txt files
% of each stimulus. Data should be saved as .txt files with each column
% being adifferent presentation of the stimulus and each row being a sample.
% .Txt file can contain header data before the traces begin, but the number
% of numeric headers must be specified.  This script will batch through all
% files matching a specified naming structure in a specifies folder and
% create and save a .tiff image of the summary figure.

% Created by EHazlett, last edited 2017-12-01
close all
clear all

folderold = cd;
%% User editted info
cd('C:\Data Processing\Processing\sub\1195\'); % Look for files in this folder
Files = dir('1195_112817_3558_1_All_trace.txt'); % Find txt files containing this phrase to batch through

sweeplength = 1000; %ms in sweep
background = 50; %ms from beginning of sweep to calculate rmp
smoothBin = 3; % How many reps do you want to smooth over to find dropreps?
badrmp = -30; % Drop reps above this threshold
deltaV = 5; % change in V for figuring out which reps to drop
headers = 5; % number of rows containing numeric data in ascii file before the traces start
IRTest = 0; %Does this have IR test data?

%% Batch through all files in the folder
for ii = 1:length(Files)
    %% Import data
    filename = Files(ii).name;
    traces = importdata(filename);
    Reps.stim = strrep(traces.textdata(1,:),' ','');
    Reps.pulsepolarity = traces.data(1,:);
    Reps.pulsevoltage = traces.data(2,:);
    Reps.timestamp = traces.data(3,:);
    if IRTest == 1
        Reps.preIR = traces.data(4,:);
        Reps.postIR = traces.data(5,:);
    else
        Reps.preIR = NaN(size(traces.data(3,:)),'double'); % Add in NaN for IR data
        Reps.postIR = NaN(size(traces.data(3,:)),'double');
    end
    Reps.trace = (traces.data(headers+1:end,:))/10; %DW increases mV by 10
    clear traces
    
    %% RMP distribution
    [samples, reps] = size(Reps.trace);
    Reps.background = Reps.trace(1:floor(background*samples/1000),:); % Chunk of trace for calc bg
    
    for pp = 1:reps
        repper = Reps.background(:,pp);
        Reps.rmp(pp) = round(...
            mean2(repper( ...
            repper < mean2(repper) + std2(repper) )) ...
            ,1); % round to 1 decimal place
        clear repper
    end
    
%     figure
%     histogram(Reps.trace)
%     title(strrep(filename,'.txt',''),'Interpreter','none','Fontsize',12)
    
    %% Select reps to be dropped because of rmp
    % Smooth rmp over several reps
    rmpSmooth = zeros(1,reps-smoothBin);
    for pp = smoothBin+1:reps
        rmpSmooth(pp-smoothBin) = mean(Reps.rmp(pp-smoothBin:pp));
        if pp+smoothBin >reps
            continue
        end
        Reps.rmpChange(pp) = mean(Reps.rmp(pp:pp+smoothBin)) - mean(Reps.rmp(pp-smoothBin:pp)); % pos means rmp is increasing
    end
    
    % Find where rmp isn't stable (e.g., break in, lose cell)
    Reps.rmpDecreasing = find(Reps.rmpChange < -deltaV);
    Reps.rmpDecreasing = Reps.rmpDecreasing(lt(Reps.rmpDecreasing,reps/2));
    
    Reps.rmpIncreasing = find(Reps.rmpChange > deltaV);
    Reps.rmpIncreasing = Reps.rmpIncreasing(gt(Reps.rmpIncreasing,reps/2)); %Only include increasing reps from second half of recording
    
    % Check if rmp is not sufficiently hyperpolarized
    % Use rmpSmooth, but add missing reps from beginning and end that
    % were lost due to the smooth bin
    repper = gt(...
        [repmat(mean(Reps.rmp(1:floor(smoothBin/2))),1,floor(smoothBin/2)), ...
        rmpSmooth, ...
        repmat(mean(Reps.rmp(reps-ceil(smoothBin/2)+1:end)),1,ceil(smoothBin/2))] ...
        , badrmp);
    Reps.badrmpBeginning = find(repper);
    Reps.badrmpBeginning = Reps.badrmpBeginning(lt(Reps.badrmpBeginning,reps/2));
    
    Reps.badrmpEnding = find(repper);
    Reps.badrmpEnding = Reps.badrmpEnding(gt(Reps.badrmpEnding,reps/2));
    clear repper
    
    %% Logicals of reps to drop based on rmp
    Reps.todrop = false(1,reps);
    
    % If rmp starts increasing too much, then you're losing the cell.
    % Mark rest of reps to end to be dropped.
    if ~isempty(Reps.rmpIncreasing)
        Reps.todrop(min(Reps.rmpIncreasing):end) = true;
    end
    
    % If rmp decreases, then you're breaking into a cell
    % Mark reps before this to be dropped.
    if ~isempty(Reps.rmpDecreasing)
        Reps.todrop(1:max(Reps.rmpDecreasing)) = true;
    end
    
    % If rmp drifts to not sufficiently polarized
    % Makr reps to be dropped before/ after then.
    if ~isempty(Reps.badrmpEnding)
        Reps.todrop(min(Reps.badrmpEnding):end) = true;
    end
    if ~isempty(Reps.badrmpBeginning)
        Reps.todrop(1:max(Reps.badrmpBeginning)) = true;
    end
    
    %% Drop Reps based on rmp
    Reps.stim(strcmp(Reps.stim,'StimName')) = [];
    Reps.stim(:,Reps.todrop) = [];
    Reps.pulsepolarity(:,Reps.todrop) = [];
    Reps.pulsevoltage(:,Reps.todrop) = [];
    Reps.timestamp(:,Reps.todrop) = [];
    Reps.preIR(:,Reps.todrop) = [];
    Reps.postIR(:,Reps.todrop) = [];
    Reps.trace(:,Reps.todrop) = [];
    
    %% Spike detect
    rmp = mean2(Reps.rmp); % calculate global rmp now that reps have been dropped
    rmpSD = std2(Reps.rmp);
    
    Reps.spikeIndex = cell(1,reps);
    Reps.spikeHeight = cell(1,reps);
    threshold = min(rmp + 30, -10); % Set spike threshold to the lesser of rmp+30 and =10mV
    for i = 1:reps
        if any(Reps.trace(:,i)>threshold)
            [Reps.spikeHeight{i}, Reps.spikeIndex{i}] = findpeaks(Reps.trace(:,i), ...
                'MinPeakHeight',threshold, ...
                'MinPeakDistance', ceil(1.5/(1000/samples)) ...
                ); %Find spike peaks that are <10 SD above rmp and with a hold time of ~1s (assuming a 1s sweep)
        end
    end
    
    %% Find changes in spike height
    Reps.meanHeight = cellfun(@mean,Reps.spikeHeight) - rmp;
    
    spikeSmooth = zeros(1,reps-smoothBin);
    for pp = smoothBin+1:reps
        spikeSmooth(pp-smoothBin) = mean(Reps.meanHeight(pp-smoothBin:pp));
        if pp+smoothBin >reps
            continue
        end
        Reps.spikeChange(pp) = mean(Reps.meanHeight(pp:pp+smoothBin)) - mean(Reps.meanHeight(pp-smoothBin:pp)); % pos means rmp is increasing
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
    

    
    %% Save txt files for each stim
    stimList = unique(Reps.stim);
    for pp = 1:length(stimList)
        grouper = strcmp(stimList{pp}, Reps.stim); % drop reps for this stim
        if isempty(Reps.stim(grouper)) %Don't do anything if all reps are dropped
            continue
        end
        stimreps = Reps.pulsepolarity(:,grouper);
        stimreps = [stimreps; Reps.pulsevoltage(:,grouper)];
        stimreps = [stimreps; Reps.timestamp(:,grouper)];
        stimreps = [stimreps; Reps.preIR(:,grouper)];
        stimreps = [stimreps; Reps.postIR(:,grouper)];
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
                
                % Add depolarization/ hyperpolarization to txtfile name
                if stimrepsX(2,1) ~= 0
                    if stimrepsX(1,1) > 0
                        txtname = strrep(filename, 'All_trace', [stimList{pp},'_trace', '_', num2str(count), '_depol']);
                    elseif stimrepsX(1,1) < 0
                        txtname = strrep(filename, 'All_trace', [stimList{pp},'_trace', '_', num2str(count), '_hyper']);
                    end
                else
                    txtname = strrep(filename, 'All_trace', [stimList{pp},'_trace', '_', num2str(count)]);
                end
                
                fid = fopen(txtname, 'w');
                fprintf(fid,[repmat('%s \t',1,length(texterX)-1), '%s'],texterX{:});
                fprintf(fid,'\n');
                fclose(fid);
                dlmwrite(txtname, stimrepsX, '-append', 'delimiter', '\t','roffset', 0);
                clear stimrepsX* texterX*
                count = count + 1;
            end
        else % save one text file if there's only one test
            % Add depolarization/ hyperpolarization to txtfile name
            if stimreps(2,1) ~= 0
                if stimreps(1,1) > 0
                    txtname = strrep(filename, 'All_trace', [stimList{pp},'_trace', '_', num2str(count), '_depol']);
                elseif stimreps(1,1) < 0
                    txtname = strrep(filename, 'All_trace', [stimList{pp},'_trace', '_', num2str(count), '_hyper']);
                end
            else
                txtname = strrep(filename, 'All_trace', [stimList{pp},'_trace', '_', num2str(count)]);
            end
            
            fid = fopen(txtname, 'w');
            fprintf(fid,[repmat('%s \t',1,length(texter)-1), '%s'],texter{:});
            fprintf(fid,'\n');
            fclose(fid);
            dlmwrite(txtname, stimreps, '-append', 'delimiter', '\t','roffset', 0);
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