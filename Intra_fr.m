% prompt = {'Folder containing files','Files to use'};
% dlg_title = 'Input';
% num_lines = [1, 100];
% defAns = {'C:\Users\emily\Desktop\bigone\good\', '*.mat'};
% options.Resize='on';
% answer = inputdlg(prompt, dlg_title, num_lines, defAns, options);
%
% folderOld = cd;
% folderNew = answer{1};
% cd(folderNew);
% Files = dir(answer{2});
%
% button = questdlg('Start at begining of folder?', ...
%     'Start place', ...
%     'Yes', 'No', 'Yes');
% switch button
%     case 'Yes'
%         ii = 1;
%     case 'No'
%         prompt = {'Which file number to start at?'};
%         dlg_title = 'Input';
%         num_lines = [1, 100];
%         defAns = {['Must be <=', num2str(length(Files))]};
%         options.Resize='on';
%         answer = inputdlg(prompt, dlg_title, num_lines, defAns, options);
%         if isempty(answer)
%             warndlg('Must enter value')
%             return
%         elseif strcmp(answer, ['Must be <=', num2str(length(Files))]) == 1
%             warndlg('Must enter value')
%             return
%         end
%         ii = str2double(answer{1});
% end

folderOld = cd;
folderNew = 'C:\Users\emily\Desktop\bigone\good\';
cd(folderNew);
Files = dir('*.mat');
ii = 1;

% p = figure('units','normalized','outerposition',[0 0 1 1]);
for i = ii:length(Files)
    filename = Files(i).name;
    load(filename);
    
    button = questdlg(filename, ...
        'Anesthesia', ...
        'Ketamine', 'Urethane', 'None', 'Urethane');
    IntraData.anesthesia = button;
    
    save(filename, 'IntraData')
    
    %     scat(i, 1) = IntraData.rmp;
    %     scat(i, 2) = IntraData.FR;
    %     clear IntraData
    %
    %     threshold = (IntraData.rmp + 40)*10;
    %     [spikeheight, spikeindex] = findpeaks(IntraData.trace(IntraData.selectstart:IntraData.selectend), ...
    %         'MinPeakHeight',threshold, ...
    %         'MinPeakDistance', 150); % 150 samples = min distance of 3 ms at Fs=50000
    %     spikeTimestampFull = spikeindex*20+IntraData.timepoints(IntraData.selectstart);
    %
    %     instafr = (60*50000)./diff(spikeTimestampFull); %instantaneous FR in Hz
    %     fr = length(spikeindex)/((IntraData.selectend - IntraData.selectstart)/50000);
    %
    %     subplot (2, 1, 1)
    %     histogram(instafr, 'binwidth', 0.5);
    %
    %     subplot(2, 1, 2)
    %     plot(IntraData.timepoints(IntraData.selectstart:IntraData.selectend), ...
    %         IntraData.trace(IntraData.selectstart:IntraData.selectend))
    %     hold on
    %     line('XData', [IntraData.timepoints(IntraData.selectstart), IntraData.timepoints(IntraData.selectend)], ...
    %         'YData', [IntraData.rmp*10 IntraData.rmp*10], 'linewidth', 2, 'color', 'r')
    %     plot(spikeTimestampFull(2:end), instafr)
    %     scatter(IntraData.markertime,repmat(100, 1, length(IntraData.markertime)))
    %     scatter(spikeTimestampFull, spikeheight)
    %     xlim([IntraData.timepoints(IntraData.selectstart) IntraData.timepoints(IntraData.selectend)])
    %     ylim([-1000 200])
    %     title([filename, ' mode FR = ', num2str(fr)], 'interpreter', 'none')
    %     figure(p)
    %     hold off
    %     keyboard
    %
    %     indicies = get(gca, 'XLim');
    %     indicies (1) = floor(max(1, indicies(1)));
    %     indicies (2) = ceil(indicies(2));
    %     IntraData.FRselect(1) = max(1, find(IntraData.timepoints>indicies(1), 1, 'first') - 1);
    %     d = find(IntraData.timepoints>indicies(2), 1, 'first');
    %     if isempty(d)
    %         IntraData.FRselect(2) = length(IntraData.timepoints);
    %     else
    %         IntraData.FRselect(2) = d-1;
    %     end
    %
    %     clear spikeheight spikeindex
    %     [spikeheight, spikeindex] = findpeaks(IntraData.trace(IntraData.FRselect(1):IntraData.FRselect(2)), ...
    %         'MinPeakHeight',threshold, ...
    %         'MinPeakDistance', 150);
    %     spikeTimestamp = spikeindex*20+IntraData.timepoints(IntraData.FRselect(1));
    %     fr = length(spikeindex)/((IntraData.FRselect(2) - IntraData.FRselect(1))/50000);
    %
    %     subplot(2, 1, 2)
    %     plot(IntraData.timepoints(IntraData.FRselect(1):IntraData.FRselect(2)), ...
    %         IntraData.trace(IntraData.FRselect(1):IntraData.FRselect(2)))
    %     hold on
    %     plot(spikeTimestampFull(2:end), instafr)
    %     scatter(IntraData.markertime,repmat(100, 1, length(IntraData.markertime)))
    %     scatter(spikeTimestamp, spikeheight)
    %     xlim([IntraData.timepoints(IntraData.FRselect(1)) IntraData.timepoints(IntraData.FRselect(2))])
    %     ylim([-1000 200])
    %     title([filename, ' mode FR = ', num2str(fr)], 'interpreter', 'none')
    %     figure(p)
    %     hold off
    %     keyboard
    %
    %     IntraData.FR = fr;
    %     IntraData.instaFR = instafr;
    %     save(filename, 'IntraData')
    %     clear IntraData fr spike*
end