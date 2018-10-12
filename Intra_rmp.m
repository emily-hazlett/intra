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
folderNew = 'C:\Users\emily\Desktop\bigone\maybe\';
cd(folderNew);
Files = dir('*.mat');
ii = 1;

p = figure('units','normalized','outerposition',[0 0 1 1]);
for i = ii:length(Files)
    filename = Files(i).name;
    load(filename);
    
    subplot(2, 1, 1)
    h = histogram(IntraData.trace(IntraData.traceselect(1):IntraData.traceselect(2)));
    hold on
    [~, m] = max(movmean(h.Values,20));
    valuers = h.BinLimits(1):h.BinWidth:h.BinLimits(2)-h.BinWidth;
    rmp(1) = round((valuers(m) - h.BinWidth/2)/10, 1);
    line('XData', [rmp*10 rmp*10], 'YData', [0 max(get(gca, 'ylim'))], ...
        'LineWidth', 3, 'Color', 'r')
    title([filename, ' rmp = ', num2str(rmp)], 'interpreter', 'none')
    hold off
    clear h valuers m
    
    subplot(2, 1, 2)
    plot(IntraData.timepoints(IntraData.traceselect(1):IntraData.traceselect(2)), ...
        IntraData.trace(IntraData.traceselect(1):IntraData.traceselect(2)))
    hold on
    line('XData', [IntraData.timepoints(IntraData.traceselect(1)), IntraData.timepoints(IntraData.traceselect(2))], ...
        'YData', [rmp*10 rmp*10], 'LineWidth', 2, 'Color', 'r')
    if isfield(IntraData, 'markertime')
        selecter = IntraData.markertime(IntraData.markertime > IntraData.timepoints(IntraData.traceselect(1)));
        selecter = selecter(selecter < IntraData.timepoints(IntraData.traceselect(2)));
        if isempty(selecter) == 0
            IntraData.markerselect (1) = find(IntraData.markertime==selecter(1));
            IntraData.markerselect (2) = find(IntraData.markertime==selecter(end));
            
            scatter(IntraData.markertime(IntraData.markerselect(1):IntraData.markerselect(2)), ...
                repmat(100, 1, length(IntraData.markertime(IntraData.markerselect(1):IntraData.markerselect(2)))))
            plot(IntraData.markertime(IntraData.markerselect(1):IntraData.markerselect(2)), ...
                IntraData.pulsedir(IntraData.markerselect(1):IntraData.markerselect(2)) .* IntraData.pulseamp(IntraData.markerselect(1):IntraData.markerselect(2)) -800)
        else
            IntraData.markerselect = [];
        end
        clear selecter
    end
    
    axis tight
    ylim([-1000 200])
    figure(p)
    hold off
    keyboard
    
    IntraData.rmp = rmp;
    save(filename, 'IntraData')
    clear IntraData
end