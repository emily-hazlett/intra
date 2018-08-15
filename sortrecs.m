prompt = {'Folder containing files','Files to use'};
dlg_title = 'Input';
num_lines = [1, 100];
defAns = {'C:\Users\emily\Desktop\importer\', 'testerview.txt'};
options.Resize='on';
answer = inputdlg(prompt, dlg_title, num_lines, defAns, options);

folderOld = cd;
folderNew = answer{1};
cd(folderNew);
Files = dir(answer{2});

button = questdlg('Start at begining of folder?', ...
    'Start place', ...
    'Yes', 'No', 'Yes');
switch button
    case 'Yes'
        ii = 1;
    case 'No'
        prompt = {'Which file number to start at?'};
        dlg_title = 'Input';
        num_lines = [1, 100];
        defAns = {['Must be <=', num2str(length(Files))]};
        options.Resize='on';
        answer = inputdlg(prompt, dlg_title, num_lines, defAns, options);
        if isempty(answer)
            warndlg('Must enter value')
            return
        elseif strcmp(answer, ['Must be <=', num2str(length(Files))]) == 1
            warndlg('Must enter value')
            return
        end
        ii = answer{1};
end

decent = cell(1,1);
decentN = 0;
shit = cell(1,1);
shitN = 0;
eh = cell(1,1);
ehN = 0;

p = figure('units','normalized','outerposition',[0 0 1 1]);

for i = ii:length(Files)
    filename = Files(i).name;
    sweeps = importdata(filename);
    rmp = mean(sweeps,2);
    
    plot(rmp)
    figure(p)
    
    button = questdlg('File remotely decent?', ...
        'Sorter', ...
        'Yes', 'Hell no', 'Eh?', 'Eh?');
    switch button
        case 'Yes'
            decentN= decentN + 1;
            decent{decentN,1} = [folderNew, filename];
        case 'Hell no'
            shitN = shitN + 1;
            shit{shitN,1} = [folderNew, filename];
        case 'Eh?'
            ehN = ehN + 1;
            eh{ehN,1} = [folderNew, filename];
    end
    
    if isempty(button) == 1
        warndlg(['File location stopped at is ', num2str(i)])
        return
    end
end

