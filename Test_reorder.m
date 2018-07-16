% Test Reorder
% This function reorders the images in directory


tic;
cd(directory);

% Copy the files over with a new name.
inputFiles = dir(fullfile('*.tif'));
fileNames = { inputFiles.name };

file = fileNames(1);
ind_z = regexp(file,'z\d\d');
ind_m = regexp(file,'m\d\d\d');

file = file{1};
ind_z = ind_z{1};
ind_m = ind_m{1};

toc;
time = 0;
for k = 1:length(fileNames)
    if mod(k,100) == 1
        tic;
    end
	filename = fileNames{k};
    %newFileName = [filename(1:27),filename(32:size(filename, 2))];
    %    filename(ind_z:(ind_m-1)),filename((ind_m+4):size(file, 2))];
    newFileName = [filename(1:(ind_z-1)),filename(ind_m:(ind_m+3)), ...
        filename(ind_z:(ind_m-1)),filename((ind_m+4):size(filename, 2))];
    movefile(filename, newFileName);
    if mod(k,100) == 0
        newtime = toc;
        time = time + newtime;
        disp([num2str(k) ' iterations done in ' num2str(time) ' seconds.']);
        
    end
end

