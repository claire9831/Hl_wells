% Created by Rui Yan, April 20th 2018
% This function reorders the images in directory

function reorder_image(directory)

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

for k = 1:length(fileNames)
	filename = fileNames{k};
    newFileName = [filename(1:(ind_z-1)),filename(ind_m:(ind_m+3)), ...
        filename(ind_z:(ind_m-1)),filename((ind_m+4):size(filename, 2))];
    movefile(filename, newFileName);
end

end