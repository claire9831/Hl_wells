% Instead of renaming files, returns an iteration vector of convenient
% order

function itvec = reorderVecNew(directory,numDepths)
tic;
cd(directory);

% Copy the files over with a new name.
inputFiles = dir(fullfile('*.tif'));
fileNames = { inputFiles.name };

file = fileNames(1);
ind_z = regexp(file,'z\d\d');
ind_m = regexp(file,'m\d\d\d');

ind_z = ind_z{1};
ind_m = ind_m{1};

itvec = zeros([length(fileNames), 1],'uint16');
time = 0;
for k = 1:length(fileNames)
    if mod(k,100) == 1
        tic;
    end
    filename = fileNames{k};
    mNum = filename(ind_m+1:ind_m+3);
    mNum = uint16(str2double(mNum));
    zNum = filename(ind_z+1:ind_z+2);
    zNum = uint16(str2double(zNum));
    itvec((mNum - 1) * numDepths + zNum) = k;
    (mNum - 1) * numDepths + zNum
    if mod(k,100) == 0
        newtime = toc;
        time = time + newtime;
        disp([num2str(k) ' iterations done in ' num2str(time) ' seconds.']);
    end
end
end