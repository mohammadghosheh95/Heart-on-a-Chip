%% Setup selecting UI

% Select a .csv file to load
[setupFile, setupPath] = uigetfile(".csv",'MultiSelect', 'off');

% If canceled, return
if (setupFile == 0)
    return;
end

% Find a complete file path
setupFilePath = fullfile(setupPath, setupFile);

% Read the setup table from the file
setupTable = readtable(setupFilePath,'ReadRowNames',true,'ReadVariableNames',true);

% Convert to array
setupArray = table2array(setupTable);

% Set the rows and cols variables
numCols = size (setupTable.Properties.VariableNames, 2);
numRows = size (setupTable.Properties.RowNames, 1);
totalNum = numCols * numRows;

% Flip every other row (for the snake pattern)
%setupArray(2:2:end, :) = fliplr(setupArray(2:2:end, :));
% Convert to linear
setupArray = setupArray';
setupArray = setupArray(:);

% Select only the name of the drug from the setup array
drugNames = cellfun (@(x) extractBefore (x, strfind (x, ' ')), setupArray, 'UniformOutput', false);
% % % Reshape into 1x(rows * cols)
% % drugNames = reshape (drugNames', 1, numRows * numCols);
% Find groups
[drugG, drugID] = findgroups (drugNames);
% % Reshape back into the matrix
% drugG = reshape (drugG, numCols, numRows)';

% Select only the concentration of the drug from the setup array
drugConcentrations = cell2mat(cellfun (@(x) str2double(extractBetween (x, strfind (x, ' ')+1, strfind (x, 'M')-2)), setupArray, 'UniformOutput', false));

%% Folder selecting UI

% Select a folder to load
path = uigetdir;

% If canceled, return
if (path == 0)
    return;
end

%% Loading the logfiles in the directory

% Find all files in the path using dir()
files = dir (path);

% Convert the files struct to a table for sorting
filesTable = struct2table(files);

% Remove '.', '..', and 'DS_STORE' from the list
filesTable([1,2,3],:) = [];
% TO FIX: DS_STORE

% Sort the table by date of the files
filesTable = sortrows(filesTable, 'date');

% Define an array with the index of each measurment (to be used for
% navigating the following matrices)
ind = ones (totalNum, 1);

% Define an empty array with the amplitude fft frequencies (each row represents
% a different organoid, each column represents a different measurment)
ampFFTFreq = NaN (totalNum, 875);

% Define an empty array with the oxygen fft frequencies (each row represents
% a different organoid, each column represents a different measurment)
o2FFTFreq = NaN (totalNum, 875);

% Define an empty array with the mean signal/noise (each row represents
% a different organoid, each column represents a different measurment)
meanSns = NaN (totalNum, 875);

% Define an empty array with the calculated amplitude frequencies from the fft of the
% calculations. Each row represents a different organoid, each column a
% different measurment, and all of the 3rd dimension represents the
% different frequencies for the amplitude FFT.
ampFs = NaN (totalNum, 875, 45);

% Define an empty array with the calculated amplitude P1s from the fft of the
% calculations. Each row represents a different organoid, each column a
% different measurment, and all of the 3rd dimension represents the
% different P1s for the amplitude FFT.
ampP1s = NaN (totalNum, 875, 45);

% Define an empty array with the calculated o2 frequencies from the fft of the
% calculations. Each row represents a different organoid, each column a
% different measurment, and all of the 3rd dimension represents the
% different frequencies for the O2 FFT.
o2Fs = NaN (totalNum, 875, 45);

% Define an empty array with the calculated O2 P1s from the fft of the
% calculations. Each row represents a different organoid, each column a
% different measurment, and all of the 3rd dimension represents the
% different P1s for the O2 FFT.
o2P1s = NaN (totalNum, 875, 45);
% Define an empty array with the calculated O2 P1s from the fft of the
% calculations. Each row represents a different organoid, each column a
% different measurment, and all of the 3rd dimension represents the
% different P1s for the O2 FFT.
ampP10s = NaN (totalNum, 875);

% Define an empty array with the raw amplitude measurments from the
% logfiles. Each row represents a different organoid, each column a
% different measurment, and all of the 3rd dimension represents the
% different amplitude measurments.
amps = NaN (totalNum, 875, 120);

% Define an empty array with the raw O2 measurments from the
% logfiles. Each row represents a different organoid, each column a
% different logfile, and all of the 3rd dimension represents the
% different O2 measurments.
o2s = NaN(totalNum, 875, 120);

% Define an empty array with the raw time measurments from the
% logfiles. Each row represents a different organoid, each column a
% different logfile, and all of the 3rd dimension represents the
% different time measurments.
times = NaN(totalNum, 875, 120);

% Define an empty array with the start times of each measurment. Each row
% represents a different organoid, each column a different logfile, and the
% third dimension is reserved for the date format (1x6)
startTimes = NaN(totalNum, 875, 6);

% Iterate through each logfile in the available files
for i = 1:size(filesTable, 1)
    
    % Find the fullfile path to the logfile
    filePath = fullfile(filesTable.folder(i), filesTable.name(i));
    filePath = filePath{:};
   
    % Read the contents of the file
    [startTime, pos, time, a0, a1, sn, o2] = readFile(filePath);
    
    % If the file has an invalid position, skip
    if (pos < 0)
        continue;
    end
    
    % Save the raw amplitude, oxygen and time values
    amps(pos, ind(pos), 1:size(a0, 1)) = a0;
    o2s(pos, ind(pos), 1:size(o2, 1)) = o2;
    times(pos, ind(pos), 1:size(time, 1)) = time;
 
    % Calculate the FFT for the amplitude
    [ampFFTFreq(pos, ind(pos)), f, p1, p10] = calculateFFT(time, a0, 0.1);
    % Save the calculated FFT Fs and P1s
    ampFs(pos, ind(pos), 1:size(f, 2)) = f;
    ampP1s(pos, ind(pos), 1:size(p1, 1)) = p1;
    ampP10s(pos, ind(pos)) = p10;
    
    % Calculate the FFT for the oxygen
    [o2FFTFreq(pos, ind(pos)), f, p1] = calculateFFT(time, o2, 0.1); 
    % Save the calculated FFT Fs and P1s
    o2Fs(pos, ind(pos), 1:size(f, 2)) = f;
    o2P1s(pos, ind(pos), 1:size(p1, 1)) = p1;
    
    % Calculate and store the mean signal/noise
    meanSns(pos, ind(pos)) = mean(sn);
    
    % Store the start time
    startTimes(pos, ind(pos), :) = startTime;
    
    % Add one to the index
    ind(pos) = ind(pos) + 1;

    % Log the progress
    if (mod(i, 100) == 0)
        i
    end
    
end


%% Loading the file's data

% Reads the log file using textscan
% Arguments:
% - filename: the complete filepath for the logfile (obtained using
% fullfile(path, file)
% Returns:
% - startTime: the start time of the experiment as a MATLAB Date (using datevec)
% - pos: the linear numerical position of this well according to the logfile
% - time: the array of the time points in seconds for the measurments
% - a0: the array of amplitude 0 values for the measurments
% - a1: the array of amplitude 0 values for the measurments
% - sn: the signal-to-noise values for the measurments
function [startTime, pos, time, a0, a1, sn, o2] = readFile (filename)

    % Open the file at the filePath
    fileID = fopen(filename);
    
    % Read the start time
    startTime = textscan(fileID, '%q%q', 1, 'HeaderLines', 1);
    % Convert start time to string
    startTime = strcat (string (startTime (1)), " - ", string (startTime (2)));
    % Convert start time to date
    startTime = datevec(startTime,'dd.mm.yyyy - HH:MM:SS');

    % Read the position
    pos = textscan(fileID, '%*q%d%*q', 1, 'HeaderLines', 6);
    pos = pos{1};
    
    % Define the format for the rest of the file
    formatSpec = '%f%f%*f%*f%*f%*f%*f%f%f%*f%*f%f%*f%*q%*f%*f%*q%*d%*q%*d%*d';

    % Continue reading through the text file with this format, automatically
    % stopping when reaching the end
    fileData = textscan(fileID, formatSpec);%, 'HeaderLines', 5);

    % Close the file
    fclose(fileID);

    % Extract the arrays from the cell-array
    time = fileData{1};
    o2 = fileData{2};
    a0 = fileData{3};
    a1 = fileData{4};
    sn = fileData{5};

end

% Calculates the fft
% Arguments:
% - time: the time array for the measurments
% - amp: the amplitude values for the measurments
% - measurmentTime: the length of each measurment in seconds
function [freq, f, P1, p10] = calculateFFT(~, amp, measurmentTime)

    % Calculate the number of measurments
    L = size(amp, 1);  

    %figure
    %plot(time,amp);
    
    % Calculate the fft of the amplitude
    Y = fft(amp);

    % Extract the fft results
    P2 = abs(Y/L);
    P1 = P2(1:floor(L/2+1));
    P1(2:end-1) = 2*P1(2:end-1);
    
    % Calculate an array of the possible frequencies
    f = (1/measurmentTime)*(0:(L/2))/L;
    
    % Remove f=0 from f and P1 (noise)
    p10 = P1(f==0);
    P1 = P1 (f > 0.5 & f < 2);
    f = f (f > 0.5 & f < 2);
    
    %figure
    %plot(f,P1) 
    
    % Find the maximum frequency 
    [val, i] = max(P1);
    freq = f(i);
    
    if (isempty(freq))
        freq = NaN;
        return;
    end

    halfPeak = P1 (f == freq/2);
    
    if (halfPeak/val > 0.7)

        freq = freq/2;

    end

end
