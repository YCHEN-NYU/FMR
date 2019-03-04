%% ========================================================================
% Fit FMR spectral data measured from microwave experiments
% Author: Yizhang Chen, NYU
% Contact: yc1224@nyu.edu
% Fit FMR spectral data with single Lorentzian model: [S12_real, S12_imag]
% simultaneously

%% ========================================================================
clear;
close all;

% Indices of field(T), S12_real, S12_imag
indexH = 2;
indexS12real = 7;
indexS12imag = 8;

% field in Oe, unit conversion
Oe_to_Tesla = 10000;

%%% read files in the current folder
folder = pwd;

% read all data file ended with *MHz.dat in the current folder
input_folder = '../../data/';
% input_folder = pwd;
cd '../data';
files = dir('*MHz.dat');

% sort out the files in natural order
[fileList, index] = sort_natural_order({files.name});

% go to the destination data folder, define the output file
% outputloc=[folder '/' outputname];
% Open or create new file for writing. Append data to the end of the file.
% output_file = fopen(outputloc,'a+');


% starting and ending indices of data files
startIndex = 1;
endIndex = numel(fileList); % number of files inside the current folder
fileIndex = startIndex;

while fileIndex <= endIndex

% Get frequency(GHz) from file names (read out number between _ and MHz)
    tmp = strsplit(char(fileList{fileIndex}),'K_');
    tmp = tmp(end);% Take '*MHz.dat' from data files 'XXX_*MHz.dat'
    tmp = strsplit(char(tmp),'MHz'); % Get '*' by splitting '*MHz.dat'
    frequency = str2double(char(tmp(1)))/1000;% convert MHz to GHz

    tmp = strsplit(char(fileList{fileIndex}),'K_');
    tmp = strsplit(tmp{1, 1}, 'CoFeB');
    temperature = tmp{1, 2};
 
    % load data
    fileloc = [input_folder char(fileList(fileIndex))];
    fprintf('%s\n',char(fileList(fileIndex)));
    rawdata = importdata(fileloc);
    data = rawdata(2:end-1,:);% get rid of the 1st data point.

    % Construct complex data y with S12_real and S12_imag
    x = data(:,indexH)/Oe_to_Tesla; %in Oe
    y = complex(data(:,indexS12real),data(:,indexS12imag)); 

    % fit [S12_real, S12_imag] - H, extract the Hres, FWHM information
    [params, confidence_intervals, fig] = fit(x,y,frequency);

    % Ask user to move on or fit current data again
    message = input('Fit next one? (y/n)', 's');

    % if input character is 'y', save fitting parameters & figure
    if strcmpi(message, 'y')
        % Save [Hres, FWHM, FWHM_LyowerBound, FWHM_UpperBound] to outputfile
        Hres = abs(params(4)*Oe_to_Tesla); % resonant (Oe)
        FWHM= abs(params(2))*Oe_to_Tesla;  % Full-width half maximum (Oe)
        FWHM_LowerBound= min(abs(confidence_intervals(2,:)))*Oe_to_Tesla; % lower bound of 95% confident interval (Oe)
        FWHM_UpperBound= max(abs(confidence_intervals(2,:)))*Oe_to_Tesla; % upper bound of 95% confident interval (Oe)
        output_filename = [temperature  'K.txt'];
        output_file = fopen(output_filename,'a+'); % specify the current temperature
        fprintf(output_file,'%10.0f, %10.2f, %10.3f, %10.3f, %10.3f\r\n', frequency, Hres, FWHM, FWHM_LowerBound, FWHM_UpperBound);

        % Save figure with frequency, temperature information
        fig.PaperPositionMode = 'auto';% set image size as auto, this makes sure the figure doesn't get disorted
        saveas(fig,strtok(string(frequency))+'GHz_'+temperature + 'K','png');

        fileIndex = fileIndex + 1; % move on to the next data file
    end
close(fig);
end

% Print out ending information & close output file
fprintf('Your fitting is done!\n')
fclose(output_file);
