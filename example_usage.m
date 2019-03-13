% Purely for identification purposes within the MultipleUnits object, and
% don't need to be set:
patient = 'patient_name';
seizure = 1;

% Directory of files with spike sorting results from UltraMegaSort2000
filepath = '~/Data/patients/patient_name';

% Load up a list of files in that directory, e.g. where files are called DataFromChannel{1-n}.mat
fls = dir([filepath filesep 'filename_channel*.mat']);

% Name of the output struct in those files (set this to whatever variable
% name you saved the UltraMegaSort2000 output as):
spikestruct = 'spikes';

% Create the MultipleUnits object to hold them (leave out unnecessary 
% parameters for your circumstances):
data = MultipleUnits('patient', patient, 'seizure', seizure);

% Loop through the files of spike sorting results:
for f = 1:length(fls)
    % Load and extract the data into spikes variable:
    temp = load([filepath filesep fls(f).name],spikestruct);
    spikes = temp.(spikestruct);
    
    % Find which clusters have been marked as good in SplitMerge or
    % splitmerge_tool (or manually):
    good = spikes.labels(spikes.labels(:,2) == 2,1);
    
    % Extract the channel number from the file name (or however else you'd
    % like to do so, can't say I recommend this way):
    channel = str2double(strrep(strrep(fls(f).name,'filename_channel',''),'.mat',''));
    
    % For each cluster marked as "good" make a SingleUnit object, and add
    % it to the MultipleUnits object:
    for g = 1:length(good)
        unit = SingleUnit(...
            'times',        spikes.spiketimes(spikes.assigns == good(g))', ...
            'waveforms',    spikes.waveforms(spikes.assigns == good(g),:), ...
            'channel',      channel ...
            );
        % Add the SingleUnit to the MultipleUnits:
        data.add_unit(unit); % add_unit auto adds a unique ID for the unit if none is there
    end
end

% We now have the full population of single units stored as SingleUnit
% objects within a parent MultipleUnits object, & can call methods such as:

% Order the single units by firing rate:
data.order_by_rate();

% Show a speedy (but ugly) raster plot in current order:
data.raster();

% Show a slow (but pretty) raster plot, highlighting unit 23, onto an axes 
% object with the handle "bigAx":
hFig = figure('Units', 'normalized', 'Position', [0.02 0.02 0.9 0.9]);  % make a big figure
bigAx = axes('Parent', hFig, 'Position', [0.05 0.05 0.9 0.9]);          % put a big axes object on it
data.beefyraster('highlight', 23, 'ax', bigAx);                         % plot the raster to that axes, highlight unit 23

% Get the 5 channels with the most units on them:
top_ch = data.top_channels(5);

% Get just the units from channel 9:
units = data.channel_units(9);

% Visually inspect the 2nd unit from this channel (show waveforms through 
% time, voltage at detection, autocorrelation, Mahalanobis distance and 
% expected chi-squared distribution, and cluster in principal component
% space):
units(2).inspect_unit();
% We don't need to extract units to do this. For example, inspect the 79th
% unit (i.e. 79th most active unit if data.order_by_rate()):
data.units(79).inspect_unit();

% Calculate and retrieve a Gaussian kernel estimate of the firing rate of 
% the first unit on this channel:
units(1).gaussian_fr(5);    % where 5 is the width of the Gaussian in the same units as units(1).times;
[gfr, tt] = units(1).retrieve_gauss_fr();

% Plot that estimate of the firing rate:
figure, plot(tt, gfr);

% Show all available methods of the MultipleUnits collection:
methods('MultipleUnits')

% And show all available methods of the individual SingleUnit units:
methods('SingleUnit')

%-------------------------------------------------------------------------%
% Other important properties for SingleUnit:
%   wideband        the mean wideband waveform
%   celltype        what cell-type it's classified as (use UnitSubclassify)
%-------------------------------------------------------------------------%
% Store any other data you want in the "extra" field in each SingleUnit.
%-------------------------------------------------------------------------%