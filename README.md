# NeuroClass

A couple of Matlab classes to make handling spike sorting data substantially easier.

## Main Files

__SingleUnit.m__: class to store individual single units, with useful methods to inspect waveforms, autocorrelations etc.

__MultipleUnits.m__: class to store populations of single units, with methods to assess population firing activities.

__Dependencies__ (minimal and easily removed): [get_zvalues.m](dependencies/get_zvalues.m) from [UltraMegaSort](https://neurophysics.ucsd.edu/software.php), and [estimateColors.m](dependencies/estimateColor.m), a method of estimating the human name for colors, to simplify interpretation of complex raster plots.


## Usage

A more thorough example of usage, starting from the raw output of UltraMegaSort is in [example_usage.m](example_usage.m).

Quickstart:

```Matlab
data = MultipleUnits('patient', 'JoeBloggs', 'seizure', 1);

unit = SingleUnit('times',spiketimes,'waveforms',waves,'channel',1);
data.addUnit(unit);
% add all desired units (usually in a loop over output files)

data.order_by_rate();
data.raster();
```

## Structure overview 

### SingleUnit object:

Selected properties:

| Property  |                                   Description                                  |
|----------:|--------------------------------------------------------------------------------|
| UID       | Unique ID for this neuron                                                      |
| waveforms | [m by n] matrix of m spikes across n datapoints                                |
| times     | Column vector of spike times in your preferred units (e.g. s, ms, datapoints)  |
| wideband  | Mean waveform from the broadband signal (used for cell-type subclassification) |
| type      | Type of neuron (e.g. FS-IN, RS, PC)                                            |
| extra     | Structure for storage of extra data not core to the object                     |

Selected methods:

| Method       |                                  Description                                  |
|-------------:|-------------------------------------------------------------------------------|
| inspect_unit | Plots an overview of this unit for visual inspection                          |
| autocorr     | Calculate the autocorrelation for this unit                                   |
| gaussian_fr  | Estimate the instantaneous firing rate of this unit with a Gaussian kernel    |
| hist_fr      | Calculate the binned firing rate of this unit                                 |
| mean_ac_lag  | Calculate the mean lag of this unit's autocorrelation (for subclassification) |
| plot_*       | Plot the output of various other calculations                                 |
| retrieve_*   | Return the values of the various other calculations                           |

### MultipleUnit object:

Need to export this. Prompt me if you want to know and it's still not here.
