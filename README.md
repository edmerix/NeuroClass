# NeuroClass

A couple of Matlab classes to make handling spike sorting data substantially easier.

## Main Files

__SingleUnit.m__: class to store individual single units, with useful methods to inspect waveforms, autocorrelations etc.

__MultipleUnits.m__: class to store populations of single units, with methods to assess population firing activities.

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
| xcorr        | Calculate the cross-correlation between this and another SingleUnit object    |
| gaussian_fr  | Estimate the instantaneous firing rate of this unit with a Gaussian kernel    |
| hist_fr      | Calculate the binned firing rate of this unit                                 |
| mean_ac_lag  | Calculate the mean lag of this unit's autocorrelation (for subclassification) |
| ISI          | Calculate the ISI for this unit (autocorrelation is better than raw ISI...)   |
| plot_*       | Plot the output of various other calculations                                 |
| retrieve_*   | Return the values of the various other calculations                           |

### MultipleUnit object:

Properties:

| Property |                             Description                            |
|---------:|--------------------------------------------------------------------|
|  patient | Patient identifier (string)                                        |
|  seizure | Seizure number (int)                                               |
|    epoch | Start and finish times of epoch [double double]                    |
|    units | Array of SingleUnit objects                                        |
|      snr | Stores the signal-to-noise ratio of this recording once calculated |
|     info | Any extra details and info for this recording                      |

Selected methods _(see [example_usage.m](example_usage.m) for input/output explanations)_:

|           Method |                                                       Description                                                      |
|-----------------:|------------------------------------------------------------------------------------------------------------------------|
|         add_unit | Add a SingleUnit object to this MultipleUnits object                                                                   |
|  all_spike_times | Return all spike times across all units within specified epoch (defaults to all)                                       |
|     beefy_raster | Make a comprehensive raster plot of these units, allowing color-coding. Not recommended: User .raster() for speed.     |
|    channel_units | Return array of all SingleUnit objects from specified channel.                                                         |
|      gaussian_fr | Calculate the Gaussian estimate of the population firing rate across all units within this MultipleUnits               |
| order_by_channel | Order the SingleUnits by channel number                                                                                |
|    order_by_rate | Order the SingleUnits by overall firing rate                                                                           |
|plot_channel_units| Plot all units from specified channel to assess separation accuracy ([see screenshot below](#screenshots))             |
|           raster | Make a basic raster plot of these units. No color-coding, not very pretty - use .beefy_raster() for comprehensive plot |
|         unit_snr | Calculate the SNR across all units in the object                                                                       |
|     top_channels | Return the specified number of channels with the most units recorded                                                   |


## Screenshots

Example output of the MultipleUnits plot_channel_units(channel) method:

![Screenshot of plot_channel_units output](Screenshots/ExampleChannelUnits.png?raw=true "plot_channel_units output figure")
