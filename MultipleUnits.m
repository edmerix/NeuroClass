classdef MultipleUnits < handle
    properties
        units       SingleUnit
        patient     string
        seizure     single
        epoch       double
        snr         double
        info        string
    end
    
    properties (SetAccess = private, Hidden = true)
        current_order = 'none';
    end
    
    methods
        % constructor
        function obj = MultipleUnits(varargin)
            allowable = fieldnames(obj);
            if mod(length(varargin),2) ~= 0
                error('Inputs must be in name, value pairs');
            end
            for v = 1:2:length(varargin)
                if find(ismember(allowable,varargin{v}))
                    obj.(varargin{v}) = varargin{v+1};
                else
                    disp([9 'Not assigning ''' varargin{v} ''': not a property of MultipleUnits class']);
                end
            end
            if isempty(obj.epoch)
                obj.epoch = [-Inf Inf];
            end
        end
        % add a SingleUnit object to this collection
        function obj = add_unit(obj,unit)
            if isempty(unit.UID)
                unit.UID = length(obj.units) + 1;
            end
            obj.units = [obj.units unit];
            if ~isempty(unit.times) && (min(unit.times) < obj.epoch(1) || max(unit.times) > obj.epoch(2))
                disp([9 'Heads up: this unit has spiketimes outside the currently set epoch for this collection']);
            end
            % need to fix the ordering if it's not currently set to none
            if ~strcmp(obj.current_order,'none')
                obj.(['order_by_' obj.current_order]);
            end
        end
        % re-order units by rate
        function order_by_rate(obj)
            tot = zeros(1,length(obj.units));
            for n = 1:length(obj.units)
                tot(n) = length(obj.units(n).times);
            end
            [~,ord] = sort(tot,'ascend');
            obj.units = obj.units(ord);
            obj.current_order = 'rate';
        end
        % re-order units by channel
        function order_by_channel(obj)
            chan = zeros(1,length(obj.units));
            for n = 1:length(obj.units)
                chan(n) = obj.units(n).channel;
            end
            [~,ord] = sort(chan,'ascend');
            obj.units = obj.units(ord);
            obj.current_order = 'channel';
        end
        % simple raster plot (speedy)
        function raster(obj,varargin)
            settings.axes = gca;
            settings.showtypes = false;
            settings.in_color = [0.64 0.08 0.18];
            for v = 1:2:length(varargin)
                settings.(varargin{v}) = varargin{v+1};
            end
            if settings.showtypes
                if length(obj.units(end).times) > 1 && isrow(obj.units(end).times)
                    nIN = length(cell2mat([obj.units(strcmpi([obj.units.type],'in')).times]));
                    nPC = length(cell2mat([obj.units(strcmpi([obj.units.type],'pc')).times]));
                else
                    nIN = length(cell2mat({obj.units(strcmpi([obj.units.type],'in')).times}'));
                    nPC = length(cell2mat({obj.units(strcmpi([obj.units.type],'pc')).times}'));
                end                
                
                nTotalSpikesIN = 3 * nIN;
                nTotalSpikesPC = 3 * nPC;
                xPointsIN = NaN(nTotalSpikesIN*3,1);
                yPointsIN = xPointsIN;
                xPointsPC = NaN(nTotalSpikesPC*3,1);
                yPointsPC = xPointsPC;
                currentIndIN = 1;
                currentIndPC = 1;
            else
                % Based on plotSpikeRaster by Jeffrey Chiou:
                nTotalSpikes = 3 * length(obj.all_spike_times);

                % We can make it a continuous plot by separating segments with NaNs
                xPoints = NaN(nTotalSpikes*3,1);
                yPoints = xPoints;
                currentInd = 1;
            end
            
            for u = 1:length(obj.units)
                nSpikes = length(obj.units(u).times);
                nanSeparator = NaN(1,nSpikes);

                trialXPoints = [obj.units(u).times'; obj.units(u).times'; nanSeparator];
                trialXPoints = trialXPoints(:);

                trialYPoints = [(u-0.5)*ones(1,nSpikes); (u+0.5)*ones(1,nSpikes); nanSeparator];
                trialYPoints = trialYPoints(:);

                % Save points and update current index
                if settings.showtypes
                    if strcmpi(obj.units(u).type,'in')
                        xPointsIN(currentIndIN:currentIndIN+nSpikes*3-1) = trialXPoints;
                        yPointsIN(currentIndIN:currentIndIN+nSpikes*3-1) = trialYPoints;
                        currentIndIN = currentIndIN + nSpikes*3;
                    else
                        xPointsPC(currentIndPC:currentIndPC+nSpikes*3-1) = trialXPoints;
                        yPointsPC(currentIndPC:currentIndPC+nSpikes*3-1) = trialYPoints;
                        currentIndPC = currentIndPC + nSpikes*3;
                    end
                else
                    xPoints(currentInd:currentInd+nSpikes*3-1) = trialXPoints;
                    yPoints(currentInd:currentInd+nSpikes*3-1) = trialYPoints;
                    currentInd = currentInd + nSpikes*3;
                end
            end
            if settings.showtypes
                plot(settings.axes,xPointsPC,yPointsPC,'k');
                hold(settings.axes,'on')
                plot(settings.axes,xPointsIN,yPointsIN,'color',settings.in_color);
            else
                plot(settings.axes,xPoints,yPoints,'k');
            end
            xlim(settings.axes,obj.epoch);
            ylim(settings.axes,[0.5 length(obj.units)+0.5])
            xlabel(settings.axes,'Time (s)')
            ylabel(settings.axes,['Neuron ranking (by ' obj.current_order ')']);
            set(settings.axes,'tickdir','out');
        end
        % full raster plot (slow, but color-able beyond just by cell type)
        function hdls = beefy_raster(obj,varargin)
            settings.base_color = 'k';
            if isequal(get(0,'DefaultAxesColor'), [0 0 0])
                settings.base_color = 'w';
            end
            settings.blackout = 0;
            settings.highlight = [];
            settings.axes = gca;
            settings.in_color = [0.64 0.08 0.18];
            settings.linewidth = 1;
            for v = 1:2:length(varargin)
                settings.(varargin{v}) = varargin{v+1};
            end
            if nargout > 0
                hdls = cell(1,length(obj.units));
            end
            for n = 1:length(obj.units)
                if length(obj.units(n).times) > 2
                    if any(strcmp(properties(obj.units(n)), 'type')) && ~settings.blackout && isempty(settings.highlight)
                        if strcmp(obj.units(n).type,'pc')
                            col = settings.base_color;
                        else
                            col = settings.in_color;
                        end
                    else
                        if isempty(settings.highlight)
                            col = settings.base_color;
                        else
                            if exist('distinguishable_colors','file')
                                cols = distinguishable_colors(length(settings.highlight));
                            else
                                cols = lines(length(settings.highlight));
                            end
                            highlighting = find(settings.highlight == obj.units(n).channel);
                            if highlighting > 0
                                col = cols(highlighting,:);
                            else
                                col = settings.base_color;
                            end
                        end
                    end
                    if nargout > 0
                        hdls{n} = line(settings.axes,[obj.units(n).times obj.units(n).times],[n-1 n],'color',col,'linewidth',settings.linewidth);
                    else
                        line(settings.axes,[obj.units(n).times obj.units(n).times],[n-1 n],'color',col,'linewidth',settings.linewidth);
                    end
                end
            end
            if min(obj.epoch) < 0 && max(obj.epoch) > 0
                line(settings.axes,[0 0],[0 n],'color','r')
            end
            
            xlim(settings.axes,obj.epoch);
            ylim(settings.axes,[0 n])
            xlabel(settings.axes,'Time (s)')
            ylabel(settings.axes,['Neuron ranking (by ' obj.current_order ')']);
            set(settings.axes,'tickdir','out');
        end
        % show top n channels with most units on them
        function n_top = top_channels(obj,n,dropNaN)
            if nargin < 2 || isempty(n)
                n = 3;
            end
            if nargin < 3 || isempty(dropNaN)
                dropNaN = false;
            end
            chans = zeros(1,length(obj.units));
            for u = 1:length(obj.units)
                chans(u) = obj.units(u).channel;
            end
            n_top = zeros(1,n);
            for nn = 1:n
                top = mode(chans);
                if nargout < 1
                    disp([9 'Channel ' num2str(top) ' has ' num2str(length(find(chans == top))) ' units']);
                end
                n_top(nn) = top;
                chans(chans == top) = [];
            end
            if dropNaN
                n_top = n_top(~isnan(n_top));
            end
        end
        % return all spike times during specified epoch
        function all_t = all_spike_times(obj,epoch)
            if nargin < 2 || isempty(epoch)
                epoch = obj.epoch;
            end
            if length(obj.units(end).times) > 1 && isrow(obj.units(end).times)
                all_t = cell2mat([obj.units.times]);
            else
                all_t = {obj.units.times};
                all_t(cellfun(@isempty,all_t)) = [];
                all_t = cell2mat(all_t');
            end
            all_t(all_t < epoch(1) | all_t > epoch(2)) = [];
        end
        % return all units from specific channel
        function units = channel_units(obj,chan)
            picked = zeros(1,length(obj.units));
            for u = 1:length(obj.units)
                if obj.units(u).channel == chan
                    picked(u) = 1;
                end
            end
            units = obj.units(picked == 1);
        end
        % SNR ratio of units:
        function unit_snr(obj)
            obj.snr = nan(1,length(obj.units));
            for u = 1:length(obj.units)
                if size(obj.units(u).waveforms,1) > 2
                    mnW = mean(obj.units(u).waveforms);
                    A = max(mnW) - min(mnW);
                    res = zeros(size(obj.units(u).waveforms));
                    for tt = 1:size(obj.units(u).waveforms,1)
                        res(tt,:) = obj.units(u).waveforms(tt,:) - mnW;
                    end
                    sd = mean(std(res));
                    obj.snr(u) = A/(2*sd);
                end
            end
            obj.snr(isinf(obj.snr)) = nan;
        end
        % Plot a specific unit from a specific channel:
        function unit_plot(obj,channel,unit)
            if nargin < 3 || isempty(unit)
                unit = 1;
            end
            chan_units = obj.units([obj.units.channel] == channel);
            if isempty(chan_units)
                disp([9 'No units from channel ' num2str(channel)])
            else
                if strcmp(unit,'max')
                    unit = length(chan_units); % this should always be the most active?
                end
                if length(chan_units) < unit
                    disp([9 'There are only ' num2str(length(chan_units)) ' on channel ' num2str(channel)])
                else
                    plot(chan_units(unit).waveforms');
                end
            end
        end
    end
    
end