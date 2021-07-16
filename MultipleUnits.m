classdef MultipleUnits < handle
    properties
        units       SingleUnit
        patient     char
        seizure     single
        epoch       double
        snr         double
        info        char
        extra       struct
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
        % re-order units by user-supplied order of UIDs
        function order_by_UID(obj,UIDorder)
            curOrder = [obj.units.UID];
            if length(UIDorder) ~= length(curOrder)
                error('User supplied order must be an array the same length as the number of units')
            end
            % note, this will be an issue if duplicate UIDs, but that would
            % be its own issue anyway:
            if length(unique([UIDorder curOrder])) ~= length(curOrder)
                error('Supplied UIDs must match the UIDs in the MulipleUnits units.UID field')
            end
            [~,userOrder] = sort(UIDorder);
            [~,myOrder] = sort(curOrder);
            newWorldOrder(userOrder) = myOrder;
            obj.units = obj.units(newWorldOrder);
            obj.current_order = 'user supplied';
        end
        % simple raster plot (speedy)
        function raster(obj,varargin)
            settings.axes = [];
            settings.showtypes = false;
            settings.base_color = [0 0 0];
            settings.in_color = [0.64 0.08 0.18];
            settings.in_maybe_color = [1 0.41 0.16];
            settings.linewidth = 1;
            for v = 1:2:length(varargin)
                settings.(varargin{v}) = varargin{v+1};
            end
            if isempty(settings.axes)
                settings.axes = gca;
            end
            if settings.showtypes
                if length(obj.units(end).times) > 1 && isrow(obj.units(end).times)
                    nIN = length(cell2mat([obj.units(strcmpi({obj.units.type},'in')).times]));
                    nPC = length(cell2mat([obj.units(strcmpi({obj.units.type},'pc')).times]));
                    nINmaybe = length(cell2mat([obj.units(strcmpi({obj.units.type},'in?')).times]));
                else
                    nIN = length(cell2mat({obj.units(strcmpi({obj.units.type},'in')).times}'));
                    nPC = length(cell2mat({obj.units(strcmpi({obj.units.type},'pc')).times}'));
                    nINmaybe = length(cell2mat({obj.units(strcmpi({obj.units.type},'in?')).times}'));
                end
                
                nTotalSpikesIN = 3 * nIN;
                nTotalSpikesPC = 3 * nPC;
                nTotalSpikesINmaybe = 3 * nINmaybe;
                xPointsIN = NaN(nTotalSpikesIN*3,1);
                yPointsIN = xPointsIN;
                xPointsPC = NaN(nTotalSpikesPC*3,1);
                yPointsPC = xPointsPC;
                xPointsINmaybe = NaN(nTotalSpikesINmaybe*3,1);
                yPointsINmaybe = xPointsINmaybe;
                currentIndIN = 1;
                currentIndPC = 1;
                currentIndINmaybe = 1;
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
                    elseif strcmpi(obj.units(u).type,'in?')
                        xPointsINmaybe(currentIndINmaybe:currentIndINmaybe+nSpikes*3-1) = trialXPoints;
                        yPointsINmaybe(currentIndINmaybe:currentIndINmaybe+nSpikes*3-1) = trialYPoints;
                        currentIndINmaybe = currentIndINmaybe + nSpikes*3;
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
                plot(settings.axes,xPointsPC,yPointsPC,'color',settings.base_color,'linewidth',settings.linewidth);
                hold(settings.axes,'on')
                plot(settings.axes,xPointsIN,yPointsIN,'color',settings.in_color,'linewidth',settings.linewidth);
                plot(settings.axes,xPointsINmaybe,yPointsINmaybe,'color',settings.in_maybe_color,'linewidth',settings.linewidth);
            else
                plot(settings.axes,xPoints,yPoints,'color',settings.base_color,'linewidth',settings.linewidth);
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
                            %{
                            if exist('distinguishable_colors','file')
                                cols = distinguishable_colors(length(settings.highlight));
                            else
                                cols = lines(length(settings.highlight));
                            end
                            %}
                            cols = obj.dstngsh_cols(length(settings.highlight));
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
        % plot all waveforms from each unit on one channel (3D separation)
        function plot_channel_units(obj,chan,darkmode)
            if nargin < 2 || isempty(chan)
                error('Need a channel number to plot units from')
            end
            if nargin < 3 || isempty(darkmode)
                darkmode = true;
            end
            
            these_units = obj.channel_units(chan);
            
            bgcol = [1 1 1];
            if darkmode
                bgcol = [0.1412 0.1529 0.1804];
            end
            %{
            if exist('distinguishable_colors','file')
                cols = distinguishable_colors(length(these_units),bgcol);
            else
                cols = lines(length(these_units));
            end
            %}
            cols = obj.dstngsh_cols(length(these_units),bgcol);
            
            hfig = figure('Position',[rand(1,1)*200+100 rand(1,1)*200+50 1200 500]);
            ax(1) = axes('Position',[0.05 0.05 0.4 0.9]);
            ax(2) = axes('Position',[0.55 0.05 0.4 0.9]);
            hold(ax(1),'all')
            hold(ax(2),'all')
            ledge = cell(1,length(these_units));
            assigns = nan(1,length(cell2mat({these_units.times}')));
            assign_offset = 0;
            for u = 1:length(these_units)
                t = -.6:1/30:((size(these_units(u).waveforms,2)-1)/30)-.6;
                zt = ones(size(these_units(u).waveforms)).*these_units(u).times;
                if isinf(max(abs(obj.epoch)))
                    warning('Epoch for this dataset does not appear to have been set, which will cause issues with data plotting in plot_channel_units')
                end
                zt = zt - obj.epoch(1);
                zt = zt / diff(obj.epoch);
                zt = zt + u;
                % "Jittering" the color code makes them easier to see, but
                % also makes it slower to plot (single waveform at a time)
                for w = 1:size(these_units(u).waveforms,1)
                    plot3(ax(1),t,zt(w,:),these_units(u).waveforms(w,:),...
                        'color',obj.jitter_color(cols(u,:),0.4));
                end
                ledge{u} = ['UID ' num2str(these_units(u).UID)];
                assigns((1:length(these_units(u).times))+assign_offset) = ones(1,length(these_units(u).times)) * these_units(u).UID;
                assign_offset = assign_offset + length(these_units(u).times);
            end
            tt = title(ax(1),['All units from channel ' num2str(chan)]);
            tt.FontSize = 13;
            xlim(ax(1),[-0.6 1])
            set(ax(1),'box','off','tickdir','out','ticklength',[0.005 0.005],...
                'linewidth',1.5,'FontSize',14,'color','none','FontName','Helvetica Neue',...
                'xcolor','k','ycolor','k','XGrid','on','YGrid','on','ZGrid','on',...
                'YTick',1:length(these_units),'YTickLabel',ledge)
            xlabel(ax(1),'Time (ms)')
            ylabel(ax(1),'Unit number')
            zlabel(ax(1),'Voltage (\muV)')
            view(ax(1),3)
            
            full_waves = cell2mat({these_units.waveforms}');
            [~,pc] = pca(full_waves);
            unq_assigns = unique(assigns,'stable');
            for u = 1:length(unq_assigns)
                plot3(ax(2),pc(assigns == unq_assigns(u),1),pc(assigns == unq_assigns(u),2),pc(assigns == unq_assigns(u),3),'.','color',cols(u,:))
            end
            set(ax(2),'box','on','tickdir','out','ticklength',[0.005 0.005],...
                'linewidth',1.5,'FontSize',14,'color','none','FontName','Helvetica Neue',...
                'xcolor','k','ycolor','k','XGrid','on','YGrid','on','ZGrid','on')
            lg = legend(ax(2),ledge);
            xlabel(ax(2),'PC 1 score')
            ylabel(ax(2),'PC 2 score')
            zlabel(ax(2),'PC 3 score')
            view(ax(2),3)
            if darkmode
                hfig.Color = bgcol;
                for a = 1:2
                    set(ax(a),'XColor',[0.6 0.6 0.6],'YColor',[0.6 0.6 0.6],'ZColor',[0.6 0.6 0.6]);
                end
                lg.Color = 'none';
                lg.TextColor = [0.6 0.6 0.6];
                tt.Color = [0.6 0.6 0.6];
            end
        end
        % show top n channels with most units on them (by channel number,
        % not electrode label. Use top_electrodes for electrode label)
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
        % show top n electrodes with most units on them (by electrode
        % label, not channel number. Use top_channels for channel number)
        function n_top = top_electrodes(obj,n)
            if nargin < 2 || isempty(n)
                n = 3;
            end
            chans = cell(1,length(obj.units));
            for u = 1:length(obj.units)
                chans{u} = obj.units(u).electrodelabel;
            end
            [unq_chans,~,which_chan] = unique(chans);
            n = min(n,length(unq_chans));
            n_top = cell(1,n);
            for nn = 1:n
                topCh = mode(which_chan);
                top = unq_chans{topCh};
                if nargout < 1
                    disp([9 'Electrode ' top ' has ' num2str(length(find(which_chan == topCh))) ' units']);
                end
                n_top{nn} = top;
                which_chan(which_chan == topCh) = [];
                which_chan(which_chan > topCh) = which_chan(which_chan > topCh) - 1;
                unq_chans(topCh) = [];
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
                if (ischar(chan) && strcmpi(obj.units(u).electrodelabel, chan)) || ...
                        (isnumeric(chan) && obj.units(u).channel == chan)
                    picked(u) = 1;
                end
            end
            units = obj.units(picked == 1);
        end
        % return all units from specific electrode label
        function units = electrode_units(obj,elec)
            unitElecs = {obj.units.electrodelabel};
            inds = contains(unitElecs,elec);
            units = obj.units(inds);
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
        % Get gaussian estimate of firing across population:
        function [fr, tt] = gaussian_fr(obj,SD,forced_timings,matchScaling)
            if nargin < 2 || isempty(SD)
                SD = 200;
            end
            if nargin < 3 || isempty(forced_timings)
                full_t = obj.all_spike_times();
                forced_timings = [floor(min(full_t)) ceil(max(full_t))];
            end
            if nargin < 4 || isempty(matchScaling)
                matchScaling = false;
            end
            
            offset = min(forced_timings);
            adjusted_timings = forced_timings - offset;
            adjusted_timings = [floor(adjusted_timings(1)) ceil(adjusted_timings(2))];
            
            tt = (adjusted_timings(1)*1e3:adjusted_timings(2)*1e3);
            full_fr = zeros(length(obj.units),length(tt));
            
            if verLessThan('matlab', '8.0.1')
                x = get(0,'CommandWindowSize');
            else
                x = matlab.desktop.commandwindow.size;
            end
            tot = x(1)-1;
            disp(['Calculating Gaussian estimate of firing rate across ' num2str(length(obj.units)) ' units'])
            disp([9 '(using a bin width of ' num2str(SD) ' ms, from ' num2str(forced_timings(1)) ' to ' num2str(forced_timings(2)) ' seconds)'])
            prntd = 0;
            for u = 1:length(obj.units)
                full_fr(u,:) = obj.units(u).gaussian_fr(SD,forced_timings,matchScaling);
                prc = (u/length(obj.units)) * tot;
                if floor(prc) > prntd
                    fprintf('|');
                    prntd = prntd + 1;
                end
            end
            fprintf('\n');
            tt = tt / 1e3; % back to seconds
            tt = tt + offset;
            fr = sum(full_fr);
        end
        % Plot firing rate changes for two epochs (and return the values)
        function [sd, rawChange, fr] = fr_changes(obj, epochA, epochB, varargin)
            % Calculate (and plot if 'plot', true) firing rate changes
            % between two epochs across the whole population.
            % Takes 2+ input arguments listing times to compare, each a 1x2
            % vector denoting start and end times for each epoch.
            % Extra input arguments are settings given as name, value pairs
            % as below:
            %  'plot':      true/false, whether or not to make a figure of
            %               changes
            %  'grid':      show or hide grid on figure ('on' or 'off')
            %  'scale':     'log' or 'linear' for the figure axes (log by
            %               default)
            %  'links':     yet to finish coding, does nothing right now.
            %  'colors':    3x3 array of rgb values for significant
            %               increases, non-significant changes, and then
            %               significant decreases.
            %  'sd_cutoff': multiple of SD away from preictal in Poisson
            %               distribution beyond which is deemed significant
            %  'scaling':   whether to use match_confidences to scale
            %               results (true/false, defaults to false)
            %  'axes':      handle to axes object to plot to (will plot to
            %               new figure if none provided)
            %
            % Returns 3 outputs:
            %   1) the standard deviation of the firing rate change, as
            %      calculated based on the SD of a Poisson distribution
            %      from each epoch's duration;
            %   2) the raw change in firing rate (where 1 = no change, < 1
            %      means epochB had a lower firing rate & > 1 means epochB
            %      had a higher firing rate; and
            %   3) the firing rates of the two epochs as an nx2 vector
            %      where n is the number of units in the object (spikes per
            %      second)
            %
            % Note that the plot will not show units with zero firing rate
            % when 'scale' is set to 'log'.
            settings.plot = true;
            settings.grid = 'on';
            settings.scale = 'log';
            settings.links = 'off';
            settings.colors = [
                0.75 0 0.15;
                0.2 0.2 0.2;
                0 0.447 0.741;
                ];
            settings.sd_cutoff = 3;
            settings.scaling = false;
            settings.axes = [];
            
            allowable = fieldnames(settings);
            if mod(length(varargin),2) ~= 0
                error('Inputs must be in name, value pairs');
            end
            for v = 1:2:length(varargin)
                if find(ismember(allowable,varargin{v}))
                    settings.(varargin{v}) = varargin{v+1};
                else
                    disp([9 'Not assigning ''' varargin{v} ''': not a setting option for MultipleUnits.fr_changes method']);
                end
            end
            
            sd = NaN(length(obj.units), 1);
            rawChange = NaN(length(obj.units), 1);
            fr = NaN(length(obj.units), 2);
            for u = 1:length(obj.units)
                [sd(u), rawChange(u), fr(u,:)] = obj.units(u).fr_change(epochA,epochB,settings.scaling);
            end
            
            if settings.plot
                warning('off','MATLAB:Axes:NegativeDataInLogAxis');
                if isempty(settings.axes) || ~isgraphics(settings.axes)
                    figure;
                    settings.axes = gca;
                end
                hold(settings.axes,'all')
                mx = max(fr(:));
                
                scale = linspace(0, mx, 1000);
                scaleAconf = settings.sd_cutoff*(sqrt(scale/range(epochA)));
                scaleBconf = settings.sd_cutoff*(sqrt(scale/range(epochB)));
                
                plot(settings.axes,scale,scale,'k','linewidth',4,'LineStyle','--')
                plot(settings.axes,scale+scaleAconf,scale-scaleBconf,':',...
                    'color',[0.4 0.4 0.4],'linewidth',4)
                plot(settings.axes,scale-scaleAconf,scale+scaleBconf,':',...
                    'color',[0.4 0.4 0.4],'linewidth',4)
                
                cols = repmat(settings.colors(2,:),[length(sd) 1]);
                cols(sd < -settings.sd_cutoff,:) = repmat(settings.colors(3,:), [length(find(sd < -settings.sd_cutoff)) 1]);
                cols(sd > settings.sd_cutoff,:) = repmat(settings.colors(1,:), [length(find(sd > settings.sd_cutoff)) 1]);
                
                scatter(settings.axes,fr(:,1),fr(:,2),120,cols,'filled');
                
                set(settings.axes,'XScale',settings.scale,'YScale',settings.scale,'FontSize',12,...
                    'XGrid',settings.grid,'YGrid',settings.grid,'TickDir','out')
                
                line(settings.axes,([mx mx]/1000)+(mx/1000/5),[mx/1000 mx+(mx/10)],'color',grey)
                line(settings.axes,[mx/1000 mx+(mx/10)],([mx mx]/1000)+(mx/1000/5),'color',grey)
                axis(settings.axes,[mx/1000 mx+(mx/10) mx/1000 mx+(mx/10)])
                xlabel(settings.axes,['Mean firing rate from ' num2str(epochA(1)) ' s  to ' num2str(epochA(2)) ' s (spikes s^{-1})'],...
                    'fontsize',14)
                ylabel(settings.axes,['Mean firing rate from ' num2str(epochB(1)) ' s  to ' num2str(epochB(2)) ' s (spikes s^{-1})'],...
                    'fontsize',14)
                title(settings.axes,'Significant firing rate changes')
                axis(settings.axes,'square')
                warning('on','MATLAB:Axes:NegativeDataInLogAxis');
            end
        end
        % Return the unit with a specific UID
        function unit = get_unit(obj,UID)
            unit = obj.units([obj.units.UID] == UID);
        end
        % Save this structure in the default/specified location
        function save(obj,varargin)
            settings.path = mfilename('fullpath');
            settings.path = [fileparts(settings.path) filesep 'Data'];
            settings.name = 'data';
            
            allowable = fieldnames(settings);
            if mod(length(varargin),2) ~= 0
                error('Inputs must be in name, value pairs');
            end
            for v = 1:2:length(varargin)
                if find(ismember(allowable,varargin{v}))
                    settings.(varargin{v}) = varargin{v+1};
                else
                    disp([9 'Not assigning ''' varargin{v} ''': not an option in PatientDB.backup()']);
                end
            end
            if ~strcmp(settings.name,'data')
                warning(['Changed the save method and now must use ''data'' as save variable, not ''' settings.name ''' '])
            end
            data = obj;
            if ~exist(settings.path,'dir')
                mkdir(settings.path);
            end
            
            savename = [settings.path filesep 'NeuroClass_' obj.patient '_s' num2str(obj.seizure) '.mat'];
            %{
            % THIS METHOD MAKES HUGE FILES (it saves the whole class object
            % & thus a load of data that it isn't dependent on if you have
            % the class definition on your path)
            mat = matfile(savename,'Writable',true);
            mat.(settings.name) = obj;
            %}
            save(savename,'data');
            disp(['Saved NeuroClass data as variable ''' settings.name ''' in:'])
            disp([9 savename]);
            clear mat
        end
        
        function calculateMetrics(obj)
            % Calculate/update the Gaussian estimations of false +ve/-ves
            % in the "metrics" field of each child SingleUnit object. To
            % update the other metrics, use the equivalent calculateMetrics
            % method within each child SingleUnit object.
            % Use the SingleUnit.metrics methods falsePositiveRate() and
            % falseNegativeRate() to return the reportable rates as
            % described in Hill et al., JNeurosci, 2011.
            % Depends on the original UltraMegaSort2000 being on the path.
            
            if ~exist('gaussian_overlap','file')
                error('Need the original gaussian_overlap() function from UltraMegaSort2000 on the path to calculate Gaussian estimates of false +ves/-ves')
            end
            % Make sure each SingleUnit has an active UnitMetrics object in
            % the metrics field:
            for u = 1:length(obj.units)
                if isempty(obj.units(u).metrics)
                    obj.units(u).metrics = UnitMetrics();
                end
            end
            chans = obj.top_channels(length(obj.units),true);
            for c = 1:length(chans)
                chanUnits = obj.channel_units(chans(c));
                if length(chanUnits) > 1
                    pairs = nchoosek(1:length(chanUnits),2);
                    fp = NaN(length(chanUnits));
                    fn = NaN(length(chanUnits));
                    for p = 1:size(pairs,1)
                        conf = gaussian_overlap(chanUnits(pairs(p,1)).waveforms,chanUnits(pairs(p,2)).waveforms);
                        fp(pairs(p,1),pairs(p,2)) = conf(1,1);
                        fp(pairs(p,2),pairs(p,1)) = conf(2,2);
                        fn(pairs(p,1),pairs(p,2)) = conf(1,2);
                        fn(pairs(p,2),pairs(p,1)) = conf(2,1);
                    end
                    for u = 1:length(chanUnits)
                        obj.units([obj.units.UID] == chanUnits(u).UID).metrics.gmFalsePos = fp(u,setdiff(1:length(chanUnits),u));
                        obj.units([obj.units.UID] == chanUnits(u).UID).metrics.gmFalseNeg = fn(u,setdiff(1:length(chanUnits),u));
                        
                        obj.units([obj.units.UID] == chanUnits(u).UID).metrics.gmUIDs = [chanUnits(setdiff(1:length(chanUnits),u)).UID];
                    end
                else
                    obj.units([obj.units.UID] == chanUnits.UID).metrics.gmFalsePos = 0;
                    obj.units([obj.units.UID] == chanUnits.UID).metrics.gmFalseNeg = 0;
                    obj.units([obj.units.UID] == chanUnits.UID).metrics.gmUIDs = [];
                end
            end
        end
    end
    
    methods (Static = true)
        function colorcode = jitter_color(colorcode,amt)
            if nargin < 1 || isempty(colorcode) || length(colorcode) ~= 3
                error('Must supply a color code in matlab format (1 x 3, 0-1 range)')
            end
            
            if nargin < 2 || isempty(amt)
                amt = 1/3;
            end
            
            transposed = false;
            
            if iscolumn(colorcode)
                colorcode = colorcode';
                transposed = true;
            end
            for a = 1:3
                jitter = rand(1,1);
                if colorcode(:,a) == 1
                    jitter = -jitter * amt;
                elseif colorcode(:,a) == 0
                    jitter = jitter * amt;
                else
                    jitter = (jitter-0.5) * amt;
                end
                colorcode(:,a) = colorcode(:,a)+jitter;
            end
            
            colorcode(colorcode > 1) = 1;
            colorcode(colorcode < 0) = 0;
            
            if transposed
                colorcode = colorcode';
            end
        end
        
        function cols = dstngsh_cols(n_colors,bg)
            % a lightweight version of distinguisable_colors by Tim Holy
            % (https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors)
            if nargin < 2 || isempty(bg)
                bg = [1 1 1];
            end
            
            x = linspace(0,1,30);
            [R,G,B] = ndgrid(x,x,x);
            rgb = [R(:) G(:) B(:)];
            
            C = makecform('srgb2lab');
            lab = applycform(rgb,C);
            bglab = applycform(bg,C);
            
            mindist2 = inf(size(rgb,1),1);
            
            cols = zeros(n_colors,3);
            lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
            for i = 1:n_colors
                dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
                dist2 = sum(dX.^2,2);  % square distance
                mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
                [~,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
                cols(i,:) = rgb(index,:);  % save for output
                lastlab = lab(index,:);  % prepare for next iteration
            end
        end
    end
    
end
