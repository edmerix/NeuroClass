classdef SingleUnit < handle
    properties
        UID                 single
        waveforms           double
        times       (:,1)   double
        channel             single
        electrodelabel      char
        patient             char
        seizure     (1,1)   single
        epoch       (1,2)   double = [-Inf Inf]
        wideband    (1,:)   double
        type                char
        typeprob    (:,1)   single
        threshold           double
        extra % used to assign anything extra by the user
    end

    properties (SetAccess = private, Hidden = true)
        gaussian_estimate
        autocorr_data
        xcorr_data
        isi
        isi_bins
    end

    methods
        % constructor
        function obj = SingleUnit(varargin)
            allowable = fieldnames(obj);
            if mod(length(varargin),2) ~= 0
                error('Inputs must be in name, value pairs');
            end
            for v = 1:2:length(varargin)
                if find(ismember(allowable,varargin{v}))
                    obj.(varargin{v}) = varargin{v+1};
                else
                    disp([9 'Not assigning ''' varargin{v} ''': not a property of SingleUnit class']);
                end
            end
            if ~isempty(obj.times)
                obj.times = obj.times(:)'; % make sure it's a column vector
            end
        end
        % create the firing rate estimate
        function [fr, tt] = gaussian_fr(obj,SD,forced_timings,matchScaling)
            % Convolve a Gaussian of supplied SD (in ms, defaults to
            % 200) over the spike times (resolution of 1 ms).
            % Can also supply "forced_timings" which will keep the min and
            % max bins the same regardless of when the earliest or latest
            % spike occurs for this unit.
            % If "matchScaling" is set to true, then output firing rate 
            % will be scaled by the probability of each spike being a true 
            % match to that unit. Tell Ed to write this better.
            if nargin < 2 || isempty(SD)
                SD = 200;
            end
            if nargin > 2 && ~isempty(forced_timings)
                offset = min(forced_timings);
            else
                forced_timings = [min(obj.times) max(obj.times)];
                offset = floor(min(obj.times));
            end
            forced_timings = forced_timings - offset;
            forced_timings = [floor(forced_timings(1)) ceil(forced_timings(2))];
            
            times_in = obj.times - offset;
            keepTimes = times_in <= forced_timings(1) | times_in > forced_timings(2);
            times_in(keepTimes) = [];
            if nargin < 4 || isempty(matchScaling)
                matchScaling = false;
            end
            
            gaussSize = [1 SD*10]; % size of Gaussian window (make sure it doesn't get clipped)
            w = fspecial('gaussian',gaussSize,SD);
            
            tt = forced_timings(1)*1e3:forced_timings(2)*1e3;
            ap_times = zeros(size(tt));
            times_in = round(times_in*1e3); % use milliseconds
            if matchScaling
                if isfield(obj.extra,'match_confidence')
                    ap_times(times_in) = obj.extra.match_confidence(~keepTimes);
                elseif isfield(obj.extra,'probabilities')
                    ap_times(times_in) = obj.extra.probabilities(~keepTimes);
                else
                    error('Need to have waveform match confidence stored in the "extra" field under either "match_confidence" or "probabilities"')
                end
            else
                ap_times(times_in) = 1;
            end
            
            fr = conv(ap_times,w,'same');
            fr = fr * 1e3;  % back to seconds
            tt = tt/1e3;    % back to seconds
            tt = tt + offset;

            gInfo.rate = fr;
            gInfo.time = tt;
            gInfo.SD = SD;
            gInfo.forced_t = forced_timings + offset;
            obj.gaussian_estimate = gInfo;
        end
        % retrieve the gaussian estimate rate for use outside the object
        % (the gaussian estimate object is hidden and private; I like clean
        % objects...)
        function [fr, tt] = retrieve_gauss_fr(obj)
            if isempty(obj.gaussian_estimate)
                disp('Gaussian estimate of firing rate has not yet been performed');
                fr = [];
                tt = [];
            else
                fr = obj.gaussian_estimate.rate;
                tt = obj.gaussian_estimate.time;
            end
        end
        % plot the gaussian firing rate estimate
        function plot_gauss_fr(obj,win,forced_timings)
            redo = false;
            if nargin < 3 || isempty(forced_timings)
                forced_timings = [min(obj.times) max(obj.times)];
            end
            if isempty(obj.gaussian_estimate)
                redo = true;
                if nargin < 2
                    error(['Gaussian estimate of firing rate has not been set yet,'...
                        ' please include a window for new estimate']);
                end
            else
                if nargin > 2
                    if min(obj.gaussian_estimate.forced_t) > min(forced_timings)...
                            || max(obj.gaussian_estimate.forced_t) < max(forced_timings)
                        redo = true;
                    end
                end
                if nargin > 1 && ~isempty(win)
                    if obj.gaussian_estimate.window ~= win
                        redo = true;
                    end
                end
            end
            if redo
                obj.gaussian_fr(win,forced_timings);
            end
            plot(obj.gaussian_estimate.time,obj.gaussian_estimate.rate,'linewidth',2);
            grid on
            xlim(forced_timings);
            xlabel('Time (s)')
            ylabel('Probability of 1 spike in window')
        end
        % plot a regular histogram firing rate
        function plot_hist_fr(obj,bins)
            if nargin < 2
                error('Need an input bin width/vector of bins');
            end
            if isscalar(bins)
                bins = min(obj.times)+(bins/2):bins:max(obj.times)-(bins/2);
            end
            binW = bins(2)-bins(1);
            x = hist(obj.times,bins);
            bar(bins,x/binW,1);
            xlabel('Time (s)')
            ylabel('Spikes s^{-1}')
            title(['Firing rate in ' num2str(binW) ' s bins'])
            if ~isempty(obj.epoch)
                xlim(obj.epoch) % in case the unit only fired irregularly, adjust xlim to the actual recording length
            end
        end
        % get histogram based firing rate:
        function [fr,bins] = hist_fr(obj,bins)
            if nargin < 2
                error('Need an input bin width/vector of bins');
            end
            if isscalar(bins)
                bins = min(obj.times):bins:max(obj.times);
            end
            d = diff(bins)/2;
            edges = [bins(1)-d(1), bins(1:end-1)+d, bins(end)+d(end)];
            fr = histcounts(obj.times,edges);
        end
        % calculate autocorrelation
        function [ac,lags] = autocorr(obj,bins,t_subset)
            if nargin < 2 || isempty(bins)
                error('Include bins (in ms) over which to calculate the autocorrelation');
            end
            if isscalar(bins)
                error('Bin input must be literal rather than a bin width so as to know min/max lags');
            end
            if nargin > 2 && ~isempty(t_subset)
                sub_times = obj.times(obj.times > min(t_subset) & obj.times < max(t_subset));
            else
                t_subset = obj.epoch;
                sub_times = obj.times;
            end
            
            d = diff(bins)/2;
            edges = [bins(1)-d(1), bins(1:end-1)+d, bins(end)+d(end)];
            
            bigMat = repmat(sub_times,1,length(sub_times));
            bigMat = bigMat - diag(bigMat)';
            bigMat = bigMat * 1e3; % convert to milliseconds

            vals = histcounts(bigMat(:),edges);
            [~,wh] = min(abs(edges));
            vals(wh) = vals(wh) - length(sub_times); % remove self from AC
            ac_data.xc = vals;
            ac_data.lags = bins;
            ac_data.time_subset = t_subset;
            obj.autocorr_data = ac_data;
            ac = ac_data.xc;
            lags = ac_data.lags;
        end
        % plot autocorrelation
        function plot_ac(obj)
            if isempty(obj.autocorr_data)
                error('Run ''autocorr'' method first');
            end
            bar(obj.autocorr_data.lags,obj.autocorr_data.xc,1);
            if min(obj.autocorr_data.lags) < 0 && max(obj.autocorr_data.lags) > 0
                hold on
                line([0 0],ylim,'color','r','linestyle','--');
            end
            xlim([min(obj.autocorr_data.lags) max(obj.autocorr_data.lags)]);
            xlabel('Lags (ms)')
            ylabel('Counts')
            title(['Autocorrelation from ' ...
                num2str(min(obj.autocorr_data.time_subset)) ' to ' ...
                num2str(max(obj.autocorr_data.time_subset)) ' seconds'])
        end
        % retrieve autocorrelation
        function [ac,lags] = retrieve_ac(obj)
            if isempty(obj.autocorr_data)
                error('Run ''autocorr'' method first');
            end
            ac = obj.autocorr_data.xc;
            lags = obj.autocorr_data.lags;
        end
        % calculate mean AC lag:
        function mean_ac = mean_ac_lag(obj,upto)
            if nargin < 2 || isempty(upto)
                upto = 100;
            end
            sub_times = obj.times*1000;
            all_t = cell(1,length(sub_times));
            for t = 1:length(sub_times)
                adj_t = sub_times - sub_times(t);
                adj_t = adj_t(adj_t > 0 & adj_t <= upto);
                all_t{t} = adj_t;
            end
            mean_ac = mean(cell2mat(all_t(:)));
        end
        % calculate ISI
        function obj = ISI(obj,bins)
             % add an extra bin to delete to avoid all values beyond that
             % value being counted in that bin
            temp_bins = [bins bins(end)+(bins(2)-bins(1))];
            ISI = diff(obj.times*1000);
            obj.isi = hist(ISI,temp_bins);
            obj.isi(end) = [];
            obj.isi_bins = bins;
        end
        % plot ISI distribution
        function plot_isi(obj)
            if isempty(obj.isi) || isempty(obj.isi_bins)
                disp('Run ''ISI'' method first');
            else
                binW = obj.isi_bins(2) - obj.isi_bins(1);
                bar(obj.isi_bins,obj.isi,1)
                xlabel('Delay (ms)')
                ylabel('Counts')
                title(['ISI distribution in ' num2str(binW) ' ms bins'])
                xlim([min(obj.isi_bins) max(obj.isi_bins)])
            end
        end
        % retrieve the ISI info for use outside the object
        function [ISI, bins] = retrieve_ISI(obj)
            if isempty(obj.isi) || isempty(obj.isi_bins)
                disp('Run ''ISI'' method first');
                ISI = [];
                bins = [];
            else
                ISI = obj.isi;
                bins = obj.isi_bins;
            end
        end
        % calculate firing rate changes
        function [sd, rawChange, rates] = fr_change(obj,epochA,epochB,scaling)
            % Takes 2 input arguments listing times to compare, each a 1x2
            % vector denoting start and end times for each epoch.
            % Allows a third input argument to use match probability
            % scaling (true/false)
            % Returns 3 outputs: 
            %   1) the standard deviation of the firing rate change, as
            %      calculated based on the SD of a Poisson distribution
            %      from each epoch's duration;
            %   2) the raw change in firing rate (where 1 = no change, < 1
            %      means epochB had a lower firing rate & > 1 means epochB
            %      had a higher firing rate; and
            %   3) the firing rates of the two epochs as a 1x2 vector
            %      (spikes per second)
            if nargin < 3 || isempty(epochA) || isempty(epochB)
                error('Need two inputs of which times to compare, each a 1x2 vector denoting start and finish times for each epoch')
            end
            if nargin < 4 || isempty(scaling)
                scaling = false;
            end
            epochA = sort(epochA);
            epochB = sort(epochB);
            if ~scaling
                rateA = (length(find(obj.times > epochA(1) & obj.times <= epochA(2))))/range(epochA);
                rateB = (length(find(obj.times > epochB(1) & obj.times <= epochB(2))))/range(epochB);
            else
                rateA = sum(obj.extra.match_confidence(obj.times > epochA(1) & obj.times <= epochA(2)))/range(epochA);
                rateB = sum(obj.extra.match_confidence(obj.times > epochB(1) & obj.times <= epochB(2)))/range(epochB);
            end
            rates = [rateA rateB];
            
            if rateA == rateB
                rawChange = 1;
                sd = 0;
                return;
            end
            
            rawChange = rateB/rateA;
            
            % calculate the distance from equal firing rate:
            intersect = (rateA+rateB)/2;
            orthDist = pdist([rateA rateB; intersect intersect]);
            
            % calculate what 1 SD firing rate changes would be in a
            % Poisson-distribution using these epoch durations:
            lowConf = [intersect+(sqrt(intersect/diff(epochA))) intersect-(sqrt(intersect/diff(epochB)))];
            highConf = [intersect-(sqrt(intersect/diff(epochA))) intersect+(sqrt(intersect/diff(epochB)))];
            % calculate what those distances would be, for 
            lowBoundDist = pdist([lowConf;(sum(lowConf))/2 (sum(lowConf))/2]);
            hiBoundDist = pdist([highConf;(sum(highConf))/2 (sum(highConf))/2]);
            
            if rawChange > 1
                sd = orthDist/hiBoundDist;
            else
                sd = -orthDist/lowBoundDist;
            end
        end
        % inspect unit
        function inspect_unit(obj)
            %TODO: replace ISI with detection metric
            detects = obj.waveforms(:,19);
            [h,x] = histcounts(detects,100);
            x = (diff(x)/2)+x(1:end-1); % convert bin edges into bin centers

            %{
            % ISI:
            bins = 0:1:100;
            temp_bins = [bins bins(end)+(bins(2)-bins(1))];
            ISI = diff(obj.times*1000);
            this_isi = hist(ISI,temp_bins);
            this_isi(end) = [];
            this_isi_bins = bins;
            totLow = length(find(ISI < 2));
            %}
            % AUTOCORR:
            bins = -100:1:100;
            bins = bins/1000; % temporarily work in seconds (easier to change bins than all times)
            d = diff(bins)/2;
            edges = [bins(1)-d(1), bins(1:end-1)+d, bins(end)+d(end)];
            maxLag = max(bins);
            minLag = min(bins);
            result = zeros(length(obj.times),length(bins));
            for t = 1:length(obj.times)
                adj_times = setdiff(obj.times,obj.times(t));
                subset = adj_times(adj_times >= obj.times(t)+minLag & adj_times <= obj.times(t)+maxLag);
                subset = subset - obj.times(t);
                result(t,:) = histcounts(subset,edges);
            end
            xc = sum(result);
            lags = bins*1000;
            % Figure
            figure('units','normalized','position',[0.02 0.02 0.96 0.96]);
            ax(1) = axes('position',[0.035 0.55 0.35 0.4]);
            t = -0.6:1/30:((length(obj.waveforms(1,:))-1)/30)-0.6;
            zt = ones(size(obj.waveforms)).*obj.times;
            plot3(ax(1),t,zt,obj.waveforms)
            view(ax(1),[0 1])
            xlim(ax(1),[-0.6 1])
            title(ax(1),[num2str(length(obj.times)) ' waveforms from channel ' num2str(obj.channel) ' (rotate for through time)'])
            xlabel(ax(1),'Time (ms)')
            ylabel(ax(1),'Time (s)')
            zlabel(ax(1),'Voltage (\muV)')
            set(ax(1),'ydir','reverse')
            ax(2) = axes('position',[0.425 0.55 0.25 0.4]);
            %{
            bar(this_isi_bins,this_isi,1)
            line([2 2],[0 max(ylim)],'color','r','linestyle','--')
            xlim([0 100])
            xlabel('Time (ms)')
            ylabel('Count')
            title(['ISI (' num2str(totLow) ' < 2 ms)'])
            %}
            bar(x,h,1)
            xlabel('Voltage (\muV)')
            ylabel('Counts')
            title('Voltage at detection (unfinished plot...)')
            ax(3) = axes('position',[0.72 0.55 0.27 0.4]);
            bar(lags,xc,1)
            line([0 0],[0 max(ylim)],'color','k','linestyle','--')
            xlim([-100 100])
            xlabel('Time (ms)')
            ylabel('Count')
            title('Autocorrelation')
            % CHI dist
            [z,dof] = zvals(obj);
            [~,x1] = hist(z,100);
            ax(4) = axes('position',[0.05 0.04 0.4 0.42]);
            hist(z,100)
            hndl = findobj(gca,'Type','patch');
            set(hndl,'FaceColor', [0 0 1])
            y =  chi2pdf(x1,dof);
            y = y * length(z) * ( x1(2)-x1(1));
            l = line(x1,y);
            set(l,'Color',[0 1 0],'LineWidth',1.5)
            title('Chi squared distribution')
            xlabel('Mahalanobis distance')
            ylabel('Count')
            [~,pc] = pca(obj.waveforms);
            ax(5) = axes('position',[0.5 0.04 0.5 0.45]);
            plot3(pc(:,1),pc(:,2),pc(:,3),'.')
            rotate3d('on')
            xlabel('PC 1'),ylabel('PC 2'),zlabel('PC 3')

            if exist('cleanupfig','file')
                for a = 1:length(ax)
                    cleanupfig(ax(a),'grid','smallfont');
                end
            end
        end

        % calculate cross-correlation with another unit
        function [xc, lags] = xcorr(obj,unit,bins,t_subset)
            if nargin < 2 || isempty(unit)
                error('Need a SingleUnit object as the first input to calculate cross-correlation relative to');
            end
            if nargin < 3 || isempty(bins)
                error('Include bins (in ms) over which to calculate the crosscorrelation');
            end
            if isscalar(bins)
                error('Bin input must be literal rather than a bin width so as to know min/max lags');
            end
            if nargin > 3 && ~isempty(t_subset)
                sub_times = obj.times(obj.times > min(t_subset) & obj.times < max(t_subset));
            else
                sub_times = obj.times;
            end
            bins = bins/1000; % temporarily work in seconds (easier to change bins than all times)
            d = diff(bins)/2;
            edges = [bins(1)-d(1), bins(1:end-1)+d, bins(end)+d(end)];
            maxLag = max(bins);
            minLag = min(bins);
            result = zeros(length(sub_times),length(bins));
            for t = 1:length(sub_times)
                subset = unit.times(unit.times >= sub_times(t)+minLag & unit.times <= sub_times(t)+maxLag);
                subset = subset - sub_times(t);
                result(t,:) = histcounts(subset,edges);
            end
            xc = sum(result);
            lags = bins*1000; % put back to milliseconds
        end
        
        function weightedMean = weightedMeanWaveform(obj)
            % Calculate the mean waveform, but weighted by each spike's
            % match confidence
            wvs = obj.waveforms;
            if isfield(obj.extra,'match_confidence')
                conf = obj.extra.match_confidence;
            elseif isfield(obj.extra,'probabilities')
                conf = obj.extra.probabilities;
            else
                error('Need to have waveform match confidence stored in the "extra" field under either "match_confidence" or "probabilities"')
            end
            weightedMean = (size(wvs,1)/sum(conf)) * mean(conf'.*wvs',2);
        end
        
        function [z, dof] = zvals(obj)
            % Tweaked/borrowed from UltraMegaSort2000 "get_zvalues.m" by
            % Hill DN, Mehta SB, & Kleinfeld D
            w = obj.waveforms;
            covar = cov(w);
            dof = round(rank(covar)/2);
            num_dims = size(w(:,:),2);
            num_spikes = size(w,1);
            last = (1:dof) + num_dims  - dof;
            [v,d] = eig(covar);            % get PCs
            for j = 1:num_dims, v(:,j) = v(:,j); end
            v = v(:,last);                 % use last r dimensions
            w = detrend(w,'constant');     % mean subtract
            w = (w*v);                     % project on to PCs

            % get Mahalanobis distance
            z = zeros([1 num_spikes]);
            dinv = inv(d(last,last));
            for j = 1:num_spikes
                z(j) = w(j,:)*dinv*w(j,:)';
            end
        end
    end
end
