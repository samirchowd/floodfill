classdef floodfill
    properties
        emg % Ephys object
        sigma % Standard deviation across each channel 
        adj % Adjacency matrix  
        strong % Logical matrix of strong threshold crossings
        weak % Logical matrix of weak threshold crossings
        sIx % Indicies of crossings in the strong matrix
        wIx % Indicies of crossings in the weak matrix 
        grid % array 
        verbose % Level of output 
    end
    
    methods
       
        function obj = floodfill(emg, sigma, adj, varargin)
            p = inputParser; 
            
            % Set up mandatory inputs
            validateObjEmg = @(x) assert(isa(x, 'Ephys'), "emg must be a valid Ephys object");
            p.addRequired('emg', validateObjEmg);
            
            validateObjSigma = @(x) assert(isnumeric(x), "sigma must be a valid numeric array");
            p.addRequired('sigma', validateObjSigma)
            
            validateObjAdj = @(x) assert(isnumeric(x) && all(all(rem(x,1) == 0) == 1), "adj must be a valid integer array");
            p.addRequired('adj', validateObjAdj)
            
            % Set up optional inputs 
            
            % Verbose 
            defaultVerbose = 2; 
            validateVerbose = @(x) assert(rem(x,1) == 0, "verbose must be a valid integer");
            p.addParameter('verbose', defaultVerbose, validateVerbose)
            
            % Strong 
            defaultStrong = [];
            validateStrong = @(x) assert(islogical(x), "strong must be a valid logical array");
            p.addParameter('strong', defaultStrong, validateStrong)
            
            % Weak
            defaultWeak = [];
            validateWeak = @(x) assert(islogical(x), "weak must be a valid logical array");
            p.addParameter('weak', defaultWeak, validateWeak)
            
            % sIx
            default_sIx = {};
            validate_sIx = @(x) assert(iscell(x), 'sIx must be a valid cell array');
            p.addParameter('sIx', default_sIx, validate_sIx)
            
            % wIx
            default_wIx = {};
            validate_wIx = @(x) assert(iscell(x), 'wIx must be a valid cell array');
            p.addParameter('wIx', default_wIx, validate_wIx)
            
            p.KeepUnmatched = 1;
            p.parse(emg, sigma, adj, varargin{:});
            
            % Adding properties
            obj.emg = p.Results.emg; 
            obj.sigma = p.Results.sigma;
            obj.adj = p.Results.adj;
            obj.strong = p.Results.strong;
            obj.weak = p.Results.weak;
            obj.sIx = p.Results.sIx;
            obj.wIx = p.Results.wIx;
            obj.verbose = p.Results.verbose;
            
            
            if obj.verbose > 0; disp([newline 'SUCCESS: floodfill instantiated']); end 
            
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Getters & Setters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = getCross(obj)
            % Recreate data in binary fashion of crossings for thresholds
            data = obj.emg.data; 
            strongc = sign(data - 4*obj.sigma) > 0;
            weakc = sign(data - 2*obj.sigma) > 0;
            obj.strong = strongc;
            obj.weak = weakc;
        end
        
        function obj = getIx(obj)
            % Find indicies of crossings
            obj.sIx = num2cell(obj.strong,1);
            obj.wIx = num2cell(obj.weak,1);
            
            obj.sIx = cellfun(@(x) find(x==1), obj.sIx, 'UniformOutput', false);
            obj.wIx = cellfun(@(x) find(x==1), obj.wIx, 'UniformOutput', false);
        end
        
        function obj = getGrid(obj)
            % Generate grid of 1's and 2's based on double threshold
            obj.grid = zeros(obj.emg.nChannels, size(obj.strong, 1));
            
            for i = 1:obj.emg.nChannels
                for j = 1:size(obj.strong,1)
                    if(obj.strong(j,i) == 1)
                        obj.grid(i, j) = 2;
                    elseif (obj.weak(j,i) == 1)
                        obj.grid(i,j) = 1;
                    end
                end
            end
        end
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%% Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function drawData(obj)
            s = stackedplot(obj.emg.data(:, :)./obj.sigma);
           
            % Setting axis properties 
            s.Title = ['EMG Data'];
            s.GridVisible = 'on';
            
            % Determining min and max values for uniform scale
            minY = min(min(obj.emg.data(:, :)));
            maxY = max(max(obj.emg.data(:, :)));
            
            for i = 1:obj.emg.nChannels
                s.AxesProperties(i).YLimits = [minY, maxY];
                s.DisplayLabels(i) = {sprintf('Channel %d', i)};
            end
       end 
        
        function drawDataWindow(obj, ix)
            s = stackedplot(obj.emg.data(ix-60:ix+60, :));
           
            % Setting axis properties 
            s.XData = s.XData + ix - 60; % Shifting x axis forward
            s.Title = ['EMG Data Centered Around ' num2str(ix)];
            s.GridVisible = 'on';
            
            % Determining min and max values for uniform scale
            minY = min(min(obj.emg.data(ix-60:ix+60, :)));
            maxY = max(max(obj.emg.data(ix-60:ix+60, :)));
            
            for i = 1:obj.emg.nChannels
                s.AxesProperties(i).YLimits = [minY, maxY];
                s.DisplayLabels(i) = {sprintf('Channel %d', i)};
            end
        end 
        
        function drawCrossWindow(obj, ix)
            % Generate a plot of threshold crossings based on the index
            timeWindow = (ix-60):(ix+60);
            window = zeros(obj.emg.nChannels,size(timeWindow,2));
            
            % Set values for color map
            for i = 1:obj.emg.nChannels
                for j = timeWindow(1):timeWindow(end)
                    if(obj.strong(j,i) == 1)
                        window(i, j - timeWindow(1) +1) = 2;
                    elseif (obj.weak(j,i) == 1)
                        window(i,j - timeWindow(1) + 1) = 1;
                    end
                end
            end
            % Generate plot
            imagesc(window)
            colormap(gray)
            hold on
            xlabel("time (in samples)")
            ylabel("channel no.")
            title("Threshold Crossing")
        end
        
        function drawSpkWindow(obj, ix)
           % Generate a plot of captured spikes based on the index
           timeWindow = (ix-60):(ix+60);
           window = obj.grid(:, timeWindow(1):timeWindow(end));
           
           imagesc(window)
           colormap(flag)
           hold on
           xlabel("time (in samples)")
           ylabel("channel no.")
           title("Detected Spikes")
           
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%% Floodfill %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function grid = floodOnChannel(obj, chanNo, ix, spkID)
            grid = floodonchan(obj.grid, chanNo, ix, spkID, obj.verbose);
        end 
        
        function grid = floodOnGrid(obj, chanNo, ix, spkID)
            grid = floodongrid(obj.grid, obj.adj, chanNo, ix, spkID, obj.verbose);
        end
 
        function grid = floodFill(obj)
           grid = flood(obj.grid, obj.adj, obj.sIx, obj.verbose);
        end
        
    end
    
end










