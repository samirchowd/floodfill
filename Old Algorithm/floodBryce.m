%%
function [spk,wav] = floodBryce(data) 
    v = 1; 
    % Loading data and instantiating objects 
    adj = load('adj');
    adj = adj.adj;
    % Setting parameters
    sigma = median(abs(data))/0.6745;
    p = 1;
    weak = 2*sigma;
    strong = 4*sigma;
    
    % Performing algorithm and extracting templates 
    spk = detectSpk(data, adj, p, weak, strong);
    plotSpk(data, spk, 19, weak, strong, v)
end
%%

%%%%%%%%%%%%%%%%%%%%%%%% Spike Detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spk = detectSpk(data, adj, p, weak, strong)
    % Function for detecting spikes on Ephys data utlizing the floodfill
    % algorithm
    
    % INPUTS: emg: Normalised Ephys object detailed in Myosort
    % adj: Matrix that indicated adjacencies between channels
    % p: Parameter to tune values of voltage weighting (See Rossant 2016)
    % weak: 1 x # of Channels vector containing weak threshold values
    % strong: 1 x # of Channels vector containing strong threshold values
    
    % OUTPUTS: A 2x1 cell. The first cell is a 1xN vector where N is the # 
    % of spikes containing the centers of each spike. The second cell is a 
    % 1xN cell 
    
    % Finding indicies of positive threshold crossings 
    ix = diff(sign(data - strong)) > 0; 
    ix = num2cell(ix, 1);
    ix = cellfun(@(x) find(x==1), ix, 'UniformOutput', false);
    
    % Instantiating binary array and final spks
    res_bin = false(size(data));
    spk = cell(1,2);
    
    % Iterating through all the positive threshold crossings
    for i = 1:size(ix,2) % Through each channel
       for j = 1:size(ix{i}, 1) % Through each crossing 
           t = ix{i}(j); % Storing current index position
           if(~res_bin(t+1,i))
               % Insantiate Spk Weighting and Spk Location Matricies 
               spk_wt = [];
               spk_loc = [];
               
               % Call floodfill on the point shifted 1 to the right (the
               % actual strong crossing).
               [spk_wt, spk_loc, res_bin] = floodfill(data, adj, t+1, i, weak, strong, spk_wt, spk_loc, res_bin);
               
               % Append spk_wt and spk_loc to the spk matrix
               spk{1}(end+1) = sum((spk_wt(:,1).*spk_wt(:,2)).^p)/sum(spk_wt(:,2)).^p;
               spk{2}{end+1} = spk_loc;
            
           end
       end
    end
    
end

function [spk_wt, spk_loc, res_bin] = floodfill(data, adj, t, c, weak, strong, spk_wt, spk_loc, res_bin)
    % Function that performs the floodfill algorithm on a given time point
    
    % INPUTS: spk_wt: vector of voltage weighted values 
    % spk_loc: values containing time-channel pairs in 2xN fashion 
    % res_bin: matrix indicating which time-channel pairs have been visited
    
    % OUTPUTS: populated spk_wt vector, spk_loc vector, and updated res_bin

    
    % Base cases
    % Looking for time/channel values out of bounds
    if(t < 0 || c < 0 || t > size(data,1) || c > size(data,2))
       return; 
    end
    
    % Ensuring data crosses weak threshold or hasn't been checked before
    if(data(t,c) <= weak(c) || res_bin(t,c))
        return; 
    end
    
    % Calulating psi value and appending it and the (t,c) to matricies 
    psiVal = psiG(data, t, c, weak(c), strong(c));
    spk_wt(end+1, :) = [t psiVal];
    spk_loc(end+1, :) = [t,c];
    
    res_bin(t,c) = 1; % Indicating time-channel pair has been visited
    
    % Recursive Section 
    % Iterating through adjacent channels
    for i = 1:size(adj,1)
       if(adj(i,1) == c)
          chanNo = adj(i,2);
          if(~res_bin(t,chanNo))
            [spk_wt, spk_loc, res_bin] = floodfill(data, adj, t, chanNo, weak, strong, spk_wt, spk_loc, res_bin);
          end
       end
    end
    
    % Checking left and right time steps 
    if(~res_bin(t+1,c))
        [spk_wt, spk_loc, res_bin] = floodfill(data, adj, t+1, c, weak, strong, spk_wt, spk_loc, res_bin);
    end
    if(~res_bin(t-1,c))
        [spk_wt, spk_loc, res_bin] = floodfill(data, adj, t-1, c, weak, strong, spk_wt, spk_loc, res_bin);
    end
end

function x = psiG(data,t,c,weak,strong)
    % Function to calculate voltage weighted values (Rossant 2016) 
    % INPUTS: Data, time, channel, weak and strong thresholds
    % OUTPUTS: The weighted voltage value for that time-chan pair 
    
    x = min(((-data(t,c) - weak) / (strong - weak)), 1);
end

%%%%%%%%%%%%%%%%%%%%%%%% Template Extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wav = getWaves(data, spk, wavedur)
    % Function that returns a 3-D matrix containing the waveforms of the
    % spks
    % OUTPUTS: 3d matrix, spike width x channel no x spike count 
    
    % Initiating variables 
    spkCenters = spk{1}; 
    wav = [];
    
    for i = 1:size(spkCenters,2) % Iterating through all the spk centers
        channels = unique(spk{2}{i}(:,2)); % Grabbing unique channels
        center = round(spkCenters(i));  
        wav = cat(3, wav, data(center-wavedur:center+wavedur,:));
        % Deviating from Rossant
        % Replace all values outside of the actual spike with scaled noise.
        % Only values in the spike range should be the actual value
        % (whether or not it crosses the threshold)
        % For multichannel spikes, use time range for each channel (not
        % total time range). 
        
        % After noise has been inserted, obtain pca for each channel across
        % all waveform (take first three components). Cluster that 
        
        for j = 1:size(wav, 2)
           if(~ismember(j, channels))
               % Calculate the stddev on chan for timewindow
               % Scale randn within that stddev
               % Replace channel with that size of noise 
               dataOnChan = data(center-wavedur:center+wavedur, j);
               avg = mean(dataOnChan);
               stdev = std(dataOnChan);
               noise = (stdev.*randn(wavedur*2 + 1,1) + avg);
               wav(:,j,i) = noise./10;
               % Is there a way to add structured noise (gaussian is eh) 
               % Fourier Analysis on noise channels, can we use that
               % fourier analysis to generate replacement noise?
               % Discuss it later, think about only noise segments 
           end
        end
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%% Data Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = plotData(data, weak, strong, v)
    % Function used to plot the data as well as weak and strong values
    % INPUTS: data matrix, weak and strong threshold values
    % OUTPUT: A graph of data on each channel on a singular plot
    
    f = axes();
    hold on;
    
    Y_PAD = 0; 
    Y_TICK = 1/2;
    
    nChannels = size(data,2);
    
    % Normalizing the waveforms 
    yAbsLim = max(abs(data),[],1);
    ws = bsxfun(@rdivide, data, yAbsLim);
    
    % Normalizing strong and weak values
    stn = bsxfun(@rdivide, strong, yAbsLim);
    wen = bsxfun(@rdivide, weak, yAbsLim);
    
    % Applying a verticle offset to each channel and thresholds
    yPos = flipud([0;cumsum((1+Y_PAD)*2*ones(nChannels-1,1))]);
    ws = ws + yPos';
    stn = stn + yPos';
    wen = wen + yPos';
    
    % y-tick position and values
    yTickVal = Y_TICK;
    yTickPos = [yPos-yTickVal, yPos, yPos+fliplr(yTickVal)];
    yTick = round((yTickPos - yPos).*yAbsLim(:));
    
    yTickPos = reshape(fliplr(yTickPos)',(1+2*length(yTickVal))*nChannels,1);
    yTick = reshape(fliplr(yTick)',(1+2*length(yTickVal))*nChannels,1);
    
    yTickPos = flipud(yTickPos);
    yTickLabel = cellfun(@num2str,flipud(num2cell(yTick)),'uni',false);
    
    % Plot
    for ii = 1:size(ws,2)
        plot(f,ws(:,ii),'k')
        if v > 0
            plot([0,size(ws,1)],[stn(ii),stn(ii)],'--r')
            plot([0,size(ws,1)],[wen(ii),wen(ii)],'--g')
        end
    end

    set(f,'yTick',yTickPos)
    set(f,'yTickLabel',yTickLabel)

    yl = [-1 1+2*(nChannels-1)];
    set(f,'yLim',yl)
    %set(f,'xLim',[min(emg.time) max(emg.time)]); % Converts x-axis to time
end 

function f = plotSpks(data, spk, weak, strong)
    % Function used to plot the detected spks
    % INPUTS: Emg (Ephys) object, spk cell, weak and strong thresholds
    % OUTPUT: A graph containing the data and the spks detected 
    
    nChannels = size(data,2);
    
    % Create axis holding the data
    f = plotData(data, weak, strong);
    Y_PAD = 0;
    
    % Getting shifted data
    % Normalizing the waveforms 
    yAbsLim = max(abs(data),[],1);
    ws = bsxfun(@rdivide, data, yAbsLim);
    
    % Applying a verticle offset to each channel and thresholds
    yPos = flipud([0;cumsum((1+Y_PAD)*2*ones(nChannels-1,1))]);
    ws = ws + yPos';
    
    % TODO: Change the for loop to work on the entire data set rather than
    % a for loop.
    
    % Plotting the spk centers and points captured in each spk
    for i = 1:size(spk{1},2) % Iterating through each spk
       % Plotting the voltage weighted spike centers
       spkCenter = round(spk{1}(i));
       chan = spk{2}{i}(1,2);
       scatter(spkCenter, ws(spkCenter, chan), 500, 'x');
       
       % Plotting the points captured by each spk
       points = [spk{2}{i}(:,1) , spk{2}{i}(:,2)];
       scatter(points(:,1), ws(sub2ind(size(ws), points(:,1), points(:,2))), 100, '.')
       
    end
    

end

function f = plotWave(wav, spkNo, strong, weak, spk, wavdur)
    % Function plots the waveform of a given spike across all channels
    wav = wav(:,:,spkNo);
    
    f = axes();
    hold on;
    
    Y_PAD = 0.5; 
    Y_TICK = 1/2;
    
    nChannels = size(wav,2);
    
    % Normalizing the waveforms 
    yAbsLim = max(abs(wav),[],1);
    yAbsLim =  reshape(ones(nChannels,1)*max(max(abs(wav))),1, nChannels);
    ws = bsxfun(@rdivide, wav, yAbsLim);
    
    % Normalizing strong and weak values
    stn = bsxfun(@rdivide, strong, yAbsLim);
    wen = bsxfun(@rdivide, weak, yAbsLim);
    
    % Applying a verticle offset to each channel and thresholds
    yPos = flipud([0;cumsum((1+Y_PAD)*2*ones(nChannels-1,1))]);
    ws = ws + yPos';
    stn = stn + yPos';
    wen = wen + yPos';

    % y-tick position and values
    yTickVal = Y_TICK;
    yTickPos = [yPos-yTickVal, yPos, yPos+fliplr(yTickVal)];
    yTick = round((yTickPos - yPos).*yAbsLim(:));
    
    yTickPos = reshape(fliplr(yTickPos)',(1+2*length(yTickVal))*nChannels,1);
    yTick = reshape(fliplr(yTick)',(1+2*length(yTickVal))*nChannels,1);
    
    yTickPos = flipud(yTickPos);
    yTickLabel = cellfun(@num2str,flipud(num2cell(yTick)),'uni',false);
    
    % Plot
    
    for ii = 1:size(ws,2)
        plot(f,ws(:,ii),'k')
        plot([0,size(ws,1)],[stn(ii),stn(ii)],'--r')
        plot([0,size(ws,1)],[wen(ii),wen(ii)],'--g')
    end

    set(f,'yTick',yTickPos)
    set(f,'yTickLabel',yTickLabel)

    yl = [-1 1+2*(nChannels-1)];
    set(f,'yLim',yl)
    %set(f,'xLim',[min(emg.time) max(emg.time)]); % Convert x-axis to time
    
    % Plotting the spk center
    spkCenter = round(spk{1}(spkNo));
    center = spkCenter - spkCenter + wavdur;
    chan = spk{2}{spkNo}(1,2);
    scatter(center, ws(center, chan), 500, 'x');
       
    % Plotting the points
    points = [spk{2}{spkNo}(:,1) , spk{2}{spkNo}(:,2)];
    points(:,1) = points(:,1) - spkCenter + wavdur + 1;
    scatter(points(:,1), ws(sub2ind(size(ws), points(:,1), points(:,2))), 100, '.')
    
end

function f = plotSpkIso(data, spk, numSpk, weak, strong)
    
    nChannels = size(data,2);

    % Extracting x values out of spike 
    xVals = min(spk{2}{numSpk}(:,1)):max(spk{2}{numSpk}(:,1)); 
    
    % Trimming data to xValues
    data = data(xVals, :); 
    
    % Plotting data
    f = plotData(data, weak, strong);
    hold on 
    Y_PAD = 0; 
    
    % Getting shifted data 
    % Normalizing the waveforms 
    yAbsLim = max(abs(data),[],1);
    ws = bsxfun(@rdivide, data, yAbsLim);
    
    % Applying a verticle offset to each channel and thresholds
    yPos = flipud([0;cumsum((1+Y_PAD)*2*ones(nChannels-1,1))]);
    ws = ws + yPos';
    
    % Plotting the voltage weighted spike centers
    spkCenter = round((xVals(end)-xVals(1))/2);
    chan = spk{2}{numSpk}(1,2);
    scatter(spkCenter, ws(spkCenter, chan), 500, 'x');
       
    % Plotting the points captured by each spk
    points = [spk{2}{numSpk}(:,1) , spk{2}{numSpk}(:,2)];
    scatter(points(:,1)-min(xVals)+1, ws(sub2ind(size(ws), points(:,1)-min(xVals)+1, points(:,2))), 100, '.') % Subtract min of xVals to get adjusted xvalues
    
end

function f = plotSpk(data, spk, numSpk, weak, strong, v)
    
    nChannels = size(data,2);

    % Extracting x values out of spike 
    xVals = min(spk{2}{numSpk}(:,1)):max(spk{2}{numSpk}(:,1));  
    
    % Plotting data
    f = plotData(data, weak, strong, v);
    hold on 
    Y_PAD = 0; 
    
    % Getting shifted data 
    % Normalizing the waveforms 
    yAbsLim = max(abs(data),[],1);
    ws = bsxfun(@rdivide, data, yAbsLim);
    
    % Applying a verticle offset to each channel and thresholds
    yPos = flipud([0;cumsum((1+Y_PAD)*2*ones(nChannels-1,1))]);
    ws = ws + yPos';
    
    % Plotting the voltage weighted spike centers
    spkCenter = round((xVals(end)+xVals(1))/2);
    chan = spk{2}{numSpk}(1,2);
    scatter(spkCenter, ws(spkCenter, chan), 1000, 'rx');
       
    % Plotting the points captured by each spk
    points = [spk{2}{numSpk}(:,1) , spk{2}{numSpk}(:,2)];
    scatter(points(:,1), ws(sub2ind(size(ws), points(:,1), points(:,2))), 100, 'b.') % Subtract min of xVals to get adjusted xvalues
    
end

function f = plotEveryNthSpk(data, spk, n, weak, strong)
    nChannels = size(data,2);
    
    % Plotting data
    f = plotData(data, weak, strong);
    hold on 
    Y_PAD = 0; 
    
    % Getting shifted data 
    % Normalizing the waveforms 
    yAbsLim = max(abs(data),[],1);
    ws = bsxfun(@rdivide, data, yAbsLim);
    
    % Applying a verticle offset to each channel and thresholds
    yPos = flipud([0;cumsum((1+Y_PAD)*2*ones(nChannels-1,1))]);
    ws = ws + yPos';
    
    % Plotting the spk centers and points captured in each spk
    for i = 1:size(spk{1},2) % Iterating through each spk
       % Plotting the voltage weighted spike centers
       if(mod(i,n) == 0)
        spkCenter = round(spk{1}(i));
        chan = spk{2}{i}(1,2);
        scatter(spkCenter, ws(spkCenter, chan), 500, 'x');
       
        % Plotting the points captured by each spk
        points = [spk{2}{i}(:,1) , spk{2}{i}(:,2)];
        scatter(points(:,1), ws(sub2ind(size(ws), points(:,1), points(:,2))), 100, '.')
       end
    end
end
















