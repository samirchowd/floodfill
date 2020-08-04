% Instantiate EMG data
% Set weighting factor for spk times >> p (See Rossant 2016) 
% Set strong and weak threshhold crossing parameters 
%%
function spk = floodBryce
    emg = load('emg');
    sigma = load('sigma');
    adj = load('adj');
    emg = emg.emg;
    sigma = sigma.sigma;
    adj = adj.adj;
    p = 1;
    weak = 2*sigma;
    strong = 4*sigma;
    %f = plotData(emg);
    spk = detectSpk(emg, adj, p, weak, strong)
end
%%
% Find positive strong threshold crossings 
% >> find(diff(sign(x-threshold)) > 0) >> ix   

% Instantiate a binary false valued  matrix, same size as the EMG data >> res_bin 
% Instantiate spk cell array >> spk = {{[1 x N spikes]} {N arrays of data coordinates}} 

% Iterate through all positive threshold crossings sequentially through
% time through each channel

% At each strong threshold crossing: 

% If ~res_bin(t,c)
% >> Set res_bin(t,c)= 1  
% >> Instantiate spk time weighting matrix >> spk_wt
% >> Instantiate spk location matrix >> spk_loc
% >> Calculate psi(t,c) --> See Rossant 2016
% >> spk_wt = [t psi(t,c)] 
% >> [spk_wt, spk_loc, res_bin] = floodfll(data, t, c, p, spk_wt, spk_loc, res_bin)
% >> spk{1}(end+1) = sum((spk_wt(1,:).*spk_wt(2,:)).^p)/sum(spk_wt(2,:)).^p
% >> spk{2}{end+1} = spk_loc

function spk = detectSpk(emg, adj, p, weak, strong)
    % Function for detecting spikes on Ephys data
    
    % Finding indicies of positive threshold crossings 
    ix = diff(sign(emg.data - strong)) > 0; 
    ix = num2cell(ix, 1);
    ix = cellfun(@(x) find(x==1), ix, 'UniformOutput', false);
    
    % Instantiating binary array and final spks
    res_bin = false(size(emg.data));
    spk = cell(1,2);
    
    % Iterating through all the positive threshold crossings
    for i = 1:size(ix,2) % Through each channel
       for j = 1:size(ix{i}, 1) % Through each crossing 
           t = ix{i}(j); % Storing current index position
           if(~res_bin(t,i))
               % res_bin(t,i) = 1;
               psiVal  = psiG(emg.data, t, i, weak, strong);
               spk_wt = [t psiVal];
               spk_loc = [t i];
               %spk_wt = [];
               %spk_loc = [];
               [spk_wt, spk_loc, res_bin] = floodfill(emg.data, adj, t, i, weak, strong, spk_wt, spk_loc, res_bin);
               spk{1}(end+1) = sum((spk_wt(:,1).*spk_wt(:,2)).^p)/sum(spk_wt(:,2)).^p;
               spk{2}{end+1} = spk_loc;
            
           end
       end
    end
    
end

% FUNCTION 
% FLOOD FILL 

% Iterate through adjacent channels to see if it is greater than weak threshold crossing 
% If greater than weak threshold crossing && ~res_bin(t,c)  
% >> Calculate psi(t,c) --> See Rossant 2016
% >> spk_wt(end+1, :) = [t psi(t,c)] --> See Rossant 2016 
% >> spk_loc(end+1, :) = [t, c]
% Set res_bin element to 1

% Check next adjacent channel 
% >> for adj(1,i) == c
% >> >> if ~res_bin(t,i)
% >> >> >> Call function on (t, i, spk_wt)
% Check adjacent times 
% >> If ~res_bin(t+1, c)
% >> >> Call function on (t+1, c, spk_wt)
% >> if ~res_bin(t-1, c)
% >> >> Call function on (t-1, c, spk_wt) 

function [spk_wt, spk_loc, res_bin] = floodfill(data, adj, t, c, weak, strong, spk_wt, spk_loc, res_bin)
    % Function that performs the floodfill algorithm on a given time point
    
    % Base cases
    if(t < 0 || c < 0 || t > size(data,1) || c > size(data,2))
       return; 
    end
    
    if(data(t,c) <= weak(c) || res_bin(t,c))
        return; 
    end
    
    % If the point crosses weak and has not been checked before
    %if(data(t,c) > weak(c) && ~res_bin(t,c))
    psiVal = psiG(data, t, c, weak(c), strong(c));
    spk_wt(end+1, :) = [t psiVal];
    spk_loc(end+1, :) = [t,c];
    
    res_bin(t,c) = 1; 
    
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
    x = min(((-data(t,c) - weak) / (strong - weak)), 1);
end

function f = plotData(emg, weak, strong)
    f = axes();
    hold on;
    
    Y_PAD = 0; 
    Y_TICK = 1/2;
    
    % Normalizing the waveforms 
    yAbsLim = max(abs(emg.data),[],1);
    ws = bsxfun(@rdivide, emg.data, yAbsLim);
    stn = bsxfun(@rdivide, strong, yAbsLim);
    wen = bsxfun(@rdivide, weak, yAbsLim);
    
    % Applying a verticle offset to each channel
    yPos = flipud([0;cumsum((1+Y_PAD)*2*ones(emg.nChannels-1,1))]);
    ws = ws + yPos';
    stn = stn + yPos';
    wen = wen + yPos';
    
    % y-tick position and values
    yTickVal = Y_TICK;
    yTickPos = [yPos-yTickVal, yPos, yPos+fliplr(yTickVal)];
    yTick = round((yTickPos - yPos).*yAbsLim(:));
    
    yTickPos = reshape(fliplr(yTickPos)',(1+2*length(yTickVal))*emg.nChannels,1);
    yTick = reshape(fliplr(yTick)',(1+2*length(yTickVal))*emg.nChannels,1);
    
    yTickPos = flipud(yTickPos);
    yTickLabel = cellfun(@num2str,flipud(num2cell(yTick)),'uni',false);
    
    % Plot
    
    for ii = 1:size(ws,2)
        plot(f,ws(:,ii))
        plot([0,size(ws,1)],[stn(ii),stn(ii)])
        plot([0,size(ws,1)],[wen(ii),wen(ii)])
    end
    drawnow

    set(f,'yTick',yTickPos)
    set(f,'yTickLabel',yTickLabel)

    yl = [-1 1+2*(emg.nChannels-1)];
    set(f,'yLim',yl)
    %set(f,'xLim',[min(emg.time) max(emg.time)]);
    
end 

function f = plotSpks(data, spk, weak, strong)
    
    % Create axis holding the data
    f = plotData(data, weak, strong);
    

end



