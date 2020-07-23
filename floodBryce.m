% Instantiate EMG data
% Set weighting factor for spk times >> p (See Rossant 2016) 
% Set strong and weak threshhold crossing parameters 

% Find positive strong threshold crossings 
% >> find(diff(sign(x-threshold)) > 0)   

% Instantiate a binary false valued  matrix, same size as the EMG data >> res_bin 
% Instantiate spk cell array >> spk = {{[1 x N spikes]} {N arrays of data coordinates}} 

% Iterate through all positive threshold crossings sequentially through
% time through each channel

% At each strong threshold crossing: 

% If ~res_bin(t,c)
% Set res_bin(t,c)= 1  
% Instantiate spk time weighting matrix >> spk_wt
% Instantiate spk location matrix >> spk_loc
% >> Calculate psi(t,c) --> See Rossant 2016
% >> spk_wt = [t psi(t,c)] 
% >> [spk_wt, spk_loc, res_bin] = floodfll(data, t, c, p, spk_wt, spk_loc, res_bin)
% spk{1}(end+1) = sum((spk_wt(1,:).*spk_wt(2,:)).^p)/sum(spk_wt(2,:)).^p
% spk{2}{end+1} = spk_loc


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

