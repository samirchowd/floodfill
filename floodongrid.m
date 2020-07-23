function [grid, adj] = floodongrid(grid, adj, chanNo, ix, spkID, v)

    if v > 0; disp(['Checking ' num2str(ix)]); end
    
    % Base Cases
    
    % Return if the spike or channels are out of bounds
    if(~ismember(chanNo, 1:size(grid,1)) || ~ismember(ix, 1:size(grid,2)))
        if v > 0; disp(['Case 1 satisfied']); end
        return 
    end

    % Return if the value is not a spike or has already been given
    % a spike ID
    
    if(grid(chanNo, ix) ~= 1 && grid(chanNo, ix) ~= 2) 
        if v > 0; disp(['Case 2 satisfied']); end
        return 
    end 
    
    % Set the value of the grid to the spike ID
    if v > 0; disp(['Setting ' num2str(ix) ' >> ' num2str(spkID)]); end
    grid(chanNo, ix) = spkID;
    
    % Recurse on adjacent channels up and down* (According to adj
    % matrix)
    for i = 1:size(adj,1)
       if(adj(i,1) == chanNo)
          if v > 0; disp(['Checking from channel ' num2str(chanNo) ' to ' num2str(adj(i,2))]); end
          [grid, adj] = floodongrid(grid, adj, adj(i,2), ix, spkID, v);
       end
    end
    
    % Recurse on the timesteps to the left and right
    if v > 0; disp(['Checking right of ' num2str(ix)]); end
    [grid, adj] = floodongrid(grid, adj, chanNo, ix+1, spkID, v);
    
    if v > 0; disp(['Checking left of ' num2str(ix)]); end
    [grid, adj] = floodongrid(grid, adj, chanNo, ix-1, spkID, v);
    
    if v > 0; disp(['Returning Control from ' num2str(ix) ]); end

end

