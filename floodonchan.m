function grid = floodonchan(grid, chanNo, ix, spkID, v)
    
    if v > 0; disp(['Checking ' num2str(ix)]); end
    % Base Cases
    if(~ismember(chanNo, 1:size(grid,1)) || ~ismember(ix, 1:size(grid, 2)))
        if v > 0; disp(['Case 1 satisfied']); end
        return 
    end

    if(grid(chanNo, ix) ~= 1 && grid(chanNo, ix) ~= 2) 
        if v > 0; disp(['Case 2 satisfied']); end
        return 
    end 
    
    if v > 0; disp(['Setting ' num2str(ix) ' >> ' num2str(spkID)]); end
    grid(chanNo, ix) = spkID; 

    if v > 0; disp(['Checking right of ' num2str(ix)]); end
    grid = floodonchan(grid, chanNo, ix+1, spkID, v);

    if v > 0; disp(['Checking left of ' num2str(ix)]); end
    grid = floodonchan(grid, chanNo, ix-1, spkID, v);

    if v > 0; disp(['Returning Control from ' num2str(ix) ]); end


end 

