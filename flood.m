 function grid = flood(grid, adj, sIx, verbose)
   spkID = 3;
   for i = 1:size(sIx, 2)
      for j = 1:size(sIx{i},1)
         grid = floodongrid(grid, adj, i, sIx{i}(j), spkID, verbose);
         spkID = spkID + 1;
      end
   end
end
