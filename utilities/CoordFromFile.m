function timeseries = CoordFromFile(RunName,nbp,nsnap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(RunName,'r') ;
timeseries = fscanf(fid,'%f') ;

timeseries = reshape(timeseries , [ nsnap , 24*nbp-18 ]) ;

end

