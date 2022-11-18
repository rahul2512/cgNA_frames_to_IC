function olig = base3(timeseries,Run,nbp,phosphate_flag,debug)

%-------------------------------------------------------
% cgDNA function: olig = base3(uvwp,fileName,nbp,phosphate_flag,debug)
%-------------------------------------------------------
% Compute oligo-based statistics.  
%  
%
% TODO: Complete documentation.
%
% Input:
%    uvwp            struct with timeseries of cgDNA+ variables             
%    nbp             number of basepair of the oligomer 
%    phosphate_flag  flag to cgDNA+ coord (flag=1) or cgDNA coord (flag=0)             
%    debug           flag to debug execution  
%
% Output:
%    olig           structure containing cgDNA statistics for the oligomer
%                    (see Note 1).

% Note 1:
%   'olig' are struct with fields:
%    - nbp       number of base pairs of the oligomer 
%    - seq       sequence of the oligomer
%    - shape     average oligomer shape : 
%                [Size: N x 1 (see note 2)]
%    - s1b       raw stiffness of the oligomer
%                [Size: N x N (see note 2)]
%    - stiff_me  banded stiffness computed as the maximum entropy fit
%                [Size: N x N (see note 2) ]
%    - nsnap     number of accepted snapshot after filtering
%
% If you find this code useful, please cite:
% TODO: add reference
%
%-------------------------------------------------------

tic;

%% 1) Compute/load the hydrogen bonds filter 
path_to_filter = [ Run.Nbr '/anl/HB/par_hbonds.mat' ] ;

% Check if hydrogen bond filter file exist otherwise comupte it
if exist(path_to_filter, 'file') == 0
  % Extract hydrogen bond filter and save it in ./anl/HB/par_hbonds.mat
  hb(1,Run.Nbr) ;
  
end

% Load the hb filter
filter = load(path_to_filter) ;
filter = find(filter.iall) ;
% Number of accepted snapshots 
nhb = size(filter);

%% 2) Compute mean and covariance of the filtered timeseries
% Filter the timeseries
y = timeseries(filter,:);
% Compute average
mu =  mean(y, 1);
% Compute covariance
cv = cov(y);

% Compute the stiffness matrix using the Maximum Entropy Fit.
% cornerset gives the set of the position of the cornsers of the stencil
[~, stiff_me] = completeMaxEntropy(cv,cornerset(nbp,phosphate_flag), 1);

%% Finalize the output
olig.nbp      = nbp ;
olig.seq      = fgetl(fopen([ './seq.' Run.Name '.txt' ]));
olig.shape    = mu ;
olig.s1b      = inv(cv) ;
olig.stiff_me = stiff_me ;
olig.nsnap    = nhb;

if debug
  fprintf('Done. Total time %6.1f min\n',toc/60);
end

end



