function [av , cv, stiff_me, nhb] = base3(timeseries,path_to_filter, sten,debug)

%-------------------------------------------------------
% cgDNA function: average, covariance = base3(timeseries,filter,debug)
%-------------------------------------------------------
% Calculates first and second moments of coordinates (eta, w, u, v),
% optionally filtering the trajectories using the provided filter.
%
% TODO: Complete documentation.
%
% Input:
%    timeseries      unraveled vector of coordinates at each snapshot
%                    (see base2 and base2p, unravel).
%    filter          vector of 0 or 1, of same length as shape, or
%                    empty vector if no filtering required.
%    debug           flag to debug execution.
%
% Output:
%    average         average coordinate vector in
%                    non-dimensional Curves+ form [size N x 1].
%    covariance      coordinate covariance matrix in
%                    non-dimensional Curves+ form [size N x N].
%    where N is the total number of coordinates.
%
% If you find this code useful, please cite:
% TODO: add reference
%
%-------------------------------------------------------

tic;


if isempty(path_to_filter)
  hb = 1:size(timeseries,1);
else
  filter = load(path_to_filter) ;
  
  iall = filter.iall;
  hb = find(iall) ;
end
nhb = size(hb);

%% filter snapshots
y = timeseries(hb,:);
%% average
av =  mean(y, 1);
%% covariance
cv = cov(y);

% Compute the stiffness matrix using the Tricky Inversion (Maximum Entropy Fit)
% sten = stencil, is the set of the position of the cornsers of the stencil
[~, stiff_me] = completeMaxEntropy(cv,sten, 1);


if debug
  fprintf('Done. Total time %6.1f min\n',toc/60);
end





