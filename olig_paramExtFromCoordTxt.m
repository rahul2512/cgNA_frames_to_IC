function olig = olig_paramExtFromCoordTxt( RunName, RunNbr , nfra, nsnap, nbp )
%function olig = olig_paramExtFromCoordTxt( RunName, RunNbr , nfra, nsnap, nbp )

%-------------------------------------------------------
% cgDNA+ function: olig_uvw, olig_uvwp = param_Ext(nbrSeq,ruName,nsnap,nbp)
%-------------------------------------------------------
% Computes statistics of intenral coordinates from MD timeseries. The
% internal coordinates should be computed before running this scripts by
% using the function coordinates_extraction.m
%
% TODO: Complete documentation.
%
% Input:
%    RunName         main file name for .fra and .pfra
%    RunNbr          Number of the run
%    nfra            number of .fra and .pfra files to process
%    nsnap           number of snapshot in a .fra and .pfra
%    nbp             number of base pairs of the oligomer
%
% Output:
%    olig            structure containing cgDNA statistics for the oligomer
%                    (see Note 1).
%
%
% Note 1:
%
%   'olig' is a struct with fields:
%    - nbp       number of base pairs of the oligomer
%    - seq       sequence of the oligomer
%    - shape     average oligomer shape :
%                [Size: N x 1]
%    - s1b       raw stiffness of the oligomer
%                [Size: N x N]
%    - stiff_me  banded stiffness computed as the maximum entropy fit
%                [Size: N x N]
%    - nsnap     number of accepted snapshot after filtering
%
% Note 2:
%
%
% If you find this code useful, please cite:
% TODO: add reference
%
%-------------------------------------------------------

warning('off')

Run.Name = [ RunName '_' RunNbr ] ;
Run.Nbr = RunNbr ;

%%
if size(nfra,2) == 1
    nfra_id = 1 : nfra ;
    fprintf('Extracting coordinates for file numbers %i -> %i for sequence %s \n', 1 , nfra, Run.Name) ;
else
    nfra_id = nfra(1) : nfra(2) ;
    fprintf('Extracting coordinates for file numbers %i -> %i for sequence %s \n', nfra(1) , nfra(2), Run.Name) ;
end

%%
S = zeros(24*nbp-18,24*nbp-18,length(nfra_id)) ;
m = zeros(24*nbp-18,length(nfra_id)) ;

path_to_filter = [ Run.Nbr '/anl/HB/par_hbonds.mat' ] ;

%% Check if hydrogen bond filter file exist otherwise comupte it
if exist(path_to_filter, 'file') == 0
    % Extract hydrogen bond filter and save it in ./anl/HB/par_hbonds.mat
    hb(1,Run.Nbr) ;
    
end

%% Load the hb filter
filter = load(path_to_filter) ;

filter = reshape(filter.iall, [ nsnap length(filter.iall)/nsnap ]) ;

filter = filter(:,nfra_id); 

nhb = 0 ;
%% Loop over coord files
for fra = nfra_id
    
    filter_tmp = find(filter(:,fra)); 
    % Number of accepted snapshots
    nhb_tmp = size(filter_tmp,1);

    if nhb_tmp ~= 0
        
        filename = [ './' Run.Nbr '/Coord/' Run.Name '.' num2str(fra) '.coord'  ] ;
        
        timeseries_tmp = CoordFromFile(filename,nbp,nsnap) ;
        timeseries_tmp = timeseries_tmp(filter_tmp,:) ;
        
        mtmp = mean(timeseries_tmp,1)' ;
        
        S(:,:,fra) = cov(timeseries_tmp)*(nhb_tmp-1) + nhb_tmp*(mtmp*mtmp') ;
        m(:,fra)  = mtmp*nhb_tmp ;
        KK(:,:,fra) = inv(cov(timeseries_tmp));
	
    end
    
    nhb = nhb + nhb_tmp ;
      
end % end loop over fra files

%% Compute final statistics
Stot = sum(S,3)/(nhb-1) ;
mu   = sum(m,2)/nhb ;
cv   = Stot - nhb/(nhb-1) * (mu*mu') ;
RunNbr
nhb
[~, stiff_me] = completeMaxEntropy(cv,cornerset(nbp,1), 1);

%% Create the output structure
olig.nbp      = nbp ;
olig.seq      = fgetl(fopen([ './' Run.Nbr '/seq.' Run.Name '.txt' ]));
olig.shape    = mu ;
olig.s1b      = inv(cv) ;
olig.stiff_me = stiff_me ;
olig.nsnap    = nhb ;

% for seq 9 trial
%for i = 1: 100;
%eigen_min(i)= min(eig(KK(:,:,i)));
%eigen_max(i)= max(eig(KK(:,:,i)));
%end;

%olig.KK    = KK;
olig.m   = m;
%olig.eigen_min   = eigen_min;
%olig.eigen_max   = eigen_max;
end
