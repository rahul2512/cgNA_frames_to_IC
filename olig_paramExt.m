function olig = olig_paramExt( Run, nfra, nsnap, nbp )

%-------------------------------------------------------
% cgDNA+ function: olig_uvw, olig_uvwp = param_Ext(nbrSeq,ruName,nsnap,nbp)
%-------------------------------------------------------
% Computes time series of cgDNA+ coordinates from base and phosphate frames 
% extracted from .fra and .pfra file (Curves+ output). it then computes
% oligomer base statistics (average, covariance, maximum entropy fit).
% 
% TODO: Complete documentation.
%
% Input:
%    Run        	structure containing the two field :
%						- Name : full name of the run
%						- Nbr  : number of the folder containing the data
%    nfra            number of .fra and .pfra files to process.
%    nsnap           number of snapshot in a .fra and .pfra.
%    nbp             number of base pairs of the oligomer.
%
%                   
%    olig            structure containing cgDNA+ statistics for the oligomer
%                    (see Note 1).
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
% If you find this code useful, please cite:
% TODO: add reference
%
%-------------------------------------------------------

warning('off')                                                             %#ok<WNOFF>

%% 1) cgDNA+ coordinates extraction from .fra and .pfra files

for fra =1:nfra
  
  file.bframes = [ './Frames/' Run.Name '.' num2str(fra) '.fra'  ] ;
  file.pframes = [ './Frames/' Run.Name '.' num2str(fra) '.pfra' ] ;
  
  % Computes cgDNA and phosphate coordinates from .fra and .pfra
  [ uvw, pho, ~ , ~ ] = base2p(nsnap, nbp, file);
  
  % Merge the cgDNA coordinate with the phosphate coordinates
  uvwp(fra) = mergecoord(uvw, pho, nsnap, nbp) ;
  
end % end loop over frames

uvwp = mergeframes(uvwp) ;

save([ Run.Name '_uvwp_nj_5scale.mat' ], '-struct', 'uvwp' ) ;

%% 2) Compute statistics
timeseries = unravel(uvwp,1,0) ;
%
% Run should be a structure containing the name and the location of
% the HB data and should contain the name and the location of text file 
% with he sequence of the current run. 
% See lines 41,46, and 70 in base3 and also look in hb.m !!
% To be adapted for Cheatham data and maybe it is time to do something more
% general.
%
% Run.XXX = YYY 
% etc
olig = base3(timeseries,Run,nbp,1,0);


end

