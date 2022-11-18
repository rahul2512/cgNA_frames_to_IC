function olig = olig_paramExtFromCoordTxt( RunName, RunNbr , nfra, nsnap, nbp,index )
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
hb_count(1,Run.Nbr,index) ;    
