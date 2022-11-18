function uvwp  = mergecoord( uvw, pho, nsnap, nbp )

%-------------------------------------------------------
% cgDNA+ function: uvwp  = mergecoord( uvw, pho, nsnap, nbp )
%-------------------------------------------------------
% Merge uvw coordinates (intra inter basepair variables) with phopshate
% coordinate to obtain uvwp coordinates
% 
% TODO: Complete documentation.
%
% Input:
%    uvw             struct variable with timeseries of cgDNA coords
%                    (see Note 1)      
%    pho             struct variable with timeseries of phosphate coords
%                    (see Note 1) 
%    nsnap           number of snapshots
%    nbp             number of basepair of the oligomer
%
% Output:
%    uvwp            struct variable with timeseries of cgDNA+ coords
%    
%
% Note 1:
% 
%   'uvw' is a single-entry struct array, with each field encoding
%   the timeseries of length nsnap of a 3-dimensional variable at each
%   level of a nbp oligomer, i.e. the name and size of each field is:        
%   - eta    [ nsnap x nbp x 3 ] 
%   - w      [ nsnap x nbp x 3 ]
%   - u      [ nsnap-1 x nbp x 3 ]
%   - v      [ nsnap-1 x nbp x 3 ]
%   'pho' is a single-entry struct array, with each field encoding
%   the timeseries of length nsnap of a 3-dimensional variable at each
%   level of a nbp oligomer, i.e. the name and size of each field is:   
%   - etapW  [ nsnap-1 x nbp x 3 ]
%   - wpW    [ nsnap-1 x nbp x 3 ]
%   - etapC  [ nsnap-1 x nbp x 3 ]
%   - wpC    [ nsnap-1 x nbp x 3 ]
%
% If you find this code useful, please cite:
% TODO: add reference
% 
%-------------------------------------------------------

uvwp = struct('eta'   , zeros(nsnap,nbp,3), ... 
              'w'     , zeros(nsnap,nbp,3), ...
              'etapW' , zeros(nsnap,nbp-1,3), ...
              'wpW'   , zeros(nsnap,nbp-1,3), ...
              'u'     , zeros(nsnap,nbp-1,3), ...
              'v'     , zeros(nsnap,nbp-1,3), ...
              'etapC' , zeros(nsnap,nbp-1,3), ...
              'wpC'   , zeros(nsnap,nbp-1,3) ) ;

 uvwp.eta = uvw.eta ;
 uvwp.w   = uvw.w ;
 uvwp.u   = uvw.u ; 
 uvwp.v   = uvw.v ;
 
 uvwp.etapW = pho.etapW ;
 uvwp.wpW   = pho.wpW ;
 uvwp.etapC = pho.etapC ;
 uvwp.wpC   = pho.wpC ;
            
            
end

