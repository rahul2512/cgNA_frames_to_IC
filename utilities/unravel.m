function timeseries = unravel(shape,flag_coord,debug)

%-------------------------------------------------------
% cgDNA+ function: timeseries = unravel(shape,field_coord,debug)
%-------------------------------------------------------
% Combines the timeseries contained in the fields of shape
% according to the model chosed by field_coord.
%
% TODO: Complete documentation.
%
% Input:
%    shape           structure with coordinates for each snapshot
%                    (see Note 1).
%    flag_coord      integer number corrseponding to type of coordinates. 
%                    1: cgDNA+ coordinates
%                    0: cgDNA coordinates
%    debug           flag to debug execution.
%
% Output:
%    timeseries      full, ordered coordinate timeseries in
%                    non-dimensional Curves+ form [size nsnap x N].
%    
%    nsnap in the number of snapshots and N is the total number of
%    coordinates. N is given by the chosen model: 12*nbp-6 for cgDNA, 
%    24*nbp-18 for cgDNA+. 
%
% Note 1:
% 
%   'shape' is a single-entry struct array, with each field encoding
%   the timeseries of length nsnap of a 3-dimensional variable at each
%   level of a nbp oligomer, i.e. the name and size of each field are:        
%   - eta    [ nsnap x nbp x 3 ] 
%   - w      [ nsnap x nbp x 3 ]
%   - etapW  [ nsnap-1 x nbp x 3 ]
%   - wpW    [ nsnap-1 x nbp x 3 ]
%   - u      [ nsnap-1 x nbp x 3 ]
%   - v      [ nsnap-1 x nbp x 3 ]
%   - etapC  [ nsnap-1 x nbp x 3 ]
%   - wpC    [ nsnap-1 x nbp x 3 ]
%
% If you find this code useful, please cite:
% TODO: add reference
% 
%-------------------------------------------------------
tic;

if flag_coord == 1 
  field_order = {'eta','w','etapC','wpC','u','v','etapW','wpW'} ;
elseif flag_coord == 0 
  field_order = {'eta','w','u','v'} ;
  
else
  error('No field_order for the input flag_coord')
  
end

nf = numel(field_order);
nflast = nf;

%% checks:
%% 1) all field names are valid
%% 2) sizes are compatible
for field = 1:nf
    assert(isfield(shape, field_order{field}))
    sizes = size(shape.(field_order{field}));
    % only 3D vectors: nsnap x nbp x 3
    assert(numel(sizes) == 3);
    tnsnap = sizes(1);
    tnbp   = sizes(2);
    if field == 1
        nsnap = tnsnap;
        nbp   = tnbp;
    else
        % check timeseries are compatible
        assert(tnsnap == nsnap);
        % accept fewer variables only at last position
        if tnbp == nbp-1
            nflast = nflast - 1;
        else 
            assert(tnbp == nbp);
        end
    end
end

nv     = nf*3;
nvlast = nflast*3;
ndim = nv*(nbp-1)+nvlast;
timeseries = zeros(nsnap,ndim);

for k = 1:nbp
    for f = 1:nf
        if k == nbp && f > nflast
            continue;
        end
        tmp = shape.(field_order{f});
        % [k,f,size(tmp),size(timeseries) nv*(k-1)+3*(f-1)+(1:3)]
        timeseries(:,nv*(k-1)+3*(f-1)+(1:3)) = squeeze(tmp(:,k,:));
    end
end

if debug 
    fprintf('Done. Total time %6.1f min\n',toc/60);
end

end
