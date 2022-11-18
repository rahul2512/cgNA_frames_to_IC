function uvwp_out = mergeframes(uvwp_in)
%-------------------------------------------------------
% cgDNA+ function: uvwp_out = mergeframes(uvwp_in)
%-------------------------------------------------------
% Merge cgDNA+ coordinates extracted from multiple .fra and .pfra files
% stor
% 
% TODO: Complete documentation.
%
% Input:
%    uvwp_in         array of struct variable with timeseries of cgDNA+ coordinates
%                    (see Note 1)      
%
%
% Output:
%    uvwp_out        struct variable with timeseries of cgDNA+ coords
%         
%
% Note 1:
% 
%   'uvwp_in' is an array of struct, with dimension equal the number of frames
%   files processed, with each field encoding the timeseries of length nsnap 
%   of a 3-dimensional variable at each level of a nbp oligomer, i.e. 
%   the name and size of each field is:        
%   - eta    [ nsnap x nbp x 3 ] 
%   - w      [ nsnap x nbp x 3 ]
%   - etapW  [ nsnap-1 x nbp x 3 ]
%   - wpW    [ nsnap-1 x nbp x 3 ]
%   - u      [ nsnap-1 x nbp x 3 ]
%   - v      [ nsnap-1 x nbp x 3 ]
%   - etapC  [ nsnap-1 x nbp x 3 ]
%   - wpC    [ nsnap-1 x nbp x 3 ]  
%
%
%
% If you find this code useful, please cite:
% TODO: add reference
% 
%-------------------------------------------------------
field_name = {'eta','w','etapW','wpW','u','v','etapC','wpC'} ;

nfields = length(field_name) ;

uvwp_out = struct(uvwp_in(1)) ;

for j = 1:nfields
  
  uvwp_out.(field_name{j}) = cat( 1 , uvwp_in.(field_name{j}) ) ;
  
end

end

