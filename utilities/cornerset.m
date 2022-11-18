function sten = cornerset(nbp,phosphate_flag)

%-------------------------------------------------------
% cgDNA+ function: sten = cornerset(nbp,flag_phosphate)
%-------------------------------------------------------
% Computes the stencil of the banded stiffness matrix for the cgDNA+ model
% or the cgDNA model. The stencil is given as a index set of the corner of
% the blocks of the stiffness matrix.
%
%
% Input:
%    nbp             number of base pairs
%    phosphate_flag  flag to cgDNA+ coord (flag=1) or cgDNA coord (flag=0)
%
% Output:
%    sten            matrix containing the index set of the corners of the
%                    blocks of the stiffness matrix. [Size (nbp-1)x2]
%
%
% Note 1:
%  TO DO: rewrtie this part
%   The cgDNA stiffness matrix are banded matrix with 18x18 block on the
%   diagonal which overlap are 6x6 blocks while the cgDNA+ stiffness matrix
%   have 42x42 blocks on the diagonal wiht 18x18 blocks as overlap but the
%   first and the last blocks are 30x30.
%
% If you find this code useful, please cite:
% TODO: add reference
%
%-------------------------------------------------------

sten = zeros(nbp-1,2) ;

% Compute stencil for cgDNA+ stiffness
if phosphate_flag
  sten(1,:) = [ 1 , 36 ] ;
  for i = 1 : nbp-3 ;
    sten(i+1,:) = [ 19 + 24*(i-1) , 36 + 24*i ] ;
    
  end
  sten(nbp-1,:) = [ 24 + sten(nbp-2,1) , 18 + sten(nbp-2,2) ] ;
  
% Compute stencil for cgDNA stiffness
else
  for i = 1 : nbp-1 ;
    sten(i,:) = [ 1 + 12*(i-1) , 18 + 12*(i-1) ] ;
    
  end
  
end

end

