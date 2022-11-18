function olig = main1(runName, nseq,  nfra, nsnap, nbp, list_seqName )
%-------------------------------------------------------
% cgDNA+ function: main(runName, nseq,  nfra, nsnap, nbp, list_seqName )
%-------------------------------------------------------
% main script running the parameter exatraction and the computation of the
% statistics of md trajectories. The main input files are the base frames
% file (.fra) and the phosphate frames files (.pfra) computed with Curves+
% (lcvmm version with single precision orthogonal frames). The main outcome  
% of the analysis are two structures containing oligomer based statistics for
% cgDNA and cgDNA+ coordinates. 
%
% This script calls : param_Ext.m 
%
% Input:
%    runName         common name of the simulations.
%    nseq            number of sequences to analyse.
%    nfra            number of .fra and .pfra to analyse.
%    nsnap           number of snapshot in a single .fra and .pfra file.
%    nbp             number of base pairs of the oligomer.
%    list_seqName    list of the sequence names to analyse. If no list is
%                    given the name of the sequences are assumed to be the 
%                    numbers from 1 to nseq.
%
% Output:
%   no output. 
%    
% Note 1:
% The set of script runned using main.m are suitable for a specific folder 
% organization set up of result of AMBER md simulation.   
% See Folder_Organization.txt for more information.
%
%
% If you find this code useful, please cite:
% TODO: add reference
%
%-------------------------------------------------------

if nargin < 6
  
  for i = 1:nseq
    list_seqName{i} = num2str(i) ; 
  end
  
end

parfor j = 1:nseq
  
  fileName = [ runName '_' num2str(j) ] ;
  
  Run.Name = fileName;
  Run.Nbr = num2str(j);
  
  cd (list_seqName{j})
    olig(j) = olig_paramExt( fileName, nfra, nsnap, nbp ) ;         
  cd ../
  fprintf('Computation for sequence %i done \n', j)
  
end

save( [ 'olig_' runName '.mat' ] , 'olig'  ) ;

end
