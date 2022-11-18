function main2(RunName, nseq,  nfra, nsnap, nbp, list_seqName_i )
%-------------------------------------------------------
% cgDNA+ function: main(runName, nseq,  nfra, nsnap, nbp, list_seqName )
%-------------------------------------------------------
% main script running the parameter exatraction and the computation of the
% statistics of md trajectories. The main input files are the base frames
% file (.fra) and the phosphate frames files (.pfra) computed with Curves+
% (lcvmm version with single precision orhigonal frames). The main outcome  
% of the analysis are two structures containing oligomer based statistics for
% cgDNA and cgDNA+ coordinates. 
%
% This script calls : param_Ext.m 
%
% Input:
%    RunName         common name of the simulations.
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
    list_seqName_s{i} = num2str(i) ; 
  end
 
else 
  count = 1;
  for i = list_seqName_i
    list_seqName_s{count} = num2str(i) ; 	
    count = count + 1 ;	  
  end
 
end

parfor j = 1:length(list_seqName_s)
 
  olig(j) = olig_paramExtFromCoordTxt( RunName,list_seqName_s{j}, nfra, nsnap, nbp ) ;  %#ok<PFOUS>
  fprintf('Computation for sequence %s done \n', list_seqName_s{j})
  
end

save( [ 'olig_uvwp' RunName '.mat' ] , 'olig'  ) ;

end
