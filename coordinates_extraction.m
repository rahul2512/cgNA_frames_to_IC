function coordinates_extraction( RunName, RunNbr , nfra, nsnap, nbp , path_input, path_output )
%-------------------------------------------------------
% cgDNA+ function: coordinates_extraction( RunName, RunNbr , nfra, nsnap, nbp )
%-------------------------------------------------------
% Computes time series of cgDNA+ coordinates from base and phosphate frames 
% extracted from .fra and .pfra file (Curves+ output) and save the
% coordnates in a text file. 
% 
% TODO: Complete documentation.
%
% Input:
%    RunName         main file name for .fra and .pfra.
%    RunNbr          number of the run.
%    nfra            number of snapshot in a .fra and .pfra.
%    nsnap           total number of snapshot in the .fra and .pfra file.
%    nbp             number of base pairs of the oligomer.
%
%
% If you find this code useful, please cite:
% TODO: add reference
%
%-------------------------------------------------------

warning('off')                                                            

Run.Name = [ RunName '_' RunNbr ] ;
Run.Nbr = RunNbr ;

if nargin < 6 
    path_input = ['./' Run.Nbr '/Frames/'] ; % default value for standard folder organization
    path_output= ['./' Run.Nbr '/Coord/'];   % default value for standard folder organization
end

if size(nfra,2) == 1 
    
    nfra_id = 1 : nfra ;
    fprintf('Extracting frames for file numbers %i -> %i for sequence %s \n', 1 , nfra, Run.Name) ;
else
    
    nfra_id = nfra(1) : nfra(2) ;
    fprintf('Extracting frames for file numbers %i -> %i for sequence %s \n', nfra(1) , nfra(2), Run.Name) ;
    
end

%% 1) cgDNA+ coordinates extraction from .fra and .pfra files
for fra = nfra_id
  file.bframes = [ path_input Run.Name '.' num2str(fra) '.fra'  ] ;
  file.pframes = [ path_input Run.Name '.' num2str(fra) '.pfra' ] ;
  
  % Computes cgDNA and phosphate coordinates from .fra and .pfra
  [ uvw, pho, ~ , ~ ] = base2p(nsnap, nbp, file,0,0,0);
  
  % Merge the cgDNA coordinate with the phosphate coordinates
  uvwp = mergecoord(uvw, pho, nsnap, nbp) ;                           %#ok<*AGROW>
  
  timeseries = unravel(uvwp,1,0) ;

  fid = fopen( [ path_output Run.Name '.' num2str(fra) '.coord' ] , 'w' ) ;
  
  fprintf(fid, [ repmat( '%10.6f ' , [1 size(timeseries,1)] ) '\n' ] , timeseries ) ; 
    
  fclose(fid) ; 
  
  fprintf( 'Coordinates file %i created for the Run %s \n ' , fra , Run.Name ) ;
  
end % end loop over frames


end

