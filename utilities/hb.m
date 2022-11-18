function hb(debug,SeqNbr)

%-------------------------------------------------------
% cgDNA+ function: hb(debug,SeqNbr)
%-------------------------------------------------------
% Collect and process all the hydrogen bond data to compute the hydrogen
% bond filter based on distance and angle
%  
%
% Input:
%    debug           flag to debug execution
%    SeqNbr          number of the sequence to process            
%
%
% Output:
% No output
%                    
% If you find this code useful, please cite:
% TODO: add reference
%
%-------------------------------------------------------

if nargin < 2
  path_to_HB = './anl/HB' ;
else
    if isnumeric(SeqNbr)
        SeqNbr = num2str(SeqNbr) ;
    end
  path_to_HB = [ SeqNbr '/anl/HB' ] ;
end

dist_files  = dir( [ path_to_HB '/all_hd*' ] ) ;
angle_files = dir( [ path_to_HB '/all_ht*' ] ) ;

nbrfiles = length(dist_files) ;

if nbrfiles ~= length(angle_files) 
  error('Error: the number of hd files must be the same of number of ht files ')
end

dis = dlmread( [ path_to_HB '/' dist_files(1).name  ] , '\t', 0, 0 ) ;
ang = dlmread( [ path_to_HB '/' angle_files(1).name ] , '\t', 0, 0 ) ;

nsnap = length(ang)*nbrfiles ;

if debug
  fprintf('loading and merging hd and ht files ...\n')
end

for i = 2 : nbrfiles
  currdist=dlmread( [ path_to_HB '/all_hd_' num2str(i) '.out' ] , '\t', 0, 0 ) ;
  n = length(currdist) ;
  dis(end+1:end+n,:) = currdist(:,:) ;
  
  currang=dlmread( [ path_to_HB '/all_ht_' num2str(i) '.out' ] , '\t', 0, 0 ) ;
  m = length(currang) ;
  
  ang(end+1:end+m,:) = currang(:,:) ;
  
end

dis = dis(:,1:end) ;
ang = ang(:,1:end) ;

if debug
  fprintf('filtering h-bonds... \n')
end

dislim = 4.0 ;
anglim = 120 ;

if debug
  fprintf('computing logical matrices ...\n')
end

idis = (dis <= dislim);
iang = (ang >= anglim);
ida = idis & iang;
idall = all(idis')';
iaall = all(iang')';
iall = idall & iaall; % iall(n)==true iff there's no H-bond broken in snapshot n
niall = ~iall; % niall(n)==true iff there's at least one H-bond broken in snapshot n
fhbtot = sum(iall)/nsnap*100;

if debug
  fprintf('Fraction of snapshots with no broken H-bond: %5.1f percent\n',fhbtot)
  fprintf('saving H-bonds logical matrices...\n')
end

filename = 'par_hbonds' ;
save([ path_to_HB '/' filename ], 'idis', 'idall', 'iang', 'ida', 'iaall', 'iall','niall','fhbtot');% fhb ;

end

