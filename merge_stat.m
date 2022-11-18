function data = merge_stat(name,nbr_files,nbr_seq)

addpath(genpath('/Users/lcvmm2/Dropbox/cgDNAolig_ver2'))

nhb = zeros(nbr_seq,1);

for i = 1:nbr_files
    
    load([ name '.' num2str(i) '.mat']);
    data = olig; 
    
    for s = 1:nbr_seq
        
        m_tmp = olig(s).shape ;
        cv_tmp = inv(olig(s).s1b);
        nhb_tmp = olig(s).nsnap(1);
        
        S(:,:,i,s) = cv_tmp*(nhb_tmp-1) + nhb_tmp*(m_tmp*m_tmp') ;
        m(:,i,s)  = m_tmp*nhb_tmp;
        
        nhb(s) = nhb(s) + nhb_tmp ;
        
    end
    
    
end % end loop over fra files



for s = 1:nbr_seq
    Stot = sum(S(:,:,:,s),3)/(nhb(s)-1) ;
    mu   = sum(m(:,:,s),2)/nhb(s) ;
    cv   = Stot - nhb(s)/(nhb(s)-1) * (mu*mu') ;
    
    data(s).shape = mu ; 
    data(s).nsnap = nhb(s);
    data(s).s1b = inv(cv);
    [~, stiff_me] = completeMaxEntropy(cv,cornerset(data(s).nbp,1), 1);
    data(s).stiff_me = stiff_me;
    
end


end
