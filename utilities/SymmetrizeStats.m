function olig = SymmetrizeStats(olig)


P = diag([-1 1 1, -1 1 1]);

for j=1:length(olig)
    w = olig(j).shape;     % shape
    C = inv(olig(j).s1b);  % covariance
    S = C + w*w';               % second moment
    nbp=olig(j).nbp;
    
    N=24*nbp-18;
    E=zeros(N);
    for i=1:nbp
        b = 24*(i-1);
        E(b+1:b+6,N-b-5:N-b) = P;
        if i<nbp
            E(b+7:b+12,N-b-11:N-b-6) = eye(6);
            E(b+13:b+18,N-b-17:N-b-12) = P;
            E(b+19:b+24,N-b-23:N-b-18) = eye(6);
        end
    end
    
    

    
    wsym = 0.5*(w + E*w);
    Ssym = 0.5*(S + E*S*E);
    Csym = Ssym - wsym*wsym';
    
    olig(j).shape_sym = wsym;
    olig(j).s1b_sym = inv(Csym);
    
    [~, stiff_me] = completeMaxEntropy(Csym,cornerset(nbp,1), 1);
    olig(j).stiff_me_sym = stiff_me;
end

end
