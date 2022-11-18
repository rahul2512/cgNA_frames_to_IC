function mkmer = group_kmer(oligomers, grouptype, debug)

%-------------------------------------------------------
% cgDNA function: mkmer = base4(oligomers,grouptype,debug)
%-------------------------------------------------------
% Groups the provided oligomer-level models by k-mer
% sequence. In particular, subsets of the oligomer-level
% ground-state configuration and stiffness are grouped according to
% their k-mer sequence.
%
% TODO: Complete documentation.
%
% Input:
%    oligomers       structure containing oligomer-level models
%                    for each oligomer (see Note 1).
%    grouptype       type of grouping to perform (see Note 2).
%    debug           flag to debug execution
%
% Output:
%    mkmer           structure with ground-state configuration and
%                    stiffness for each dinucleotide sequence
%                    (see Note 3).
%
% Note 1:
% 
%   'oligomers' is a struct array which must have 3 fields:
%    - w        ground-state configuration [Size: N x 1]
%    - stiff    ground-state stiffness [Size: N x N]
%    - sequence sequence of the oligomer [Size: nbp x 1]
%    where N is 12*nbp - 6.
% 
% Note 2:
%   Currently 2 types of groupings are supported:
%   - 'ab'      18x18 blocks including (intra, inter, intra) are
%               grouped by the dinucleotide sequence.
%   - 'abc'     6x6 blocks of intras are grouped according to the
%               trinucleotide sequence, including the previous and
%               next basepair.
%
% Note 3:
%
%   'mkmer' is a struct array which must have 3 fields:
%    - m    average stiffness matrix [Size: 18x18]
%    - c    average weighted shape vector [Size: 18x1]
%    - w    average shape vector [Size: 18x1]
% 
% If you find this code useful, please cite:
% TODO: add reference
% 
%-------------------------------------------------------

%% initialize
offb = 0; % distance of output block from start of sequence that
          % determines grouping. (e.g. nv for trimer-dependent
          % intra-parameters).  Can be < 0.
roff = 0; % how far ahead does the symmetric start
bInv = 0; % bInv = 1: invert before selecting; 0: use stiffness
expb = 0; % >0 : how far to expand BACK before inverting (must bInv); 
expf = 0; % >0 : how far to expand FRWD before inverting (must bInv);
exi = 2;  % skip basepairs at beginning of oligomers
exf = 2;  % skip basepairs at end of oligomers
P = diag([-1 1 1 -1 1 1]);
nv = 12;
nover=6;
nl =  6;
linP = {P P};

%% configure
bigPfun = @symm_rel;
%%% DIMER
% slen = 2;
% bsiz = nv+nover;
%%% MONOMER from TRIMER
% slen = 3;
% bsiz = nover;
% offb = nv;
% expb = nv;
% expf = nv;
% exi=1; % allow second: offb > 0
%%% 30x30 MONOMER (xAx)
% bsiz=2*nv+nover;
% slen=1;
% offb = -nv;
%%% 18x18 DIMER (AB)
% bsiz=nv+nover;
% slen=2;
% offb = 0;
%%% 42x42 DIMER (xABx)
% bsiz=3*nv+nover;
% slen=2;
% offb = -nv;
%%% 30x30 TRIMER (ABC)
% slen=3;
%%% P 42x42 MONOMER (pAp)
% offb = -(nv-nover);
% roff = 0;
% bsiz=2*nv-nover;
% slen=1;
%%% P 30x30 DIMER (AB)
% bsiz=nv+nover;
% slen = 2;
% offb = 0;
% roff = 1;
%%% P 6x6 TRIMER (ABC)
% bsiz=nover;
% slen = 3;
% offb = nv;
% roff = 2;
%%% group types
if strcmpi(grouptype, 'ab') == 1
    slen = 2;
    bsiz = nv+nover;
    roff = 1;
elseif strcmpi(grouptype, 'abc') == 1
    slen = 3;
    bsiz = nover;
    offb = nv;
    exi=1; % allow second: offb > 0
else
    fprintf('Unknown group type %s...Exiting.\n', grouptype);
    return;
end
bigP = bigPfun(linP, P, bsiz, nv, nover);

%% define k-mers
if slen == 1
    Su = ['A'; 'G']; % 2 -> 4
elseif slen == 2
    Su = ['AT'; 'GC'; 'TA'; 'CG';    'GT'; 'TG'; 'AG'; 'GA'; 'AA'; 'GG';]; ...
    % 10 -> 16
elseif slen == 3
    Su = ['AAC'; 'GAT'; 'AAT'; 'GAC'; 'TAA'; 'CAG'; 'CAA'; 'TAG'; 'AAA'; ...
          'AAG'; 'GAA'; 'GAG'; 'CAC'; 'CAT'; 'TAC'; 'TAT'; ...
          'AGC'; 'GGT'; 'AGT'; 'GGC'; 'TGA'; 'CGG'; 'CGA'; 'TGG'; 'AGA'; ...
          'AGG'; 'GGA'; 'GGG'; 'CGC'; 'CGT'; 'TGC'; 'TGT']; % 32 -> 64
else
    fprintf('Sequence dependence of %d not implemented...\n',slen);
    return
end
% complete set with complements
S=Su;
for i=1:size(Su,1)
    c = wcc(Su(i,:),-1);
    if fi(Su, c) < 0
        S=[S; c];
    end
end

%% initialize (2)
for i = 1:size(S,1)
    mkmer(i).S = S(i,:);
    mkmer(i).n = 0;
end 

tic;
if debug
    fprintf('base4: using %d oligomers\n', numel(oligomers));
end
for ioli = 1:numel(oligomers)
    seq = oligomers(ioli).sequence;
    nbp = numel(seq);
    stiff = oligomers(ioli).stiff;
    w = oligomers(ioli).w;

    if bInv
        cov = inv(stiff);
    end
    Kw = stiff*w;
    
    for i = 1+exi:(nbp-exf)   %-- loop over basepairs, removing (or not) the ends     
        k = (i-1)*nv +1;   
        idx = k+offb:k+offb+bsiz-1;
        %--- if idx extends over stiff, pad with zeros
        % for offb<0 && exi<?, or offb+bsiz>nover and exf>?
        pad = [0 0];
        if idx(1) < 1                     % left-pad 
            i0 = find(idx>0,1);
            pad(1) = i0-1;
            idx = idx(i0:end);
        end
        if idx(end) > size(stiff,1)
            tsiz = size(stiff,1);
            iF = find(idx<tsiz,1,'last') + 1;
            pad(2) = numel(idx)-iF;
            idx = idx(1:iF);
        end

        if debug
            fprintf('%3d) %s[%2d]: [%3d:%3d], padding: %2d, %2d\n',...
                    iseq, seq, i, idx([1 end]), pad);
        end
        
        if bInv == 1
            bigidx = idx(1)-expb:idx(end)+expf;
            M = inv(cov(bigidx,bigidx));
            M = M(expb+1:expb+bsiz, expb+1:expb+bsiz);
        else
            M = stiff(idx,idx);
        end
        M = tril(M) + tril(M)' - diag(diag(M));  %in case there were small errors inverting 
        z = Kw(idx);
        c = z;
        ww = w(idx);
        
        %% pad
        if any(pad>0)
            tsiz = size(M,1)+sum(pad);
            MS = zeros(tsiz);
            zs = zeros(tsiz, 1);
            ws = zeros(tsiz, 1);
            I = [1+pad(1):tsiz-pad(2)];
            MS(I,I) = M; M = MS;
            zs(I) = z;   z = zs; c = z;
            ws(I) = ww;  ww= ws;
        end
        
        %% changing reference strand, e.g. transforming a GT block to an AC block 
        MT = bigP*M*bigP;
        MT = triu(MT) + triu(MT)' - diag(diag(MT));
        zT = bigP*z;
        cT  = zT; 
        wwT  = bigP*ww;
        
        %% only 1 match possible in S
        mkmer = madd(mkmer, fsi(mkmer, seq(i:i+slen-1)), M, c, ww, seq, i, 0, i);
        i2 = numel(seq) - i + 1 - roff;
        mkmer = madd(mkmer, fsi(mkmer, wcc(seq(i:i+slen-1),-1)), MT, cT, wwT, seq, i, 1, i2);
        
   end % i basepairs
end % ioli oligomers

%% compute averages
for l = 1:numel(mkmer)
    mkmer(l).m =   reshape(mean(mkmer(l).all,1),bsiz,bsiz).*sten;
    mkmer(l).c =   reshape(mean(mkmer(l).call,1),bsiz,1).*stenv;
    mkmer(l).w =   reshape(mean(mkmer(l).wall,1),bsiz,1).*stenv;
end
end

function bigP = symm_rel(linP, P, bsiz, nv, nover)
% 
% the default bigP with nover
%
    tlinP = repmat(linP, 1, (bsiz-nover)/size(P,1)/numel(linP));
    totP = {tlinP{:} P};     % XXX bsiz must have an extra 'nover'
    bigP = blkantidiag(totP{:});
end

