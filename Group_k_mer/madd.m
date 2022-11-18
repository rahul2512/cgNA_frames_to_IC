function mdimer = madd(mdimer, l, M, c, w, seq, i, w2c, j)
    n = mdimer(l).n + 1;
    mdimer(l).n  = n;
    mdimer(l).all(n,:,:) = M;
    mdimer(l).call(n,:) = c;
    mdimer(l).wall(n,:) = w;
    mdimer(l).w2c(n) = w2c;
    mdimer(l).origin.seq{n} = seq;
    mdimer(l).origin.lev(n) = i;
    mdimer(l).origin.len(n) = numel(seq);
    mdimer(l).origin.lev2(n)= j;
end
