function nots = wcc(s,varargin)
%
% nots = wcc(s)
% 
% Watson-Crick complementary sequence
%
    W=['ACGT_'];
    C=['TGCA_'];
    nots=s;
    % assert numel(s) == max(size(s))
    
    dir = 0;
    if nargin > 1
        dir = varargin{1};
    end
        
    for i=1:numel(s)
        nots(i) = C(find(W==s(i)));
    end
    
    if dir < 0
        nots = nots(end:-1:1);
    end
end

