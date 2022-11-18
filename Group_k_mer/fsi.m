function i = fsi(struc, s, varargin)
    % fieldname = 'S';
    % if nargin > 2
    %     fieldname = varargin{1};
    % end
    % i = fi({getfield(struc,fieldname)}, s);
    i = fi({struc.S}, s);
end
