function M = blkantidiag(varargin)
  T = blkdiag(varargin{:});
  l = size(T,1)/nargin;
  M = blkproc(fliplr(T),[l l],@fliplr);
end
