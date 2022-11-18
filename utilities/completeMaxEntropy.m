function [ B, BInv ] = completeMaxEntropy( C, cornerSet, inverseOnly, ...
    zeroMargin )
% File: completeMaxEntropy.m
%
% Copyright (C) 2013 Jaroslaw Glowacki <jarek.glowacki@gmail.com> 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% MAXENTFORINVSTENCIL Computes the matrix B that is the maximum entropy
%    completion of the input matrix C and/or its inverse. C is assumed to
%    have non-zero elements only within a sparsity pattern defined by the
%    corner set given as the second input argument. The values of elements
%    of C outside the sparsity pattern are ignored.
%
%    The inverse of B is computed and returned as a second output argument,
%    if a second ouput argument is requested.
%
%    If a third input argument evaluates to true ONLY the inverse
%    computation is performed (the first output argument will be the
%    ORIGINAL matrix C, not the maximum entropy completion)
%
%    If a fourth input argument is given and is positive, at each step s
%    the SPD condition is checked for the input matrix and every leading
%    diagonal sub-block B_s. The value of the argument is used as the
%    tolerance margin for 0. Numbers with absolute value lower than the
%    margin will be treated as 0. This fourth argument is intended mainly
%    for testing purposes.
%
%    Eamples:
%    * compute ONLY the maximum entropy completion
%       B = completeMaxEntropy(C, s)
%    * computes the completion AND its inverse
%       [B, BInv] = completeMaxEntropy(C, s)
%    * computes ONLY the inverse of the maximum entropy completion directly
%       [~, BInv] = completeMaxEntropy(C, s, true)
%
%    Further eamples:
%    * compute the completion ONLY, no SPD checks
%       B = completeMaxEntropy(C, s)
%       B = completeMaxEntropy(C, s, false)
%       B = completeMaxEntropy(C, s, false, -1.0)
%    * computes the completion ONLY, and asks minimum eigenval of sub-block
%    after ech step to be > 0.001
%       B = completeMaxEntropy(C, s, false, 0.001)
%    * computes the completion and the inverse, without SPD checks
%       [B, BInv] = completeMaxEntropy(C, s)
%       [B, BInv] = completeMaxEntropy(C, s, false)
%       [B, BInv] = completeMaxEntropy(C, s, false, -1.0)
%    * computes the completion AND the inverse, without SPD checks,
%    value of B is ignored though computed
%       [~, BInv] = completeMaxEntropy(C, s)
%       [~, BInv] = completeMaxEntropy(C, s, false)
%    * computes ONLY the inverse without SPD checks, the returned B = C
%       [B, BInv] = completeMaxEntropy(C, s, true)
%    * computes ONLY the inverse without SPD checks, the returned value of
%    the first argument (= C) is ignored
%       [~, BInv] = completeMaxEntropy(C, s, true)
%
%    Note that
%       B = completeMaxEntropy(C, s, true)
%    returns immediatelly the value B = C (because the value of the inverse
%    is the second output argument, which is ignored in this particular
%    call and because of the third argument that is true no comutations are
%    performed here)

    % If the fourth argument is not given is is assumend negative
    % which means that the SPD test is not to be performed
    if nargin < 4
        zeroMargin = -1.0;
        % If the third argument is not given it is assumed to be false
        % ans so the "forward" entropy completion is computed
        if nargin < 3
            inverseOnly = false;
        end
    end
    
    [k, colsS] = size(cornerSet);
    [n, colsC] = size(C);
      
    % Check input data
    assert(n == colsC, 'Input matrix not square');
    assert(k <= n && colsS == 2, ['The corner set has to be a k x 2 ' ...
        'matrix with k less or equal to the dim of the input matrix']); 
    assert(cornerSet(1, 1) == 1 && cornerSet(k, 2) == n, ...
        ['The stencil data has to start with (i, j) = (1, *) ' ...
        'and end with (i, j) = (*, n) with n - dim of the input matrix']);
    
    for s = 1:k - 1
        cornerSetCorr = cornerSet(s, 1) < cornerSet(s + 1, 1) ...
           && cornerSet(s + 1, 1) <= cornerSet(s, 2) ...
           && cornerSet(s, 2) < cornerSet(s + 1, 2);
        assert(cornerSetCorr, ...
            'Corner set conditions not satisfied for corner %d', s);
    end
    
    % If SPD test requested
    if zeroMargin > 0.0
        symError = norm(C - C');
        assert(symError < zeroMargin, ['Input matrix not symmetric\n', ...
            '(||C - C''|| = %f > %f)'], symError, zeroMargin);
        minEigv = min(eig(C));
        assert(minEigv > zeroMargin, ['Input matrix not pos. def.\n', ...
            '(min eigval = %f < %f)'], minEigv, zeroMargin);
    end

    % The sequence of sub-blocks B_[s] will be stored directly in the
    % result matrix B initiated with the original content of C
    B = C;
    
    % If the second OUTPUT argument requested compute the inverse
    if(nargout > 1)
        BInv = zeros(n);
        % Start with the inverse of the first of the overlapping blocks
        % inv(D[1])
        j1 = cornerSet(1,2);
        BInv(1:j1, 1:j1) = inv(C(1:j1, 1:j1));
    elseif inverseOnly
        warning(['Only inverse required, but the second output ' ...
            'argumentis ignored.', sprintf('\n'), ...
            'No computations will be done']);
    end
    
    % Iterate through the corners from the top left one
    for s = 2:k
        % j_(s-1)
        jsm1 = cornerSet(s - 1, 2);
        % i_s
        is = cornerSet(s, 1);
        % j_s
        js = cornerSet(s, 2);
        
        % D[s]_1,1 (the overlap)
        Ds_11 = B(is:jsm1, is:jsm1);
        nDs_11 = jsm1 - is  + 1;

        % If the second OUTPUT argument requested compute the inverse
        if(nargout > 1)
            % The entire D_[s] (the overlapping block)
            Ds = B(is:js, is:js);
            
            % Do explicit Cholesky factorization as we need both the D_[s]
            % one and the part corresponding to D_[s]_1,1 (the overlap)
            L_Ds = chol(Ds, 'lower');
            LT_Ds = L_Ds';
            
            % Matlab should recognize z triangular matrix
            % Forward substitution for D_[s] and D_[s]_11 is the same
            Y_Ds = L_Ds \ eye(js - is + 1);
            DsInv = LT_Ds \ Y_Ds;
            Ds_11Inv = LT_Ds(1:nDs_11, 1:nDs_11) \ Y_Ds(1:nDs_11, 1:nDs_11);
            
            BInv(is:js, is:js) = BInv(is:js, is:js) + DsInv;
            BInv(is:jsm1, is:jsm1) = BInv(is:jsm1, is:jsm1) - Ds_11Inv; 
        end
        
        % If the "forward" maximum entropy completion is also requested
        if ~inverseOnly
            % D[s]_1,2
            Ds_12 = B(is:jsm1, jsm1 + 1:js);
        
            % B[s-1]_1,2
            Bsm1_12 = B(1:is - 1, is:jsm1);
            
            % Compute Cholesky decomposition of D_[s]_1,1 (the overlap)
            % or use the already computed value if it exists
            % (if the inverse computations have been done the partial results
            % can be reused here)
            if exist('L_Ds')
                L_Ds_11 = L_Ds(1:nDs_11, 1:nDs_11);
            else
                L_Ds_11 = chol(Ds_11, 'lower');
            end
            
            % Omega_s
            % Matlab should recognize triangular matrices
            Omega = L_Ds_11' \ (L_Ds_11 \ Ds_12);
        
            % Compute the elements of the completion
            % The block that replaces entries outside of the sparcity
            % pattern (the ones within E_s) in of B_[s] above the diagonal
            B(1:is - 1, jsm1 + 1:js) = Bsm1_12 * Omega;
            % Replace entries outside of the sparcity pattern
            % (the ones within E_s) in of B_[s] below the diagonal
            % with transpose of the above
            B(jsm1 + 1:js, 1:is - 1) = B(1:is - 1, jsm1 + 1:js)';
        end
       
        % See if the sub-block B_[s] that is currently ready is SPD
        % (only if such a check is requested)
        if zeroMargin > 0.0
            symError = norm(B(1:js, 1:js) - B(1:js, 1:js)');
            assert(symError < zeroMargin, ...
                ['The leading diagonal subblock after step %d ', ...
                'not symmetric\n(||B_[%d] - B_[%d]''|| = %f > %f)'], ...
                s, s, s, symError, zeroMargin);
            minEigv = min(eig(B(1:js, 1:js)));
            assert(minEigv > zeroMargin, ...
                ['The leading diagonal subblock after step %d ', ...
                'not pos. def.\n(min eigval = %f < %f)'], ...
                s, minEigv, zeroMargin);
        end
    end
end
