function i = fi(A, s)
%
%
%
    if ~iscell(A)
        A = cellstr(A);
    end
    if numel(A) == 1
        fprintf('sfsi Warning: only %d possible sequences!\n', numel(A));
    end
    i=-1;
    for j = 1:numel(A)
        if(strcmpi(A{j},s))
            i=j; 
            return;
        end  
    end    
end
