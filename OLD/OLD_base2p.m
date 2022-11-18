function pho = base2p(nsnap,nbp,file,jframes,debug)

%-------------------------------------------------------
% cgDNA function: shape = base2p(nsnap,nbp,bframes,basepair_flag,junction_flag,debug)
%-------------------------------------------------------
% Calculates phosphate variables (eta_p, w_p) for each strand
% (Watson, Crick).
%
% TODO: Complete documentation (e.g. describe embedding).
%
% Input:
%    nsnap           number of snapshots
%    nbp             number of base pairs
%    pframes         name of file with phosphate frames
%    jframes         structure with junction frames (see Note 2)
%    debug           flag to debug execution
%
% Output:
%    shape           structure with coordinates for each snapshot
%                    (see Note 1).
%
% Note 1:
% 
%   'shape' is a single-entry struct array with fields:
%    - etapW  phosphate rotational coords for the Watson strand
%             [Size: nsnap x nbp-1 x 3]
%    - wpW    phosphate translational coords for the Watson strand
%             [Size: nsnap x nbp-1 x 3]
%    - etapC  phosphate rotational coords for the Crick strand
%             [Size: nsnap x nbp-1 x 3]
%    - wpC    phosphate translational coords for the Crick strand
%             [Size: nsnap x nbp-1 x 3]
%
% Note 2:
%
%   In 'jframes', the reference point and frame are
%   encoded in a 3 x 4 matrix, where the last column contains the
%   reference point: 'jframes' has size [nsnap x nbp-1 x 3 x 4].
% 
% If you find this code useful, please cite:
% TODO: add reference
% 
%-------------------------------------------------------

tic;
% open the file with phosphate frames
pfra_ID=fopen(file.pframes,'r');

%% initialize
op = zeros(3,2,nbp-1); 
ap = zeros(3,3,2,nbp-1);
pho.etapW = zeros(nsnap,nbp-1,3); 
pho.wpW = zeros(nsnap,nbp-1,3);
pho.etapC = zeros(nsnap,nbp-1,3); 
pho.wpC = zeros(nsnap,nbp-1,3);

%%% buffering
bufsize = 100;
tep = zeros(bufsize,nbp-1,2,3);
twp = zeros(bufsize,nbp-1,2,3);

% loop over snapshots
partime = 0;
tns = 1;
for ns=1:nsnap
    %% progress and debug
    if debug && mod(ns,1000)==0
        t = toc;
        fprintf('%d snaps, total time%6.1f min\n',ns,t/60)
        delta = t-partime;
        partime = t;
        % detail = 0;
        % if delta > 120
        %     detail = 1;
        % end
    end
    
    %% 1) read phosphate frames from file
    for k=1:2*(nbp-1)
        descr=fscanf(pfra_ID,'%d',2);
        if descr(2) == 1
            descr(1) = descr(1)-1;
        end
        ap(:,:,descr(2),descr(1))=reshape(fscanf(pfra_ID,'%f',9), [3 3])/1000.0;
        op(:,descr(2),descr(1))=fscanf(pfra_ID,'%f',3)/1000.0;
    end
    %if detail; fprintf('|>>|%2d|%6.4f|\n',0,toc); end

    %% 2) calculate phosphate frames and variables (etap, wp, for each strand)
    for k=1:nbp-1
        for l=[1 2]
            %if detail; fprintf('|>>|%2d|%2d|%6.4f|\n',k,l,toc); end
            ob = squeeze(jframes(ns,k,:,4));
            oa = op(:,l,k);
            rb = squeeze(jframes(ns,k,:,1:3));
            ra = squeeze(ap(:,:,l,k)); 
            if l==2
                rb(:,2) = -rb(:,2);
                rb(:,3) = -rb(:,3);
            end
            %if detail; fprintf('|>>|%2d|%2d|%6.4f|\n',k,l,toc); end
            twp(tns,k,l,:) = rb'*(oa-ob);
            tep(tns,k,l,:) = pcay(rb'*ra) * 5;
        end % l strands
    end % k levels
    
    if tns == bufsize
        %if detail; fprintf('|>>|FF|%6.4f|\n',toc); end
        I = ns-bufsize+1:ns;
        pho.etapW(I,:,:) = squeeze(tep(:,:,1,:));
        pho.wpW(I,:,:) = squeeze(twp(:,:,1,:));
        pho.etapC(I,:,:) = squeeze(tep(:,:,2,:));
        pho.wpC(I,:,:) = squeeze(twp(:,:,2,:));
        tns = 1;
        %if detail; fprintf('|>>|FF|%6.4f|\n',toc); end
    else
        tns = tns + 1;
    end
    % if detail
    %     detail=0;
    % end
end % ns snapshots
fclose(pfra_ID);
if debug 
    fprintf('Done. Total time %6.1f min\n',toc/60);
end

end

function cay = pcay(m)
angle = real(acos(0.5*(trace(m)-1)));
if ((pi - angle) < 0.05) 
    %% angle is close to pi
    A = m+eye(3);
    [eigvect, eigval] = eig(A);
    [~, order] = sort(diag(eigval), 'descend'); 
    cay = eigvect(:,order(1));
    cay = cay/norm(cay);
    cay = cay*2*tan(angle); 
else  
    anti = m - m';
    cay = [anti(3,2) anti(1,3) anti(2,1)];
    cay = cay*(2/(1+trace(m)));
end
end
