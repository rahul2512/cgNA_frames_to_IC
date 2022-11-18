function [uvw, pho, basepair, jframes] = base2p(nsnap,nbp,file,basepair_flag,junction_flag,debug)

%-------------------------------------------------------
% cgDNA function: shape = base2(nsnap,nbp,bframes,basepair_flag,junction_flag,debug)
%-------------------------------------------------------
% Calculates intra and inter basepair variables (eta, w, u, v) and base
% phosphate variables (etapW,wpW,etapC,wpC).
% Optionally returns the basepair frame vector and the junction frame vector.
%
% TODO: Complete documentation (e.g. describe embedding).
%
% Input:
%    nsnap           number of snapshots
%    nbp             number of base pairs
%    bframes         name of file with base frames
%    basepair_flag   flag to return basepair frames
%    junction_flag   flag to return junction frames
%    debug           flag to debug execution
%
% Output:
%    uvw             structure with base coordinates for each snapshot
%                    (see Note 1).
%    pho             structure with phosphate coordinates for each snapshot
%                    (see Note 1).
%    basepair        matrix with reference point and frame for each
%                    basepair, for each snapshot (see Note 2).
%    jframes         matrix with reference point and frame for each
%                    junction, for each snapshot (see Note 2).
%
% Note 1:
%
%   'uvw' is a struct array with fields:
%    - eta   intra-basepair rotational coords
%            (Buckle,Propeller,Opening) [Size: nsnap x nbp x 3]
%    - w     intra-basepair translational coords
%            (Shear,Stretch,Stagger) [Size: nsnap x nbp x 3]
%    - u     list of inter-basepair rotational coords
%            (Tilt,Roll,Twist) [Size: nsnap x nbp-1 x 3]
%    - v     list of inter-basepair translational coords
%            (Shift,Slide,Rise) [Size: nsnap x nbp-1 x 3]
%   'pho' is a struct array with fields:
%    - etapW list of base-Watson phosphate rotational coords
%            [Size: nsnap x nbp-1 x 3]
%    - wpW   list of base-Watson phosphate translational coords
%            [Size: nsnap x nbp-1 x 3]
%    - etapC list of base-Crick phosphate rotational coords
%            [Size: nsnap x nbp-1 x 3]
%    - wpC   list of base-Crick phosphate translational coords
%            [Size: nsnap x nbp-1 x 3]
%
% Note 2:
%
%   In 'basepair' and 'junction', the reference point and frame are
%   encoded in a 3 x 4 matrix, where the last column contains the
%   reference point.
%   Therefore, 'basepair' has size [nsnap x nbp x 3 x 4]; while
%   'jframes' has size [nsnap x nbp-1 x 3 x 4].
%
% If you find this code useful, please cite:
% TODO: add reference
%
%-------------------------------------------------------

if nargin < 4
  basepair_flag = 0 ;
  junction_flag = 0 ;
  debug = 0 ;
  
end


tic;
% open the file with base frames
fra_ID=fopen(file.bframes,'r');
pfra_ID=fopen(file.pframes,'r');

%% initialization for uvw coordinates
o = zeros(3,2,nbp);
a = zeros(3,3,2,nbp);
oab = zeros(3,nbp);
qab = zeros(3,3,nbp);
uvw.eta = zeros(nsnap,nbp,3);
uvw.w = zeros(nsnap,nbp,3);
uvw.u = zeros(nsnap,nbp-1,3);
uvw.v = zeros(nsnap,nbp-1,3);

%% initialization for pho coordinates
op = zeros(3,2,nbp-1);
ap = zeros(3,3,2,nbp-1);
pho.etapW = zeros(nsnap,nbp-1,3);
pho.wpW = zeros(nsnap,nbp-1,3);
pho.etapC = zeros(nsnap,nbp-1,3);
pho.wpC = zeros(nsnap,nbp-1,3);

basepair = 1;
jframes = 1;


if basepair_flag
  basepair = zeros(nsnap,nbp,3,4);
end
if junction_flag
  jframes = zeros(nsnap,nbp-1,3,4);
end

%%% buffering
bufsize = 1;
tet = zeros(bufsize,nbp-1,3);
tww = zeros(bufsize,nbp-1,3);
tuu = zeros(bufsize,nbp-1,3);
tvv = zeros(bufsize,nbp-1,3);

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
  end
  
  %%
  % 1.1) read base frames from file
  for k=1:2*nbp
    descr=fscanf(fra_ID,'%d',2);
    a(:,:,descr(1),descr(2))=reshape(fscanf(fra_ID,'%f',9), [3 3]);
    o(:,descr(1),descr(2))=fscanf(fra_ID,'%f',3);
  end
  
  % 1.2) read phosphate frames from file
  for k=1:2*(nbp-1)
    descr=fscanf(pfra_ID,'%d',2);
    if descr(2) == 1
      descr(1) = descr(1)-1;
    end
    ap(:,:,descr(2),descr(1))=reshape(fscanf(pfra_ID,'%f',9), [3 3]);
    op(:,descr(2),descr(1))=fscanf(pfra_ID,'%f',3);
  end
  
  
  %% 2) calculate bpframes and intra variables (eta, w - from b to a)
  for k=1:nbp
    
    oa = o(:,1,k);
    ob = o(:,2,k);
    oab(:,k) = 0.5*(oa+ob);
    ra = a(:,:,1,k);
    rb = a(:,:,2,k);
    rb(:,2) = -rb(:,2);
    rb(:,3) = -rb(:,3);
    aa=rb'*ra;
    [p, eta] = xcay(aa) ;
    rab=rb*p;
    qab(:,:,k)=rab;
    
    tww(tns,k,:) = rab'*(oa-ob);
    tet(tns,k,:) = eta * 5;
  end % k levels
  
  if basepair_flag
    basepair(ns, :, :, 1:3) = permute(qab, [3 1 2]);
    basepair(ns, :, :, 4) = oab';
  end
  
  %% 3) calculate jframes, inter variables (u, v) and phosphate variables (etapW,wpW,etapC,wpC)
  for k=1:nbp-1
    
    % 3.1) calculate jframes, inter variables (u, v)
    r1 = oab(:,k);
    r2 = oab(:,k+1);
    q1 = qab(:,:,k);
    q2 = qab(:,:,k+1);
    lam = q1'*q2;
    [pp, u] = xcay(lam);
    q12 = q1*pp;
    
    if junction_flag
      jframes(ns,k,:,1:3) = q12;
      jframes(ns,k,:,4) = 0.5*(r1+r2);
    end
    
    tuu(tns,k,:) = u * 5;
    tvv(tns,k,:) = q12'*(r2-r1);
    
    % 3.2) calculate phosphate variables (etapW,wpW,etapC,wpC)
    % 3.2.1) compute the Watson phosphate variables (etapW,wpW)
    oa = op(:,1,k);
    ob = o(:,1,k+1);
    ra = ap(:,:,1,k);
    rb = a(:,:,1,k+1);
    
    twp(tns,k,1,:) = rb'*(oa-ob);
    tep(tns,k,1,:) = pcay(rb'*ra)*5;
  
    % 3.2.2) compute the Crick phosphate variables (etapC,wpC)
    oa = op(:,2,k);
    ob = o(:,2,k);
    ra = ap(:,:,2,k);
    rb = a(:,:,2,k);
    
    twp(tns,k,2,:) = rb'*(oa-ob);
    tep(tns,k,2,:) = pcay(rb'*ra)*5;
    
  end % k levels
  
  if tns == bufsize
    I = ns-bufsize+1:ns;
    
    uvw.eta(I,:,:) = tet;
    uvw.w(I,:,:) = tww;
    uvw.u(I,:,:) = tuu;
    uvw.v(I,:,:) = tvv;
    
    pho.etapW(I,:,:) = squeeze(tep(:,:,1,:));
    pho.wpW(I,:,:) = squeeze(twp(:,:,1,:));
    pho.etapC(I,:,:) = squeeze(tep(:,:,2,:));
    pho.wpC(I,:,:) = squeeze(twp(:,:,2,:));
    
    tns = 1;
    
  else
    tns = tns + 1;
  end
  
end % ns snapshots
fclose(fra_ID);
fclose(pfra_ID);

if debug
  fprintf('Done. Total time %6.1f min\n',toc/60);
end

end

function [p, cay] = xcay(m)
angle = real(acos(0.5*(trace(m)-1)));
if ((pi - angle) < 0.05)
  %% angle is close to pi
  A = m+eye(3);
  [eigvect, eigval] = eig(A);
  [~, order] = sort(diag(eigval), 'descend');
  cay = eigvect(:,order(1));
  cay = cay/norm(cay);
  X = [0 -cay(3)  cay(2);
    cay(3) 0 -cay(1);
    -cay(2)  cay(1) 0];
  angle2 = 0.5*angle;
  p = cos(angle2) * eye(3) + (1-cos(angle2))*(cay*cay') + sin(angle2)*X;
  cay = cay*2*tan(angle2);
else
  p = sqrtm(m);
  anti = m - m';
  cay = [anti(3,2) anti(1,3) anti(2,1)];
  cay = cay*(2/(1+trace(m)));
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

