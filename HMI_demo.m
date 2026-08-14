function results = HMI_demo(matFile)
% INPUT
%   matFile : MAT file containing one scalar structure named "model".
%
% EXAMPLES
%   results = HMI_demo('mesh_cube.mat');
%   results = HMI_demo('mesh_sphere.mat');
%   results = HMI_demo('mesh_torus.mat');
%   results = HMI_demo('mesh_spiky_virus.mat');
%   results = HMI_demo('mesh_nanostar.mat');
%
% A user can generate an arbitrary conforming tetrahedral mesh
%
% REQUIRED MODEL DATA CONTRACT
%   model.schemaVersion
%   model.name
%   model.nodes                    [nNode x 3]
%   model.elements.tet4            [nTet x 4]
%   model.region.solidMask         [nTet x 1] logical
%   model.region.fluidMask         [nTet x 1] logical
%   model.region.portIndex         [nTet x 1], 0 outside ports
%   model.boundary.outerFaces      [nFace x 3]
%   model.scale.a                  characteristic length
%
% OPTIONAL
%   model.ports(k).elements
%   model.ports(k).nominalCenter
%   model.meta
%
% NUMERICAL METHOD
%   - Conforming 3D P2/P1 mixed finite-element formulation.
%   - P2 velocity obtained by edge enrichment of the input P1 tetrahedral mesh.
%   - Separate solid and fluid pressures with one pressure gauge.
%   - Kelvin--Voigt solid, Stokes fluid, and localized force/torque ports.

close all;
clc;

matFile = char(matFile);
if exist(matFile,'file') ~= 2
    error('HMIDemo:Input','Input MAT file does not exist:\n%s',matFile);
end

S = load(matFile,'model');
if ~isfield(S,'model')
    error('HMIDemo:Input','The MAT file must contain a structure named model.');
end
model = S.model;
demo_validate_model(model);

model.elements.tet4 = demo_orient_tet4_local( ...
    model.nodes,model.elements.tet4);

[inputDir,stem,~] = fileparts(matFile);
if isempty(inputDir)
inputDir = pwd;
end
outputDir = fullfile(inputDir,[stem '_HMI_results']);
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

% =====================================================================
% Demonstration material parameters and frequency window
% =====================================================================
cfg = struct();
cfg.a = model.scale.a;
cfg.rhoS = 1.0;
cfg.G = 1.0;
cfg.LambdaF = 0.50;
cfg.LambdaS = 0.10;
cfg.tRef = cfg.a*sqrt(cfg.rhoS/cfg.G);
cfg.etaF = cfg.LambdaF*cfg.G*cfg.tRef;
cfg.etaS = cfg.LambdaS*cfg.G*cfg.tRef;
cfg.omegaHat = [5e-5 1e-4 2e-4 4e-4];
cfg.rankRelTol = 1e-9;
cfg.pngResolution = 300;

fprintf('Input mesh: %s\n',matFile);
fprintf('Geometry name: %s\n',char(string(model.name)));
fprintf('Output folder: %s\n\n',outputDir);

% =====================================================================
% 1. Build P2 velocity topology and geometry data
% =====================================================================
fprintf('[1/7] Building quadratic velocity topology...\n');
p2 = demo_enrich_p2(model.nodes,model.elements.tet4,model.boundary.outerFaces);

fprintf('  P1 nodes               %d\n',size(model.nodes,1));
fprintf('  P2 velocity nodes      %d\n',size(p2.nodes,1));
fprintf('  tetrahedra             %d\n',size(model.elements.tet4,1));
fprintf('  fixed outer P2 nodes   %d\n',nnz(p2.fixedNode));

% =====================================================================
% 2. Solid connectivity and port resultants
% =====================================================================
fprintf('[2/7] Detecting solid components and port resultants...\n');
geom = demo_geometry_data(model);

fprintf('  solid components C     %d\n',geom.nComponent);
fprintf('  ports P                %d\n',geom.nPort);
fprintf('  load coordinates n     %d\n',6*geom.nPort);

% =====================================================================
% 3. Assemble mixed FEM operators and generalized loads
% =====================================================================
fprintf('[3/7] Assembling 3D P2/P1 mixed FEM operators...\n');
fem = demo_assemble_fem(model,p2,geom,cfg);

fprintf('  free velocity dofs     %d\n',fem.nV);
fprintf('  pressure dofs          %d\n',fem.nP);
fprintf('  mixed+gauge size       %d\n',fem.systemSize);

% =====================================================================
% 4. Observable and hidden dimensions
% =====================================================================
fprintf('[4/7] Computing static, harmonic, and hidden dimensions...\n');
rankData = demo_rank_structure(fem,geom,cfg);

fprintf('  q0 static observable   %d\n',rankData.q0);
fprintf('  qw harmonic observable %d\n',rankData.qw);
fprintf('  qh hidden              %d\n',rankData.qh);
fprintf('  projected-load rank    %d\n',rankData.qw);
fprintf('  rigid-resultant rank   %d\n',rankData.q0);

% =====================================================================
% 5. Frequency solves
% =====================================================================
fprintf('[5/7] Solving harmonic mobility at %d frequencies...\n',numel(cfg.omegaHat));
nW = numel(cfg.omegaHat);
nLoad = size(fem.Fphysical,2);
mobility = complex(zeros(nLoad,nLoad,nW));
singularValues = zeros(nLoad,nW);
linearResidual = zeros(nW,1);
divergenceResidual = zeros(nW,1);
gaugeResidual = zeros(nW,1);
reciprocityResidual = zeros(nW,1);
passivityMinEigenvalue = zeros(nW,1);
solveSeconds = zeros(nW,1);

for iw = 1:nW
    fprintf('  omega_hat = %.4e\n',cfg.omegaHat(iw));
    tSolve = tic;
    sol = demo_dynamic_solve(fem,cfg,cfg.omegaHat(iw));
    solveSeconds(iw) = toc(tSolve);

    mobility(:,:,iw) = sol.Mhat;
    singularValues(:,iw) = svd(sol.Mhat);
    linearResidual(iw) = sol.linearResidual;
    divergenceResidual(iw) = sol.divergenceResidual;
    gaugeResidual(iw) = sol.gaugeResidual;
    reciprocityResidual(iw) = sol.reciprocityResidual;
    passivityMinEigenvalue(iw) = sol.passivityMinEigenvalue;

    fprintf('    rank Mhat           %d\n', ...
        demo_numerical_rank(sol.Mhat,cfg.rankRelTol));
    fprintf('    linear residual     %.3e\n',linearResidual(iw));
    fprintf('    divergence residual %.3e\n',divergenceResidual(iw));
    fprintf('    reciprocity         %.3e\n',reciprocityResidual(iw));
    fprintf('    solve time          %.1f s\n',solveSeconds(iw));
end

% =====================================================================
% 6. Hidden activation and normalized opening
% =====================================================================
fprintf('[6/7] Extracting hidden activation and opening ratios...\n');
qh = rankData.qh;
activation = zeros(qh,1);
activationMatrix = zeros(qh);
hiddenSigma = zeros(qh,nW);
collapseRatio = zeros(qh,nW);
fittedOrder = zeros(qh,1);
projectedActivationByFrequency = zeros(qh,nW);

if qh > 0
    QH = rankData.Qhidden;
    Mh0 = QH.'*mobility(:,:,1)*QH;
    A0 = real(Mh0/(1i*cfg.omegaHat(1)));
    A0 = 0.5*(A0+A0.');
    activationMatrix = A0;
    activation = sort(real(eig(A0)),'descend');
    if any(~isfinite(activation)) || any(activation<=0)
        error('HMIDemo:ActivationPositivity', ...
            ['The lowest-frequency hidden activation estimate is not ', ...
             'positive definite. Refine the mesh or lower the demo ', ...
             'frequency window before interpreting the result.']);
    end

    for iw=1:nW
        Mh = QH.'*mobility(:,:,iw)*QH;
        Aw = real(Mh/(1i*cfg.omegaHat(iw)));
        Aw = 0.5*(Aw+Aw.');
        projectedActivationByFrequency(:,iw) = ...
            sort(real(eig(Aw)),'descend');

        sv = singularValues(:,iw);
        hiddenSigma(:,iw) = sv(rankData.q0+(1:qh));
        collapseRatio(:,iw) = hiddenSigma(:,iw) ./ ...
            (abs(cfg.omegaHat(iw))*activation);
    end

    for j=1:qh
        p = polyfit(log(cfg.omegaHat(:)),log(hiddenSigma(j,:).'),1);
        fittedOrder(j)=p(1);
    end

    fprintf('  activation eigenvalues\n    ');
    fprintf('%.8g ',activation);
    fprintf('\n');
    fprintf('  fitted opening orders\n    ');
    fprintf('%.6f ',fittedOrder);
    fprintf('\n');
    fprintf('  max |R-1| at lowest omega  %.3e\n', ...
        max(abs(collapseRatio(:,1)-1)));
else
    QH = zeros(nLoad,0);
    fprintf('  No hidden directions were detected for this loading family.\n');
end

% =====================================================================
% 7. Save text, MAT, and three separate figures
% =====================================================================
fprintf('[7/7] Writing results and separate figures...\n');
figureFiles = demo_write_figures(outputDir,model,geom,rankData,activation, ...
    cfg.omegaHat,collapseRatio,cfg);

textFile = fullfile(outputDir,'HMI_results.txt');
demo_write_text_results(textFile,matFile,model,p2,geom,fem,cfg,rankData, ...
    activation,projectedActivationByFrequency,hiddenSigma,collapseRatio, ...
    fittedOrder,singularValues,linearResidual,divergenceResidual, ...
    gaugeResidual,reciprocityResidual,passivityMinEigenvalue,solveSeconds, ...
    figureFiles);

results = struct();
results.inputFile = matFile;
results.outputDir = outputDir;
results.config = cfg;
results.q0 = rankData.q0;
results.qw = rankData.qw;
results.qOmega = rankData.qw;
results.qh = rankData.qh;
results.B = rankData.B;
results.Qhidden = QH;
results.activationMatrix = activationMatrix;
results.activation = activation;
results.omegaHat = cfg.omegaHat;
results.mobility = mobility;
results.singularValues = singularValues;
results.hiddenSigma = hiddenSigma;
results.collapseRatio = collapseRatio;
results.projectedActivationByFrequency = projectedActivationByFrequency;
results.fittedOrder = fittedOrder;
results.linearResidual = linearResidual;
results.divergenceResidual = divergenceResidual;
results.gaugeResidual = gaugeResidual;
results.reciprocityResidual = reciprocityResidual;
results.passivityMinEigenvalue = passivityMinEigenvalue;
results.solveSeconds = solveSeconds;
results.textFile = textFile;
results.figureFiles = figureFiles;

save(fullfile(outputDir,'hidden_mechanical_information_results.mat'), ...
    'results','-v7.3');

fprintf('\nCompleted.\n');
fprintf('Text summary:\n  %s\n',textFile);
fprintf('Figures:\n');
fprintf('  %s\n',figureFiles.geometryPng);
fprintf('  %s\n',figureFiles.activationPng);
fprintf('  %s\n',figureFiles.openingPng);
end

% =====================================================================
% Input validation
% =====================================================================
function demo_validate_model(model)
required = {'schemaVersion','name','nodes','elements','region','boundary','scale'};
for k=1:numel(required)
    if ~isfield(model,required{k})
        error('HMIDemo:ModelField','Missing model.%s.',required{k});
    end
end
if ~(isscalar(model.schemaVersion) && model.schemaVersion==1)
    error('HMIDemo:Schema', ...
        'This solver currently supports model.schemaVersion = 1.');
end
if ~isfield(model.elements,'tet4')
    error('HMIDemo:ModelField','Missing model.elements.tet4.');
end
for f={'solidMask','fluidMask','portIndex'}
    if ~isfield(model.region,f{1})
        error('HMIDemo:ModelField','Missing model.region.%s.',f{1});
    end
end
if ~isfield(model.boundary,'outerFaces')
    error('HMIDemo:ModelField','Missing model.boundary.outerFaces.');
end
if ~isfield(model.scale,'a') || ~(isscalar(model.scale.a) && model.scale.a>0)
    error('HMIDemo:ModelField','model.scale.a must be positive.');
end
X=model.nodes;
T=model.elements.tet4;
if ~(isnumeric(X)&&ismatrix(X)&&size(X,2)==3&&all(isfinite(X(:))))
    error('HMIDemo:Nodes','model.nodes must be finite N-by-3 data.');
end
if ~(isnumeric(T)&&ismatrix(T)&&size(T,2)==4&&all(T(:)>=1)&& ...
        all(T(:)<=size(X,1))&&all(abs(T(:)-round(T(:)))<1e-12))
    error('HMIDemo:Tets','model.elements.tet4 is invalid.');
end
nTet=size(T,1);
sm=logical(model.region.solidMask(:));
fm=logical(model.region.fluidMask(:));
pi=model.region.portIndex(:);
if numel(sm)~=nTet || numel(fm)~=nTet || numel(pi)~=nTet
    error('HMIDemo:RegionLength','Region arrays must match the tetrahedron count.');
end
if any(sm==fm) || any(~(sm|fm))
    error('HMIDemo:RegionPartition','solidMask and fluidMask must partition all tetrahedra.');
end
if ~any(sm) || ~any(fm)
    error('HMIDemo:RegionPartition', ...
        'Both solid and exterior-fluid tetrahedra are required.');
end
if any(pi<0) || any(abs(pi-round(pi))>1e-12)
    error('HMIDemo:PortIndex','portIndex must contain nonnegative integers.');
end
if any(pi>0 & ~sm)
    error('HMIDemo:PortIndex','Every port tetrahedron must be solid.');
end
if max(pi)<1
    error('HMIDemo:Ports','At least one loading port is required.');
end
F=model.boundary.outerFaces;
if ~(isnumeric(F)&&size(F,2)==3&&all(F(:)>=1)&&all(F(:)<=size(X,1)))
    error('HMIDemo:OuterFaces','model.boundary.outerFaces is invalid.');
end
end

function T=demo_orient_tet4_local(X,T)
x1=X(T(:,1),:);
x2=X(T(:,2),:);
x3=X(T(:,3),:);
x4=X(T(:,4),:);
detJ=dot(x2-x1,cross(x3-x1,x4-x1,2),2);
if any(abs(detJ)<=1e-14*max(1,max(abs(detJ))))
    error('HMIDemo:DegenerateTet','Degenerate tetrahedron detected.');
end
neg=detJ<0;
T(neg,[2 3])=T(neg,[3 2]);
end

% =====================================================================
% P2 edge enrichment
% =====================================================================
function p2=demo_enrich_p2(X,T,outerFaces)
n0=size(X,1);
edges=[T(:,[1 2]);
T(:,[2 3]);
T(:,[1 3]);
T(:,[1 4]);
T(:,[2 4]);
T(:,[3 4])];
edges=sort(edges,2);
[edgeTable,~,ic]=unique(edges,'rows');
nE=size(edgeTable,1);
edgeNode=(1:nE).'+n0;
Xmid=0.5*(X(edgeTable(:,1),:)+X(edgeTable(:,2),:));
X2=[X;
Xmid];

nTet=size(T,1);
ic=reshape(ic,nTet,6);
T10=[T, n0+ic];

outerEdges=[outerFaces(:,[1 2]);
outerFaces(:,[2 3]);
outerFaces(:,[1 3])];
outerEdges=sort(outerEdges,2);
[tf,loc]=ismember(outerEdges,edgeTable,'rows');
outerMid=n0+loc(tf);
fixedNode=false(size(X2,1),1);
fixedNode(unique(outerFaces(:)))=true;
fixedNode(outerMid)=true;

p2=struct('nodes',X2,'tet10',T10,'edgeTable',edgeTable, ...
    'fixedNode',fixedNode,'nOriginalNode',n0);
end

% =====================================================================
% Geometry: connected solid components and port moments
% =====================================================================
function geom=demo_geometry_data(model)
X=model.nodes;
T=model.elements.tet4;
solidMask=logical(model.region.solidMask(:));
portIndex=model.region.portIndex(:);
solidIds=find(solidMask);
Ts=T(solidMask,:);

% Connectivity through shared triangular faces.
faces=[Ts(:,[1 2 3]);
Ts(:,[1 2 4]);
Ts(:,[1 3 4]);
Ts(:,[2 3 4])];
faces=sort(faces,2);
owners=repmat((1:size(Ts,1)).',4,1);
[~,~,fc]=unique(faces,'rows');
count=accumarray(fc,1);
ownerMin=accumarray(fc,owners,[],@min,0);
ownerMax=accumarray(fc,owners,[],@max,0);
shared=find(count==2);
A=sparse(ownerMin(shared),ownerMax(shared),1,size(Ts,1),size(Ts,1));
A=spones(A+A'+speye(size(A,1)));
compLocal=demo_sparse_components(A);
nComp=max(compLocal);
compFull=zeros(size(T,1),1);
compFull(solidIds)=compLocal;

% Component volumes and centroids.
compVolume=zeros(nComp,1);
compCentroid=zeros(nComp,3);
vol=demo_tet_volumes(X,T);
tetCent=(X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:)+X(T(:,4),:))/4;
for c=1:nComp
    idx=find(compFull==c);
    compVolume(c)=sum(vol(idx));
    compCentroid(c,:)=sum(tetCent(idx,:).*vol(idx),1)/compVolume(c);
end

nPort=max(portIndex);
portVolume=zeros(nPort,1);
portCentroid=zeros(nPort,3);
portJ=zeros(3,3,nPort);
[Lq,Wq]=demo_tet_quadrature4();
for p=1:nPort
    idx=find(portIndex==p);
    if isempty(idx);
    error('HMIDemo:Port','Port %d is empty.',p);
    end
    vp=vol(idx);
    cp=tetCent(idx,:);
    portVolume(p)=sum(vp);
    portCentroid(p,:)=sum(cp.*vp,1)/portVolume(p);

    Jp=zeros(3);
    for ee=idx(:).'
        Xe=X(T(ee,:),:);
        for iq=1:numel(Wq)
            xq=Lq(iq,:)*Xe;
            r=xq-portCentroid(p,:);
            Jp=Jp+vol(ee)*Wq(iq)*((r*r.')*eye(3)-r.'*r);
        end
    end
    Jp=0.5*(Jp+Jp.');
    if rcond(Jp)<1e-10
        error('HMIDemo:PortInertia','Port %d inertia tensor is singular.',p);
    end
    portJ(:,:,p)=Jp;
end

portComponent=zeros(nPort,1);
for p=1:nPort
    labels=unique(compFull(portIndex==p));
    labels(labels==0)=[];
    if numel(labels)~=1
        error('HMIDemo:PortComponent', ...
            'Port %d must belong to exactly one solid component.',p);
    end
    portComponent(p)=labels;
end

geom=struct('nComponent',nComp,'componentOfTet',compFull, ...
    'componentVolume',compVolume,'componentCentroid',compCentroid, ...
    'nPort',nPort,'portVolume',portVolume,'portCentroid',portCentroid, ...
    'portJ',portJ,'portComponent',portComponent);
end

function labels=demo_sparse_components(A)
n=size(A,1);
labels=zeros(n,1);
c=0;
for i=1:n
    if labels(i)~=0
    continue;
    end
    c=c+1;
    stack=i;
    labels(i)=c;
    while ~isempty(stack)
        v=stack(end);
        stack(end)=[];
        nb=find(A(v,:));
        nb=nb(labels(nb)==0);
        labels(nb)=c;
        stack=[stack;
        nb(:)];
    end
end
end

function vol=demo_tet_volumes(X,T)
x1=X(T(:,1),:);
x2=X(T(:,2),:);
x3=X(T(:,3),:);
x4=X(T(:,4),:);
detJ=dot(x2-x1,cross(x3-x1,x4-x1,2),2);
vol=abs(detJ)/6;
if any(vol<=1e-14*max(vol))
    error('HMIDemo:DegenerateTet','Degenerate tetrahedron detected.');
end
end

% =====================================================================
% FEM assembly
% =====================================================================
function fem=demo_assemble_fem(model,p2,geom,cfg)
X=model.nodes;
T=model.elements.tet4;
T10=p2.tet10;
solidMask=logical(model.region.solidMask(:));
fluidMask=logical(model.region.fluidMask(:));
portIndex=model.region.portIndex(:);
nTet=size(T,1);
nVnode=size(p2.nodes,1);
nVfull=3*nVnode;

% Separate P1 pressures in the solid and fluid.
solidPnodes=unique(T(solidMask,:));
fluidPnodes=unique(T(fluidMask,:));
mapS=zeros(size(X,1),1);
mapF=zeros(size(X,1),1);
mapS(solidPnodes)=1:numel(solidPnodes);
mapF(fluidPnodes)=1:numel(fluidPnodes);
nPs=numel(solidPnodes);
nPf=numel(fluidPnodes);
nP=nPs+nPf;

% Upper-triangle local indices for symmetric stiffness.
[ku,kj]=find(triu(true(30)));
nKu=numel(ku);
[mu,mj]=demo_mass_upper_pairs();
nMu=numel(mu);

nSolid=nnz(solidMask);
nFluid=nnz(fluidMask);
IKs=zeros(nKu*nSolid,1);
JKs=IKs;
VKs=IKs;
IKf=zeros(nKu*nFluid,1);
JKf=IKf;
VKf=IKf;
IM=zeros(nMu*nSolid,1);
JM=IM;
VM=IM;
IC=zeros(120*nTet,1);
JC=IC;
VC=IC;
ptrS=0;
ptrF=0;
ptrM=0;
ptrC=0;

nLoad=6*geom.nPort;
Ffull=zeros(nVfull,nLoad);
[Lq,Wq]=demo_tet_quadrature11();

dLamRef=[-1 -1 -1;
1 0 0;
0 1 0;
0 0 1];

for e=1:nTet
    vtx=T(e,:);
    Xe=X(vtx,:);
    J=[Xe(2,:)-Xe(1,:);
    Xe(3,:)-Xe(1,:);
    Xe(4,:)-Xe(1,:)].';
    detJ=det(J);
    if detJ<=0
        error('HMIDemo:Orientation','Tetrahedron %d is not positively oriented.',e);
    end
    vol=detJ/6;
    gradLam=dLamRef/J;

    Kloc=zeros(30);
    M10=zeros(10);
    Cloc=zeros(4,30);
    floc=zeros(30,nLoad);

    for iq=1:numel(Wq)
        L=Lq(iq,:);
        w=vol*Wq(iq);
        [N,gN]=demo_p2_shapes(L,gradLam);

        for i=1:10
            gi=gN(i,:);
            for j=1:10
                gj=gN(j,:);
                dij=dot(gi,gj);
                for alpha=1:3
                    ia=3*(i-1)+alpha;
                    for beta=1:3
                        jb=3*(j-1)+beta;
                        Kloc(ia,jb)=Kloc(ia,jb)+w*( ...
                            (alpha==beta)*dij + gi(beta)*gj(alpha));
                    end
                end
                M10(i,j)=M10(i,j)+w*N(i)*N(j);
            end
            for p=1:4
                for alpha=1:3
                    Cloc(p,3*(i-1)+alpha)=Cloc(p,3*(i-1)+alpha)+ ...
                        w*L(p)*gN(i,alpha);
                end
            end
        end

        pidx=portIndex(e);
        if pidx>0
            xq=L*Xe;
            base=6*(pidx-1);
            Vp=geom.portVolume(pidx);
            cp=geom.portCentroid(pidx,:);
            Jp=geom.portJ(:,:,pidx);

            profiles=zeros(3,6);
            profiles(:,1:3)=eye(3)/Vp;
            r=xq-cp;
            for alpha=1:3
                ealpha=zeros(3,1);
                ealpha(alpha)=1;
                omegaVec=Jp\ealpha;
                profiles(:,3+alpha)=cross(omegaVec.',r).';
            end

            for i=1:10
                rows=3*(i-1)+(1:3);
                floc(rows,base+(1:6))=floc(rows,base+(1:6))+ ...
                    w*N(i)*profiles;
            end
        end
    end

    ldof=reshape([3*T10(e,:)-2;
    3*T10(e,:)-1;
    3*T10(e,:)],1,[]);
    Kvals=Kloc(sub2ind([30 30],ku,kj));

    if solidMask(e)
        idx=ptrS+(1:nKu);
        ptrS=ptrS+nKu;
        IKs(idx)=ldof(ku);
        JKs(idx)=ldof(kj);
        VKs(idx)=Kvals;

        Mloc=cfg.rhoS*kron(M10,eye(3));
        mvals=Mloc(sub2ind([30 30],mu,mj));
        idx=ptrM+(1:nMu);
        ptrM=ptrM+nMu;
        IM(idx)=ldof(mu);
        JM(idx)=ldof(mj);
        VM(idx)=mvals;

        prow=mapS(vtx);
    else
        idx=ptrF+(1:nKu);
        ptrF=ptrF+nKu;
        IKf(idx)=ldof(ku);
        JKf(idx)=ldof(kj);
        VKf(idx)=Kvals;
        prow=nPs+mapF(vtx);
    end

    [pr,vc]=ndgrid(prow,ldof);
    vals=Cloc;
    idx=ptrC+(1:120);
    ptrC=ptrC+120;
    IC(idx)=pr(:);
    JC(idx)=vc(:);
    VC(idx)=vals(:);

    Ffull(ldof,:)=Ffull(ldof,:)+floc;
end

portResultantRaw=zeros(6,6,geom.nPort);
portResultantNormalized=zeros(6,6,geom.nPort);
portResultantCondition=zeros(geom.nPort,1);
for p=1:geom.nPort
    cols=6*(p-1)+(1:6);
    Rlocal=demo_rigid_coefficients_node_major( ...
        p2.nodes,geom.portCentroid(p,:));
    Braw=Rlocal.'*Ffull(:,cols);
    portResultantCondition(p)=cond(Braw);
    if ~isfinite(portResultantCondition(p)) || ...
            portResultantCondition(p)>1e8
        error('HMIDemo:PortResultant', ...
            'Port %d discrete resultant matrix is ill-conditioned.',p);
    end
    Ffull(:,cols)=Ffull(:,cols)/Braw;
    portResultantRaw(:,:,p)=Braw;
    portResultantNormalized(:,:,p)=Rlocal.'*Ffull(:,cols);
end
maxNormalizedPortResultantError=0;
for p=1:geom.nPort
    maxNormalizedPortResultantError=max( ...
        maxNormalizedPortResultantError, ...
        norm(portResultantNormalized(:,:,p)-eye(6),'fro')/sqrt(6));
end

KsU=sparse(IKs,JKs,VKs,nVfull,nVfull);
Ks=KsU+KsU.'-spdiags(diag(KsU),0,nVfull,nVfull);

KfU=sparse(IKf,JKf,VKf,nVfull,nVfull);
Kf=KfU+KfU.'-spdiags(diag(KfU),0,nVfull,nVfull);

MU=sparse(IM,JM,VM,nVfull,nVfull);
Ms=MU+MU.'-spdiags(diag(MU),0,nVfull,nVfull);

C=sparse(IC(1:ptrC),JC(1:ptrC),VC(1:ptrC),nP,nVfull);

freeNode=find(~p2.fixedNode);
freeDof=reshape([3*freeNode-2,3*freeNode-1,3*freeNode].',[],1);
freeDof=sort(freeDof);

Ks=Ks(freeDof,freeDof);
Kf=Kf(freeDof,freeDof);
Ms=Ms(freeDof,freeDof);
C=C(:,freeDof);
Fphysical=Ffull(freeDof,:);
pressureNullResidual=norm(C.'*ones(nP,1))/max(norm(C,'fro'),eps);
if pressureNullResidual>1e-8
    error('HMIDemo:PressureGauge', ...
        ['The expected single global pressure-null mode is not resolved. ', ...
         'Check the conforming solid-fluid mesh and outer boundary.']);
end

% Reference force/velocity scaling.
a=cfg.a;
Fref=cfg.G*a^2;
Tref=cfg.G*a^3;
Vref=a/cfg.tRef;
Oref=1/cfg.tRef;
DG=zeros(1,nLoad);
DV=zeros(1,nLoad);
for p=1:geom.nPort
    cols=6*(p-1)+(1:6);
    DG(cols)=[Fref Fref Fref Tref Tref Tref];
    DV(cols)=[Vref Vref Vref Oref Oref Oref];
end

fem=struct();
fem.Ks=Ks;
fem.Kf=Kf;
fem.Ms=Ms;
fem.C=C;
fem.Fphysical=Fphysical;
fem.DG=spdiags(DG(:),0,nLoad,nLoad);
fem.DV=spdiags(DV(:),0,nLoad,nLoad);
fem.freeDof=freeDof;
fem.nV=numel(freeDof);
fem.nP=nP;
fem.systemSize=fem.nV+nP+1;
fem.p2=p2;
fem.portResultantRaw=portResultantRaw;
fem.portResultantNormalized=portResultantNormalized;
fem.portResultantCondition=portResultantCondition;
fem.maxNormalizedPortResultantError=maxNormalizedPortResultantError;
fem.pressureNullResidual=pressureNullResidual;
end

function R=demo_rigid_coefficients_node_major(X,c)
n=size(X,1);
r=X-c;
R=zeros(3*n,6);
for i=1:n
    rows=3*(i-1)+(1:3);
    R(rows,1:3)=eye(3);
    R(rows,4:6)=-demo_skew(r(i,:));
end
end

function [ui,uj]=demo_mass_upper_pairs()
ui=[];
uj=[];
for i=1:10
    for j=i:10
        for a=1:3
            ui(end+1,1)=3*(i-1)+a;
            uj(end+1,1)=3*(j-1)+a;
        end
    end
end
end

function [N,gN]=demo_p2_shapes(L,gradLam)
N=zeros(10,1);
gN=zeros(10,3);
for i=1:4
    N(i)=L(i)*(2*L(i)-1);
    gN(i,:)=(4*L(i)-1)*gradLam(i,:);
end
pairs=[1 2;
2 3;
1 3;
1 4;
2 4;
3 4];
for k=1:6
    i=pairs(k,1);
    j=pairs(k,2);
    N(4+k)=4*L(i)*L(j);
    gN(4+k,:)=4*(L(j)*gradLam(i,:)+L(i)*gradLam(j,:));
end
end

% =====================================================================
% Rank structure
% =====================================================================
function r=demo_rank_structure(fem,geom,cfg)
nLoad=size(fem.Fphysical,2);

% Rigid resultant map in nondimensional generalized coordinates.
B=zeros(6*geom.nComponent,nLoad);
for p=1:geom.nPort
    c=geom.portComponent(p);
    rr=(geom.portCentroid(p,:)-geom.componentCentroid(c,:))/cfg.a;
    rows=6*(c-1)+(1:6);
    cols=6*(p-1)+(1:6);
    B(rows(1:3),cols(1:3))=eye(3);
    B(rows(4:6),cols(1:3))=demo_skew(rr);
    B(rows(4:6),cols(4:6))=eye(3);
end

% Discrete incompressibility projection of the scaled loading columns.
Fscaled=fem.Fphysical*fem.DG;
C=fem.C;
nV=fem.nV;
nP=fem.nP;
g=ones(nP,1);
g=g/norm(g);
sys=[speye(nV),C',sparse(nV,1);
...
     C,sparse(nP,nP),g;
     ...
     sparse(1,nV),g',0];
rhs=[Fscaled;
zeros(nP,nLoad);
zeros(1,nLoad)];
Uproj=decomposition(sys,'lu')\rhs;
Uproj=Uproj(1:nV,:);
Gp=real(0.5*(Fscaled.'*Uproj+(Fscaled.'*Uproj).'));
[V,D]=eig(full(Gp),'vector');
[d,ord]=sort(real(D),'descend');
V=V(:,ord);
tol=max(abs(d))*cfg.rankRelTol;
qw=sum(d>max(tol,eps));
if qw==0
error('HMIDemo:ProjectedLoad','Projected loading rank is zero.');
end
Uact=V(:,1:qw);

Bact=B*Uact;
q0=demo_numerical_rank(Bact,cfg.rankRelTol);
[~,~,Vb]=svd(Bact);
Qact=Vb(:,q0+1:end);
Qhidden=Uact*Qact;
qh=size(Qhidden,2);

r=struct('B',B,'projectedGram',Gp,'activeBasis',Uact, ...
    'Qhidden',Qhidden,'q0',q0,'qw',qw,'qh',qh);
end

function S=demo_skew(r)
S=[0 -r(3) r(2);
r(3) 0 -r(1);
-r(2) r(1) 0];
end

function r=demo_numerical_rank(A,relTol)
s=svd(full(A));
if isempty(s) || max(s)==0
r=0;
return;
end
r=sum(s>relTol*max(s));
end

% =====================================================================
% Dynamic mixed solve
% =====================================================================
function sol=demo_dynamic_solve(fem,cfg,omegaHat)
s=1i*omegaHat/cfg.tRef;
A=cfg.etaF*fem.Kf + (cfg.etaS+cfg.G/s)*fem.Ks + s*fem.Ms;
C=fem.C;
nV=fem.nV;
nP=fem.nP;
g=ones(nP,1);
g=g/norm(g);
system=[A,C',sparse(nV,1);
...
        C,sparse(nP,nP),g;
        ...
        sparse(1,nV),g',0];
RHSvel=fem.Fphysical*fem.DG;
rhs=[RHSvel;
zeros(nP,size(RHSvel,2));
zeros(1,size(RHSvel,2))];
X=decomposition(system,'lu')\rhs;
U=X(1:nV,:);
P=X(nV+(1:nP),:);
Y=fem.Fphysical.'*U;
Mhat=fem.DV\Y;
res=system*X-rhs;
sol=struct();
sol.Mhat=Mhat;
sol.linearResidual=norm(res,'fro')/max(norm(rhs,'fro'),eps);
sol.divergenceResidual=norm(C*U,'fro')/max(norm(RHSvel,'fro'),eps);
sol.gaugeResidual=norm(g.'*P,'fro')/max(norm(P,'fro'),eps);
sol.reciprocityResidual=norm(Mhat-Mhat.','fro')/max(norm(Mhat,'fro'),eps);
H=0.5*(Mhat+Mhat');
sol.passivityMinEigenvalue=min(real(eig(full(H))));
end

% =====================================================================
% Figures
% =====================================================================
function files=demo_write_figures(outdir,model,geom,rankData,lambda,omega,R,cfg)
files=struct();

% Figure 1: geometry mesh only.
f1=figure('Color','w','Name','Geometry mesh','Renderer','painters');
ax=axes(f1);
hold(ax,'on');
axis(ax,'equal');
axis(ax,'off');
view(ax,34,22);
F=demo_free_boundary_faces(model.elements.tet4(model.region.solidMask,:));
patch(ax,'Faces',F,'Vertices',model.nodes,'FaceColor',[0.92 0.92 0.92], ...
    'EdgeColor',[0.25 0.25 0.25],'LineWidth',0.35);
title(ax,sprintf( ...
    '%s,  $C=%d$,   $P=%d$,   $q_0=%d$,   $q_{\\omega}=%d$,   $q_h=%d$', ...
    char(string(model.name)),geom.nComponent,geom.nPort, ...
    rankData.q0,rankData.qw,rankData.qh), ...
    'Interpreter','latex');
files.geometryPng=fullfile(outdir,'01_geometry_mesh.png');
files.geometryFig=fullfile(outdir,'01_geometry_mesh.fig');
files.geometryEps=fullfile(outdir,'01_geometry_mesh.eps');
exportgraphics(f1,files.geometryPng,'Resolution',cfg.pngResolution,'BackgroundColor','white');
print(f1,files.geometryEps,'-depsc2','-painters');
savefig(f1,files.geometryFig,'compact');

% Figure 2: activation spectrum.
f2=figure('Color','w','Name','Hidden activation spectrum','Renderer','painters');
ax=axes(f2);
box(ax,'on');
grid(ax,'on');
hold(ax,'on');
if isempty(lambda)
    text(ax,0.5,0.5,'No hidden directions','Units','normalized', ...
        'HorizontalAlignment','center');
        axis(ax,'off');
else
    semilogy(ax,1:numel(lambda),lambda,'o-','LineWidth',1.4,'MarkerSize',6, ...
        'MarkerFaceColor','w');
    xlabel(ax,'hidden branch index $j$','Interpreter','latex');
    ylabel(ax,'$\widehat{\lambda}_j$','Interpreter','latex');
    xlim(ax,[0.7 numel(lambda)+0.3]);
    xticks(ax,1:numel(lambda));
end
files.activationPng=fullfile(outdir,'02_activation_spectrum.png');
files.activationFig=fullfile(outdir,'02_activation_spectrum.fig');
files.activationEps=fullfile(outdir,'02_activation_spectrum.eps');
exportgraphics(f2,files.activationPng,'Resolution',cfg.pngResolution,'BackgroundColor','white');
print(f2,files.activationEps,'-depsc2','-painters');
savefig(f2,files.activationFig,'compact');

% Figure 3: normalized opening.
f3=figure('Color','w','Name','Normalized hidden opening','Renderer','painters');
ax=axes(f3);
box(ax,'on');
grid(ax,'on');
hold(ax,'on');
if isempty(R)
    text(ax,0.5,0.5,'No hidden directions','Units','normalized', ...
        'HorizontalAlignment','center');
        axis(ax,'off');
else
    for j=1:size(R,1)
        semilogx(ax,omega,R(j,:),'o-','LineWidth',1.2,'MarkerSize',5, ...
            'DisplayName',sprintf('$j=%d$',j))
    end
    yline(ax,1,'--','HandleVisibility','off');
    xlabel(ax,'$\widehat{\omega}$','Interpreter','latex');
    ylabel(ax, ...
    '$\widehat{\sigma}_{q_0+j}/(|\widehat{\omega}|\,\widehat{\lambda}_j)$', ...
    'Interpreter','latex');
    legend(ax,'Location','best','Interpreter','latex');
end
files.openingPng=fullfile(outdir,'03_normalized_opening.png');
files.openingFig=fullfile(outdir,'03_normalized_opening.fig');
files.openingEps=fullfile(outdir,'03_normalized_opening.eps');
exportgraphics(f3,files.openingPng,'Resolution',cfg.pngResolution,'BackgroundColor','white');
print(f3,files.openingEps,'-depsc2','-painters');
savefig(f3,files.openingFig,'compact');
end

function faces=demo_free_boundary_faces(T)
allF=[T(:,[1 2 3]);
T(:,[1 2 4]);
T(:,[1 3 4]);
T(:,[2 3 4])];
key=sort(allF,2);
[~,~,ic]=unique(key,'rows');
count=accumarray(ic,1);
faces=allF(count(ic)==1,:);
end

% =====================================================================
% Text output
% =====================================================================
function demo_write_text_results(filename,matFile,model,p2,geom,fem,cfg,r, ...
    lambda,lambdaW,hiddenSigma,R,order,sv,lin,div,gauge,recip,pass,solveT,files)
fid=fopen(filename,'w');
assert(fid>0,'Cannot create %s.',filename);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'RESULTS\n');
fprintf(fid,'Input MAT                         %s\n',matFile);
fprintf(fid,'Geometry                          %s\n',char(string(model.name)));
fprintf(fid,'P1 nodes                          %d\n',size(model.nodes,1));
fprintf(fid,'P2 velocity nodes                 %d\n',size(p2.nodes,1));
fprintf(fid,'Tetrahedra                        %d\n',size(model.elements.tet4,1));
fprintf(fid,'Solid tetrahedra                  %d\n',nnz(model.region.solidMask));
fprintf(fid,'Fluid tetrahedra                  %d\n',nnz(model.region.fluidMask));
fprintf(fid,'Solid components C                %d\n',geom.nComponent);
fprintf(fid,'Ports P                           %d\n',geom.nPort);
fprintf(fid,'Load coordinates n                %d\n',6*geom.nPort);
fprintf(fid,'Free velocity dofs                %d\n',fem.nV);
fprintf(fid,'Pressure dofs                     %d\n',fem.nP);
fprintf(fid,'Mixed+gauge system size           %d\n',fem.systemSize);
fprintf(fid,'max normalized port resultant err %.12e\n\n', ...
    fem.maxNormalizedPortResultantError);
fprintf(fid,'DEMONSTRATION PARAMETERS\n');
fprintf(fid,'a                                 %.12e\n',cfg.a);
fprintf(fid,'rho_s                             %.12e\n',cfg.rhoS);
fprintf(fid,'G                                 %.12e\n',cfg.G);
fprintf(fid,'Lambda_f                          %.12e\n',cfg.LambdaF);
fprintf(fid,'Lambda_s                          %.12e\n',cfg.LambdaS);
fprintf(fid,'frequencies                       ');
fprintf(fid,'%.6e ',cfg.omegaHat);
fprintf(fid,'\n\n');
fprintf(fid,'OBSERVABLE DIMENSIONS\n');
fprintf(fid,'q0 static observable              %d\n',r.q0);
fprintf(fid,'q_omega harmonic observable       %d\n',r.qw);
fprintf(fid,'qh hidden                         %d\n\n',r.qh);

if ~isempty(lambda)
    fprintf(fid,'HIDDEN ACTIVATION\n');
    for j=1:numel(lambda)
        fprintf(fid,'lambda_%d                          %.12e\n',j,lambda(j));
    end
    fprintf(fid,'\nPROJECTED ACTIVATION BY FREQUENCY\n');
    for iw=1:numel(cfg.omegaHat)
        fprintf(fid,'omega %.6e  ',cfg.omegaHat(iw));
        fprintf(fid,'%.10e ',lambdaW(:,iw));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\nHIDDEN OPENING\n');
    for j=1:size(hiddenSigma,1)
        fprintf(fid,'branch %d fitted order             %.8f\n',j,order(j));
        for iw=1:numel(cfg.omegaHat)
            fprintf(fid,'  omega %.6e sigma %.10e R %.10e\n', ...
                cfg.omegaHat(iw),hiddenSigma(j,iw),R(j,iw));
        end
    end
    fprintf(fid,'\n');
end

fprintf(fid,'NUMERICAL DIAGNOSTICS\n');
for iw=1:numel(cfg.omegaHat)
    fprintf(fid,'omega %.6e\n',cfg.omegaHat(iw));
    rankMob=sum(sv(:,iw)>cfg.rankRelTol*max(sv(:,iw)));
    fprintf(fid,'  mobility rank                   %d\n',rankMob);
    fprintf(fid,'  linear residual                 %.12e\n',lin(iw));
    fprintf(fid,'  divergence residual             %.12e\n',div(iw));
    fprintf(fid,'  gauge residual                  %.12e\n',gauge(iw));
    fprintf(fid,'  reciprocity residual            %.12e\n',recip(iw));
    fprintf(fid,'  minimum passivity eigenvalue    %.12e\n',pass(iw));
    fprintf(fid,'  solve seconds                   %.3f\n',solveT(iw));
end
fprintf(fid,'\nFIGURE FILES\n');
fprintf(fid,'%s\n%s\n%s\n',files.geometryPng,files.activationPng,files.openingPng);
fprintf(fid,'\nNOTE\n');
fprintf(fid,['This lightweight solver and the example meshes are intended to ', ...
    'reproduce the mechanism rather than the publication numerical values.\n']);
end

% =====================================================================
% Tetrahedral quadrature
% =====================================================================
function [L,w]=demo_tet_quadrature4()
a=0.5854101966249685;
b=0.1381966011250105;
L=[a b b b;
b a b b;
b b a b;
b b b a];
w=0.25*ones(4,1);
end

function [L,w]=demo_tet_quadrature11()
% Degree-4 symmetric tetrahedral rule. Weights are normalized to sum to 1,
% so an affine tetrahedron integral is volume * sum(w_q f_q).
a=0.785714285714286;
b=0.0714285714285714;
c=0.399403576166799;
d=0.100596423833201;
L=zeros(11,4);
w=zeros(11,1);
L(1,:)=[0.25 0.25 0.25 0.25];
w(1)=-0.0789333333333333;
for i=1:4
    q=b*ones(1,4);
    q(i)=a;
    L(1+i,:)=q;
    w(1+i)=0.0457333333333333;
end
pairs=nchoosek(1:4,2);
for k=1:6
    q=d*ones(1,4);
    q(pairs(k,:))=c;
    L(5+k,:)=q;
    w(5+k)=0.149333333333333;
end
if abs(sum(w)-1)>1e-12
    error('HMIDemo:Quadrature','Tetrahedral quadrature weights do not sum to one.');
end
end
