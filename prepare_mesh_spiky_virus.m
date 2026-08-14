function outputFile = prepare_mesh_spiky_virus(gmshExecutable)
% Usage
%   outputFile = prepare_mesh_spiky_virus();
%   outputFile = prepare_mesh_spiky_virus('path\to\gmsh.exe');
%
% Requirements
%   - Gmsh command-line executable.

if nargin < 1
    gmshExecutable = '';
end

close all;
clc;

% =====================================================================
% USER-EDITABLE GEOMETRY
% =====================================================================
a = 1.0;

coreRadius = 0.70*a;

spikeBaseRadius = 0.20*a;
spikeBaseCenterRadius = 0.60*a;
spikeLength = 0.46*a;

outerRadius = 1.80*a;

portRadius = 0.15*a;
portCenters = [ ...
    -0.25*a 0 0; ...
     0.25*a 0 0];

hPort = 0.18*a;
hSolid = 0.24*a;
hFluid = 0.65*a;

outputFile = fullfile(pwd,'mesh_spiky_virus.mat');
workDir = fullfile(pwd,'mesh_spiky_virus_work');
% =====================================================================
% END USER-EDITABLE GEOMETRY
% =====================================================================

s = 1/sqrt(3);

directions = [ ...
     1  0  0; ...
    -1  0  0; ...
     0  1  0; ...
     0 -1  0; ...
     0  0  1; ...
     0  0 -1; ...
     s  s  s; ...
     s  s -s; ...
     s -s  s; ...
     s -s -s; ...
    -s  s  s; ...
    -s  s -s; ...
    -s -s  s; ...
    -s -s -s];

nSpike = size(directions,1);

spikeBaseCenters = spikeBaseCenterRadius*directions;
spikeVectors = spikeLength*directions;

if coreRadius >= outerRadius
    error('HMIDemo:Geometry', ...
        'outerRadius must exceed coreRadius.');
end

tipRadius = spikeBaseCenterRadius+spikeLength;

if tipRadius >= outerRadius
    error('HMIDemo:Geometry', ...
        'outerRadius must enclose all spike tips.');
end

if sqrt(spikeBaseCenterRadius^2+spikeBaseRadius^2) >= coreRadius
    error('HMIDemo:Geometry', ...
        'Each spike base must lie strictly inside the core.');
end

for p=1:size(portCenters,1)
    if norm(portCenters(p,:))+portRadius >= coreRadius
        error('HMIDemo:Geometry', ...
            'Port %d is not strictly inside the core.',p);
    end
end

for p=1:size(portCenters,1)
    for q=p+1:size(portCenters,1)
        if norm(portCenters(p,:)-portCenters(q,:)) <= 2*portRadius
            error('HMIDemo:Geometry', ...
                'Ports %d and %d overlap.',p,q);
        end
    end
end

if exist(workDir,'dir') ~= 7
    mkdir(workDir);
end

[gmshExecutable,gmshVersion] = demo_resolve_gmsh(gmshExecutable);

geoFile = fullfile(workDir,'mesh_spiky_virus.geo');
mshFile = fullfile(workDir,'mesh_spiky_virus.msh');

demo_write_spiky_geo( ...
    geoFile, ...
    coreRadius, ...
    spikeBaseCenters,spikeVectors,spikeBaseRadius, ...
    outerRadius, ...
    portCenters,portRadius, ...
    hPort,hSolid,hFluid);

cmd = sprintf('"%s" "%s" -3 -format msh2 -o "%s" -v 2', ...
    gmshExecutable,geoFile,mshFile);

fprintf('Running Gmsh...\n%s\n',cmd);

[status,out] = system(cmd);

if status ~= 0 || exist(mshFile,'file') ~= 2
    error('HMIDemo:GmshFailure', ...
        'Gmsh failed with status %d.\n%s',status,out);
end

raw = demo_read_msh2_3d(mshFile);

X = raw.nodes;
T = demo_orient_tet4(X,raw.tet4);

centroid = ( ...
    X(T(:,1),:) + ...
    X(T(:,2),:) + ...
    X(T(:,3),:) + ...
    X(T(:,4),:))/4;

solidMask = sqrt(sum(centroid.^2,2)) <= coreRadius*(1+1e-10);

for k=1:nSpike
    solidMask = solidMask | demo_points_in_cone( ...
        centroid, ...
        spikeBaseCenters(k,:), ...
        spikeVectors(k,:), ...
        spikeBaseRadius);
end

fluidMask = ~solidMask;

portIndex = zeros(size(T,1),1);

for p=1:size(portCenters,1)
    rp = sqrt(sum((centroid-portCenters(p,:)).^2,2));
    hit = rp <= portRadius*(1+1e-10);

    if any(hit & ~solidMask)
        error('HMIDemo:PortOutsideSolid', ...
            'Port %d contains elements outside the solid.',p);
    end

    if any(hit & portIndex>0)
        error('HMIDemo:PortOverlap', ...
            'Port %d overlaps a previous port.',p);
    end

    portIndex(hit) = p;
end

if any(portIndex>0 & ~solidMask)
    error('HMIDemo:PortClassification', ...
        'A port element is not solid.');
end

if ~any(solidMask) || ~any(fluidMask)
    error('HMIDemo:RegionClassification', ...
        'Solid or fluid classification is empty.');
end

for p=1:size(portCenters,1)
    if ~any(portIndex==p)
        error('HMIDemo:EmptyPort', ...
            'Port %d contains no tetrahedra.',p);
    end
end

tri = raw.tri3;

rv = [ ...
    sqrt(sum(X(tri(:,1),:).^2,2)), ...
    sqrt(sum(X(tri(:,2),:).^2,2)), ...
    sqrt(sum(X(tri(:,3),:).^2,2))];

outerFaceMask = max(abs(rv-outerRadius),[],2) <= 1e-6*outerRadius;
outerFaces = tri(outerFaceMask,:);

if isempty(outerFaces)
    error('HMIDemo:OuterBoundary', ...
        'Could not identify the outer boundary faces.');
end

% =====================================================================
% MODEL DATA CONTRACT
% =====================================================================
model = struct();

model.schemaVersion = 1;
model.name = 'spiky virus';
model.nodes = X;
model.elements = struct('tet4',T);

model.region = struct();
model.region.solidMask = logical(solidMask(:));
model.region.fluidMask = logical(fluidMask(:));
model.region.portIndex = portIndex(:);

model.boundary = struct();
model.boundary.outerFaces = outerFaces;
model.boundary.outerNodes = unique(outerFaces(:));

model.ports = repmat( ...
    struct('elements',[],'nominalCenter',[]), ...
    size(portCenters,1),1);

for p=1:size(portCenters,1)
    model.ports(p).elements = find(portIndex==p);
    model.ports(p).nominalCenter = portCenters(p,:);
end

model.scale = struct('a',a);

model.meta = struct();
model.meta.geometry = 'spiky virus';
model.meta.coreRadius = coreRadius;
model.meta.spikeBaseRadius = spikeBaseRadius;
model.meta.spikeBaseCenterRadius = spikeBaseCenterRadius;
model.meta.spikeLength = spikeLength;
model.meta.spikeDirections = directions;
model.meta.outerRadius = outerRadius;
model.meta.gmshExecutable = gmshExecutable;
model.meta.gmshVersion = gmshVersion;
model.meta.geoFile = geoFile;
model.meta.mshFile = mshFile;

save(outputFile,'model','-v7.3');

fprintf('\nSpiky-virus mesh prepared.\n');
fprintf('  core radius/a      %.8g\n',coreRadius/a);
fprintf('  spike base/a       %.8g\n',spikeBaseRadius/a);
fprintf('  spike length/a     %.8g\n',spikeLength/a);
fprintf('  spikes             %d\n',nSpike);
fprintf('  nodes              %d\n',size(X,1));
fprintf('  tetrahedra         %d\n',size(T,1));
fprintf('  solid tetrahedra   %d\n',nnz(solidMask));
fprintf('  fluid tetrahedra   %d\n',nnz(fluidMask));
fprintf('  port tetrahedra    ');

for p=1:numel(model.ports)
    fprintf('%d ',numel(model.ports(p).elements));
end

fprintf('\n');
fprintf('  outer faces        %d\n',size(outerFaces,1));
fprintf('  output MAT         %s\n',outputFile);
end

% =====================================================================
% Gmsh geometry
% =====================================================================
function demo_write_spiky_geo( ...
    filename, ...
    coreRadius, ...
    spikeBaseCenters,spikeVectors,spikeBaseRadius, ...
    outerRadius, ...
    portCenters,portRadius, ...
    hPort,hSolid,hFluid)

fid = fopen(filename,'w');
assert(fid>0,'Cannot create %s.',filename);
cleanup = onCleanup(@()fclose(fid));

fprintf(fid,'SetFactory("OpenCASCADE");\n');
fprintf(fid,'General.Terminal = 1;\n');
fprintf(fid,'Mesh.MshFileVersion = 2.2;\n');
fprintf(fid,'Mesh.Binary = 0;\n');
fprintf(fid,'Mesh.SaveAll = 1;\n');
fprintf(fid,'Mesh.ElementOrder = 1;\n');
fprintf(fid,'Mesh.Optimize = 1;\n');
fprintf(fid,'Mesh.MeshSizeFromPoints = 0;\n');
fprintf(fid,'Mesh.MeshSizeFromCurvature = 0;\n');
fprintf(fid,'Mesh.MeshSizeExtendFromBoundary = 0;\n');
fprintf(fid,'Mesh.MeshSizeMin = %.17g;\n',hPort);
fprintf(fid,'Mesh.MeshSizeMax = %.17g;\n',hFluid);

fprintf(fid,'Sphere(1) = {0,0,0,%.17g};\n',outerRadius);
fprintf(fid,'Sphere(2) = {0,0,0,%.17g};\n',coreRadius);

firstSpikeTag = 3;

for k=1:size(spikeBaseCenters,1)
    c = spikeBaseCenters(k,:);
    v = spikeVectors(k,:);
    tag = firstSpikeTag+k-1;

    fprintf(fid, ...
        'Cone(%d) = {%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,0};\n', ...
        tag, ...
        c(1),c(2),c(3), ...
        v(1),v(2),v(3), ...
        spikeBaseRadius);
end

firstPortTag = firstSpikeTag+size(spikeBaseCenters,1);

for p=1:size(portCenters,1)
    c = portCenters(p,:);
    tag = firstPortTag+p-1;

    fprintf(fid, ...
        'Sphere(%d) = {%.17g,%.17g,%.17g,%.17g};\n', ...
        tag,c(1),c(2),c(3),portRadius);
end

solidTags = 2:(firstPortTag+size(portCenters,1)-1);
tagList = sprintf('%d,',solidTags);
tagList(end) = [];

fprintf(fid,['fr[] = BooleanFragments{ Volume{1}; Delete; }', ...
    '{ Volume{%s}; Delete; };\n'],tagList);

fprintf(fid,'Coherence;\n');

extent = coreRadius+norm(spikeVectors(1,:))+0.10;

fprintf(fid,'Field[1] = Box;\n');
fprintf(fid,'Field[1].VIn = %.17g;\n',hSolid);
fprintf(fid,'Field[1].VOut = %.17g;\n',hFluid);
fprintf(fid,'Field[1].XMin = %.17g;\n',-extent);
fprintf(fid,'Field[1].XMax = %.17g;\n', extent);
fprintf(fid,'Field[1].YMin = %.17g;\n',-extent);
fprintf(fid,'Field[1].YMax = %.17g;\n', extent);
fprintf(fid,'Field[1].ZMin = %.17g;\n',-extent);
fprintf(fid,'Field[1].ZMax = %.17g;\n', extent);

fieldTags = 1;

for p=1:size(portCenters,1)
    c = portCenters(p,:);
    pad = 2.0*portRadius;
    tag = 1+p;

    fieldTags(end+1) = tag;

    fprintf(fid,'Field[%d] = Box;\n',tag);
    fprintf(fid,'Field[%d].VIn = %.17g;\n',tag,hPort);
    fprintf(fid,'Field[%d].VOut = %.17g;\n',tag,hSolid);
    fprintf(fid,'Field[%d].XMin = %.17g;\n',tag,c(1)-pad);
    fprintf(fid,'Field[%d].XMax = %.17g;\n',tag,c(1)+pad);
    fprintf(fid,'Field[%d].YMin = %.17g;\n',tag,c(2)-pad);
    fprintf(fid,'Field[%d].YMax = %.17g;\n',tag,c(2)+pad);
    fprintf(fid,'Field[%d].ZMin = %.17g;\n',tag,c(3)-pad);
    fprintf(fid,'Field[%d].ZMax = %.17g;\n',tag,c(3)+pad);
end

minTag = 20;
flist = sprintf('%d,',fieldTags);
flist(end) = [];

fprintf(fid,'Field[%d] = Min;\n',minTag);
fprintf(fid,'Field[%d].FieldsList = {%s};\n',minTag,flist);
fprintf(fid,'Background Field = %d;\n',minTag);
end

% =====================================================================
% Geometry helpers
% =====================================================================
function inside = demo_points_in_cone(X,base,vec,baseRadius)

L2 = dot(vec,vec);

if L2 <= 0
    error('HMIDemo:Geometry', ...
        'Cone vector must be nonzero.');
end

D = X-base;
t = (D*vec.')/L2;

proj = base+t.*vec;
radialDistance = sqrt(sum((X-proj).^2,2));

localRadius = baseRadius*(1-t);

inside = t >= -1e-10 & ...
    t <= 1+1e-10 & ...
    radialDistance <= localRadius.*(1+1e-10);
end

% =====================================================================
% Gmsh discovery and MSH2 parser
% =====================================================================
function [exe,versionText] = demo_resolve_gmsh(requested)

candidates = {};

if ~isempty(strtrim(char(requested)))
    candidates{end+1} = char(requested);
end

envExe = strtrim(getenv('GMSH_EXE'));

if ~isempty(envExe)
    candidates{end+1} = envExe;
end

candidates = [ ...
    candidates, ...
    {'gmsh',fullfile(pwd,'gmsh.exe'),'C:\Program Files\gmsh\gmsh.exe'}];

for k=1:numel(candidates)
    cmd = sprintf('"%s" -version',candidates{k});
    [status,out] = system(cmd);

    if status==0
        exe = candidates{k};
        versionText = strtrim(out);

        if isempty(versionText)
            versionText = 'unknown';
        end

        return;
    end
end

error('HMIDemo:GmshNotFound', ...
    ['Gmsh was not found. Supply its full path, define GMSH_EXE, ', ...
     'or place gmsh on the system PATH.']);
end

function raw = demo_read_msh2_3d(filename)

fid = fopen(filename,'r');
assert(fid>0,'Cannot open %s.',filename);
cleanup = onCleanup(@()fclose(fid));

nodeTags = [];
nodes = [];
tetTags = zeros(0,4);
triTags = zeros(0,3);

while ~feof(fid)
    line = strtrim(fgetl(fid));

    if strcmp(line,'$Nodes')
        n = str2double(strtrim(fgetl(fid)));
        nodeTags = zeros(n,1);
        nodes = zeros(n,3);

        for i=1:n
            v = sscanf(fgetl(fid),'%f').';
            nodeTags(i) = v(1);
            nodes(i,:) = v(2:4);
        end

        assert( ...
            strcmp(strtrim(fgetl(fid)),'$EndNodes'), ...
            'Malformed $Nodes section.');

    elseif strcmp(line,'$Elements')
        n = str2double(strtrim(fgetl(fid)));
        tetTags = zeros(n,4);
        triTags = zeros(n,3);
        nt = 0;
        nf = 0;

        for i=1:n
            v = sscanf(fgetl(fid),'%f').';

            if numel(v)<3
                continue;
            end

            typ = v(2);
            nTags = v(3);
            first = 4+nTags;

            if typ==4
                nt = nt+1;
                tetTags(nt,:) = v(first:first+3);

            elseif typ==2
                nf = nf+1;
                triTags(nf,:) = v(first:first+2);
            end
        end

        tetTags = tetTags(1:nt,:);
        triTags = triTags(1:nf,:);

        assert( ...
            strcmp(strtrim(fgetl(fid)),'$EndElements'), ...
            'Malformed $Elements section.');
    end
end

if isempty(nodes) || isempty(tetTags)
    error('HMIDemo:MeshFormat', ...
        'No tetrahedra found.');
end

maxTag = max(nodeTags);
map = zeros(maxTag,1);
map(nodeTags) = 1:numel(nodeTags);

raw = struct();
raw.nodes = nodes;
raw.tet4 = reshape(map(tetTags(:)),size(tetTags));
raw.tri3 = reshape(map(triTags(:)),size(triTags));
end

function T = demo_orient_tet4(X,T)

x1 = X(T(:,1),:);
x2 = X(T(:,2),:);
x3 = X(T(:,3),:);
x4 = X(T(:,4),:);

detJ = dot(x2-x1,cross(x3-x1,x4-x1,2),2);
neg = detJ<0;

T(neg,[2 3]) = T(neg,[3 2]);
end
