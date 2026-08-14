function outputFile = prepare_mesh_cube(gmshExecutable)
% Usage
%   outputFile = prepare_mesh_cube();
%   outputFile = prepare_mesh_cube('path\to\gmsh.exe');
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

side = (4*pi/3)^(1/3)*a;
halfSide = side/2;
outerRadius = 1.8*a;

portRadius = 0.16*a;
portCenters = [ ...
    -0.40*a 0 0; ...
     0.40*a 0 0];

hPort = 0.20*a;
hSolid = 0.34*a;
hFluid = 0.70*a;

outputFile = fullfile(pwd,'mesh_cube.mat');
workDir = fullfile(pwd,'mesh_cube_work');
% =====================================================================
% END USER-EDITABLE GEOMETRY
% =====================================================================

if outerRadius <= sqrt(3)*halfSide
    error('HMIDemo:Geometry', ...
        'outerRadius must enclose the entire cube.');
end

for p=1:size(portCenters,1)
    if any(abs(portCenters(p,:))+portRadius >= halfSide)
        error('HMIDemo:Geometry', ...
            'Port %d is not strictly inside the cube.',p);
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

geoFile = fullfile(workDir,'mesh_cube.geo');
mshFile = fullfile(workDir,'mesh_cube.msh');

demo_write_cube_geo( ...
    geoFile,a,side,outerRadius,portCenters,portRadius, ...
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

solidMask = all(abs(centroid) <= halfSide*(1+1e-10),2);
fluidMask = ~solidMask;

portIndex = zeros(size(T,1),1);

for p=1:size(portCenters,1)
    rp = sqrt(sum((centroid-portCenters(p,:)).^2,2));
    hit = rp <= portRadius*(1+1e-10);

    if any(hit & ~solidMask)
        error('HMIDemo:PortOutsideSolid', ...
            'Port %d contains elements outside the cube.',p);
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
model.name = 'cube';
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
model.meta.geometry = 'cube';
model.meta.side = side;
model.meta.outerRadius = outerRadius;
model.meta.gmshExecutable = gmshExecutable;
model.meta.gmshVersion = gmshVersion;
model.meta.geoFile = geoFile;
model.meta.mshFile = mshFile;

save(outputFile,'model','-v7.3');

fprintf('\nCube mesh prepared.\n');
fprintf('  side/a             %.8g\n',side/a);
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
function demo_write_cube_geo( ...
    filename,a,side,outerRadius,portCenters, ...
    portRadius,hPort,hSolid,hFluid)

fid = fopen(filename,'w');
assert(fid>0,'Cannot create %s.',filename);
cleanup = onCleanup(@()fclose(fid));

half = side/2;

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

fprintf(fid, ...
    'Box(2) = {%.17g,%.17g,%.17g,%.17g,%.17g,%.17g};\n', ...
    -half,-half,-half,side,side,side);

for p=1:size(portCenters,1)
    c = portCenters(p,:);

    fprintf(fid, ...
        'Sphere(%d) = {%.17g,%.17g,%.17g,%.17g};\n', ...
        2+p,c(1),c(2),c(3),portRadius);
end

portTags = 3:(2+size(portCenters,1));
plist = sprintf('%d,',portTags);
plist(end) = [];

fprintf(fid,['fr[] = BooleanFragments{ Volume{1}; Delete; }', ...
    '{ Volume{2,%s}; Delete; };\n'],plist);

fprintf(fid,'Coherence;\n');

padSolid = 0.10*a;

fprintf(fid,'Field[1] = Box;\n');
fprintf(fid,'Field[1].VIn = %.17g;\n',hSolid);
fprintf(fid,'Field[1].VOut = %.17g;\n',hFluid);
fprintf(fid,'Field[1].XMin = %.17g;\n',-half-padSolid);
fprintf(fid,'Field[1].XMax = %.17g;\n', half+padSolid);
fprintf(fid,'Field[1].YMin = %.17g;\n',-half-padSolid);
fprintf(fid,'Field[1].YMax = %.17g;\n', half+padSolid);
fprintf(fid,'Field[1].ZMin = %.17g;\n',-half-padSolid);
fprintf(fid,'Field[1].ZMax = %.17g;\n', half+padSolid);

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
