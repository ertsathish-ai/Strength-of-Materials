%% VM Mechanical Engineering Virtual Laboratory - SFD and BMD
% Euler-Bernoulli beam finite-element solver for:
% simply supported, cantilever, overhanging, and fixed beams.
% Units: m, kN, kN.m, GPa, m^4. Loads are entered as positive magnitudes.
% Point/UDL/UVL act downward. Applied moments act clockwise.
clear; clc; close all;

%% USER INPUTS
beamType = "ss";       % "ss", "cant", "over", or "fixed"
L = 6;                 % total beam length (m)
E_GPa = 200;           % Young's modulus (GPa)
I = 1.0e-4;            % second moment of area (m^4)

% Overhanging beam support coordinates (used only for beamType="over")
xA = 1.0; xB = 5.0;

% Load table columns: {type, magnitude, xStart, xEnd}
% type: "point", "udl", "uvl", "moment"
% For point/moment, xEnd is ignored. UVL rises linearly from zero at xStart
% to the entered magnitude at xEnd.
loads = {
    "point", 20, 3, 6
    % "udl", 5, 1, 5
    % "uvl", 8, 0, 4
    % "moment", 12, 3, 3
};

%% SUPPORTS
switch beamType
    case "ss", supports = [0 L];
    case "cant", supports = 0;
    case "fixed", supports = [0 L];
    case "over", supports = [xA xB];
    otherwise, error("Unknown beamType");
end
if L<=0 || E_GPa<=0 || I<=0, error("L, E and I must be positive."); end
E = E_GPa*1e6; % GPa to kN/m^2

%% VALIDATE LOADS AND CREATE MESH
nLoad = size(loads,1);
loadData = struct('type',{},'mag',{},'x1',{},'x2',{});
coords = [0 L supports];
for k=1:nLoad
    typ = string(loads{k,1}); mag = loads{k,2};
    x1 = loads{k,3}; x2 = loads{k,4};
    if mag<0 || x1<0 || x1>L, error("Invalid load row %d",k); end
    if any(typ==["udl","uvl"]) && (x2<=x1 || x2>L)
        error("Invalid distributed-load limits in row %d",k);
    end
    loadData(k)=struct('type',typ,'mag',mag,'x1',x1,'x2',x2);
    coords=[coords x1];
    if any(typ==["udl","uvl"]), coords=[coords x2]; end
end
coords=unique(coords);
x=[];
for k=1:numel(coords)-1
    x=[x coords(k)];
    ne=max(1,ceil((coords(k+1)-coords(k))/(L/80)));
    temp=linspace(coords(k),coords(k+1),ne+1); x=[x temp(2:end)];
end
x=unique([x L]); n=length(x); ndof=2*n;
K=zeros(ndof); F=zeros(ndof,1); EI=E*I;

%% ASSEMBLE BEAM STIFFNESS AND CONSISTENT DISTRIBUTED LOADS
for e=1:n-1
    le=x(e+1)-x(e); kfac=EI/le^3;
    ke=kfac*[12 6*le -12 6*le; 6*le 4*le^2 -6*le 2*le^2; ...
             -12 -6*le 12 -6*le; 6*le 2*le^2 -6*le 4*le^2];
    id=[2*e-1 2*e 2*e+1 2*e+2];
    K(id,id)=K(id,id)+ke;
    xm=mean(x(e:e+1)); q=0;
    for j=1:nLoad
        typ=loadData(j).type; a=loadData(j).x1; b=loadData(j).x2;
        if typ=="udl" && xm>=a && xm<=b
            q=q-loadData(j).mag;
        elseif typ=="uvl" && xm>=a && xm<=b
            q=q-loadData(j).mag*(xm-a)/(b-a);
        end
    end
    fe=[q*le/2; q*le^2/12; q*le/2; -q*le^2/12];
    F(id)=F(id)+fe;
end

%% APPLY POINT LOADS AND APPLIED COUPLES
for j=1:nLoad
    typ=loadData(j).type; a=loadData(j).x1;
    [~,ii]=min(abs(x-a));
    if typ=="point", F(2*ii-1)=F(2*ii-1)-loadData(j).mag; end
    if typ=="moment", F(2*ii)=F(2*ii)-loadData(j).mag; end
end

%% BOUNDARY CONDITIONS
fixedDof=[];
switch beamType
    case "cant"
        fixedDof=[1 2];
    case "fixed"
        fixedDof=[1 2 ndof-1 ndof];
    case {"ss","over"}
        for j=1:numel(supports)
            [~,ii]=min(abs(x-supports(j)));
            fixedDof=[fixedDof 2*ii-1]; % vertical displacement restrained
        end
end
freeDof=setdiff(1:ndof,fixedDof);
U=zeros(ndof,1);
U(freeDof)=K(freeDof,freeDof)\F(freeDof);
R=K*U-F;

%% REPORT SUPPORT REACTIONS
fprintf('\nSUPPORT REACTIONS\n');
for j=1:numel(supports)
    [~,ii]=min(abs(x-supports(j)));
    fprintf('x = %.3f m: V = %.4f kN, M = %.4f kN.m\n', ...
        supports(j),R(2*ii-1),R(2*ii));
end

%% SECTION CUTS: SHEAR AND SAGGING-POSITIVE MOMENT
xs=linspace(0,L,1001); V=zeros(size(xs)); BM=zeros(size(xs));
for i=1:numel(xs)
    z=xs(i);
    for j=1:numel(supports)
        if supports(j)<=z+1e-10
            [~,ii]=min(abs(x-supports(j)));
            V(i)=V(i)+R(2*ii-1);
            BM(i)=BM(i)+R(2*ii-1)*(z-supports(j))-R(2*ii);
        end
    end
    for j=1:nLoad
        typ=loadData(j).type; a=loadData(j).x1; b=loadData(j).x2; P=loadData(j).mag;
        if typ=="point" && a<=z
            V(i)=V(i)-P; BM(i)=BM(i)-P*(z-a);
        elseif typ=="moment" && a<=z
            BM(i)=BM(i)+P;
        elseif typ=="udl" || typ=="uvl"
            zz=min(z,b);
            if zz>a
                xx=linspace(a,zz,101);
                if typ=="udl", q=-P*ones(size(xx));
                else, q=-P*(xx-a)/(b-a); end
                V(i)=V(i)+trapz(xx,q);
                BM(i)=BM(i)+trapz(xx,q.*(z-xx));
            end
        end
    end
end

%% PLOT BEAM, SFD AND BMD
figure('Name','Beam loading','Color','w'); hold on; plot([0 L],[0 0],'k','LineWidth',4);
for j=1:numel(supports)
    plot(supports(j),-0.12,'kv','MarkerFaceColor',[.7 .85 .95],'MarkerSize',10);
    [~,ii]=min(abs(x-supports(j)));
    quiver(supports(j),0,0,R(2*ii-1)*.025,0,'Color',[0 .5 .2],'LineWidth',1.8,'MaxHeadSize',1.5);
end
title('Beam, supports and reactions'); xlabel('Position (m)'); axis tight; grid on;

figure('Name','Shear Force Diagram','Color','w');
area(xs,V,'FaceColor',[.75 .86 .96],'EdgeColor','none'); hold on;
plot(xs,V,'b','LineWidth',1.8); yline(0,'k-'); grid on;
xlabel('Position along beam (m)'); ylabel('Shear force (kN)'); title('Shear Force Diagram');

figure('Name','Bending Moment Diagram','Color','w');
area(xs,BM,'FaceColor',[.77 .91 .80],'EdgeColor','none'); hold on;
plot(xs,BM,'Color',[.05 .48 .22],'LineWidth',1.8); yline(0,'k-'); grid on;
xlabel('Position along beam (m)'); ylabel('Bending moment (kN.m)'); title('Bending Moment Diagram');

fprintf('\nExtrema from sampled section cuts:\n');
fprintf('Max shear = %.4f kN; Min shear = %.4f kN\n',max(V),min(V));
fprintf('Max moment = %.4f kN.m; Min moment = %.4f kN.m\n',max(BM),min(BM));
