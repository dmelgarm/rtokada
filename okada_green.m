function [G]=okada_green(e,n,depth,strike,dip,W,L)

%This is modified from Francois Beauducel's okada85 function available from
%the mathworks file exchange

%OKADA85 Surface deformation due to a finite rectangular source.
%	[uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(E,N,...
%	   DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN,NU)
%	computes Okada's 1985 solution for displacements, tilts and strains in a
%	geographic referential (East,North,Up). Coordinates (0,0) and depth
%	correspond to the fault centroid.
%
%	   E,N    : Coordinates of observation points (relative to fault centroid)
%	   DEPTH  : Depth of the fault centroid (DEP > 0)
%	   STRIKE : Strike-angle from North (in degrees)
%	   DIP    : Dip-angle (in degrees, must be scalar)
%	   LENGTH : Fault length in strike direction (LEN > 0)
%	   WIDTH  : Fault width in dip direction (WIDTH > 0)
%	   RAKE   : Slip-angle direction on fault plane (in degrees)
%	   SLIP   : Dislocation in rake direction (length unit)
%	   OPEN   : Dislocation in tensile component (same unit as SLIP)
%	   NU     : Poisson's ratio (optional, default is 0.25)
%
%	returns the following variables (same matrix size as E and N):
%	   uN,uE,uZ        : Displacements (unit of SLIP and OPEN)
%	   uZE,uZN         : Tilts (in radian)
%	   uNN,uNE,uEN,uEE : Strains (unit of SLIP and OPEN)/(Unit of N,E,..,WIDTH)
%	                     POSITIVE = COMPRESSION
% Default values for optional input arguments
nu = 0.25;	% isotropic Poisson's ratio

% Assigns input arguments
strike = strike.*pi./180;	% converting STRIKE in radian
dip = dip.*pi./180;	% converting DIP in radian ('delta' in Okada's equations)

G=[];
for j = 1:length(e(1,:))
    GR=[];
    for i = 1:length(e(:,1))
        
        % Converts fault coordinates (E,N,DEPTH) relative to centroid
        % into Okada's reference system (X,Y,D)
        d = depth(i,j) + sin(dip(i)).*W(i)/2;	% d is fault's top edge
        ec = e(i,j) + cos(strike(i))*cos(dip(i)).*W(i)/2;
        nc = n(i,j) - sin(strike(i))*cos(dip(i)).*W(i)/2;
        x = cos(strike(i)).*nc + sin(strike(i)).*ec + L(i)/2;
        y = sin(strike(i)).*nc - cos(strike(i)).*ec + cos(dip(i)).*W(i);
        
        % Variable substitution (independent from xi and eta)
        p = y.*cos(dip(i)) + d.*sin(dip(i));
        q = y.*sin(dip(i)) - d.*cos(dip(i));
        
        % Displacements
        g1 = -1/(2*pi) * chinnery(@ux_ss,x,p,L(i),W(i),q,dip(i),nu); % strike-slip, x
        g2 = -1/(2*pi) * chinnery(@ux_ds,x,p,L(i),W(i),q,dip(i),nu); % dip-slip, x
        
        g3 = -1/(2*pi) * chinnery(@uy_ss,x,p,L(i),W(i),q,dip(i),nu);% strike-slip, y
        g4 = -1/(2*pi) * chinnery(@uy_ds,x,p,L(i),W(i),q,dip(i),nu); % dip-slip, y
        
        g5 = -1/(2*pi) * chinnery(@uz_ss,x,p,L(i),W(i),q,dip(i),nu); % strike-slip, z
        g6 = -1/(2*pi) * chinnery(@uz_ds,x,p,L(i),W(i),q,dip(i),nu); % dip-slip, z
        
        % Rotation from Okada's axes to geographic
        g1n = sin(strike(i)).*g1 - cos(strike(i)).*g3;
        g3n = cos(strike(i)).*g1 + sin(strike(i)).*g3;
        
        g2n = sin(strike(i)).*g2 - cos(strike(i)).*g4;
        g4n = cos(strike(i)).*g2 + sin(strike(i)).*g4;
        
        GA = [g1n g2n;g3n g4n;g5 g6];
        GR = [GR GA];
        
    end
    G = [G;GR];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Notes for I... and K... subfunctions:
%
%	1. original formulas use Lame's constant as mu/(mu+lambda) which
%	   depends only on the Poisson's ratio = 1 - 2*nu
%	2. tests for cos(dip) == 0 are made with "cos(dip) > eps"
%	   because cos(90*pi/180) is not zero but = 6.1232e-17 (!)
%	   NOTE: don't use cosd and sind because of incompatibility
%	   with Matlab v6 and earlier...


% =================================================================
% Chinnery's notation [equation (24) p. 1143]

% -----------------------------------------------------------------
function u=chinnery(f,x,p,L,W,q,dip,nu)
u = feval(f,x,p,q,dip,nu) ...
    - feval(f,x,p-W,q,dip,nu) ...
    - feval(f,x-L,p,q,dip,nu) ...
    + feval(f,x-L,p-W,q,dip,nu);


% =================================================================
% Displacement subfunctions

% strike-slip displacement subfunctions [equation (25) p. 1144]

% -----------------------------------------------------------------
function u=ux_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = xi.*q./(R.*(R + eta)) ...
    + atan(xi.*eta./(q.*R)) ...
    + I1(xi,eta,q,dip,nu,R).*sin(dip);

% -----------------------------------------------------------------
function u=uy_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = (eta*cos(dip) + q*sin(dip)).*q./(R.*(R + eta)) ...
    + q.*cos(dip)./(R + eta) ...
    + I2(eta,q,dip,nu,R).*sin(dip);

% -----------------------------------------------------------------
function u=uz_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sin(dip) - q*cos(dip);
u = (eta*sin(dip) - q*cos(dip)).*q./(R.*(R + eta)) ...
    + q.*sin(dip)./(R + eta) ...
    + I4(db,eta,q,dip,nu,R).*sin(dip);

% dip-slip displacement subfunctions [equation (26) p. 1144]
% -----------------------------------------------------------------
function u=ux_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = q./R ...
    - I3(eta,q,dip,nu,R).*sin(dip).*cos(dip);

% -----------------------------------------------------------------
function u=uy_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = (eta*cos(dip) + q*sin(dip)).*q./(R.*(R + xi)) ...
    + cos(dip).*atan(xi.*eta./(q.*R)) ...
    - I1(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);

% -----------------------------------------------------------------
function u=uz_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta*sin(dip) - q*cos(dip);
u = db.*q./(R.*(R + xi)) ...
    + sin(dip).*atan(xi.*eta./(q.*R)) ...
    - I5(xi,eta,q,dip,nu,R,db).*sin(dip).*cos(dip);


% I... displacement subfunctions [equations (28) (29) p. 1144-1145]
% -----------------------------------------------------------------
function I=I1(xi,eta,q,dip,nu,R)
db = eta*sin(dip) - q*cos(dip);
if cos(dip) > eps
    I = (1 - 2*nu) * (-xi./(cos(dip).*(R + db))) - sin(dip)./cos(dip) .*I5(xi,eta,q,dip,nu,R,db);
else
    I = -(1 - 2*nu)/2 * xi.*q./(R + db).^2;
end

% -----------------------------------------------------------------
function I=I2(eta,q,dip,nu,R)
I = (1 - 2*nu) * (-log(R + eta)) - I3(eta,q,dip,nu,R);

% -----------------------------------------------------------------
function I=I3(eta,q,dip,nu,R)
yb = eta*cos(dip) + q*sin(dip);
db = eta*sin(dip) - q*cos(dip);
if cos(dip) > eps
    I = (1 - 2*nu) * (yb./(cos(dip)*(R + db)) - log(R + eta)) ...
        + sin(dip)./cos(dip) * I4(db,eta,q,dip,nu,R);
else
    I = (1 - 2*nu)/2 * (eta./(R + db) + yb.*q./(R + db).^2 - log(R + eta));
end

% -----------------------------------------------------------------
function I=I4(db,eta,q,dip,nu,R)
if cos(dip) > eps
    I = (1 - 2*nu) * 1/cos(dip) * (log(R + db) - sin(dip)*log(R + eta));
else
    I = -(1 - 2*nu) * q./(R + db);
end

% -----------------------------------------------------------------
function I=I5(xi,eta,q,dip,nu,R,db)
X = sqrt(xi.^2 + q.^2);
if cos(dip) > eps
    I = (1 - 2*nu) * 2./cos(dip) ...
        .* atan((eta.*(X + q.*cos(dip)) + X.*(R + X).*sin(dip))./(xi.*(R + X).*cos(dip)));
else
    I = -(1 - 2*nu) * xi.*sin(dip)./(R + db);
end


