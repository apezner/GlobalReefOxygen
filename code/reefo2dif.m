%-------------------------------------
% function yp = reefo2dif(t,y)
% provides difequations for reef oxygen assessment
%-------------------------------------
function yp = reefo2dif(t,y)

% Note that surface variable F is denoted by y
global oceanO2 reefvol tau resp

%Residence time (tau) and input/output flow
%tau=1;              %hours
Flow=reefvol/tau;   %L/h

%Input (Fi) and output (Fo) fluxes
Fi=oceanO2*1e-06*Flow;          %Oxygen input mol/h
Fo=y(1)./reefvol*Flow;          %

%Photosynthesis
Pf=0.04;            %Max rate: 40 mmol/m2/h
PS=Pf*(sin(pi*t/12).^1.2);

if PS>0    
    Fp=PS;
else
    Fp=0;
end

%Respiration
Fr=resp;                  %12 mmol/m2/hr in base scenario

% differential equations
yp = [Fi-Fo+Fp-Fr];
return;
  