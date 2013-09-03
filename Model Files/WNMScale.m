 % script created by Richard Balson 14/10/2012

% description
% ~~~~~~~~~~~
% This function describes the Wendling neural mass model with 8 differential equations and its
% observation function, the simplified version

% last edit
% ~~~~~~~~~

% Yout function defined correctly

% next edit
% ~~~~~~~~~

% Beginning of function
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [Xout Yout] = WNM(Xin,dt,CoV, gain, tcon,C,Scaling)
Xout = Xin;
zs1 =  Xin(2,:)-C(4)*Xin(3,:)-Xin(4,:);
[Xout(1,:) Xout(5,:)] = PSPkernel([Xin(1,:); Xin(5,:)],dt,gain(1,:),tcon(1),sigmoid(zs1));
zs2 = C(1)*Xin(1,:);
[Xout(2,:) Xout(6,:)] = PSPkernel([Xin(2,:); Xin(6,:)],dt,gain(1,:),tcon(1),CoV +C(2)*sigmoid(zs2));
zs3 = C(3)*Xin(1,:);
[Xout(3,:) Xout(7,:)] = PSPkernel([Xin(3,:); Xin(7,:)],dt,gain(2,:),tcon(2),sigmoid(zs3));
zs4 = C(5)*Xin(1,:)-C(6)*Xin(3,:);
[Xout(4,:) Xout(8,:)] = PSPkernel([Xin(4,:); Xin(8,:)],dt,gain(3,:),tcon(3),C(7)*sigmoid(zs4));


Xout(9,:) = Scaling(1,:).*(Xout(2,:)-C(4)*Xout(3,:)-Xout(4,:))+Scaling(2,:); % Observation function

Yout = Xout(9,:);
end

% Parameter specifiction
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Units
% ~~~~~~~~

% Equations
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% End of function description    


