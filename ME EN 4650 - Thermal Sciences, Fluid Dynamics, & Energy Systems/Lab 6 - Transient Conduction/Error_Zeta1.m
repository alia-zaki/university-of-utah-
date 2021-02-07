function sigma_z = Error_Zeta1(z,Fo,theta)
%*************************************************************************
%function sigma_z = Error_Zeta1(z,Fo,theta)
%
% INPUTS
%    z            best-fit parameter zeta1 (1x1, scalar)
%    Fo           Fourier number for the curve-fit data (Nx1, vector)
%    theta        theta^* values for the curve-fit data (Nx1, vector)
%
% OUTPUT
%    sigma_z      error in z (1x1, scalar)
%
% 
% This function calculates the uncertainty in the curve fit parameter
% zeta1 obtained from a nonlinear least squares curve fit to the 
% one-term approximation of the exact solution to the unsteady heat
% conduction problem for a sphere with convective boundary conditions.
%
% The error in the curve-fit parameter z is given by
%
%    sigma_z = sqrt(inv(X'*X)*s^2)
%
% where s^2 is the mean square error of data relative to the model, and
% X is the Jacobian, which is simply the derivative of the total curve-
% fit error with respect to z. The total curve-fit error is the difference
% between the data and the model,
% 
%    e = theta - (4*(sin(z)-z*cos(z))/(2*z-sin(2*z))*exp(-z^2*Fo))
%
% Therefore, X=de/dz. The derivative was calculated using Matlab's symbolic
% toolbox. The mean square error is
% 
%    s^2 = sum(e^2)/(N-1),
%
% where N is the total number of data points and also equal to the length
% of Fo. 
%
% M Metzger
% 2-2020
%*************************************************************************

% check that Fo is a column vector
[r,c]=size(Fo);
if (c~=1), Fo=Fo'; end

% check that theta is a column vector
[r,c]=size(theta);
if (c~=1), theta=theta'; end

% number of data points
N=length(Fo);

% Jacobian - from Matlab's symbolic toolbx
X= (exp(-Fo*z^2)*(4*sin(z) - 4*z*cos(z))*(2*cos(2*z) - 2))/...
   (2*z - sin(2*z))^2 + (4*z*exp(-Fo*z^2)*sin(z))/(2*z - sin(2*z)) - ...
   (2*Fo*z.*exp(-Fo*z^2)*(4*sin(z) - 4*z*cos(z)))/(2*z - sin(2*z));

% transpose of Jacobian
XT=X';

% X'*X: since X is a column array, this operation will produce a scalar
XTX=X'*X;

% inv(X'*X): since X is a column array, inv(X'*X) is simply 1/(X'*x)
inv_XTX = 1/XTX;

% curve fit error: difference between data and model
e = theta - (4*(sin(z)-z*cos(z))/(2*z-sin(2*z))*exp(-z^2*Fo));

% mean squared error
s2=sum(e.^2)/(N-1);

% error in z
sigma_z = sqrt(inv_XTX*s2);

end

