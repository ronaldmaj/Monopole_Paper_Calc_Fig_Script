function [V2, Integrand, w_vec] = plasma_noiseTEIono58_4(l,N,T,a,...
    fmin,fmax,n,dipole,z_max,z_min,kap)

% Alterations to TEIono58 to include the option of using a kappa
% distribution instead of a Maxwellian. This kappa code is adapted
% from the work of Elias Odelstad: 
% https://github.com/eliasodelstad/irfuproj_JUICE_noise

% Code for calculating the voltage power spectrum similar to that done in
% the Meyer-Vernet/Perche 1989 Paper "Tool Kit for Antennae and Thermal
% Noise Near the Plasma Frequency". Produces a plot of V^2/sqrt(T_e) at the
% end, and outputs this voltage spectrum and the frequency range.

% Input: 
% 

% Calculates the voltage power spectrum V_omega^2 and voltage spectrum V 
% taking in as input the length of the antenna "l". This code differs from 
% the function plasma_noiseCout1.m in that parameters are specifically set 
% for the ionosphere (T = 2000K, LD = 2e-3), processes to increase speed 
% of calculation are incorporated and the spectrum for voltage is also 
% determined. Details are given below:

% Physical constants
kB = 1.3806488e-23;     % Boltzmann constant
me = 9.109383e-31;      % Mass of electron
e = 1.60217657e-19;     % Electron charge
eps0 = 8.85418782e-12;  % Vacuum Permittivity 

% Plasma parameters
wp = sqrt((N*e^2)/(me*eps0));

vT = sqrt(2*kB*T/me);   % Thermal velocity of electrons (maxwellian)
v0 = @(kap) vT*sqrt((2*kap - 3)/(2*kap));
LD = sqrt((kB*T*eps0)/(N*e^2));
LDk = @(kap) (v0(kap)./wp).*(kap/(2*kap-1));
u = l/LDk(kap);
%wp = vT/(sqrt(2)*LD); % Plasma frequency

% Functions in integrand
if dipole == 1
    F1 = @(x) ((sinint(x))-(sinint(2*x)./2)-((2*((sin(x/2)).^4))./x))./x;
    J = @(x) besselj(0,x);
elseif dipole == 0
    % Assuming wire monopole
    F1 = @(x) 0.125*(x.^2 + 2*cos(x) - 2)./(x.^2);
    J = @(x) besselj(0,x);
else
    F1 = @(x) 0.25*(1 - sin(x)./x);
    J = @(x) (sin(x).^2)./(x.^2);
end

% w/wpC dependent calculations
f = logspace(fmin, fmax, n);
w_vec = f*(2*pi); % Determine omega from omega/omega_p vector
w_wp_vec = w_vec./wp;
% Vector for recording voltage power spectrum
sz = size(w_wp_vec);
Integrand = zeros(sz);
V2 = zeros(sz);

for i = 1:length(w_wp_vec)
    
    % Display calculation 
    disp([int2str(i) ' out of ' int2str(length(w_vec)) ' complete.'])
    w = w_vec(i); % Set specific omega value 
    r = w_wp_vec(i);
    
    % Make z and p same size matrices
    P = @(z,p) p'*ones(1,length(z));
    Z = @(z,p) ones(length(p),1)*z;
    % Implement summation term
    S = @(z,p) factorial(kap + P(z,p))./(factorial(P(z,p)).*(2*1i)...
        .^(kap+1+P(z,p)).*(Z(z,p)+1i).^(kap+1-P(z,p)));
    % Dielectric permittivity
    
    % Epsilon_L
    epsL = @(z) 1 + ((z.^2)./(r^2)).*...
    (2*kap - 1 + ((-2).^(kap+1))./prod(1:2:(2*kap-3)).*...
    (1i*z).*sum(S(z,0:kap)));
    % Integral calculation
    C = (2^(kap+3)/(pi^2*eps0))*...
        (factorial(kap)/(prod(1:2:2*kap-3)*sqrt(kap)))*(me*v0(kap)/r^2);
    Int = @(z) z.*F1(r*u./(z*sqrt(2*kap-1))).*((1+z.^2).^kap.*abs(epsL(z)).^2).^(-1);            
    Integrand(i) = integral(Int,z_min,z_max); 
    V2(i) = C*Integrand(i);
    
end



end