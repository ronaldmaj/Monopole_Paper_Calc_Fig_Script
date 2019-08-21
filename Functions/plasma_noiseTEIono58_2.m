function [V2, Integrand, w_vec] = plasma_noiseTEIono58_2(l,N,T,a,...
    fmin,fmax,n,dipole,z_max,z_min)

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
vT = sqrt(2*kB*T/me);   % Thermal velocity of electrons
LD = sqrt((kB*T*eps0)/(N*e^2));
wp = vT/(sqrt(2)*LD); % Plasma frequency

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

Zdash = @(x) -2*(1+x.*plasmaZ(x));

% w/wpC dependent calculations
f = logspace(fmin, fmax, n);
w_vec = f*(2*pi); % Determine omega from omega/omega_p vector
w_div_wp_vec = w_vec./wp; % Set omega/omega_p vector

% Vector for recording voltage power spectrum
sz = size(w_div_wp_vec);
Integrand = zeros(sz);


for i = 1:length(w_div_wp_vec)
    
   disp([int2str(i) ' out of ' int2str(length(w_div_wp_vec)) ' complete.'])

    alpha = w_div_wp_vec(i);
    w = w_vec(i); % Set specific omega value
    k = @(x) w./(x.*vT); 

% Permittivities
    epsL = @(x) -x.^2/alpha^2.*Zdash(x);
    DL = @(x) 1 + epsL(x);

% Integral calculation
    Int = @(z) (F1(k(z)*l).*J(k(z)*a).^2./(z.^2)).*abs(imag(1./(DL(z))));
    zmin = z_min;
    zmax = z_max;                 
    Integrand(i) = integral(Int,zmin,zmax); 

    
end

V2 = Integrand*(4*kB*T*(4/((pi^2)*eps0*vT))); % Convert from unitless to
                                            % units of V^2/Hz

%figure; plot(f,V2,'k')
%set(gca,'XScale','log','YScale','log')
%xlabel('Frequency (Hz)'), ylabel('Noise Spectral Density V_\omega^2 (V^2 Hz^{-1})')

end