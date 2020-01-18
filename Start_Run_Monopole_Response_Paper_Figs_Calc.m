%% Calculations and Figures for the Monopole vs Dipole Response Paper
%
% This script will run to produce all the figures and calculations found
% within Maj R., Cairns I.H., and Martinovic M.M., (2019) "Impedance and 
% voltage power spectra of a monopole antenna in a warm plasma - derivation
% and application to CubeSats". 
%
% The script in its current state loads many of the results of calculations
% to increase the run time of the code (aim being to produce the figures as
% quickly as possible). To run the calculations from scratch, please
% comment out "load('Monopole_paper_QTN_data2.mat')" and 
% "load('Monopole_paper_ShotNoise_Z_data_updated2.mat')" and uncomment the 
% associated code following these.
%
% The code used for running the kappa distribution is adapted from the work
% of Elias Odelstad in his work "Noise sources in the electric field 
% antenna on the ESA JUICE satellite" (2013). His code can be found at 
% https://github.com/eliasodelstad/irfuproj_JUICE_noise

clc
clear all

addpath('Data')
addpath('Functions')

% Physical constants
kB = 1.3806488e-23;     % Boltzmann constant
me = 9.109383e-31;      % Mass of electron
e = 1.60217657e-19;     % Electron charge
eps0 = 8.85418782e-12;  % Vacuum Permittivity

% Load data
load('Ne300D.mat');
load('Te300D.mat');

% Averages
n_e = mean(mean(Ne300D));
T_e = mean(mean(Te300D));

% Antenna parameters
a = 2e-4; 
l = 0.3;

% Plasma parameters
LD = sqrt((kB*T_e*eps0)/(n_e*e^2));
vT = sqrt(2*kB*T_e/me);
wp = vT/(sqrt(2)*LD);
fp = wp/(2*pi);
v0 = @(kap) vT*sqrt((2*kap - 3)/2*kap);
LDk = @(kap) (v0(kap)./wp).*(kap/(2*kap-1));

% Frequency range
f = logspace(5,8,400);
w_vec = 2*pi*f;

% Choose specific frequency to work with:
w = w_vec(100);
% Ratio of frequency omega to plasma frequency:
w_wp_ex = w/wp;

% Define faddeeva
W = @(x) faddeeva1(x);

% Definition of z variable
z = @(x) w./(x.*vT);
z_kap = @(x,kap) w./(x.*v0(kap)*sqrt(kap));

% D_L or eps_L - dielectric permittivity tensor
epsL = @(k) (1 + ((1 + (1i)*sqrt(pi)*z(k).*W(z(k)))./((k.^2)*(LD^2))));

% For a kappa distribution, we need definitions for r:
r = @(w) w/wp; 
% Odelstad's implementation for kappa distribution:
P = @(z,p) p'*ones(1,length(z));
Z = @(z,p) ones(length(p),1)*z;
% Implement summation term
S = @(k,kap,p) factorial(kap + P(z_kap(k,kap),p))./...
    (factorial(P(z_kap(k,kap),p)).*(2*1i).^(kap+1+P(z_kap(k,kap),p)).*...
    (Z(z_kap(k,kap),p)+1i).^(kap+1-P(z_kap(k,kap),p)));
% Dielectric permittivity
% Epsilon_L
epsL_kap = @(k, kap, w) 1 + ((z_kap(k,kap).^2)/(r(w).^2)).*...
(2*kap - 1 + ((-2).^(kap+1))./prod(1:2:(2*kap-3)).*...
(1i*z(k)).*sum(S(k,kap,0:kap)));

% Antenna response functions
% Dipole wires:
F1 = @(x) ((sinint(x))-(sinint(2*x)./2)-((2*((sin(x/2)).^4))./x))./x;
% Dipole spheres:
F2 = @(x) 0.25*(1 - sin(x)./x);
% Monopole:
F3 = @(x) 0.125*(x.^2 + 2*cos(x) - 2)./(x.^2);
% New monopole function - Fcos:
%Fcos = @(x) (1/16)* (2*pi*cos(x) - 2*pi*sin(x)./x +...
%    (3*pi^3*atanh(2*x./pi) - 6*pi^2*x - 8*x.^3)./(24*x));

% A k-range for plotting:
k = logspace(-4,4,5000);

%% Figure 1 - Plot the antenna response functions against k
Fig1=figure('pos',[5 30 700 750]);
P1 = plot(k,F1(k),k,F3(k),k,F2(k),'LineWidth',1);
set(gca,'XScale','log','YScale','linear','FontSize',15);
ylim([-0.06 0.31])
xlabel('x','FontSize',16,'FontWeight','bold');
ylabel('F(x)','FontSize',16,'FontWeight','bold');
Le1 = legend('F_1(x) - Dipole (wire)',...
    'F_{m1}(x) - Monopole',...
    'F_{s1}(x) - Dipole (spheres)',...
    'location','best');
%saveas(Fig1,'Fig1-Di_vs_Mono_Response_Fn5','epsc')

%% Figure 2a - Plot the integrand for F1(x) and F3(x) - Maxwellian

k2 = logspace(-2,10,5000);
Integ_Di = (F1(k2*l).*(besselj(0,k2*a).^2))./(epsL(k2));
Integ_Mono = (F3(k2*l).*(besselj(0,k2*a).^2))./(epsL(k2));

Fig2a=figure('pos',[20 -150 1400 750]);

subplot(1,2,1);
P2a = plot(k2,real(Integ_Di),k2,real(Integ_Mono),k2,imag(Integ_Di),...
k2,imag(Integ_Mono),'LineWidth',1);
[~, Idx_Max] = max(Integ_Di);
%vertline2a = xline(k2(Idx_Max),'LineWidth',1.5,'LineStyle',':','Color','b');
P2a(3).LineStyle = '--';
P2a(4).LineStyle = '--';
ax2a_1 = gca;
ax2a_1.XScale = 'log';
ax2a_1.YScale = 'linear';
ax2a_1.XLim = [1e0 1e7];
ax2a_1.YLim = [-6e-3 10e-3];
ax2a_1.FontSize = 15;
ax2a_1.YAxis.Exponent = 0;
box off
ax2a_2 = axes('Position',ax2a_1.Position,...
'XAxisLocation','top',...
'YAxisLocation','right',...
'XScale','log',...
'Color','None');
hold on
P2a_2 = plot(...
    k2*LD,real(Integ_Di),...
    k2*LD,real(Integ_Mono),...
    k2*LD,imag(Integ_Di),...
    k2*LD,imag(Integ_Mono),...
    'Parent',ax2a_2,'LineWidth',1);
%vertline2a = xline(k2(Idx_Max),'LineWidth',1);
P2a_2(3).LineStyle = '--';
P2a_2(4).LineStyle = '--';
ax2a_2.FontSize = 15;
ax2a_2.XLim = [1e0*LD 1e7*LD];
ax2a_2.YLim = [-6e-3 10e-3];
ax2a_2.YAxis.Exponent = 0;
ax2a_2.YTickLabel = [];
Le2a = legend('Dipole-wire Z Integrand - Real',...
'Monopole Z Integrand - Real',...
'Dipole-wire Z Integrand - Imag',...
'Monopole Z Integrand - Imag');
ax2a_1.XLabel.String = 'Wavenumber (k)';
ax2a_1.XLabel.FontSize = 16;
ax2a_1.XLabel.FontWeight = 'bold';
ax2a_1.YLabel.String = 'Integrand';
ax2a_1.YLabel.FontSize = 16;
ax2a_1.YLabel.FontWeight = 'bold';
ax2a_1.Title.String = 'k\lambda_D';

subplot(1,2,2);
P2b = loglog(k2,real(Integ_Di),k2,real(Integ_Mono),...
'LineWidth',1);
ax2b_1 = gca;
ax2b_1.XScale = 'log';
ax2b_1.YScale = 'log';
ax2b_1.XLim = [1e0 1e7];
ax2b_1.YLim = [1e-13 1e0];
ax2b_1.FontSize = 15;
ax2b_1.YAxis.Exponent = 0;
box off
ax2b_2 = axes('Position',ax2b_1.Position,...
'XAxisLocation','top',...
'YAxisLocation','right',...
'XScale','log',...
'YScale','log',...
'Color','None');
hold on
P2b_2 = loglog(...
    k2*LD,real(Integ_Di),...
    k2*LD,real(Integ_Mono),...
    k2*LD,((k2*LD).^-1),...
    k2*LD,((k2*LD).^-2),...
    'Parent',ax2b_2,'LineWidth',1);
P2b_2(3).LineStyle = '-.';
P2b_2(4).LineStyle = '-.';
P2b_2(3).LineWidth = 2;
P2b_2(4).LineWidth = 2;
P2b_2(3).Color = [0.4660, 0.6740, 0.1880];
P2b_2(4).Color = [0.3010, 0.7450, 0.9330];
ax2b_2.FontSize = 15;
ax2b_2.XLim = [1e0*LD 1e7*LD];
ax2b_2.YLim = [1e-13 1e0];
ax2b_2.YTickLabel = [];
Le2b = legend('Dipole-wire Z Integrand - Real',...
'Monopole Z Integrand - Real',...
'(k\lambda_D)^{-1}','(k\lambda_D)^{-2}',...
'Location','best');
ax2b_1.XLabel.String = 'Wavenumber (k)';
ax2b_1.XLabel.FontSize = 16;
ax2b_1.XLabel.FontWeight = 'bold';
ax2b_1.YLabel.String = 'Integrand';
ax2b_1.YLabel.FontSize = 16;
ax2b_1.YLabel.FontWeight = 'bold';
ax2b_1.Title.String = 'k\lambda_D';
%saveas(Fig2a,'Fig2a-Di_vs_Mono_Integrand_new5','epsc')

%% Figure 2b - Plot the integrand for F1(x) and F3(x) - kappa
kap = 4;
Integ_Di_kap = (F1(k2*l).*(besselj(0,k2*a).^2))./(epsL_kap(k2,kap,w));
Integ_Mono_kap = (F3(k2*l).*(besselj(0,k2*a).^2))./(epsL_kap(k2,kap,w));

Fig2b=figure('pos',[20 -150 1400 750]);
[~, Idx_kap] = max(Integ_Di_kap);

subplot(1,2,1);
P2a2 = plot(...
    k2,real(Integ_Di_kap),...
    k2,real(Integ_Mono_kap),...
    k2,imag(Integ_Di_kap),...
    k2,imag(Integ_Mono_kap),...
    [k2(Idx_Max) k2(Idx_Max)],[-2.5e-2 0.02],...
    [k2(Idx_kap) k2(Idx_kap)],[-2.5e-2 0.02],...
    'LineWidth',1);
P2a2(3).LineStyle = '--';
P2a2(4).LineStyle = '--';
P2a2(5).LineStyle = ':';
P2a2(6).LineStyle = ':';
P2a2(5).LineWidth = 2;
P2a2(6).LineWidth = 2;
ax2a2_1 = gca;
ax2a2_1.XScale = 'log';
ax2a2_1.YScale = 'linear';
ax2a2_1.XLim = [1e0 1e7];
ax2a2_1.YLim = [-2.5e-2 0.02];
ax2a2_1.YAxis.Exponent = 0;
ax2a2_1.FontSize = 15;
box off
ax2a2_2 = axes('Position',ax2a2_1.Position,...
'XAxisLocation','top',...
'YAxisLocation','right',...
'XScale','log',...
'Color','None');
hold on
P2a2_2 = plot(...
    k2*LDk(kap),real(Integ_Di_kap),...
    k2*LDk(kap),real(Integ_Mono_kap),...
    k2*LDk(kap),imag(Integ_Di_kap),...
    k2*LDk(kap),imag(Integ_Mono_kap),...
    [k2(Idx_Max)*LDk(kap) k2(Idx_Max)*LDk(kap)],[-2.5e-2 0.02],...
    [k2(Idx_kap)*LDk(kap) k2(Idx_kap)*LDk(kap)],[-2.5e-2 0.02],...    
    'Parent',ax2a2_2,'LineWidth',1);
P2a2_2(3).LineStyle = '--';
P2a2_2(4).LineStyle = '--';
P2a2_2(5).LineStyle = ':';
P2a2_2(6).LineStyle = ':';
P2a2_2(5).LineWidth = 2;
P2a2_2(6).LineWidth = 2;
ax2a2_2.FontSize = 15;
ax2a2_2.XLim = [1e0*LDk(kap) 1e7*LDk(kap)];
ax2a2_2.YLim = [-2.5e-2 0.02];
ax2a2_2.YAxis.Exponent = 0;
ax2a2_2.YTickLabel = [];
Le2a = legend('Dipole-wire Z Integrand - Real',...
'Monopole Z Integrand - Real',...
'Dipole-wire Z Integrand - Imag',...
'Monopole Z Integrand - Imag',...
'Location of Maxwellian Peak',...
'Location of Kappa Peak');
ax2a2_1.XLabel.String = 'Wavenumber (k)';
ax2a2_1.XLabel.FontSize = 16;
ax2a2_1.XLabel.FontWeight = 'bold';
ax2a2_1.YLabel.String = 'Integrand';
ax2a2_1.YLabel.FontSize = 16;
ax2a2_1.YLabel.FontWeight = 'bold';
ax2a2_1.Title.String = 'k\lambda_{D-\kappa}';

subplot(1,2,2);
P2b2 = loglog(...
    k2,real(Integ_Di_kap),...
    k2,real(Integ_Mono_kap),...
    k2,imag(Integ_Di_kap),...
    k2,imag(Integ_Mono_kap),'LineWidth',1);
P2b2(3).LineStyle = '--';
P2b2(4).LineStyle = '--';
ax2b2_1 = gca;
ax2b2_1.XScale = 'log';
ax2b2_1.YScale = 'log';
ax2b2_1.XLim = [1e0 1e7];
ax2b2_1.YLim = [1e-13 1e0];
ax2b2_1.FontSize = 15;
ax2b2_1.YAxis.Exponent = 0;
box off
ax2b2_2 = axes('Position',ax2b2_1.Position,...
'XAxisLocation','top',...
'YAxisLocation','right',...
'YScale','log',...
'XScale','log',...
'Color','None');
hold on
P2b2_2 = loglog(...
    k2*LDk(kap),real(Integ_Di_kap),...
    k2*LDk(kap),real(Integ_Mono_kap),...
    k2*LDk(kap),((k2*LDk(kap)).^-1),...
    k2*LDk(kap),((k2*LDk(kap)).^-2),...
    'Parent',ax2b2_2,'LineWidth',1);
P2b2_2(3).LineStyle = '-.';
P2b2_2(4).LineStyle = '-.';
P2b2_2(3).LineWidth = 2;
P2b2_2(4).LineWidth = 2;
P2b2_2(3).Color = [0.4660, 0.6740, 0.1880];
P2b2_2(4).Color = [0.3010, 0.7450, 0.9330];
ax2b2_2.FontSize = 15;
ax2b2_2.XLim = [1e0*LDk(kap) 1e7*LDk(kap)];
ax2b2_2.YLim = [1e-13 1e0];
ax2b2_2.YTickLabel = [];
ax2b2_2.YAxis.Exponent = 0;
Le2b_2 = legend('Dipole-wire Z Integrand - Real',...
'Monopole Z Integrand - Real',...
'(k\lambda_{D-\kappa})^{-1}',...
'(k\lambda_{D-\kappa})^{-2}',...
'Location','best');
ax2b2_1.XLabel.String = 'Wavenumber (k)';
ax2b2_1.XLabel.FontSize = 16;
ax2b2_1.XLabel.FontWeight = 'bold';
ax2b2_1.YLabel.String = 'Integrand';
ax2b2_1.YLabel.FontSize = 16;
ax2b2_1.YLabel.FontWeight = 'bold';
ax2b2_1.Title.String = 'k\lambda_{D-\kappa}';
%saveas(Fig2b,'Fig2b-Di_vs_Mono_Integrand_kap5_peaks','epsc')

%% Figure 2c - Plot the integrand for F1(x) and F3(x) - Maxwellian/spheres

Integ_Disph = (F2(k2*l).*(besselj(0,k2*a).^2))./(epsL(k2));

Fig2c=figure('pos',[20 -150 1400 750]);

subplot(1,2,1);
P2a3 = plot(...
    k2,real(Integ_Disph),...
    k2,real(Integ_Mono),...
    k2,imag(Integ_Disph),...
    k2,imag(Integ_Mono),...
    'LineWidth',1);
P2a3(3).LineStyle = '--';
P2a3(4).LineStyle = '--';
ax2a3_1 = gca;
ax2a3_1.XScale = 'log';
ax2a3_1.YScale = 'linear';
ax2a3_1.XLim = [1e0 1e7];
ax2a3_1.YLim = [-10e-3 0.145];
ax2a3_1.FontSize = 15;
ax2a3_1.YAxis.Exponent = 0;
box off
ax2a3_2 = axes('Position',ax2a3_1.Position,...
'XAxisLocation','top',...
'YAxisLocation','right',...
'XScale','log',...
'Color','None');
hold on
P2a3_2 = plot(...
    k2*LD,real(Integ_Disph),...
    k2*LD,real(Integ_Mono),...
    k2*LD,imag(Integ_Disph),...
    k2*LD,imag(Integ_Mono),...
    'Parent',ax2a3_2,'LineWidth',1);
P2a3_2(3).LineStyle = '--';
P2a3_2(4).LineStyle = '--';
ax2a3_2.FontSize = 15;
ax2a3_2.XLim = [1e0*LD 1e7*LD];
ax2a3_2.YLim = [-10e-3 0.145];
ax2a3_2.YAxis.Exponent = 0;
ax2a3_2.YTickLabel = [];
Le3 = legend('Dipole-spheres Z Integrand - Real',...
'Monopole Z Integrand - Real',...
'Dipole-spheres Z Integrand - Imag',...
'Monopole Z Integrand - Imag');
ax2a3_1.XLabel.String = 'Wavenumber (k)';
ax2a3_1.XLabel.FontSize = 16;
ax2a3_1.XLabel.FontWeight = 'bold';
ax2a3_1.YLabel.String = 'Integrand';
ax2a3_1.YLabel.FontSize = 16;
ax2a3_1.YLabel.FontWeight = 'bold';
ax2a3_1.Title.String = 'k\lambda_D';

subplot(1,2,2);
P2b3 = loglog(...
    k2,real(Integ_Disph),...
    k2,real(Integ_Mono),...
    'LineWidth',1);
ax2b3_1 = gca;
ax2b3_1.XScale = 'log';
ax2b3_1.YScale = 'log';
ax2b3_1.XLim = [1e0 1e7];
ax2b3_1.YLim = [1e-13 1e0];
ax2b3_1.FontSize = 15;
ax2b3_1.YAxis.Exponent = 0;
box off
ax2b3_2 = axes('Position',ax2b3_1.Position,...
'XAxisLocation','top',...
'YAxisLocation','right',...
'XScale','log',...
'YScale','log',...
'Color','None');
hold on
P2b3_2 = loglog(...
    k2*LD,real(Integ_Disph),...
    k2*LD,real(Integ_Mono),...
    k2*LD,((k2*LD).^-1),...
    'Parent',ax2b3_2,'LineWidth',1);
P2b3_2(3).LineStyle = '-.';
P2b3_2(3).LineWidth = 2;
P2b3_2(3).Color = [0.4660, 0.6740, 0.1880];
ax2b3_2.FontSize = 15;
ax2b3_2.XLim = [1e0*LD 1e7*LD];
ax2b3_2.YLim = [1e-13 1e0];
ax2b3_2.YTickLabel = [];
Le3b = legend('Dipole-spheres Z Integrand - Real',...
'Monopole Z Integrand - Real',...
'(k\lambda_D)^{-1}',...
'Location','best');
ax2b3_1.XLabel.String = 'Wavenumber (k)';
ax2b3_1.XLabel.FontSize = 16;
ax2b3_1.XLabel.FontWeight = 'bold';
ax2b3_1.YLabel.String = 'Integrand';
ax2b3_1.YLabel.FontSize = 16;
ax2b3_1.YLabel.FontWeight = 'bold';
ax2b3_1.Title.String = 'k\lambda_D';
%saveas(Fig2c,'Fig2c-Disph_vs_Mono_Integrand_new5','epsc')

%% Figure 3 - Plot QTN for wire, sphere and mono using F1(x), F2(x) & F3(x)

% Load QTN Data:
load('Monopole_paper_QTN_data2.mat') 
% If you want to run the calculations from scratch, comment out above line
% and uncomment line below:
%Run_QTN_Spectra_Calculations

% Use the f/fp vector for figures:
f = f./fp;

Fig3=figure('pos',[5 30 1400 750]);

subplot(1,2,1);
P3a = loglog(...
    f,V2QTN_Di,...
    f,V2QTN_Mono,...
    f,V2QTN_Disph,...
    'LineWidth',1);
ax3a = gca;
ax3a.FontSize = 15; ax3a.YLim = [1e-18 1e-12];
Le3a = legend('Dipole-wire QTN Spectrum - Maxwellian VDF',...
    'Monopole QTN Spectrum - Maxwellian VDF',...
    'Dipole-spheres QTN Spectrum - Maxwellian VDF');
Le3a.Location = 'southwest';
ax3a.XLabel.String = 'Normalized Frequency (f/f_p)'; 
ax3a.XLabel.FontSize = 16;
ax3a.XLabel.FontWeight = 'bold';
ax3a.YLabel.String = 'Voltage Power Spectral Density (V^2.Hz^{-1})';
ax3a.YLabel.FontSize = 16; ax3a.YLabel.FontWeight = 'bold';

subplot(1,2,2);
P3b = loglog(...
    f,V2QTN_Di_kap,...
    f,V2QTN_Mono_kap,...
    f,V2QTN_Disph_kap,...
    'LineWidth',1);
ax3b = gca;
ax3b.FontSize = 15; ax3b.YLim = [1e-18 1e-12];
Le3b = legend('Dipole-wire QTN Spectrum - \kappa VDF',...
    'Monopole QTN Spectrum - \kappa VDF',...
    'Dipole-spheres QTN Spectrum - \kappa VDF');
Le3b.Location = 'southwest';
ax3b.XLabel.String = 'Normalized Frequency (f/f_p)'; 
ax3b.XLabel.FontSize = 16;
ax3b.XLabel.FontWeight = 'bold';
ax3b.YLabel.String = 'Voltage Power Spectral Density (V^2.Hz^{-1})';
ax3b.YLabel.FontSize = 16; ax3b.YLabel.FontWeight = 'bold';
%saveas(Fig3,'Fig3-Di_vs_Mono_QTN_kap_norm5','epsc')

% Find the power law relationship between frequency and QTN for f>>f_p
x = f(301:end)'; % Take values from f/f_p = 2.62 to 14.57

% Assign the exponent coefficient from the power law fits to variables:
[~, diw_expo, diw_lb, diw_ub] = find_power_rel(x,V2QTN_Di(301:end)');

[~, mon_expo, mon_lb, mon_ub] = find_power_rel(x,V2QTN_Mono(301:end)');

[~, dis_expo, dis_lb, dis_ub] = find_power_rel(x,V2QTN_Disph(301:end)');

[~, diw_expok, diw_lbk, diw_ubk] = ...
    find_power_rel(x,V2QTN_Di_kap(301:end)');

[~, mon_expok, mon_lbk, mon_ubk] = ...
    find_power_rel(x,V2QTN_Mono_kap(301:end)');

[~, dis_expok, dis_lbk, dis_ubk] = ...
    find_power_rel(x,V2QTN_Disph_kap(301:end)');


%% Figure 3b - Plot QTN using F1(x), F2(x) & F3(x) - normalized units
% Take the spectra calculated previously and normalize it to the V2 level
% for f << f_p by dividing by the first value:
norm_spec = @(spec) spec./spec(1) ;

Fig3b=figure('pos',[5 30 1400 750]);

subplot(1,2,1);
P3a2 = loglog(...
    f,norm_spec(V2QTN_Di),...
    f,norm_spec(V2QTN_Mono),...
    f,norm_spec(V2QTN_Disph),...
    'LineWidth',1);
ax3a2 = gca;
ax3a2.FontSize = 15; ax3a2.YLim = [1e-4 1e4];
Le3a2 = legend('Dipole-wire QTN Spectrum - Maxwellian VDF',...
    'Monopole QTN Spectrum - Maxwellian VDF',...
    'Dipole-spheres QTN Spectrum - Maxwellian VDF');
Le3a2.Location = 'best';
ax3a2.XLabel.String = 'Normalized Frequency (f/f_p)'; 
ax3a2.XLabel.FontSize = 16 ;
ax3a2.XLabel.FontWeight = 'bold';
ax3a2.YLabel.String = 'Normalized Units';
ax3a2.YLabel.FontSize = 16; ax3a2.YLabel.FontWeight = 'bold';

subplot(1,2,2);
P3b2 = loglog(...
    f,norm_spec(V2QTN_Di_kap),...
    f,norm_spec(V2QTN_Mono_kap),...
    f,norm_spec(V2QTN_Disph_kap),...
    'LineWidth',1);
ax3b2 = gca;
ax3b2.FontSize = 15; ax3b2.YLim = [1e-4 1e4];
Le3b2 = legend('Dipole-wire QTN Spectrum - \kappa VDF',...
    'Monopole QTN Spectrum - \kappa VDF',...
    'Dipole-spheres QTN Spectrum - \kappa VDF');
Le3b2.Location = 'best';
ax3b2.XLabel.String = 'Normalized Frequency (f/f_p)'; 
ax3b2.XLabel.FontSize = 16 ;
ax3b2.XLabel.FontWeight = 'bold';
ax3b2.YLabel.String = 'Normalized Units';
ax3b2.YLabel.FontSize = 16; ax3b2.YLabel.FontWeight = 'bold';

%saveas(Fig3b,'Fig3b-Di_vs_Mono_QTN_fnorm5','epsc')



% Find the width of the peak (from frequency just before peak to
% corresponding value on the other side)
[~, Imax_Di] = max(V2QTN_Di);
Iwidth_Di = sum(V2QTN_Di>V2QTN_Di(Imax_Di-1));
f_width_Di = f(Imax_Di+Iwidth_Di) - f(Imax_Di);
[~, Imax_Mono] = max(V2QTN_Mono);
Iwidth_Mono = sum(V2QTN_Mono>V2QTN_Mono(Imax_Mono-1));
f_width_Mono = f(Imax_Mono+Iwidth_Mono) - f(Imax_Mono);

%% Figure 4 - Plot Shot Noise for dipole and monopole using F1(x) and F3(x)

% New frequency vector with values 1 to 150 of the previous vector (this is
% done because the shot noise expression is only valid for 
f2 = logspace(4,6,150);

% Load Shot noise and impedance data:
load('Monopole_paper_ShotNoise_Z_data_updated2.mat');
% If you want to run the calculations from scratch, comment out above line
% and uncomment line below:
%Run_Shot_Noise_and_Impedance_Calculations

% Normalize frequency vector for plotting:
f2 = f2./fp;

% Plot figure 4 - Shot noise comaparison using different ranges - Maxwell
Fig4=figure('pos',[5 30 1400 750]);

subplot(1,2,1);
P4a = loglog(...
f2(1:5:end),V2S_Di_mink_neg4(1:5:end),...
f2(3:5:end),V2S_Di_mink_neg5(3:5:end),...
f2(1:5:end),V2S_Di_mink_neg6(1:5:end),...
f2(1:5:end),V2S_Mono_mink_neg4(1:5:end),...
f2(1:5:end),V2S_Mono_mink_neg5(1:5:end),...
f2(1:5:end),V2S_Mono_mink_neg6(1:5:end),...
f2(1:5:end),V2S_Disph_mink_neg4(1:5:end),...
f2(3:5:end),V2S_Disph_mink_neg5(3:5:end),...
f2(1:5:end),V2S_Disph_mink_neg6(1:5:end),...
'LineWidth',1);
P4a(1).Marker = '*';
P4a(2).Marker = 'o';
P4a(4).Marker = '*';
P4a(5).Marker = 'o';
P4a(7).Marker = '*';
P4a(8).Marker = 'o';
ax4a = gca;
ax4a.FontSize = 15; ax4a.YLim = [1e-15 1e-3]; 
ax4a.XLabel.String = 'Normalized Frequency (f/f_p)'; 
ax4a.XLabel.FontSize = 16 ;
ax4a.XLabel.FontWeight = 'bold';
ax4a.YLabel.String = 'Voltage Power Spectral Density (V^2.Hz^{-1})';
ax4a.YLabel.FontSize = 16; ax4a.YLabel.FontWeight = 'bold';
Le4a = legend('Dipole-wire - k range: 10^{-4} - 10^{4}',...
'Dipole-wire - k range: 10^{-5} - 10^{5}',...
'Dipole-wire - k range: 10^{-6} - 10^{6}',...
'Monopole - k range: 10^{-4} - 10^{4}',...
'Monopole - k range: 10^{-5} - 10^{5}',...
'Monopole - k range: 10^{-6} - 10^{6}',...
'Dipole-sph - k range: 10^{-4} - 10^{4}',...
'Dipole-sph - k range: 10^{-5} - 10^{5}',...
'Dipole-sph - k range: 10^{-6} - 10^{6}');
Le4a.Location = 'best';
title('Distribution - Maxwellian');

% Plot next figure
subplot(1,2,2);
P4b = loglog(...
    f2(1:5:end),V2S_Di_mink_neg4_kap(1:5:end),...
    f2(3:5:end),V2S_Di_mink_neg5_kap(3:5:end),...
    f2(1:5:end),V2S_Di_mink_neg6_kap(1:5:end),...
    f2(1:5:end),V2S_Mono_mink_neg4_kap(1:5:end),...
    f2(1:5:end),V2S_Mono_mink_neg5_kap(1:5:end),...
    f2(1:5:end),V2S_Mono_mink_neg6_kap(1:5:end),...
    f2(1:5:end),V2S_Disph_mink_neg4_kap(1:5:end),...
    f2(3:5:end),V2S_Disph_mink_neg5_kap(3:5:end),...
    f2(1:5:end),V2S_Disph_mink_neg6_kap(1:5:end),...
    'LineWidth',1);
P4b(1).Marker = '*';
P4b(2).Marker = 'o';
P4b(4).Marker = '*';
P4b(5).Marker = 'o';
P4b(7).Marker = '*';
P4b(8).Marker = 'o';
ax4b = gca;
ax4b.FontSize = 15; ax4b.YLim = [1e-15 1e-3]; 
ax4b.XLabel.String = 'Normalized Frequency (f/f_p)'; 
ax4b.XLabel.FontSize = 16 ;
ax4b.XLabel.FontWeight = 'bold';
ax4b.YLabel.String = 'Voltage Power Spectral Density (V^2.Hz^{-1})';
ax4b.YLabel.FontSize = 16; ax4b.YLabel.FontWeight = 'bold';
Le4b = legend('Dipole-wire - k range: 10^{-4} - 10^{4}',...
'Dipole-wire - k range: 10^{-5} - 10^{5}',...
'Dipole-wire - k range: 10^{-6} - 10^{6}',...
'Monopole - k range: 10^{-4} - 10^{4}',...
'Monopole - k range: 10^{-5} - 10^{5}',...
'Monopole - k range: 10^{-6} - 10^{6}',...
'Dipole-sph - k range: 10^{-4} - 10^{4}',...
'Dipole-sph - k range: 10^{-5} - 10^{5}',...
'Dipole-sph - k range: 10^{-6} - 10^{6}');
Le4b.Location = 'best';
title('Distribution - \kappa = 4');
%saveas(Fig4,'Fig4-Di_vs_Mono_Shot_kap_fnorm5','epsc')

%% Figure 5 - Plot Shot Noise based on physical restriction of k
% Also need to plot out the physically restricted k examples. These
% include limiting k to 2pi/L and 2pi/LD

% Plot figure
Fig5=figure('pos',[5 30 1400 750]);
subplot(1,2,1)
P5a = loglog(...
    f2,V2S_Di_maxk_2pi_L,...
    f2,V2S_Di_maxk_2pi_LD,...
    f2,V2S_Mono_maxk_2pi_L,...
    f2,V2S_Mono_maxk_2pi_LD,...
    f2,V2S_Di_mink_neg6,...
    'LineWidth',1);
P5a(5).LineStyle = '--';
set(gca,'FontSize',15);
xlabel('Normalized Frequency (f/f_p)','FontSize',16,'FontWeight','bold');
ylabel('Voltage Power Spectral Density (V^2.Hz^{-1})','FontSize',16,...
    'FontWeight','bold');
title('Maxwellian VDF')

subplot(1,2,2)
P5b = loglog(...
    f2,V2S_Di_maxk_2pi_L_kap,...
    f2,V2S_Di_maxk_2pi_LD_kap,...
    f2,V2S_Mono_maxk_2pi_L_kap,...
    f2,V2S_Mono_maxk_2pi_LD_kap,...
    f2,V2S_Di_mink_neg6_kap,...
    'LineWidth',1);
P5b(5).LineStyle = '--';
set(gca,'FontSize',15);
Le5b = legend('Dipole-wire - k_{max}: 2\pi/L',...
    'Dipole-wire - k_{max}: 2\pi/\lambda_D',...
    'Monopole  - k_{max}: 2\pi/L',...
    'Monopole  - k_{max}: 2\pi/\lambda_D',...
    'Dipole-wire  - k_{max}: 10^{6}');
Le5b.Position = [0.40    0.81    0.1    0.1];
xlabel('Normalized Frequency (f/f_p)','FontSize',16,'FontWeight','bold');
title('Kappa VDF \kappa = 4')
%saveas(Fig5,'Fig5-Di_vs_Mono_Shot_kmax_kap_norm5','epsc')

%% Figure 6 - Plot multiples of 2pi/LD to find max k variable

Fig6=figure('pos',[5 30 1400 750]);
subplot(1,2,1)
P6a = loglog(...
    f2,V2S_Di_maxk_4pi_LD,...
    f2,V2S_Di_maxk_8pi_LD,...
    f2,V2S_Di_maxk_16pi_LD,...
    f2,V2S_Mono_maxk_4pi_LD,...
    f2,V2S_Mono_maxk_8pi_LD,...
    f2,V2S_Mono_maxk_16pi_LD,...
    f2,V2S_Mono_maxk_32pi_LD,...
    f2,V2S_Di_mink_neg6,...
    'LineWidth',1);
P6a(8).LineStyle = '--';
set(gca,'FontSize',15);
xlabel('Normalized Frequency (f/f_p)','FontSize',16,'FontWeight','bold');
ylabel('Voltage Power Spectral Density (V^2.Hz^{-1})','FontSize',16,...
    'FontWeight','bold');
title('Maxwellian VDF')

subplot(1,2,2)
P6b = loglog(...
    f2,V2S_Di_maxk_4pi_LD_kap,...
    f2,V2S_Di_maxk_8pi_LD_kap,...
    f2,V2S_Di_maxk_16pi_LD_kap,...
    f2,V2S_Mono_maxk_4pi_LD_kap,...
    f2,V2S_Mono_maxk_8pi_LD_kap,...
    f2,V2S_Mono_maxk_16pi_LD_kap,...
    f2,V2S_Mono_maxk_32pi_LD,...
    f2,V2S_Di_mink_neg6_kap,...
    'LineWidth',1);
P6b(8).LineStyle = '--';
set(gca,'FontSize',15);
Le6b = legend('Dipole-wire - k_{max}: 4\pi/\lambda_D',...
    'Dipole-wire - k_{max}: 8\pi/\lambda_D',...
    'Dipole-wire - k_{max}: 16\pi/\lambda_D',...
    'Monopole  - k_{max}: 4\pi/\lambda_D',...
    'Monopole  - k_{max}: 8\pi/\lambda_D',...
    'Monopole  - k_{max}: 16\pi/\lambda_D',...
    'Monopole  - k_{max}: 32\pi/\lambda_D',...
    'Dipole-wire  - k_{max}: 10^{6}');
Le6b.Position = [0.40    0.775    0.1    0.1];
xlabel('Normalized Frequency (f/f_p)','FontSize',16,'FontWeight','bold');
title('Kappa VDF \kappa = 4')
%saveas(Fig6,'Fig6-Di_vs_Mono_Shot_kmax_kap_npi-LD5','epsc')

%% Figure 7 - zoomed in version of Figure 6

Fig7=figure('pos',[5 30 1400 750]);
subplot(1,2,1)
P7a = loglog(...
    f2,V2S_Di_maxk_4pi_LD,...
    f2,V2S_Di_maxk_8pi_LD,...
    f2,V2S_Di_maxk_16pi_LD,...
    f2,V2S_Mono_maxk_4pi_LD,...
    f2,V2S_Mono_maxk_8pi_LD,...
    f2,V2S_Mono_maxk_16pi_LD,...
    f2,V2S_Di_mink_neg6,...
    'LineWidth',1);
P7a(7).LineStyle = '--';
set(gca,'FontSize',15,'YLim',[0.8e-11 2.5e-11],'XLim',[1.5e-3 2e-3]);
xlabel('Normalized Frequency (f/f_p)','FontSize',16,'FontWeight','bold');
ylabel('Voltage Power Spectral Density (V^2.Hz^{-1})','FontSize',16,...
    'FontWeight','bold');
title('Dipole')

subplot(1,2,2)
P7b = loglog(...
    f2,V2S_Di_maxk_4pi_LD,...
    f2,V2S_Di_maxk_8pi_LD,...
    f2,V2S_Di_maxk_16pi_LD,...
    f2,V2S_Mono_maxk_4pi_LD,...
    f2,V2S_Mono_maxk_8pi_LD,...
    f2,V2S_Mono_maxk_16pi_LD,...
    f2,V2S_Di_mink_neg6,...
    'LineWidth',1);
P7b(7).LineStyle = '--';
set(gca,'FontSize',15,'YLim',[0.45e-5 2.6e-5],'XLim',[1.5e-3 2e-3]);
Le7b = legend('Dipole-wire - k_{max}: 4\pi/\lambda_D',...
    'Dipole-wire - k_{max}: 8\pi/\lambda_D',...
    'Dipole-wire - k_{max}: 16\pi/\lambda_D',...
    'Monopole  - k_{max}: 4\pi/\lambda_D',...
    'Monopole  - k_{max}: 8\pi/\lambda_D',...
    'Monopole  - k_{max}: 16\pi/\lambda_D',...
    'Dipole-wire  - k_{max}: 10^{6}');
Le7b.Position = [0.39    0.75    0.1    0.1];
xlabel('Normalized Frequency (f/f_p)','FontSize',16,'FontWeight','bold');
title('Monopole')
%saveas(Fig7,'Fig7-Di_vs_Mono_Shot_kmax_npi-LD_zoom5','epsc')


%% Calculation of ratios from increasing kmax by multiples of 2pi/LD

% Ratio of increase from 4 to 8 to 16 pi/LD for Dipole case:
r_Di_4pi_ref = V2S_Di_maxk_4pi_LD(1)/V2S_Di_mink_neg6(1);
r_Di_8pi_ref = V2S_Di_maxk_8pi_LD(1)/V2S_Di_mink_neg6(1);
r_Di_16pi_ref = V2S_Di_maxk_16pi_LD(1)/V2S_Di_mink_neg6(1);

r_Di_16pi_ref_kap = V2S_Di_maxk_16pi_LD_kap(1)/V2S_Di_mink_neg6_kap(1);

r_Di_4pi_8pi = V2S_Di_maxk_4pi_LD(1)/V2S_Di_maxk_8pi_LD(1);
r_Di_8pi_16pi = V2S_Di_maxk_8pi_LD(1)/V2S_Di_maxk_16pi_LD(1);

r_diff4pi8pi_4pi = (V2S_Di_maxk_8pi_LD(1)-V2S_Di_maxk_4pi_LD(1))...
    /V2S_Di_maxk_4pi_LD(1);
r_diff8pi16pi_8pi = (V2S_Di_maxk_16pi_LD(1)-V2S_Di_maxk_8pi_LD(1))...
    /V2S_Di_maxk_8pi_LD(1);
r_diff16piref_16pi = (V2S_Di_mink_neg6(1)-V2S_Di_maxk_16pi_LD(1))...
    /V2S_Di_maxk_16pi_LD(1);

r_diff16piref_16pi_kap = ...
    (V2S_Di_mink_neg6_kap(1)-V2S_Di_maxk_16pi_LD_kap(1))...
    /V2S_Di_maxk_16pi_LD_kap(1);

r_diff4pi8pi_4pi_kap = ...
    (V2S_Di_maxk_8pi_LD_kap(1)-V2S_Di_maxk_4pi_LD_kap(1))...
    /V2S_Di_maxk_4pi_LD_kap(1);

r_diff8pi16pi_8pi_kap = ...
    (V2S_Di_maxk_16pi_LD_kap(1)-V2S_Di_maxk_8pi_LD_kap(1))...
    /V2S_Di_maxk_8pi_LD_kap(1);

% Ratio of increase from 4 to 8 to 16 pi/LD for Monopole case:
r_Mono_4pi_ref = V2S_Mono_maxk_4pi_LD(1)/V2S_Mono_mink_neg6_kap(1);
r_Mono_8pi_ref = V2S_Mono_maxk_8pi_LD(1)/V2S_Mono_mink_neg6_kap(1);
r_Mono_16pi_ref = V2S_Mono_maxk_16pi_LD(1)/V2S_Mono_mink_neg6_kap(1);

r_Mono_4pi_8pi = V2S_Mono_maxk_4pi_LD(1)/V2S_Mono_maxk_8pi_LD(1);
r_Mono_8pi_16pi = V2S_Mono_maxk_8pi_LD(1)/V2S_Mono_maxk_16pi_LD(1);

r_diff4pi8pi_4pi_Mono = ...
    (V2S_Mono_maxk_8pi_LD(1)-V2S_Mono_maxk_4pi_LD(1))...
    /V2S_Mono_maxk_4pi_LD(1);

r_diff8pi16pi_8pi_Mono = ...
    (V2S_Mono_maxk_16pi_LD(1)-V2S_Mono_maxk_8pi_LD(1))...
    /V2S_Mono_maxk_8pi_LD(1);

r_diff16pi32pi_16pi_Mono = ...
    (V2S_Mono_maxk_32pi_LD(1)-V2S_Mono_maxk_16pi_LD(1))...
    /V2S_Mono_maxk_16pi_LD(1);

r_diff4pi8pi_4pi_kap_Mono = ...
    (V2S_Mono_maxk_8pi_LD_kap(1)-V2S_Mono_maxk_4pi_LD_kap(1))...
    /V2S_Mono_maxk_4pi_LD_kap(1);
r_diff8pi16pi_8pi_kap_Mono = ...
    (V2S_Mono_maxk_16pi_LD_kap(1)-V2S_Mono_maxk_8pi_LD_kap(1))...
    /V2S_Mono_maxk_8pi_LD_kap(1);


%% Calculate Capacitance Comparison
% Values for impedance have been calculated in previous set of results but
% we do need to calculate the approximate dipole capacitance from the
% antenna parameters:
C_di_apprx = (pi*eps0*l)/(log(LD/a));
C_mono_apprx = 2*C_di_apprx;

% Frequency vector
f3 = logspace(4,6,150);

% Capacitance based on the impedance - Mawellian
C_di_mink_neg4 = 1./((2*pi*f3).*imag(Z_di_mink_neg4));
C_di_mink_neg5 = 1./((2*pi*f3).*imag(Z_di_mink_neg5));
C_di_mink_neg6 = 1./((2*pi*f3).*imag(Z_di_mink_neg6));
C_di_mink_neg7 = 1./((2*pi*f3).*imag(Z_di_mink_neg7));

C_mono_mink_neg4 = 1./((2*pi*f3).*imag(Z_mono_mink_neg4));
C_mono_mink_neg5 = 1./((2*pi*f3).*imag(Z_mono_mink_neg5));
C_mono_mink_neg6 = 1./((2*pi*f3).*imag(Z_mono_mink_neg6));
C_mono_mink_neg7 = 1./((2*pi*f3).*imag(Z_mono_mink_neg7));

C_disph_mink_neg4 = 1./((2*pi*f3).*imag(Z_disph_mink_neg4));
C_disph_mink_neg5 = 1./((2*pi*f3).*imag(Z_disph_mink_neg5));
C_disph_mink_neg6 = 1./((2*pi*f3).*imag(Z_disph_mink_neg6));
C_disph_mink_neg7 = 1./((2*pi*f3).*imag(Z_disph_mink_neg7));

% Kappa
C_di_mink_neg6_kap = 1./((2*pi*f3).*imag(Z_di_mink_neg6_kap));
C_mono_mink_neg4_kap = 1./((2*pi*f3).*imag(Z_mono_mink_neg4_kap));
C_mono_mink_neg5_kap = 1./((2*pi*f3).*imag(Z_mono_mink_neg5_kap));
C_mono_mink_neg6_kap = 1./((2*pi*f3).*imag(Z_mono_mink_neg6_kap));
C_mono_mink_neg7_kap = 1./((2*pi*f3).*imag(Z_mono_mink_neg7_kap));

C_disph_mink_neg4_kap = 1./((2*pi*f3).*imag(Z_disph_mink_neg4_kap));
C_disph_mink_neg5_kap = 1./((2*pi*f3).*imag(Z_disph_mink_neg5_kap));
C_disph_mink_neg6_kap = 1./((2*pi*f3).*imag(Z_disph_mink_neg6_kap));
C_disph_mink_neg7_kap = 1./((2*pi*f3).*imag(Z_disph_mink_neg7_kap));

% Capacitance values for k_max 
C_di_2piL = 1./((2*pi*f3).*imag(Z_Di_2pi_L));
C_di_2piLD = 1./((2*pi*f3).*imag(Z_Di_2pi_LD));
C_di_16piLD = 1./((2*pi*f3).*imag(Z_Di_16pi_LD));
C_di_32piLD = 1./((2*pi*f3).*imag(Z_Di_32pi_LD));

C_mono_2piL = 1./((2*pi*f3).*imag(Z_Mono_2pi_L));
C_mono_2piLD = 1./((2*pi*f3).*imag(Z_Mono_2pi_LD));
C_mono_16piLD = 1./((2*pi*f3).*imag(Z_Mono_16pi_LD));
C_mono_32piLD = 1./((2*pi*f3).*imag(Z_Mono_32pi_LD));

% Approx. Capacitance at 800 and 1500 km (approx. for f << f_p)
% 800 km
Ne800 = xlsread('Ne_800km.xlsx');
Ne800 = Ne800(2:end,2:end);
Te800 = xlsread('Te_800km.xlsx');
Te800 = Te800(2:end,2:end);
n_e800 = mean(mean(Ne800));
T_e800 = mean(mean(Te800));

% 1500 km
Ne1500 = xlsread('Ne_1500km.xlsx');
Ne1500 = Ne1500(2:end,2:end);
Te1500 = xlsread('Te_1500km.xlsx');
Te1500 = Te1500(2:end,2:end);
n_e1500 = mean(mean(Ne1500));
T_e1500 = mean(mean(Te1500));
LD800 = sqrt((kB*T_e800*eps0)/(n_e800*e^2));
vT800 = sqrt(2*kB*T_e800/me);
wp800 = vT800/(sqrt(2)*LD800);

LD1500 = sqrt((kB*T_e1500*eps0)/(n_e1500*e^2));
vT1500 = sqrt(2*kB*T_e1500/me);
wp1500 = vT1500/(sqrt(2)*LD1500);

C_di_800km_apprx = (pi*eps0*l)/(log(LD800/a));
C_di_1500km_apprx = (pi*eps0*l)/(log(LD1500/a));
C_mono_800km_apprx = 2*(pi*eps0*l)/(log(LD800/a));
C_mono_1500km_apprx = 2*(pi*eps0*l)/(log(LD1500/a));

% Use the monopole approximation as the criterion for choosing k_max, based
% on the adjustment of k to reach this aim, we have found k = 207.4 brings
% the capacitance 

% Capacitance
C_di_207_4 = 1./((2*pi*f3).*imag(Z_Di_207_4));
C_mono_207_4 = 1./((2*pi*f3).*imag(Z_Mono_207_4));

% Also at 800 and 1500 km:
C_di_800km = 1./((2*pi*f3).*imag(Z_Di_800km));
C_mono_800km = 1./((2*pi*f3).*imag(Z_Mono_800km));
C_di_1500km = 1./((2*pi*f3).*imag(Z_Di_1500km));
C_mono_1500km = 1./((2*pi*f3).*imag(Z_Mono_1500km));
