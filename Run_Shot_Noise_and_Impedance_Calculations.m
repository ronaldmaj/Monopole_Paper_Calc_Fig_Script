%% Calculations of Shot Noise/Impedance over certain frequency ranges
% Used for Figure 4

% Calculate the shot noise spectra over the new frequency range - Maxwell
[~, V2S_Di_mink_neg4, Z_di_mink_neg4, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,1e4,1e-4,0,0);

[~, V2S_Di_mink_neg5, Z_di_mink_neg5, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,1e5,1e-5,0,0);

[~, V2S_Di_mink_neg6, Z_di_mink_neg6, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,1e6,1e-6,0,0);

[~, V2S_Di_mink_neg7, Z_di_mink_neg7, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,1e7,1e-7,0,0);

[~, V2S_Mono_mink_neg4, Z_mono_mink_neg4, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,1e4,1e-4,0,0);

[~, V2S_Mono_mink_neg5, Z_mono_mink_neg5, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,1e5,1e-5,0,0);

[~, V2S_Mono_mink_neg6, Z_mono_mink_neg6, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,1e6,1e-6,0,0);

[~, V2S_Mono_mink_neg7, Z_mono_mink_neg7, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,1e7,1e-7,0,0);

[~, V2S_Disph_mink_neg4, Z_disph_mink_neg4, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,1,4,6,150,1,1e4,1e-4,0,0);

[~, V2S_Disph_mink_neg5, Z_disph_mink_neg5, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,1,4,6,150,1,1e5,1e-5,0,0);

[~, V2S_Disph_mink_neg6, Z_disph_mink_neg6, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,1,4,6,150,1,1e6,1e-6,0,0);

[~, V2S_Disph_mink_neg7, Z_disph_mink_neg7, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,1,4,6,150,1,1e6,1e-6,0,0);

% Calculate the shot noise spectra over the new frequency range - Kappa
[~, V2S_Di_mink_neg4_kap, Z_di_mink_neg4_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,1e4,1e-4,1,4);

[~, V2S_Di_mink_neg5_kap, Z_di_mink_neg5_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,1e5,1e-5,1,4);

[~, V2S_Di_mink_neg6_kap, Z_di_mink_neg6_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,1e6,1e-6,1,4);

[~, V2S_Di_mink_neg7_kap, Z_di_mink_neg7_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,1e7,1e-7,1,4);

[~, V2S_Mono_mink_neg4_kap, Z_mono_mink_neg4_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,1e4,1e-4,1,4);

[~, V2S_Mono_mink_neg5_kap, Z_mono_mink_neg5_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,1e5,1e-5,1,4);

[~, V2S_Mono_mink_neg6_kap, Z_mono_mink_neg6_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,1e6,1e-6,1,4);

[~, V2S_Mono_mink_neg7_kap, Z_mono_mink_neg7_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,1e7,1e-7,1,4);

[~, V2S_Disph_mink_neg4_kap, Z_disph_mink_neg4_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,1,4,6,150,1,1e4,1e-4,1,4);

[~, V2S_Disph_mink_neg5_kap, Z_disph_mink_neg5_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,1,4,6,150,1,1e5,1e-5,1,4);

[~, V2S_Disph_mink_neg6_kap, Z_disph_mink_neg6_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,1,4,6,150,1,1e6,1e-6,1,4);

[~, V2S_Disph_mink_neg7_kap, Z_disph_mink_neg7_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,1,4,6,150,1,1e7,1e-7,1,4);

%% Calculation based on restricted k to a physically specificed criterion
% Used for Figure 5

% Physically restricted k examples. These include limiting k to 2pi/L and 
% 2pi/LD
kmax1 = (2*pi)/l;
kmax2 = (2*pi)/LD;
kmax2k = (2*pi)/LDk(kap);

% Maxwellian
[~, V2S_Di_maxk_2pi_L, Z_Di_2pi_L, ~] = ...
    plasma_noise_mono_vs_di(l,n_e,T_e,a,0,4,6,150,1,kmax1,1e-6);
    
[~, V2S_Di_maxk_2pi_LD, Z_Di_2pi_LD, ~] = ...
    plasma_noise_mono_vs_di(l,n_e,T_e,a,0,4,6,150,1,kmax2,1e-6);
    
[~, V2S_Mono_maxk_2pi_L, Z_Mono_2pi_L, ~] = ...
    plasma_noise_mono_vs_di(l,n_e,T_e,a,2,4,6,150,0,kmax1,1e-6);
    
[~, V2S_Mono_maxk_2pi_LD, Z_Mono_2pi_LD, ~] = ...
    plasma_noise_mono_vs_di(l,n_e,T_e,a,2,4,6,150,0,kmax2,1e-6);
    

% Kappa
[~, V2S_Di_maxk_2pi_L_kap, Z_Di_2pi_L_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,kmax1,1e-6,1,4);

[~, V2S_Di_maxk_2pi_LD_kap, Z_Di_2pi_LD_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,(kmax2k),1e-6,1,4);

[~, V2S_Mono_maxk_2pi_L_kap, Z_Mono_2pi_L_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,kmax1,1e-6,1,4);

[~, V2S_Mono_maxk_2pi_LD_kap, Z_Mono_2pi_LD_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,kmax2k,1e-6,1,4);

%% Calculation based on restricted k to multiples of 2pi/LD:
% Used for Figure 6/7

% Multiples of 2pi/LD:
kmax3 = (4*pi)/LD;
kmax4 = (8*pi)/LD;
kmax5 = (16*pi)/LD;
kmax6 = (32*pi)/LD;

kmax3k = (4*pi)/LDk(kap);
kmax4k = (8*pi)/LDk(kap);
kmax5k = (16*pi)/LDk(kap);
kmax6k = (32*pi)/LDk(kap);

% Maxwellian calculations:
[~, V2S_Di_maxk_4pi_LD, Z_Di_4pi_LD, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,kmax3,1e-6,0,0);
    
[~, V2S_Mono_maxk_4pi_LD, Z_Mono_4pi_LD, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,kmax3,1e-6,0,0);
    
[~, V2S_Di_maxk_8pi_LD, Z_Di_8pi_LD, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,kmax4,1e-6,0,0);
    
[~, V2S_Mono_maxk_8pi_LD, Z_Mono_8pi_LD, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,kmax4,1e-6,0,0);
    
[~, V2S_Di_maxk_16pi_LD, Z_Di_16pi_LD, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,kmax5,1e-6,0,0);
    
[~, V2S_Mono_maxk_16pi_LD, Z_Mono_16pi_LD, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,kmax5,1e-6,0,0);
    
[~, V2S_Di_maxk_32pi_LD, Z_Di_32pi_LD, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,kmax6,1e-6,0,0);
    
[~, V2S_Mono_maxk_32pi_LD, Z_Mono_32pi_LD, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,kmax6,1e-6,0,0);
    
% Kappa calculations:
[~, V2S_Di_maxk_4pi_LD_kap, Z_Di_4pi_LD_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,kmax3k,1e-6,1,4);

[~, V2S_Mono_maxk_4pi_LD_kap, Z_Mono_4pi_LD_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,kmax3k,1e-6,1,4);

[~, V2S_Di_maxk_8pi_LD_kap, Z_Di_8pi_LD_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,kmax4k,1e-6,1,4);

[~, V2S_Mono_maxk_8pi_LD_kap, Z_Mono_8pi_LD_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,kmax4k,1e-6,1,4);

[~, V2S_Di_maxk_16pi_LD_kap, Z_Di_16pi_LD_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,kmax5k,1e-6,1,4);

[~, V2S_Mono_maxk_16pi_LD_kap, Z_Mono_16pi_LD_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,kmax5k,1e-6,1,4);

[~, V2S_Di_maxk_32pi_LD_kap, Z_Di_32pi_LD_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,0,4,6,150,1,kmax6k,1e-6,0,0);
    
[~, V2S_Mono_maxk_32pi_LD_kap, Z_Mono_32pi_LD_kap, ~] = ...
    plasma_noise_mono_vs_di_kap(l,n_e,T_e,a,2,4,6,150,0,kmax6k,1e-6,0,0);
    
%% Calculation based on restricted k based on capacitance approximation:
% Used in Table 2

% Calculate the average ionospheric conditions for 800 and 1500km altitude:
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

% Use the monopole approximation as the criterion for choosing k_max, based
% on the adjustment of k to reach this aim, we have found k = 207.4 brings
% the capacitance 
[~, ~, Z_Di_207_4, ~] = plasma_noise_mono_vs_di(l,n_e,T_e,a,...
    0,4,6,150,1,207.4,1e-6);

[~, ~, Z_Mono_207_4, ~] = plasma_noise_mono_vs_di(l,n_e,T_e,a,...
    2,4,6,150,0,207.4,1e-6);


% Also at 800 and 1500 km:
[~, ~, Z_Di_800km, ~] = plasma_noise_mono_vs_di(l,n_e800,T_e800,a,...
    0,4,6,150,1,207.4,1e-6);

[~, ~, Z_Mono_800km, ~] = plasma_noise_mono_vs_di(l,n_e800,T_e800,a,...
    2,4,6,150,0,207.4,1e-6);

[~, ~, Z_Di_1500km, ~] = plasma_noise_mono_vs_di(l,n_e1500,T_e1500,a,...
    0,4,6,150,1,207.4,1e-6);

[~, ~, Z_Mono_1500km, ~] = plasma_noise_mono_vs_di(l,n_e1500,T_e1500,a,...
    2,4,6,150,0,207.4,1e-6);