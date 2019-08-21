% Script to calculate the QTN spectra for monopole and dipoles for both
% Maxwellian and Kappa distributions

% Maxwellian
% Dipole - wire QTN
[V2QTN_Di,~] = plasma_noiseTEIono58_2(l,n_e,T_e,a,5,8,400,1,1e6,1e-6);

% Dipole - spheres QTN
[V2QTN_Disph,~] = plasma_noiseTEIono58_2(l,n_e,T_e,a,5,8,400,2,1e6,1e-6);

% Monopole QTN
[V2QTN_Mono,~] = plasma_noiseTEIono58_2(l,n_e,T_e,a,5,8,400,0,1e6,1e-6);

% Kappa versions of each of the above:
[V2QTN_Di_kap,~] = plasma_noiseTEIono58_4(l,n_e,T_e,a,5,8,400,...
    1,1e6,1e-6,kap);
[V2QTN_Disph_kap,~] = plasma_noiseTEIono58_4(l,n_e,T_e,a,5,8,400,...
    2,1e6,1e-6,kap);
[V2QTN_Mono_kap,~] = plasma_noiseTEIono58_4(l,n_e,T_e,a,5,8,400,...
    0,1e6,1e-6,kap);