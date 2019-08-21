function [V2QTN, V2S, Z_a, w_vec] = plasma_noise_mono_vs_di_kap(l,N,T,a,...
    sph,fmin,fmax,n,dipole,k_max,k_min,k_dis,kap)
    % Code for calculating the voltage power spectrum similar to that done in
    % the Meyer-Vernet/Perche 1989 Paper "Tool Kit for Antennae and Thermal
    % Noise Near the Plasma Frequency". Produces a plot of V^2 at the
    % end, and outputs this voltage spectrum and the frequency range.
    %
    % Inputs: 
    % l     - length of antenna                 (m)
    % N     - density of electrons              (m^-3)
    % T     - temperature of electrons          (K)
    % a     - radius of antenna                 (m)
    % sph   - spherical antenna?                (1 = yes, 0 = cylindrical)
    % fmin  - lower bound for frequency vector  (power of 10)
    % fmax  - upper bound for frequency vector  (power of 10)
    % n     - no. of elements in freq. vector
    % dipole- considering dipole antenna?       (1 = yes, 0 = mono)
    % k_max - maximum k value in integral       (m^-1)
    % k_min - minimum k value in integral       (m^-1)
    % k_dis - kappa VDF for electrons?         (1 = yes, 0 = no)
    % kap   - kappa value
    % 
    % Physical constants
    kB = 1.3806488e-23;     % Boltzmann constant
    me = 9.109383e-31;      % Mass of electron
    e = 1.60217657e-19;     % Electron charge
    eps0 = 8.85418782e-12;  % Vacuum Permittivity 

    % Plasma parameters
    vT = sqrt(2*kB*T/me);   % Thermal velocity of electrons
    LD = sqrt((kB*T*eps0)/(N*e^2));
    wp = vT/(sqrt(2)*LD); % Plasma frequency

    % Antenna Response Function in integrand
    if dipole == 1
        F1 = @(x) ((sinint(x))-(sinint(2*x)./2)-((2*((sin(x/2)).^4))./x))./x;
    else
        F1 = @(x) 0.125*(x.^2 + 2*cos(x) - 2)./(x.^2);
    end

    % Faddeeva function:
    W = @(x) faddeeva(x);

    % Frequency vector
    f = logspace(fmin, fmax, n); 
    w_vec = f*(2*pi); % Determine omega from f

    % Surface area based on whether spherical or cylindrical antenna
    if sph == 0
        S=2*pi*a*l;
    elseif sph == 1 
        S=4*pi*a^2;
        if dipole == 1
            F1 = @(x) 0.25*(1 - sin(x)./x);
        end
    else
        S = 0.01*6;
    end

    % Impact rate of electrons onto antenna
    N_e = ((4*pi)^(-0.5))*N*vT*S;

    % Vector for recording resistance
    V2QTN = zeros(size(w_vec));
    V2S = zeros(size(w_vec));
    Z_a = zeros(size(w_vec));

    % Progress bar initialization:
    textprogressbar('Progress: ');
    for i = 1:length(w_vec)

        % Display which frequency data point 
        %disp([int2str(i) ' out of ' int2str(length(w_vec)) ' complete.'])
        textprogressbar((i/length(w_vec))*100);

        % Permittivity
        w = w_vec(i);           % Set specific omega value
        z = @(x) w./(x.*vT); 

        if k_dis == 1
            % For a kappa distribution, we need definitions for r:
            r = @(wi) wi/wp; 
            % Make z and p same size matrices
            P = @(z,p) p'*ones(1,length(z));
            Z = @(z,p) ones(length(p),1)*z;
            % Implement summation term
            S = @(k,p) factorial(kap + P(z(k),p))./(factorial(P(z(k),p)).*...
                (2*1i).^(kap+1+P(z(k),p)).*(Z(z(k),p)+1i).^(kap+1-P(z(k),p)));
            % Dielectric permittivity
            % Epsilon_L
            epsL = @(k, kap, w) 1 + ((z(k).^2)/(r(w).^2)).*...
            (2*kap - 1 + ((-2).^(kap+1))./prod(1:2:(2*kap-3)).*...
            (1i*z(k)).*sum(S(k,0:kap)));
            % Integral function:
            Int = @(k) (F1(k*l).*(besselj(0,k*a)).^2).*((1./(epsL(k, kap, w))));
        else
            % Dielectric permittivity
            epsL = @(k) 1 + ((1+(1i)*sqrt(pi)*z(k).*W(z(k)))./((k.^2)*(LD^2)));
            % Integral calculation
            Int = @(k) (F1(k*l).*(besselj(0,k*a)).^2).*((1./(epsL(k))));
        end
        % Integral calculation
        kmin = k_min;
        kmax = k_max;
        Integrand = integral(Int,kmin,kmax); 

        % Impedance
        Z = (4i/(pi^2*eps0*w))*Integrand;
        Z_a(i) = Z;

        % Voltage power for QTN and shot noise at frequency w
        V2QTN(i) = 4*kB*T*real(Z);
        V2S(i) = (2*(e^2)*N_e)*((abs(Z))^2);

        clear Z
    end
    textprogressbar('Complete');
end