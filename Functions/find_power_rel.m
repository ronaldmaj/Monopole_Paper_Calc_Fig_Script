function [const, expo, lb, ub] = find_power_rel(X,Y)

    % Take the x and y values and find the power law fit:
    pow_fit = fit(X,Y,'power1');
    coef = coeffvalues(pow_fit);
    const = coef(1);
    expo = coef(2);
    
    % Confidence intervals:
    ci = confint(pow_fit,0.95);
    lb = ci(1,2);
    ub = ci(2,2);