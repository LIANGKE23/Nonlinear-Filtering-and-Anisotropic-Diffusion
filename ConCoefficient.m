function [ConC] = ConCoefficient(mode,gradient, K)

if mode == 1
    ConC = exp(-abs(gradient)/K/K);
elseif mode == 2
    ConC = 1/(1+abs(gradient)/K/K);
end

