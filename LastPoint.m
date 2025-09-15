function [y_last] = LastPoint(u_last, M_axis, gamma, R, p_st, T_st, MassFlow_IVL)

t = T_st/(1 + (gamma - 1)/2*M_axis^2);                       % Static temperature
p = p_st/(1 + (gamma - 1)/2*M_axis^2)^(gamma/(gamma - 1));   % Static pressure
rho = p/(R*t);                                               % Density
y_last = round((MassFlow_IVL/(pi*rho*u_last))^(1/2),4);

end


