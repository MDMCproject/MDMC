%% Exact integration over constant g(r)=g_value for all r to r_max
% see page 27 in my handwritten MDMC notes
function retval = integrate_over_r(Q, g_value)
  r_max = 21.2;
  n_atom = 1372;
  volume = 42.718286^3;
  retval = (g_value-1)*((n_atom-1)/volume)*4*pi*(sin(Q*r_max)/Q-r_max*cos(Q*r_max))/Q^2;