function matlab_cal_s_q_omega  %(filename)
%matlab_cal_s_q_omega
%
% Written for debugging MDMC code. This file aim to cal 
% s_q_omega from a s_q_time file to directly compare with that
% returned by MDMC.
%
% Input: s_q_time file 

s_q_time_info = read_s_q_time('s_q_time1.xml');
s_q_omega_info = read_s_q_omega('s_q_omega1.xml');

omega = s_q_omega_info.omega;
n_omega = length(omega);

time_length = s_q_time_info.t(2) - s_q_time_info.t(1);
n_t = length(s_q_time_info.t);

R = zeros([n_t n_omega]);

which_method = 1;  % just trying different ways of integrating over time
                   % see page 31 of my notes

if which_method == 1
  
  prefac = time_length / pi;
  
  for i_omega = 1 : n_omega
       
    R(1, i_omega) = 0.5;
    
    for i_t = 2 : n_t
      t = (i_t-1)* time_length;

      R(i_t, i_omega) = cos(omega(i_omega)*t);
    end
  end
  
  R = R * prefac;
  
else

  for i_omega = 1 : n_omega
    prefac = 1 / ( pi * omega(i_omega) );

    lower_lim = 0;  % since sin(omega*t) = 0 for t=0

    for i_t = 0 : n_t-1
      t = (i_t+0.5) * time_length;

      upper_lim = sin(omega(i_omega)*t);

      R(i_t+1, i_omega) = (upper_lim - lower_lim) * prefac;

      lower_lim = upper_lim;
    end
  end

end

s_q_omega_diff = s_q_time_info.S_d * R;
s_q_omega_self = s_q_time_info.S_s * R;

q = s_q_omega_info.q;

subplot(1,2,1)
surf(q, omega, s_q_omega_diff')
xlabel('q')
ylabel('\omega')
zlabel('S^d (q,\omega)')

subplot(1,2,2)
surf(q, omega, s_q_omega_self')
xlabel('q')
ylabel('\omega')
zlabel('S^s (q,\omega)')


end 