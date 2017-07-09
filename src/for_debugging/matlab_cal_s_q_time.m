% this function requires you to first set volume and number of atoms

function matlab_cal_s_q_time(filename_s_d, filename_s_s)
%matlab_cal_s_q_time
%
% Written for debugging MDMC code. This file aim to plot 
% s_q_time from a h_d_hist file to directly compare with that
% returned by MDMC.
%
% Input:
%  filename : filename of a normalised h_d file generated 
%             from MDMC using print_h_d_hist()

n_atom = 1372; %864; %1372;
volume = 42.718286^3; %36.615673^3; %42.718286^3;



h_d_info = read_h_d(filename_s_d);
h_s_info = read_h_s(filename_s_s);

q = [0.1:0.1:4.1];

n_r = h_d_info.n_bin;
n_q = length(q);

for i_q = 1 : n_q
      prefac = 4 * pi / q(i_q)^2;
      
      lower_lim = 0;  % since sin(Q*r) = 0 for r=0
      
      for i_r = 1 : n_r
        r = i_r * h_d_info.bin_length;
        
        upper_lim = sin(q(i_q)*r) / q(i_q) - r * cos(q(i_q)*r);
        
        R(i_r, i_q) = (upper_lim - lower_lim) * prefac;
        
        lower_lim = upper_lim ;
      end
end

% calculate histogram volume plus going toward 1 normalisation constant

r_center = h_d_info.r;
i_num = [1:n_r];
delta_r = h_d_info.bin_length;
prefac_s_d = (volume/(n_atom*(n_atom-1))) ./ ...
  ( (4*pi/3)*delta_r^3*(i_num.^3-(i_num-1).^3) );
prefac_s_s = (volume/n_atom) ./ ...
  ( (4*pi/3)*delta_r^3*(i_num.^3-(i_num-1).^3) );

% include the histogram volumes with the histogram

prefac_s_combined = (volume/n_atom^2) ./ ( (4*pi/3)*delta_r^3*(i_num.^3-(i_num-1).^3) );
h_combined = zeros(size(h_s_info.val));
for i_t = 1 : length(h_d_info.t)
  h_combined(:,i_t) =   (h_d_info.val(:,i_t)+h_s_info.val(:,i_t)) .* prefac_s_combined';
  h_d_info.val(:,i_t) = h_d_info.val(:,i_t) .* prefac_s_d';
  h_s_info.val(:,i_t) = h_s_info.val(:,i_t) .* prefac_s_s'; 
end

% R (n_r x n_q), h_d_info.val (n_r x n_q) 

s_q_time_diff = ((n_atom-1)/volume)*R' * (h_d_info.val-1);  % should have dim n_q x n_t
s_q_time_self = (1/volume)*R' * (h_s_info.val-1);  % should have dim n_q x n_t

% try also to perform the integration directly

s_q_time_combined = (n_atom/volume)*R' * (h_combined-1);

% We know S_self(Q,t) = 1 exactly when t = 0

s_q_time_self(:,1) = 1;

% Robert corrected S^s and S^d

s_q_time_robert_diff = zeros(size(s_q_time_diff));
s_q_time_robert_self = zeros(size(s_q_time_self));
for i = 1 : length(h_d_info.t)
   s_q_time_robert_diff(:,i) =  s_q_time_diff(:,i) + s_q_time_self(:,end);
   s_q_time_robert_self(:,i) =  s_q_time_self(:,i) - s_q_time_self(:,end); 
end

subplot(1,2,1)
surf(q, h_d_info.t, s_q_time_diff'+s_q_time_self')
xlabel('q')
ylabel('t')
%zlabel('S\^d (q,t)')

subplot(1,2,2)
surf(q, h_d_info.t, s_q_time_combined')
xlabel('q')
ylabel('t')
%zlabel('S\^s (q,t)')


end 