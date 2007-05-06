function s_q_omega_info = read_s_q_omega(filename)

s = xmlread(filename);

mylist = s.getElementsByTagName('SQomega');

q_temp = zeros([1 mylist.getLength()]);
omega_temp = zeros([1 mylist.getLength()]);
Ss_temp = zeros([1 mylist.getLength()]);
Sd_temp = zeros([1 mylist.getLength()]);

n_q = 0; 
last_q = -1000; % q should never be negative I believe
for i = 0 : mylist.getLength()-1
  q_temp(i+1) = str2num(mylist.item(i).getAttribute('q'));
  omega_temp(i+1) = str2num(mylist.item(i).getAttribute('omega'));
  Ss_temp(i+1) = str2num(mylist.item(i).getAttribute('S-self'));
  Sd_temp(i+1) = str2num(mylist.item(i).getAttribute('S-diff'));
  if q_temp(i+1) > last_q
    n_q = n_q + 1;
    last_q = q_temp(i+1);
  end
end 

top_element = s.getElementsByTagName('s-q-omega');

n_time = length(q_temp) / n_q;

q = q_temp(1:n_q);
omega = zeros([1 n_time]);
S_d = zeros([n_q n_time]);
S_s = zeros([n_q n_time]);

time_index = 1;
for j = 1 : n_time
  q_index = 1;
  omega(time_index) = omega_temp(n_q*(j-1)+1);
  for i = 1 : n_q   
    S_d(q_index, time_index) = Sd_temp(n_q*(j-1)+i);
    S_s(q_index, time_index) = Ss_temp(n_q*(j-1)+i);
    q_index = q_index + 1;
  end
  time_index = time_index + 1;
end


if nargout == 1
  s_q_omega_info.S_s = S_s;
  s_q_omega_info.S_d = S_d;
  s_q_omega_info.q = q;
  s_q_omega_info.omega = omega;
 
  return;
end


subplot(1,3,1)
surf(q, omega, S_d')
xlabel('q [AA\^-1]')
ylabel('\omega 1/[10\^-13 s]')
zlabel('S_d (q,\omega)')
title(char(top_element.item(0).getAttribute('title')))

subplot(1,3,2)
surf(q, omega, S_s')
xlabel('q [AA\^-1]')
ylabel('\omega 1/[10\^-13 s]')
zlabel('S_s (q,\omega)')

subplot(1,3,3)
surf(q, omega, S_s'+S_d')
xlabel('q [AA\^-1]')
ylabel('\omega 1/[10\^-13 s]')
zlabel('S_s + S_d')