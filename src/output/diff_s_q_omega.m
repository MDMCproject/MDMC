% here it is implicitely assumed that filename2 refers to the data

function s_q_omega_info = diff_s_q_omega(filename1, filename2)

s = xmlread(filename1);

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

s2 = xmlread(filename2);

mylist2 = s2.getElementsByTagName('SQomega');

q_temp2 = zeros([1 mylist2.getLength()]);
omega_temp2 = zeros([1 mylist2.getLength()]);
%Ss_temp2 = zeros([1 mylist2.getLength()]);
%Sd_temp2 = zeros([1 mylist2.getLength()]);
Ssum2 = zeros([1 mylist2.getLength()]);
error2 = zeros([1 mylist2.getLength()]);

n_q = 0; 
last_q = -1000; % q should never be negative I believe
for i = 0 : mylist2.getLength()-1
  q_temp2(i+1) = str2num(mylist2.item(i).getAttribute('q'));
  omega_temp2(i+1) = str2num(mylist2.item(i).getAttribute('omega'));
%  Ss_temp2(i+1) = str2num(mylist2.item(i).getAttribute('S-self'));
%  Sd_temp2(i+1) = str2num(mylist2.item(i).getAttribute('S-diff'));
   if strcmp(mylist2.item(i).getAttribute('S'), 'no data') == 0
     Ssum2(i+1) = str2num(mylist2.item(i).getAttribute('S'));
     error2(i+1) = str2num(mylist2.item(i).getAttribute('error'));
   else
     Ssum2(i+1) = NaN; 
     error2(i+1) = NaN;
   end
  if q_temp2(i+1) > last_q
    n_q = n_q + 1;
    last_q = q_temp2(i+1);
  end
end 

top_element2 = s2.getElementsByTagName('s-q-omega');

n_time2 = length(q_temp2) / n_q;

q2 = q_temp2(1:n_q);
omega2 = zeros([1 n_time2]);
%S_d2 = zeros([n_q n_time2]);
%S_s2 = zeros([n_q n_time2]);
S_tot2 = zeros([n_q n_time2]);
error2Plot = zeros([n_q n_time2]);

time_index = 1;
for j = 1 : n_time2
  q_index = 1;
  omega2(time_index) = omega_temp2(n_q*(j-1)+1);
  for i = 1 : n_q   
    %S_d2(q_index, time_index) = Sd_temp2(n_q*(j-1)+i);
    %S_s2(q_index, time_index) = Ss_temp2(n_q*(j-1)+i);
    S_tot2(q_index, time_index) = Ssum2(n_q*(j-1)+i);
    error2Plot(q_index, time_index) = error2(n_q*(j-1)+i);
    q_index = q_index + 1;
  end
  time_index = time_index + 1;
end


if nargout == 1
  %s_q_omega_info2.S_s = S_s2;
  %s_q_omega_info2.S_d = S_d2;
  s_q_omega_info2.S_tot = S_tot2;
  s_q_omega_info2.error = error2Plot;
  s_q_omega_info2.q = q2;
  s_q_omega_info2.omega = omega2;
 
  return;
end

subplot(3,2,1)
surf(q, omega, S_s'+S_d')
xlabel('q [AA\^-1]')
ylabel('\omega 1/[10\^-13 s]')
zlabel('S_1(q,\omega)')
colorbar

subplot(3,2,3)
surf(q, omega, S_tot2')
xlabel('q [AA\^-1]')
ylabel('\omega 1/[10\^-13 s]')
zlabel('S_2(q,\omega)')
title(char(top_element2.item(0).getAttribute('title')))
colorbar

subplot(3,2,5)
surf(q, omega, S_tot2'- (S_s'+S_d'))
xlabel('q [AA\^-1]')
ylabel('\omega 1/[10\^-13 s]')
zlabel('S_2 - S_1')
colorbar

subplot(3,2,2)

S_tot1 = S_s+S_d;

Sq1 = zeros([1 length(q)]);
for i = 1 : length(q)
    for j = 1 : length(omega)-1
        if isnan(S_tot1(i,j)) == 0 && isnan(S_tot1(i,j+1)) == 0
          Sq1(i) = Sq1(i)+(omega(j+1)-omega(j))*(S_tot1(i,j+1)+S_tot1(i,j+1))/2;
        end
    end
end

plot(q,Sq1);
xlabel('q [AA\^-1]')
ylabel('\int S_1(q) d\omega')


subplot(3,2,4)

Sq2 = zeros([1 length(q)]);
for i = 1 : length(q)
    for j = 1 : length(omega)-1
        if isnan(S_tot1(i,j)) == 0 && isnan(S_tot1(i,j+1)) == 0
          Sq2(i) = Sq2(i)+(omega(j+1)-omega(j))*(S_tot2(i,j+1)+S_tot2(i,j+1))/2;
        end
    end
end

plot(q,Sq2);
xlabel('q [AA\^-1]')
ylabel('\int S_2(q) d\omega')

subplot(3,2,6)
surf(q, omega, error2Plot')
xlabel('q [AA\^-1]')
ylabel('\omega 1/[10\^-13 s]')
zlabel('error on data')
colorbar

