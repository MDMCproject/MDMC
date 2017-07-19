% Take as input two S(q,omega) files and plot these and differences
% between these, for example to compare a MDMC generated S(q,omega) file
% with a data S(q,omega). The required XML format is:
% <s-q-omega>
%
%   <SQomega q="0.42000" omega="0.00000" S-self="25.50111" S-diff="-24.30186" S="1.23" error="0.1" />
%   ... etc
% 
% It is currently assumed that filename1 contains S-self and S-diff and
% S_1(q,omega) = S-self + S-diff. filename2 may either contain S-self and S-diff
% or S attribute. If the latter is the case then S_2(q,omega) = S

function s_q_omega_info = diff_s_q_omega(filename1, filename2)

%%%% reading info from the 1st file %%%%

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

%top_element1 = s.getElementsByTagName('s-q-omega');
%title1 = char(top_element1.item(0).getAttribute('title'))

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

S_tot1 = S_s+S_d;

if nargout == 1
  s_q_omega_info.S_s = S_s;
  s_q_omega_info.S_d = S_d;
  s_q_omega_info.q = q;
  s_q_omega_info.omega = omega;
 
  return;
end


%%%% reading info from the 2nd file %%%%

s2 = xmlread(filename2);

mylist2 = s2.getElementsByTagName('SQomega');

if ( mylist.getLength() ~= mylist2.getLength() )
   disp('Number of SQomega elements in filename1 and filename2 must be the same')
   return
end

q_temp2 = zeros([1 mylist2.getLength()]);
omega_temp2 = zeros([1 mylist2.getLength()]);
Ss_temp2 = zeros([1 mylist2.getLength()]);
Sd_temp2 = zeros([1 mylist2.getLength()]);
Stot_temp2 = zeros([1 mylist2.getLength()]);
error2 = zeros([1 mylist2.getLength()]);

% if errors in filename2

errors_in_filename2 = mylist2.item(i).hasAttribute('error');

n_q = 0; 
last_q = -1000; % q should never be negative I believe
for i = 0 : mylist2.getLength()-1
  q_temp2(i+1) = str2num(mylist2.item(i).getAttribute('q'));
  omega_temp2(i+1) = str2num(mylist2.item(i).getAttribute('omega'));
  if mylist2.item(i).hasAttribute('S')
   if strcmp(mylist2.item(i).getAttribute('S'), 'no data') == 0
     Stot_temp2(i+1) = str2num(mylist2.item(i).getAttribute('S'));
     error2(i+1) = str2num(mylist2.item(i).getAttribute('error'));     
   else
     Stot_temp2(i+1) = NaN; 
     error2(i+1) = NaN;
   end
  else
   Ss_temp2(i+1) = str2num(mylist2.item(i).getAttribute('S-self'));
   Sd_temp2(i+1) = str2num(mylist2.item(i).getAttribute('S-diff'));
   Stot_temp2(i+1) = Ss_temp2(i+1) + Sd_temp2(i+1);
   if errors_in_filename2
     error2(i+1) = str2num(mylist2.item(i).getAttribute('error'));
   end
  end
  if q_temp2(i+1) > last_q
    n_q = n_q + 1;
    last_q = q_temp2(i+1);
  end
end 

%top_element2 = s2.getElementsByTagName('s-q-omega');
%title2 = char(top_element2.item(0).getAttribute('title'))

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
    S_tot2(q_index, time_index) = Stot_temp2(n_q*(j-1)+i);
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


% calculate FOM difference between the two input S(q, omega)

val = 0;
val_with_error = 0;
for i_q = 1 : length(q)
    for i_o = 1 : length(omega)
        if ~isnan(S_tot2(i_q,i_o))
            val = val + (S_tot2(i_q,i_o)-S_tot1(i_q,i_o))^2;
            if errors_in_filename2
              val_with_error = val_with_error + (S_tot2(i_q,i_o)-S_tot1(i_q,i_o))^2 ...
              / error2Plot(i_q,i_o)^2;
            end
        end
    end
end
FOM_val_without_errors = val
if errors_in_filename2
  FOM_val_with_errors = val_with_error
end


% now time for plotting

subplot(3,2,1)
surf(q, omega, S_tot1')
xlabel('q [AA\^-1]')
ylabel('\omega 1/[10\^-13 s]')
zlabel('S_1(q,\omega)')
colorbar

subplot(3,2,3)
surf(q, omega, S_tot2')
xlabel('q [AA\^-1]')
ylabel('\omega 1/[10\^-13 s]')
zlabel('S_2(q,\omega)')
colorbar

subplot(3,2,5)
surf(q, omega, S_tot2' - S_tot1')
xlabel('q [AA\^-1]')
ylabel('\omega 1/[10\^-13 s]')
zlabel('S_2 - S_1')
colorbar

subplot(3,2,2)

Sq1 = zeros([1 length(q)]);
for i = 1 : length(q)
    for j = 1 : length(omega)-1
        if isnan(S_tot1(i,j)) == 0 && isnan(S_tot1(i,j+1)) == 0
          Sq1(i) = Sq1(i)+(omega(j+1)-omega(j))*(S_tot1(i,j)+S_tot1(i,j+1))/2;
        end
    end
end

plot(q, Sq1, 'b');
hold on;

Sq2 = zeros([1 length(q)]);
for i = 1 : length(q)
    for j = 1 : length(omega)-1
        if isnan(S_tot2(i,j)) == 0 && isnan(S_tot2(i,j+1)) == 0
          Sq2(i) = Sq2(i)+(omega(j+1)-omega(j))*(S_tot2(i,j)+S_tot2(i,j+1))/2;
        end
    end
end

plot(q, Sq2, 'r--');
legend('\int S_1(q) d\omega', '\int S_2(q) d\omega')
xlabel('q [AA\^-1]')
hold off

subplot(3,2,4)

if errors_in_filename2
  surf(q, omega, abs(S_tot2' - S_tot1')./error2Plot')
  zlabel('abs(S_2-S_1)/error')
else
  surf(q, omega, abs(S_tot2' - S_tot1'))
  zlabel('abs(S_2-S_1)')
end
xlabel('q [AA\^-1]')
ylabel('\omega 1/[10\^-13 s]')

colorbar

subplot(3,2,6)
if errors_in_filename2
  surf(q, omega, error2Plot')
  xlabel('q [AA\^-1]')
  ylabel('\omega 1/[10\^-13 s]')
  zlabel('error on data')
  colorbar
end

