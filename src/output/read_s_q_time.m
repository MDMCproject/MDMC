function read_s_q_time(filename)

s = xmlread(filename);

mylist = s.getElementsByTagName('SQt');

q_temp = zeros([1 mylist.getLength()]);
t_temp = zeros([1 mylist.getLength()]);
Ss_temp = zeros([1 mylist.getLength()]);
Sd_temp = zeros([1 mylist.getLength()]);

n_q = 0; 
last_q = -1000; % q should never be negative I believe
for i = 0 : mylist.getLength()-1
  q_temp(i+1) = str2num(mylist.item(i).getAttribute('q'));
  t_temp(i+1) = str2num(mylist.item(i).getAttribute('t'));
  Ss_temp(i+1) = str2num(mylist.item(i).getAttribute('S-self'));
  Sd_temp(i+1) = str2num(mylist.item(i).getAttribute('S-diff'));
  if q_temp(i+1) > last_q
    n_q = n_q + 1;
    last_q = q_temp(i+1);
  end
end 

top_element = s.getElementsByTagName('s-q-time');

n_time = length(q_temp) / n_q;

q = q_temp(1:n_q);
t = zeros([1 n_time]);
S_d = zeros([n_q n_time]);
S_s = zeros([n_q n_time]);

time_index = 1;
for j = 1 : n_time
  q_index = 1;
  t(time_index) = t_temp(n_q*(j-1)+1);
  for i = 1 : n_q   
    S_d(q_index, time_index) = Sd_temp(n_q*(j-1)+i);
    S_s(q_index, time_index) = Ss_temp(n_q*(j-1)+i);
    q_index = q_index + 1;
  end
  time_index = time_index + 1;
end

subplot(1,2,1)
surf(q, t, S_d')
xlabel('q [AA\^-1]')
ylabel('t [10\^-13 s]')
zlabel('S\_d (q,t)')
title(char(top_element.item(0).getAttribute('title')))

subplot(1,2,2)
surf(q, t, S_s')
xlabel('q [AA\^-1]')
ylabel('t [10\^-13 s]')
zlabel('S\_s (q,t)')