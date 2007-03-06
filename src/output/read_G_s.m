function read_G_s(filename)

s = xmlread(filename);

mylist = s.getElementsByTagName('G-s');

t_temp = zeros([1 mylist.getLength()]);
r_temp = zeros([1 mylist.getLength()]);
g_temp = zeros([1 mylist.getLength()]);

for i = 0 : mylist.getLength()-1
  r_temp(i+1) = str2num(mylist.item(i).getAttribute('r'));
  t_temp(i+1) = str2num(mylist.item(i).getAttribute('t'));
  g_temp(i+1) = str2num(mylist.item(i).getAttribute('G'));
end 

top_element = s.getElementsByTagName('G_s-space-time-pair-correlation-function');

n_bin = ceil(max(r_temp) / str2num(top_element.item(0).getAttribute('bin-length')));

if mod(length(r_temp), n_bin) 
  error('n_bin calculated incorrectly');
end

n_time = length(r_temp) / n_bin;

r = r_temp(1:n_bin);
t = zeros([1 n_time]);
G_s = zeros([n_bin n_time]);


time_index = 1;
for j = 1 : n_time
  bin_index = 1;
  t(time_index) = t_temp(n_bin*(j-1)+1);
  for i = 1 : n_bin   
    G_s(bin_index, time_index) = g_temp(n_bin*(j-1)+i);
    bin_index = bin_index + 1;
  end
  time_index = time_index + 1;
end


subplot(2,1,1)
surf(r, t, G_s')
xlabel('r [AA]')
ylabel('t [10\^-13 s]')
zlabel('G\_s (r,t)')
title(char(top_element.item(0).getAttribute('title')))
subplot(2,1,2)
t_lower_cut = floor(length(t)/5)*4;
t_new = t(t_lower_cut:end);
G_s_new = G_s(:,t_lower_cut:end);
surf(r, t_new, G_s_new')
xlabel('r [AA]')
ylabel('t [10\^-13 s]')
zlabel('G\_s (r,t)')