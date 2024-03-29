% Used to plot MDMC generated h_d functions

function h_d_info = read_h_d(filename)

s = xmlread(filename);

mylist = s.getElementsByTagName('h-d');

t_temp = zeros([1 mylist.getLength()]);
r_temp = zeros([1 mylist.getLength()]);
h_temp = zeros([1 mylist.getLength()]);

for i = 0 : mylist.getLength()-1
  r_temp(i+1) = str2num(mylist.item(i).getAttribute('r'));
  t_temp(i+1) = str2num(mylist.item(i).getAttribute('t'));
  h_temp(i+1) = str2num(mylist.item(i).getAttribute('h'));
end 

top_element = s.getElementsByTagName('normalised-h-d-histogram');

n_bin = ceil(max(r_temp) / str2num(top_element.item(0).getAttribute('bin-length')));

if mod(length(r_temp), n_bin) 
  error('n_bin calculated incorrectly');
end

n_time = length(r_temp) / n_bin;

r = r_temp(1:n_bin);
t = zeros([1 n_time]);
h_d = zeros([n_bin n_time]);


time_index = 1;
for j = 1 : n_time
  bin_index = 1;
  t(time_index) = t_temp(n_bin*(j-1)+1);
  for i = 1 : n_bin   
    h_d(bin_index, time_index) = h_temp(n_bin*(j-1)+i);
    bin_index = bin_index + 1;
  end
  time_index = time_index + 1;
end

surf(r, t, h_d')
xlabel('r [AA]')
ylabel('t [10\^-13 s]')
zlabel('h\_d (r,t)')
title(char(top_element.item(0).getAttribute('title')))

h_d_info.r = r;
h_d_info.t = t;
h_d_info.val = h_d;
h_d_info.bin_length = str2num(top_element.item(0).getAttribute('bin-length'));
h_d_info.n_bin = n_bin;
