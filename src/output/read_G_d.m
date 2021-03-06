% Used to plot MDMC generated g_d space-time pair correlation functions. 
% As of this writing MDMC outputs the "normalised" version of this function
% which has the property that (g^norm)_d(r,t) -> 1 when t->infinity and
% r->infinity. Note (g^norm)_d(r,t) = N * g_d(r,t) / (N-1).
%
% Note the space-time correlation function is related to the space-time
% pair correlation function by G_d = rho * g_d, where rho is N/V.

function read_g_d(filename)

s = xmlread(filename);

mylist = s.getElementsByTagName('G-d');

t_temp = zeros([1 mylist.getLength()]);
r_temp = zeros([1 mylist.getLength()]);
g_temp = zeros([1 mylist.getLength()]);

for i = 0 : mylist.getLength()-1
  r_temp(i+1) = str2num(mylist.item(i).getAttribute('r'));
  t_temp(i+1) = str2num(mylist.item(i).getAttribute('t'));
  g_temp(i+1) = str2num(mylist.item(i).getAttribute('G'));
end 

top_element = s.getElementsByTagName('G_d-space-time-pair-correlation-function');

n_bin = ceil(max(r_temp) / str2num(top_element.item(0).getAttribute('bin-length')));

if mod(length(r_temp), n_bin) 
  error('n_bin calculated incorrectly');
end

n_time = length(r_temp) / n_bin;

r = r_temp(1:n_bin);
t = zeros([1 n_time]);
G_d = zeros([n_bin n_time]);


time_index = 1;
for j = 1 : n_time
  bin_index = 1;
  t(time_index) = t_temp(n_bin*(j-1)+1);
  for i = 1 : n_bin   
    G_d(bin_index, time_index) = g_temp(n_bin*(j-1)+i);
    bin_index = bin_index + 1;
  end
  time_index = time_index + 1;
end

surf(r, t, G_d')
xlabel('r [AA]')
ylabel('t [10\^-13 s]')
zlabel('"normalised" g_d (r,t)')
title(char(top_element.item(0).getAttribute('title')))