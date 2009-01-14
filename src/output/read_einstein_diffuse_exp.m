function read_einstein_diffuse_exp(filename, format)

s = xmlread(filename);

mylist = s.getElementsByTagName('D-of-t');

t = zeros([1 mylist.getLength()]);
D = zeros([1 mylist.getLength()]);

for i = 0 : mylist.getLength()-1
  mylist.item(i).getTagName();
  t(i+1) = str2num(mylist.item(i).getAttribute('t'));
  D(i+1) = str2num(mylist.item(i).getAttribute('D'));
end 


top_element = s.getElementsByTagName('einstein-diffuse-exp');

n_atom = str2num(top_element.item(0).getAttribute('n-atom'));

density = str2double(top_element.item(0).getAttribute('density'));


subplot(2,1,1)
plot(t, D, format)
xlabel('t[10\^-13 s]')
ylabel('Einstein diffuse func')
title_element = s.getElementsByTagName('einstein-diffuse-exp');
title(char(title_element.item(0).getAttribute('title')))

subplot(2,1,2)

% to convert from diffusion constant as defined in code and on pages
% 19b and 19bb in notes then the mean value of squared differences between
% atoms at t=0 and the atoms at some later t is

meanSquareDist = D.*t*6;

plot(t, meanSquareDist, format)
xlabel('t[10\^-13 s]')
ylabel('<(r(t)-r(0)^2> [(0.1nm)^2]')
title(['Natom=' num2str(n_atom) '  box-length=' num2str( (n_atom/density)^(1/3) )])