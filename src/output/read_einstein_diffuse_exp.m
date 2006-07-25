function read_einstein_diffuse_exp(filename, format)

s = xmlread(filename);

mylist = s.getElementsByTagName('D-of-t');

r = zeros([1 mylist.getLength()]);
g = zeros([1 mylist.getLength()]);

for i = 0 : mylist.getLength()-1
  mylist.item(i).getTagName();
  r(i+1) = str2num(mylist.item(i).getAttribute('t'));
  g(i+1) = str2num(mylist.item(i).getAttribute('D'));
end 

plot(r, g, format)
xlabel('t[10\^-13 s]')
ylabel('Einstein diffuse func')

title_element = s.getElementsByTagName('einstein-diffuse-exp');
title(char(title_element.item(0).getAttribute('title')))