% Used to plot MDMC generated rdf functions

function read_rdf(filename, format)

s = xmlread(filename);

mylist = s.getElementsByTagName('g-of-r');

r = zeros([1 mylist.getLength()]);
g = zeros([1 mylist.getLength()]);

for i = 0 : mylist.getLength()-1
  mylist.item(i).getTagName();
  r(i+1) = str2num(mylist.item(i).getAttribute('r'));
  g(i+1) = str2num(mylist.item(i).getAttribute('g'));
end 

plot(r, g, format)
xlabel('r[0.1nm]')
ylabel('g(r)')