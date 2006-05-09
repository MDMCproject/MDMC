function plot_many_rdf(how_many)



% read in target dataset

s = xmlread('rdf_density0_02125temp85.xml');

mylist = s.getElementsByTagName('g-of-r');

  g_target = zeros([1 mylist.getLength()]);
  r_target = zeros([1 mylist.getLength()]);

%g = zeros([1 mylist.getLength()]);

for i = 0 : mylist.getLength()-1
  mylist.item(i).getTagName();
  r_target(i+1) = str2num(mylist.item(i).getAttribute('r'));
  g_target(i+1) = str2num(mylist.item(i).getAttribute('g'));
end 



for ii = 1 : how_many
  s = xmlread(['rdf' num2str(ii) '.xml']);

  mylist = s.getElementsByTagName('g-of-r');

  if ii == 1
    g = zeros([1 mylist.getLength()]);
    r = zeros([1 mylist.getLength()]);
    index = ones([1 mylist.getLength()]);
  end
  
  %g = zeros([1 mylist.getLength()]);

  for i = 0 : mylist.getLength()-1
    mylist.item(i).getTagName();
    r(i+1) = str2num(mylist.item(i).getAttribute('r'));
    g(i+1) = str2num(mylist.item(i).getAttribute('g'));
  end 
  
  plot3(index, r, g)
  hold on
  plot3(index, r_target, g_target, 'r')
  
  index = index + 1;
  
  % print out FOM
  
  up_to = 10/0.05;
  ii = ii
  FOM = sum((g(1:up_to)-g_target(1:up_to)).^2)/up_to
end
  
ylabel('r[0.1nm]')
zlabel('g(r)')


hold off