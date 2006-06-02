function plot_many_rdf(how_many)



% read in target dataset

s = xmlread('xml_copy_of_robert_Ar.xml');

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

  if length(r) > length(r_target)    
    plot3(index, r(1:length(r_target)), g(1:length(r_target)))
    hold on
    plot3(index, r_target, g_target, 'r')
  else
    plot3(index, r, g)
    hold on
    plot3(index, r_target(1:length(r)), g_target(1:length(r)), 'r')    
  end
  
  index = index + 1;
  
  display('print out FOM')
  display('BE AWARE THAT FOR THIS TO WORK THE OBS DATA AND CAL DATA MUST') 
  display('BE STORED ON THE SAME BINNING SETUP')
  
  up_to = 10.0/0.1;   % = r-max element of rdf-fom / bin-length from rdf data file
  ii = ii
  FOM = sum((g(1:up_to)-g_target(1:up_to)).^2)/up_to
end
  
ylabel('r[0.1nm]')
zlabel('g(r)')


hold off