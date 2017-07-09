% This was created to create a movie of the particular case of how the 
% atoms move during a creation of the data (MD-moves) for calculation of
% a time correction container.
%
% To use it is currently extremely primitive. You are required to uncomment 
% MDMC, which you can find by searching for 'movie'.
%
% This script may be modified to be used for similar purposes.

function make_movie(filename_in, n_frames)
filename = 'some text';
filename = filename_in;
bin_size = 0.1;
%bin_size = bin_size_in;
c_level = 2.0;
%c_level = c_level_in;


% read in raw structure and draw as a reference 




for ii = 1 : n_frames
  % read in data

  s = xmlread([filename num2str(ii) '.xml']);

  mylist = s.getElementsByTagName('atom');

  atom_pos = zeros([mylist.getLength() 3]);

  for i = 0 : mylist.getLength()-1
    atom_pos(i+1, 1) = str2num(mylist.item(i).getAttribute('x3'));
    atom_pos(i+1, 2) = str2num(mylist.item(i).getAttribute('y3'));
    atom_pos(i+1, 3) = str2num(mylist.item(i).getAttribute('z3'));
  end 

  box_element = s.getElementsByTagName('box-edges');

  a = str2num(box_element.item(0).getAttribute('x'));
  b = str2num(box_element.item(0).getAttribute('y'));
  c = str2num(box_element.item(0).getAttribute('z'));

  %if mod(length(r_temp), n_bin) 
  %  error('n_bin calculated incorrectly');
  %end

  %hgload('raw_structure.fig');

  alpha = pi/2;
  beta = pi/2;
  gamme = pi/2;

  x_min = 0;
  x_max = a/2;
  y_min = 0;
  y_max = b/2;
  z_min = 0;
  z_max = c/2;

  daspect([1 1 1])
  view(3)
  camlight; lighting phong
  axis([0 a/2 0 b/2 0 c/2])
  daspect([1 1 1])
  %alpha(p,.5)
  hold on;
  box on;

%    for iz = 0 : 2
%      for i = 0 : 2
%        plot3([0 a/2], [a/24+i*a/6 a/24+i*a/6], [a/24+iz*a/6 a/24+iz*a/6], 'linewidth', 1);
%      end
%    end
%      
%    for iz = 0 : 2
%      for i = 0 : 2
%        plot3([a/24+i*a/6 a/24+i*a/6], [0 a/2], [a/24+iz*a/6 a/24+iz*a/6], 'linewidth', 1);
%      end
%    end   
%  
%    for iz = 0 : 2
%      for i = 0 : 2
%        plot3([a/24+i*a/6 a/24+i*a/6], [a/24+iz*a/6 a/24+iz*a/6], [0 a/2], 'linewidth', 1);
%      end
%    end    
   
   
   scale = 0.3;
   [spX,spY,spZ] = sphere;
   for i = 1:length(atom_pos(:,1))
     if atom_pos(i,1) > x_min && atom_pos(i,1) < x_max
       if atom_pos(i,2) > y_min && atom_pos(i,2) < y_max
         if atom_pos(i,3) > z_min && atom_pos(i,3) < z_max
           h = surf(spX*scale+atom_pos(i,1),spY*scale+atom_pos(i,2),spZ*scale+atom_pos(i,3));
           set(h, 'edgecolor', 'none', 'facecolor', 'green', 'facelighting', 'phong');   
         end
       end
     end
   end
 
   %hgsave('raw_structure.fig'); return;
   
   set(gca, 'units', 'pixels')
   M(ii) = getframe(gcf, get(gca, 'position'));
   %mov = addframe(mov,M(ii));
   hold off;
   cla; %delete(gcf); %cla;
   
end

save movie.mat M
%mov = close(mov);

% drawing a unit cell box

% s2 = (cos(gamma) - cos(alpha)*cos(beta))/sin(alpha);
% s3 = cos(beta);
% s1 = sqrt(1 - s2*s2 - s3*s3);
% s4 = sin(alpha);
% s5 = cos(alpha);
% 
% UnitCell1 = [s1  s2  s3] * a;
% UnitCell2 = [0.0  s4  s5] * b;
% UnitCell3 = [0.0 0.0 1.0] * c;
% 
% plot3([0 UnitCell1(1)], [0 UnitCell1(2)], [0 UnitCell1(3)], 'linewidth', 2);
% plot3([0 UnitCell2(1)], [0 UnitCell2(2)], [0 UnitCell2(3)], 'linewidth', 2);
% plot3([0 UnitCell3(1)], [0 UnitCell3(2)], [0 UnitCell3(3)], 'linewidth', 2);
% plot3([UnitCell3(1) UnitCell3(1)+UnitCell1(1)], [UnitCell3(2) UnitCell3(2)+UnitCell1(2)], ...
%   [UnitCell3(3) UnitCell3(3)+UnitCell1(3)], 'linewidth', 2);
% plot3([UnitCell3(1) UnitCell3(1)+UnitCell2(1)], [UnitCell3(2) UnitCell3(2)+UnitCell2(2)], ...
%   [UnitCell3(3) UnitCell3(3)+UnitCell2(3)], 'linewidth', 2);
% plot3([UnitCell2(1) UnitCell2(1)+UnitCell1(1)], [UnitCell2(2) UnitCell2(2)+UnitCell1(2)], ...
%   [UnitCell2(3) UnitCell2(3)+UnitCell1(3)], 'linewidth', 2);
% plot3([UnitCell2(1) UnitCell2(1)+UnitCell3(1)], [UnitCell2(2) UnitCell2(2)+UnitCell3(2)], ...
%   [UnitCell2(3) UnitCell2(3)+UnitCell3(3)], 'linewidth', 2);
% plot3([UnitCell1(1) UnitCell1(1)+UnitCell2(1)], [UnitCell1(2) UnitCell1(2)+UnitCell2(2)], ...
%   [UnitCell1(3) UnitCell1(3)+UnitCell2(3)], 'linewidth', 2);
% plot3([UnitCell1(1) UnitCell1(1)+UnitCell3(1)], [UnitCell1(2) UnitCell1(2)+UnitCell3(2)], ...
%   [UnitCell1(3) UnitCell1(3)+UnitCell3(3)], 'linewidth', 2);
% 
% 
% % drawing the molecule as ball and stick
% 
%  dummy = size(connect_list);
%  for i = 1:dummy(1)
%    plot3([atom_pos(connect_list(i,1),1) atom_pos(connect_list(i,2),1)], ...
%      [atom_pos(connect_list(i,1),2) atom_pos(connect_list(i,2),2)], ...
%      [atom_pos(connect_list(i,1),3) atom_pos(connect_list(i,2),3)], ...
%      'y', 'linewidth', 2)
%  end


hold off;
xlabel('x')
ylabel('y')
zlabel('z')