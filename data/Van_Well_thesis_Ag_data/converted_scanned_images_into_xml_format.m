% read in 1st column

fid = fopen('table_V_column1.txt','r');

num_omega_1 = 12;
S_Q_Omega_1 = zeros([1 num_omega_1]);
S_Q_Omega_1 = NaN;
sigma_S_Q_Omega_1 = zeros([1 num_omega_1]);
sigma_S_Q_Omega_1 = NaN;

for i = 1 : num_omega_1
  dummy = fgetl(fid);
  
  fisse = findstr(dummy, ')');
  fisse_begin = findstr(dummy, '(');
  
  S_Q_Omega_1(1, i) = str2double(dummy(1:6));
  sigma_S_Q_Omega_1(1, i) = str2double(dummy((fisse_begin+1):(fisse-1)))*0.0001;
end
  
fclose(fid);


% read in 2st column

fid = fopen('table_V_column2.txt','r');

num_omega_2 = 29;
S_Q_Omega_2 = zeros([1 num_omega_2]);
S_Q_Omega_2 = NaN;
sigma_S_Q_Omega_2 = zeros([1 num_omega_2]);
sigma_S_Q_Omega_2 = NaN;

for i = 1 : num_omega_2
  dummy = fgetl(fid);
  
  fisse = findstr(dummy, ')');
  fisse_begin = findstr(dummy, '(');
  
  S_Q_Omega_2(1, i) = str2num(dummy(1:6));
  sigma_S_Q_Omega_2(1, i) = str2num(dummy((fisse_begin+1):(fisse-1)))*0.0001;
end
  
fclose(fid);


% read in column 3 to 8

fid = fopen('table_V_column3_to_10.txt','r');

num_omega_3 = 37;
S_Q_Omega_3 = zeros([8 num_omega_3]);
S_Q_Omega_3 = NaN;
sigma_S_Q_Omega_3 = zeros([8 num_omega_3]);
sigma_S_Q_Omega_3 = NaN;

i_omega = 0;
while 1
  dummy = fgetl(fid);
  
  if dummy == -1
    break;
  end
  
  i_q = 1;
  i_omega = i_omega + 1;
  
  fisse = findstr(dummy, ')');
  fisse_begin = findstr(dummy, '(');
  
  S_Q_Omega_3(i_q, i_omega) = str2num(dummy(1:6));
  sigma_S_Q_Omega_3(i_q, i_omega) = str2num(dummy((fisse_begin(1)+1):(fisse(1)-1)))*0.0001;
  
  for i_q = 2 : 8
    %kus = str2num(dummy((fisse(i_q-1)+2):(fisse(i_q-1)+2+5)))
    S_Q_Omega_3(i_q, i_omega) = str2num(dummy((fisse(i_q-1)+2):(fisse(i_q-1)+2+5)));
    sigma_S_Q_Omega_3(i_q, i_omega) = str2num(dummy((fisse_begin(i_q)+1):(fisse(i_q)-1)))*0.0001;
  end
  
end

fclose(fid);


% put everything in one matrix

S_Q_Omega = zeros([10 num_omega_3]);
%S_Q_Omega = NaN;
sigma_S_Q_Omega = zeros([10 num_omega_3]);
%sigma_S_Q_Omega = NaN;

S_Q_Omega(1,1:num_omega_1) = S_Q_Omega_1;
S_Q_Omega(1,num_omega_1+1:end) = NaN;
S_Q_Omega(2,1:num_omega_2) = S_Q_Omega_2;
S_Q_Omega(2,num_omega_2+1:end) = NaN;
S_Q_Omega(3:10,1:num_omega_3) = S_Q_Omega_3;

sigma_S_Q_Omega(1,1:num_omega_1) = sigma_S_Q_Omega_1;
sigma_S_Q_Omega(1,num_omega_1+1:end) = NaN;
sigma_S_Q_Omega(2,1:num_omega_2) = sigma_S_Q_Omega_2;
sigma_S_Q_Omega(2,num_omega_2+1:end) = NaN;
sigma_S_Q_Omega(3:10,1:num_omega_3) = sigma_S_Q_Omega_3;

% To convert van Well S(q,omega) from having
% unit of ps to having unit of 0.1ps=10^-13s
S_Q_Omega = 10*S_Q_Omega;
sigma_S_Q_Omega = 10*sigma_S_Q_Omega;

% Create the XML representation of Table V

% q values as read from Table V converted to inverse AA
q = [4.2 8.4 12.6 16.2 19.2 22.2 25.8 30.0 34.2 39.0]/10;
% omega values as read from Table V in units of ps
omega = [0.0:0.2:2.0];
omega = [omega 2.5:0.5:10.0];
omega = [omega 11:1:20];
omega = omega / 10.0;  % to convert from 1/ps to 1/[10^-13 s]

docNode = com.mathworks.xml.XMLUtils.createDocument('s-q-omega')
docRootNode = docNode.getDocumentElement;
docRootNode.setAttribute('n-omega-points' , num2str(length(omega)));
docRootNode.setAttribute('n-q-points' , num2str(length(q)));
docRootNode.setAttribute('q-unit' , 'AA^-1');
docRootNode.setAttribute('omega-unit' , '1/[10^-13 s]');
docRootNode.setAttribute('S-unit' , '0.1ps');
%docRootNode.setAttribute('title' , ['Copy of data from A. A. Van Well thesis.' ...
%    ' Table V page 41-2. Ag S(Q,Omega) data.']);
thisElement = docNode.createElement('description');
thisElement.appendChild... 
    (docNode.createTextNode(['Copy of data from A. A. Van Well thesis.' ...
    ' Table V page 41-2. Ag S(Q,Omega) data.']));
docRootNode.appendChild(thisElement); 
for j = 1 : length(omega)
  for i = 1 : length(q)
    omegaElement = docNode.createElement('SQomega');
    omegaElement.setAttribute('q', num2str(q(i)));
    omegaElement.setAttribute('omega', num2str(omega(j)));
    if isnan(S_Q_Omega(i,j))
      omegaElement.setAttribute('S', 'no data');     
    else
      omegaElement.setAttribute('S', num2str(S_Q_Omega(i,j)));
      omegaElement.setAttribute('error', num2str(sigma_S_Q_Omega(i,j), '%6.4f'));
    end
    docRootNode.appendChild(omegaElement);
  end
  %docRootNode.appendChild(thisElement);
end
 
% Save the sample XML document.
xmlFileName = 'Well_s_q_omega_Ag_data.xml';
xmlwrite(xmlFileName,docNode);
%edit(xmlFileName);
mesh(q,omega, S_Q_Omega')
title('Van Well Ag S(q,\omega )')
xlabel('q [AA^{-1}]')
ylabel('\omega  [1/(0.1ps)]')
zlabel('S(q,\omega) [0.1ps]')