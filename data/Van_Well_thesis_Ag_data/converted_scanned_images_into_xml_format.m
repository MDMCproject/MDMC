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

% Convert all zero errors into upper min error value
sigma_S_Q_Omega( find(sigma_S_Q_Omega == 0) ) = 0.001;

% Create the XML representation of Table V

% q values as read from Table V converted to inverse AA
q = [4.2 8.4 12.6 16.2 19.2 22.2 25.8 30.0 34.2 39.0]/10;
% omega values as read from Table V in units of ps
omega = [0.0:0.2:2.0];
omega = [omega 2.5:0.5:10.0];
omega = [omega 11:1:20];
omega = omega / 10.0;  % to convert from 1/ps to 1/[10^-13 s]


% create symmetried van Well data (i.e. as van Well's presents his data)

docNode = com.mathworks.xml.XMLUtils.createDocument('s-q-omega');
docRootNode = docNode.getDocumentElement;
docRootNode.setAttribute('n-omega-points' , num2str(length(omega)));
docRootNode.setAttribute('n-q-points' , num2str(length(q)));
docRootNode.setAttribute('q-unit' , 'AA^-1');
docRootNode.setAttribute('omega-unit' , '1/[10^-13 s]');
docRootNode.setAttribute('S-unit' , '0.1ps');
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
end
 
% Save symmetries data
xmlFileName = 'Well_s_q_omega_Ag_data.xml';
xmlwrite(xmlFileName,docNode);
clear docNode;

% create unsymmetried van Well data

kB_T = 10.3408; % at t=120K
omega_convertion = 6.582; % i.e. omega = 1*10^13 s^-1   translate to   hbar*omega = 6.582 meV

docNode_un = com.mathworks.xml.XMLUtils.createDocument('s-q-omega');
docRootNode_un = docNode_un.getDocumentElement;
docRootNode_un.setAttribute('n-omega-points' , num2str(length(omega)));
docRootNode_un.setAttribute('n-q-points' , num2str(length(q)));
docRootNode_un.setAttribute('q-unit' , 'AA^-1');
docRootNode_un.setAttribute('omega-unit' , '1/[10^-13 s]');
docRootNode_un.setAttribute('S-unit' , '0.1ps');
thisElement_un = docNode_un.createElement('description');
thisElement_un.appendChild... 
    (docNode_un.createTextNode(['Un-symmetried copy of data from A. A. Van Well thesis.' ...
    ' Table V page 41-2. Ag S(Q,Omega) data.']));
docRootNode_un.appendChild(thisElement_un); 
for j = 1 : length(omega)
  for i = 1 : length(q)
    omegaElement_un = docNode_un.createElement('SQomega');
    omegaElement_un.setAttribute('q', num2str(q(i)));
    omegaElement_un.setAttribute('omega', num2str(omega(j)));
    if isnan(S_Q_Omega(i,j))
      omegaElement_un.setAttribute('S', 'no data');     
    else
      omegaElement_un.setAttribute('S', num2str(S_Q_Omega(i,j)*exp(0.5*omega(j)*omega_convertion/kB_T)));
      omegaElement_un.setAttribute('error', num2str(sigma_S_Q_Omega(i,j)*exp(0.5*omega(j)*omega_convertion/kB_T), '%6.4f'));
    end
    docRootNode_un.appendChild(omegaElement_un);
  end
end
 
% Save un-symmetries data
xmlFileName_un = 'Well_s_q_omega_Ag_data_unsymmetried.xml';
xmlwrite(xmlFileName_un,docNode_un);



subplot(1,3,1)
surf(q,omega, S_Q_Omega')
title('Van Well Ag S(q,\omega )')
xlabel('q [AA^{-1}]')
ylabel('\omega  [10^{13} s^{-1}]')
zlabel('S(q,\omega) [10^{-13} s]')

subplot(1,3,2)
S_Q_Omega_un = zeros(size(S_Q_Omega));
for i = 1 : length(q)
   S_Q_Omega_un(i,:) = S_Q_Omega(i,:) .* exp(0.5*omega*omega_convertion/kB_T);
end
surf(q,omega, S_Q_Omega_un')
title('un-symmetried S(q,\omega )')
xlabel('q [AA^{-1}]')
ylabel('\omega  [10^{13} s^{-1}]')
zlabel('[10^{-13} s]')

subplot(1,3,3)
surf(q,omega, S_Q_Omega_un' - S_Q_Omega')
title('unsymmetried S(q,\omega ) - S(q,\omega )')
xlabel('q [AA^{-1}]')
ylabel('\omega  [10^{13} s^{-1}]')
zlabel('[10^{-13} s]')