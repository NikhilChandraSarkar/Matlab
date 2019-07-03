%% Import data file
clear all;
%data=xlsread('training_data.xlsx');
data=xlsread('combined_data.xlsx');


%% plot the vehicle trajectories between certain space time region

t1=3;
t2=1808.5;

d1=0;
d2=245;
lat1=0;
lat2=12;

%data2consider=data(data(:,3)>= t1 & data(:,3) <= t2 & data(:,6)>= d1 & data(:,6)<=d2,:);
data2consider=data(data(:,3)>= t1 & data(:,3) <= t2 & data(:,6)>= d1 & data(:,6)<=d2 & data(:,9)>=lat1 & data(:,9)<=lat2,:);
c =unique(data2consider(:,1)); % Unique vehicle ID
cc=hsv(size(c,1)); % Colourmap
cc=sortrows(cc,1);
figure; 

for v=1:1:size(c,1)
      cv=c(v,1);
      vehDataBIN=ismember(data2consider(:,1),cv);
      cvehData=data2consider(vehDataBIN,:); % Unique vehicle data extraction
     %plot(cvehData(:,3),cvehData(:,6),'color',cc(v,:));
     plot(cvehData(:,3),cvehData(:,6),'Marker','.','MarkerSize',8,'Color','r','linestyle','-.');
      %plot(cvehData(:,3),cvehData(:,6),'color',[random('norm',0.5,0.1),random('norma',0.5,0.1), random('norma',0.4,0.1)], 'Marker','.','MarkerSize',8);
      %plot3(cvehData(:,6),cvehData(:,9),cvehData(:,3),'color',[random('norm',0.5,0.1),random('norma',0.5,0.1), random('norma',0.4,0.1)],'Marker','.','MarkerSize',8);
hold on;
end
%  xlabel('Longitudinal distance (m)')
%  ylabel('Lateral distance (m)')
%axis([0 60 0 245])
xlabel('Time (se)')
ylabel('Observed Longitudinal Distance (meter)')
% zlabel('Time (second)')
hold off

    axis([1732.5 1792.5 50 150])
    xticks([1732.5 1792.5])
    yticks([50 60 70 80 90 100 110 120 130 140 150])





%% Simualted Trajectory Plotting
data2=load('ReplicationNo1.mat');

%% plot the vehicle trajectories between certain space time region

t1=3;
t2=1808.5;

d1=0;
d2=245;
lat1=0;
lat2=12;

%data2consider=data(data(:,3)>= t1 & data(:,3) <= t2 & data(:,6)>= d1 & data(:,6)<=d2,:);
data2considerS=data2.SimulatedData(data2.SimulatedData(:,1)>= t1 & data2.SimulatedData(:,1) <= t2 & data2.SimulatedData(:,2)>= d1 & data2.SimulatedData(:,2)<=d2 & data2.SimulatedData(:,3)>=lat1 & data2.SimulatedData(:,3)<=lat2,:);
c2 =unique(data2considerS(:,8)); % Unique vehicle ID
cc2=hsv(size(c2,1)); % Colourmap
cc2=sortrows(cc2,1);
figure; 

for v=1:1:size(c2,1)
      cv2=c2(v,1);
      vehDataBINS=ismember(data2considerS(:,8),cv2);
      cvehDataS=data2considerS(vehDataBINS,:); % Unique vehicle data extraction
     
      
     % plot(cvehDataS(:,1),cvehDataS(:,2),'color',cc2(v,:));
      
      plot(cvehDataS(:,1),cvehDataS(:,2),'Marker','.','MarkerSize',8,'Color','b','linestyle','-.');
      
      %plot(cvehData(:,3),cvehData(:,6),'color',[random('norm',0.5,0.1),random('norma',0.5,0.1), random('norma',0.4,0.1)], 'Marker','.','MarkerSize',8);
      %plot3(cvehData(:,6),cvehData(:,9),cvehData(:,3),'color',[random('norm',0.5,0.1),random('norma',0.5,0.1), random('norma',0.4,0.1)],'Marker','.','MarkerSize',8);
hold on;
end
%  xlabel('Longitudinal distance (m)')
%  ylabel('Lateral distance (m)')
xlabel('Time (s)')
ylabel('Simulated Longitudinal distance (m)')
% zlabel('Time (sec)')
hold off
  
    axis([1732.5 1792.5 50 150])
    xticks([1732.5 1792.5])
    yticks([50 60 70 80 90 100 110 120 130 140 150])


