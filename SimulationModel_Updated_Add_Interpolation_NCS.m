 %% Simulation Model for Area-Based Traffic Flow for Chennai dataset
 % Model is develeoped by SARKAR NIKHIL CHANDRA; Date: 1/10/2018
 % Domain: QUT
 
 clear all
 
 % 
 
 %% Empirical Macroscopic properties for empirical dataset
TableUpdate_E=[];
EmpiricalData=[];
SimulationData=xlsread('/Applications/NikhilChandraSarkar/WholeSimulationDataForOthersVehicles.xlsx');

SimulationDataM=xlsread('/Applications/NikhilChandraSarkar/WholeSimulationDataForMotorcycle.xlsx');
EmpiricalDataOthers=SimulationData(1:end,1:12);
EmpiricalDataCar=SimulationDataM(1:end,1:12);
EmpiricalData=[EmpiricalDataOthers;EmpiricalDataCar];

% EmpiricalData=xlsread('combined_data.xlsx'); % whole data
% EmpiricalData=xlsread('validation_data.xlsx'); % Validation data
 To_E=52.5; % starting time
 dt_E=60; % time interval (second)
 Tend_E=1808.5-dt_E; % Ending time;
 d1_E=50;    % starting distance
 d2_E=195;  % ending distance 
 
for t=To_E:dt_E:Tend_E
t1_E=t;
t2_E=t+dt_E;
dt_E=t2_E-t1_E;
Time_E =t2_E;
area_E=(t2_E-t1_E)*(d2_E-d1_E); % rectanguar area (time interval X length of the road segment)
data2region_E=EmpiricalData(EmpiricalData(:,3)>= t1_E & EmpiricalData(:,3) <= t2_E & EmpiricalData(:,6)>= d1_E & EmpiricalData(:,6)<=d2_E,:);
VehID_E = data2region_E(:,1);
UniqueVehID_E = unique(data2region_E(:,1)); % Unique vehicle ID
N_E =numel(UniqueVehID_E); 
TT_E=[];
TD_E=[];

for i=min(VehID_E):max(VehID_E)
data2eachVeh_E=data2region_E(data2region_E(:,1)==i,:); % Individual vehicle data
TimeEnd_E = max(data2eachVeh_E(:,3));
TimeStart_E = min(data2eachVeh_E(:,3));
TravelTime_E = TimeEnd_E -TimeStart_E;
DistEnd_E = max(data2eachVeh_E(:,6));
DistStart_E = min(data2eachVeh_E(:,6));
TravelDist_E = DistEnd_E-DistStart_E;
TT_E=[TT_E;TravelTime_E]; % Travel time for each vehicle
TD_E=[TD_E;TravelDist_E]; % Travel distance for each vehicle
end
% % TT_E_Int=[];
% % TD_E_Int=[];
% % 
% % TT_E_Int=[TT_E_Int;TT_E];
% % TD_E_Int=[TD_E_Int;TD_E];
% % 
% % % Interpolation adding from here
% % 
% % TT_E_Int=linspace(min(TT_E),max(TT_E),numel(TT_E));
% % 
% % TD_E_Int=linspace(min(TD_E),max(TD_E),numel(TD_E));
% % 
% % TT_E_New =min(TT_E_Int):0.01:max(TT_E_Int);
% % 
% % TD_E_New=interp1(TT_E_Int,TD_E_Int,TT_E_New,'linear');


% %   TTT_E = sum(TT_E_New);
% %   TDT_E= sum(TD_E_New);

   TTT_E = sum(TT_E);
   TDT_E= sum(TD_E);

density_E = TTT_E/area_E; %(veh/m)
Density_E=(density_E*1000); %(Veh/km)
flow_E = TDT_E/area_E;
Flow_E =flow_E*3600; % (Veh/h)
Table_E =[];
Table_E=[Time_E TTT_E TDT_E N_E Density_E Flow_E]; % time, total time taken, total distance taken, No. of vehicles, Density, Flow
TableUpdate_E=[TableUpdate_E;Table_E];
end
Speed_E =TableUpdate_E(:,6)./TableUpdate_E(:,5);





% %%Plotting macroscopic properties
% 
% figure
% plot(TableUpdate_E(:,5),TableUpdate_E(:,6),'.','MarkerSize', 20) % 2D plot
% axis([0 300 0 10000])
% %plot3(TableUpdate(:,5),TableUpdate(:,6),TableUpdate(:,1),'.','MarkerSize', 20)% 3D plot
%   xlabel('Density (veh/km)')
%   ylabel('Flow (veh/h)')
% %zlabel('Time (s)')
% 
% figure
% plot(TableUpdate_E(:,2),TableUpdate_E(:,3),'o')
% xlabel('Total time taken (TTT) (s)')
% ylabel('Total distance travel (TDT) (m)')
% 
% figure
% plot(TableUpdate_E(:,5),Speed_E,'.','MarkerSize', 20)
% axis([0 300 0 100])
% xlabel('Density (veh/km)')
% ylabel('Speed (km/h)')
% 
% figure
% plot(TableUpdate_E(:,6),Speed_E,'.','MarkerSize', 20)
% axis([0 10000 0 100])
% xlabel('Flow (veh/h)')
% ylabel('Speed (km/h)')
% 
% figure
%  plot(TableUpdate_E(:,1),TableUpdate_E(:,5),'-r.','MarkerSize', 12)% time vs density graph
%   xlabel('Time (s)')
%   ylabel('Density (veh/km)')
%   figure
%   plot(TableUpdate_E(:,1),TableUpdate_E(:,6),'-r.','MarkerSize', 12)% time vs flow graph
%   xlabel('Time (s)')
%   ylabel('Flow (veh/h)')
%   figure
%   plot(TableUpdate_E(:,1),Speed_E,'-r.','MarkerSize', 12)% time vs speed graph
%   xlabel('Time (s)')
%   ylabel('Speed (km/h)')
 
%% MIDM Development
 
%SimulationData=xlsread('CalibrationData_MIDM.xlsx');  % Simualtion Car data for MIDM
 %SimulationData=xlsread('ValidationDataForOthersVehicles.xlsx');
 
 SimulationData=xlsread('/Applications/NikhilChandraSarkar/WholeSimulationDataForOthersVehicles.xlsx'); 
 Eparameters=xlsread('/Applications/NikhilChandraSarkar/Eparameters_Updated.xlsx'); % Calibrated data for MIDM
 VId_min=min(SimulationData(:,1));
 VId_max=max(SimulationData(:,1));
 F_idm_dataset=[];
 TT=[];
 TD=[];
 
  for k= VId_min:VId_max
 
 data=SimulationData(SimulationData(:,1)==k,:); % Extract specific car data from Simulation data set
 
 if isempty(data)
     DataEmpty=data;
 else
     data=data;
     k
 
 % Random selection of a raw from estimated parameters
    idx = randperm(size(Eparameters,1),1); % Index of selected row
    B = Eparameters(idx,:); % Randomly selected a row from a matrix of estimated parameters
    v0 = (5/18)*B(1,3);
    a_max = B(1,4);
    b = B(1,5); % gamrnd(1.0,std(Eparameters(:,5)));
    t0=B(1,6); % Safety time headway (sec)
    s0=B(1,7); % Jam distance (m)
    s1=B(1,8); % Non-linear Jam distance (m)
    delta=B(1,9); % Acceleration exponent
 
 % Initial values of Modified IDM (MIDM)input:
 
 delT=0.5; % (sec)% Time increment
 tstart=data(1,3); % Starting time
 tend=data(end,3); % Ending time
 trange=data(:,3)>=tstart & data(:,3)<=tend; % Time duration
 length= size(data(trange,:),1);
 
 dist(1,1)=sqrt(data(1,6).^2+data(1,9).^2); % resultant position of subject vehicle
 F_veh=[data(trange,3) data(trange,6) data(trange,9) data(trange,7) data(trange,10) NaN(length,1) ]; % time, longdist,latdist,longspeed, latspeed
 F_veh(1,6)=sqrt(data(1,6).^2+data(1,9).^2); % resultant position for real subject vehicle
 
 F_idm = [data(trange,3) NaN(length,7)];

 F_idm(1,2)=F_veh(1,2); % initial longitudinal position of simulated follower
 F_idm(1,3)=F_veh(1,3); % initial lateral position of simulated follower

F_idm(1,4)=data(1,7);
F_idm(1,5)=data(1,10);

% Initial resultant speed for subject vehicle

%Column Indext for Predicted Alternative
CID_Pred_Alt=27;
%Column Index For Observed Alternative
OID_Obs_Alt=14;

if data(1,CID_Pred_Alt)==1
v_R=sqrt(data(1,7).^2+data(1,10).^2 );
v_S(1,1)=v_R.*cosd(atand(data(1,10)./data(1,7))-data(1,13));
elseif data(1,CID_Pred_Alt)==2
v_R=sqrt(data(1,7).^2+data(1,10).^2 );
v_S(1,1)=v_R.*cosd(atand(data(1,10)./data(1,7))-(data(1,13)));
elseif data(1,CID_Pred_Alt)==3
v_R=sqrt(data(1,7).^2+data(1,10).^2 );
v_S(1,1)=v_R.*cosd(atand(data(1,10)./data(1,7))-(data(1,13)));
end
F_idm(1,6)=v_S(1,1); %initial speed of subject vehicle
F_idm(1,7)=sqrt(data(1,6).^2+data(1,9).^2); % initial resultant position of simulated follower
F_idm(1,8)=k;% Vehicle ID
 
 % Iteration perform based on MIDM
 spacing=[];
 rspeed=[];
 second_term=[];
 s_star=[];
 ds=[];
 S=[];
 dx=[];
 dy=[];
 
 THETA1=0; % Direction for centre alternative
 THETA2=5; % Direction for right alternative
 THETA3=-5; % Direction for left alternative
 LATDist=11.2; % Lateral distance
 LONGDist=245; % Longitudinal distance
 
 acc_x(1,1)=data(1,8);
 acc_y(1,1)=data(1,11);
 v_x(1,1)=data(1,7);
 v_y(1,1)=data(1,10);
 
 acc(1,1)=sqrt(data(1,8).^2+data(1,11).^2);
 
 for i=1:length-1
    
    dist(i+1,1)=sqrt(data(i+1,6).^2+data(i+1,9).^2);
   
    F_veh(i+1,6)=dist(i+1,1);
    
     if data(i,CID_Pred_Alt)==1
         spacing(i,1)=data(i,15); % spaceing
         
         rspeed(i,1)=data(i,16); % relative speed
         
         v=F_idm(i,6).*cosd(atand(F_idm(i,5)./F_idm(i,4))-THETA1); % correction
         
         second_term(i,1) =v*t0+v*rspeed(i,1)./(2*sqrt(a_max*b));
         s_star(i,1)=s0+max(0,second_term(i,1)); % desired minimum space
         
       
         extra_term = s1*sqrt(v/v0);
         s_star(i,1)=s_star(i,1)+extra_term; % modification is added here
         
         acc(i+1,1)=a_max*(1-(v/v0).^delta-(s_star(i,1)/spacing(i,1)).^2);  % acceleration of follower from IDM 
          acc_x(i+1,1)=acc(i+1,1).*cosd(THETA1); % Model case
          acc_y(i+1,1)=acc(i+1,1).*sind(THETA1); % Model case
         
         %Derive speed and distance of subject vehicle froem acceleration of IDM
        
         F_v(i+1,1)= max(0,abs(v+acc(i+1,1)*delT)); % Correction: speed of follower from IDM 
         
         
          v_x(i+1,1)=F_v(i+1,1).*cosd(THETA1); % Model case 
          v_y(i+1,1)=F_v(i+1,1).*sind(THETA1); % Model case 
         ds(i+1,1) = max(0, abs(v*delT+0.5*acc(i+1,1)*delT^2)); % distane of follower from IDM
          dx(i+1,1)=ds(i+1,1)*cosd(THETA1);    % Model case      
          dy(i+1,1)=ds(i+1,1)*sind(THETA1);    % Model case

         DX(i+1,1)=min(F_idm(i,2)+dx(i+1,1),LONGDist); % longitudinal distance for SV in MIDM
     
         DY(i+1,1)=min(F_idm(i,3)+dy(i+1,1),LATDist); %lateral distance for SV in MIDM
         S(1,1)=F_veh(1,2);
        
         S(i+1,1)=F_idm(i,7)+ds(i+1,1); % Total distanec of subject vehicle in  IDM model
         
         
         F_idm(i+1,2)=DX(i+1,1);
         F_idm(i+1,3)=DY(i+1,1);
         
          F_idm(i+1,4)=F_v(i+1,1)*cosd(THETA1); % Model case
          F_idm(i+1,5)=F_v(i+1,1)*sind(THETA1); % Model case
          
          F_idm(i+1,6)=F_v(i+1,1);
          F_idm(i+1,7)=S(i+1,1);
         F_idm(i+1,8)=k; % Vehicle ID
         
     elseif data(i,CID_Pred_Alt)==2
         
         spacing(i,1)=data(i,18);
        
         rspeed(i,1)=data(i,19); % relative speed
      
        v=F_idm(i,6).*cosd(atand(F_idm(i,5)./F_idm(i,4))-THETA2);
         second_term(i,1) =v*t0+v*rspeed(i,1)./(2*sqrt(a_max*b));
         s_star(i,1)=s0+max(0,second_term(i,1));
         
         extra_term = s1*sqrt(v/v0);
         s_star(i,1)=s_star(i,1)+extra_term; % modification is added here
         
         acc(i+1,1)=a_max*(1-(v/v0).^delta-(s_star(i,1)/spacing(i,1)).^2);  % acceleration of follower from IDM 
         
         acc_x(i+1,1)=acc(i+1,1).*cosd(THETA2);
         acc_y(i+1,1)=acc(i+1,1).*sind(THETA2);
   
         F_v(i+1,1)= max(0,abs(v+acc(i+1,1)*delT)); % Correction: speed of follower from IDM 
         v_x(i+1,1)=F_v(i+1,1).*cosd(THETA2);
         v_y(i+1,1)=F_v(i+1,1).*sind(THETA2);
        
         ds(i+1,1) = max(0, abs(v*delT+0.5*acc(i+1,1)*delT^2)); % distane of follower from IDM
         
         dx(i+1,1)=ds(i+1,1)*cosd(THETA2);
         dy(i+1,1)=ds(i+1,1)*sind(THETA2);

         S(1,1)=F_veh(1,2);
         S(i+1,1)=F_idm(i,7)+ds(i+1,1); % Total distanec of subject vehicle in  IDM model
         
        
         DX(i+1,1)=min(F_idm(i,2)+dx(i+1,1),LONGDist); % longitudinal distance for SV in MIDM
     
         DY(i+1,1)=min(F_idm(i,3)+dy(i+1,1),LATDist); %lateral distance for SV in MIDM
         
         
          F_idm(i+1,2)=DX(i+1,1);
          F_idm(i+1,3)=DY(i+1,1);
      
          F_idm(i+1,4)=F_v(i+1,1)*cosd(THETA2); % Model case
          F_idm(i+1,5)=F_v(i+1,1)*sind(THETA2); % Model case
         
          F_idm(i+1,6)=F_v(i+1,1);
          F_idm(i+1,7)=S(i+1,1); 
          F_idm(i+1,8)=k; % Vehicle ID
        
     elseif data(i,CID_Pred_Alt)==3
        
         spacing(i,1)=data(i,21);
      
         rspeed(i,1)=data(i,22); % relative speed
     
        v=F_idm(i,6).*cosd(atand(F_idm(i,5)./F_idm(i,4))-THETA3);
         second_term(i,1)=v*t0+v*rspeed(i,1)/(2*(a_max*b)^0.5);
         s_star(i,1)=s0+max(0,second_term(i,1));
       
         
         extra_term = s1*sqrt(v/v0);
         s_star(i,1)=s_star(i,1)+extra_term; % modification is added here
         
         acc(i+1,1)=a_max*(1-(v/v0).^delta-(s_star(i,1)/spacing(i,1)).^2);  % acceleration of follower from IDM 
         acc_x(i+1,1)=acc(i+1,1).*cosd(THETA3);
         acc_y(i+1,1)=acc(i+1,1).*sind(THETA3);
   
         F_v(i+1,1)= max(0,abs(v+acc(i+1,1)*delT)); % Correction: speed of follower from IDM

         v_x(i+1,1)=F_v(i+1,1).*cosd(THETA3);
         v_y(i+1,1)=F_v(i+1,1).*sind(THETA3);
        
         ds(i+1,1) = max(0, abs(v*delT+0.5*acc(i+1,1)*delT^2)); % distane of follower from IDM
           dx(i+1,1)=ds(i+1,1)*cosd(THETA3);
            dy(i+1,1)=ds(i+1,1)*sind(THETA3);

         S(1,1)=F_veh(1,2);
         
         S(i+1,1)=F_idm(i,7)+ds(i+1,1); % Total distanec of subject vehicle in  IDM model
         
         
         DX(i+1,1)=min(F_idm(i,2)+dx(i+1,1),LONGDist); % longitudinal distance for SV in MIDM
     
         DY(i+1,1)=min(F_idm(i,3)+dy(i+1,1),LATDist); %lateral distance for SV in MIDM
         
        
         F_idm(i+1,2)=DX(i+1,1);
         F_idm(i+1,3)=DY(i+1,1);
         
          F_idm(i+1,4)=F_v(i+1,1)*cosd(THETA3);
          F_idm(i+1,5)=F_v(i+1,1)*sind(THETA3);
         
         F_idm(i+1,6)=F_v(i+1,1);
        
         F_idm(i+1,7)=S(i+1,1);
         F_idm(i+1,8)=k; % Vehicle ID
     end
      
     % Objective function for calibration of model parameters
     
    Error(i+1,1)= (F_veh(i+1,2)-F_idm(i+1,2)).^2+(F_veh(i+1,3)-F_idm(i+1,3)).^2;
 end
 RMSE=sqrt(mean(Error));


plot(F_idm(:,1),F_idm(:,2),'Marker','.','MarkerSize',8,'Color','b','linestyle','-.')  % Simulated Trajectory
hold on

%  plot3(F_idm(:,2),F_idm(:,3), F_idm(:,1),'Marker','.','MarkerSize',8,'Color','b','linestyle','-.')  % Simulated Trajectory
%  hold on
% plot3(F_veh(:,2),F_veh(:,3), F_veh(:,1),'Marker','.','MarkerSize',8,'Color','r','linestyle','-.') % Observed Trajectory
 
 end
 F_idm_dataset=[F_idm_dataset;F_idm]; % Simulated dataset build-up for car 
  end
  
  
  
% MIDM for Motorcycle
  %SimulationDataM=xlsread('FiveAlternativeDataforMotorbikeToCalibrateMIDM.xlsx'); % Five alternative base Motorbike data for calibration
%   SimulationDataM=xlsread('ValidationDataForMotorcycles.xlsx');
  SimulationDataM=xlsread('/Applications/NikhilChandraSarkar/WholeSimulationDataForMotorcycle.xlsx');
  
 EparametersM=xlsread('/Applications/NikhilChandraSarkar/Eparameters_Motorbikes.xlsx'); %  Calibrated data of MIDM for Motorcycle
 VId_minM=min(SimulationDataM(:,1));
 VId_maxM=max(SimulationDataM(:,1));
 F_idm_datasetM=[];
 TT_M=[];
 TD_M=[];
  LATDist=11.2; % Lateral distance
  LONGDist=245; % Longitudinal distance
 
  for k=VId_minM:VId_maxM
 
 dataM=SimulationDataM(SimulationDataM(:,1)==k,:); % Extract specific car data from Simulation data set
 
 if isempty(dataM)
     DataEmpty=dataM;
 else
     dataM=dataM;
     k

 %Random selection of a raw from estimated parameters
    idxM = randperm(size(EparametersM,1),1); % Index of selected row
    B_M = EparametersM(idxM,:); % Randomly selected a row from a matrix of estimated parameters
    v0_M = (5/18)*B_M(1,3);
    a_max_M = B_M(1,4);
    b_M = B_M(1,5); % gamrnd(1.0,std(Eparameters(:,5)));
    t0_M=B_M(1,6); % Safety time headway (sec)
    s0_M=B_M(1,7); % Jam distance (m)
    s1_M=B_M(1,8); % Non-linear Jam distance (m)
    delta_M=B_M(1,9); % Acceleration exponent
 
 % Initial values of Modified IDM (MIDM)input:
 
 delT_M=0.5; % (sec)% Time increment
 tstart_M=dataM(1,3); % Starting time
 tend_M=dataM(end,3); % Ending time
 trange_M=dataM(:,3)>=tstart_M & dataM(:,3)<=tend_M; % Time duration
 length_M= size(dataM(trange_M,:),1);
 
 distM(1,1)=sqrt(dataM(1,6).^2+dataM(1,9).^2); % resultant position of subject vehicle
 F_vehM=[dataM(trange_M,3) dataM(trange_M,6) dataM(trange_M,9) dataM(trange_M,7) dataM(trange_M,10) NaN(length_M,1) ]; % time, longdist,latdist,longspeed, latspeed
 F_vehM(1,6)=sqrt(dataM(1,6).^2+dataM(1,9).^2); % resultant position for real subject vehicle
 
 F_idmM = [dataM(trange_M,3) NaN(length_M,7)];

 F_idmM(1,2)=F_vehM(1,2); % initial longitudinal position of simulated follower
 F_idmM(1,3)=F_vehM(1,3); % initial lateral position of simulated follower

F_idmM(1,4)=dataM(1,7);
F_idmM(1,5)=dataM(1,10);


% Initial resultant speed for subject vehicle

%Column Indext for Predicted Alternative
CID_Pred_AltM=35;
%Column Index For Observed Alternative
OID_Obs_AltM=14;

if dataM(1,CID_Pred_AltM)==1
v_R_M=sqrt(dataM(1,7).^2+dataM(1,10).^2 );
v_S_M(1,1)=v_R_M.*cosd(atand(dataM(1,10)./dataM(1,7))-dataM(1,13));
elseif dataM(1,CID_Pred_AltM)==2
v_R_M=sqrt(dataM(1,7).^2+dataM(1,10).^2 );
v_S_M(1,1)=v_R_M.*cosd(atand(dataM(1,10)./dataM(1,7))-(dataM(1,13)));
elseif dataM(1,CID_Pred_AltM)==3
v_R_M=sqrt(dataM(1,7).^2+dataM(1,10).^2 );
v_S_M(1,1)=v_R_M.*cosd(atand(dataM(1,10)./dataM(1,7))-(dataM(1,13)));


elseif dataM(1,CID_Pred_AltM)==4
v_R_M=sqrt(dataM(1,7).^2+dataM(1,10).^2 );
v_S_M(1,1)=v_R_M.*cosd(atand(dataM(1,10)./dataM(1,7))-(dataM(1,13)));

elseif dataM(1,CID_Pred_AltM)==5
v_R_M=sqrt(dataM(1,7).^2+dataM(1,10).^2 );
v_S_M(1,1)=v_R_M.*cosd(atand(dataM(1,10)./dataM(1,7))-(dataM(1,13)));

end
F_idmM(1,6)=v_S_M(1,1); %initial speed of subject vehicle
F_idmM(1,7)=sqrt(dataM(1,6).^2+dataM(1,9).^2); % initial resultant position of simulated follower
F_idmM(1,8)=k;% Vehicle ID

 % Iteration perform based on MIDM
 spacingM=[];
 rspeedM=[];
 second_termM=[];
 s_starM=[];
 dsM=[];
 S_M=[];
 dxM=[];
 dyM=[];
 
 
 THETA1_M=-4; % Direction for alt1
 THETA2_M=-2; % Direction for alt2
 THETA3_M=0; % Direction for alt3
 THETA4_M=2; % Directionfor alt4
 THETA5_M=4; % Direction for alt5
 
 acc_xM(1,1)=dataM(1,8);
 acc_yM(1,1)=dataM(1,11);
 v_xM(1,1)=dataM(1,7);
 v_yM(1,1)=dataM(1,10);
 
 accM(1,1)=sqrt(dataM(1,8).^2+dataM(1,11).^2);
 
 for i=1:length_M-1
    
    distM(i+1,1)=sqrt(dataM(i+1,6).^2+dataM(i+1,9).^2);
   
    F_vehM(i+1,6)=distM(i+1,1);
    
     if dataM(i,CID_Pred_AltM)==3
         spacingM(i,1)=dataM(i,21); % spaceing
         
         rspeedM(i,1)=dataM(i,22); % relative speed
        
         v_M=F_idmM(i,6).*cosd(atand(F_idmM(i,5)./F_idmM(i,4))-0); % correction
         second_termM(i,1) =v_M*t0_M+v_M*rspeedM(i,1)./(2*sqrt(a_max_M*b_M));
         s_starM(i,1)=s0_M+max(0,second_termM(i,1)); % desired minimum space
         
         
         extra_term = s1_M*sqrt(v_M/v0_M);
        extra_term = 0;
         s_starM(i,1)=s_starM(i,1)+extra_term; % modification is added here
         
         accM(i+1,1)=a_max_M*(1-(v_M/v0_M).^delta_M-(s_starM(i,1)/spacingM(i,1)).^2);  % acceleration of follower from IDM 
         
          acc_xM(i+1,1)=accM(i+1,1).*cosd(0); % Model case
          acc_yM(i+1,1)=accM(i+1,1).*sind(0); % Model case
         
        % Derive speed and distance of subject vehicle froem acceleration of IDM
        
         F_v(i+1,1)= max(0,abs(v_M+accM(i+1,1)*delT_M)); % speed of follower from IDM 
         
          v_xM(i+1,1)=F_v(i+1,1).*cosd(0); % Model case 
          v_yM(i+1,1)=F_v(i+1,1).*sind(0); % Model case
         dsM(i+1,1) = max(0, abs(v_M*delT_M+0.5*accM(i+1,1)*delT_M^2)); % distane of follower from ID
          dxM(i+1,1)=dsM(i+1,1)*cosd(0);    % Model case      
          dyM(i+1,1)=dsM(i+1,1)*sind(0);    % Model case

         DX(i+1,1)=min(F_idmM(i,2)+dxM(i+1,1),LONGDist); % longitudinal distance for SV in MIDM
     
         DY(i+1,1)=min(F_idmM(i,3)+dyM(i+1,1),LATDist); %lateral distance for SV in MIDM
         
         S_M(1,1)=F_vehM(1,2);
        
         S_M(i+1,1)=F_idmM(i,7)+dsM(i+1,1); % Total distanec of subject vehicle in  IDM model
         
         
         F_idmM(i+1,2)=DX(i+1,1);
         F_idmM(i+1,3)=DY(i+1,1);
         
          F_idmM(i+1,4)=F_v(i+1,1)*cosd(0); % Model case
          F_idmM(i+1,5)=F_v(i+1,1)*sind(0); % Model case
          
          F_idmM(i+1,6)=F_v(i+1,1);
          F_idmM(i+1,7)=S_M(i+1,1);  
          F_idmM(i+1,8)=k; % Vehicle ID
     elseif dataM(i,CID_Pred_AltM)==2
         
         spacingM(i,1)=dataM(i,18);
        
         rspeedM(i,1)=dataM(i,19); % relative speed
      
        
         v_M=F_idmM(i,6).*cosd(atand(F_idmM(i,5)./F_idmM(i,4))-THETA2_M);
         second_termM(i,1) =v_M*t0_M+v_M*rspeedM(i,1)./(2*sqrt(a_max_M*b_M));
         s_starM(i,1)=s0_M+max(0,second_termM(i,1));
         
         
         extra_term = s1_M*sqrt(v_M/v0_M);
       
         s_starM(i,1)=s_starM(i,1)+extra_term; % modification is added here
         
         accM(i+1,1)=a_max_M*(1-(v_M/v0_M).^delta_M-(s_starM(i,1)/spacingM(i,1)).^2);  % acceleration of follower from IDM 
         
         acc_xM(i+1,1)=accM(i+1,1).*cosd(THETA2_M);
         acc_yM(i+1,1)=accM(i+1,1).*sind(THETA2_M);
   
         F_v(i+1,1)= max(0,abs(v_M+accM(i+1,1)*delT_M)); % speed of follower from IDM 
         v_xM(i+1,1)=F_v(i+1,1).*cosd(THETA2_M);
         v_yM(i+1,1)=F_v(i+1,1).*sind(THETA2_M);
        
         dsM(i+1,1) = max(0, abs(v_M*delT_M+0.5*accM(i+1,1)*delT_M^2)); % distane of follower from IDM
         
         dxM(i+1,1)=dsM(i+1,1)*cosd(THETA2_M);
         dyM(i+1,1)=dsM(i+1,1)*sind(THETA2_M);

         S_M(1,1)=F_vehM(1,2);
         S_M(i+1,1)=F_idmM(i,7)+dsM(i+1,1); % Total distanec of subject vehicle in  IDM model
         
        
         DX(i+1,1)=min(F_idmM(i,2)+dxM(i+1,1),LONGDist); % longitudinal distance for SV in MIDM
     
         DY(i+1,1)=min(F_idmM(i,3)+dyM(i+1,1),LATDist); %lateral distance for SV in MIDM
         
         
         
          F_idmM(i+1,2)=DX(i+1,1);
          F_idmM(i+1,3)=DY(i+1,1);
          
          F_idmM(i+1,4)=F_v(i+1,1)*cosd(THETA2_M); % Model case
          F_idmM(i+1,5)=F_v(i+1,1)*sind(THETA2_M); % Model case
         
          F_idmM(i+1,6)=F_v(i+1,1);
          F_idmM(i+1,7)=S_M(i+1,1); 
          F_idmM(i+1,8)=k; % Vehicle ID
     elseif dataM(i,CID_Pred_AltM)==1
        
         spacingM(i,1)=dataM(i,15);
      
         rspeedM(i,1)=dataM(i,16); % relative speed
         v_M=F_idmM(i,6).*cosd(atand(F_idmM(i,5)./F_idmM(i,4))-THETA1_M);
         second_termM(i,1)=v_M*t0_M+v_M*rspeedM(i,1)/(2*(a_max_M*b_M)^0.5);
         s_starM(i,1)=s0_M+max(0,second_termM(i,1));
         
         
         extra_term = s1_M*sqrt(v_M/v0_M);
         extra_term = 0;
         s_starM(i,1)=s_starM(i,1)+extra_term; % modification is added here
         
         accM(i+1,1)=a_max_M*(1-(v_M/v0_M).^delta_M-(s_starM(i,1)/spacingM(i,1)).^2);  % acceleration of follower from IDM 
         acc_xM(i+1,1)=accM(i+1,1).*cosd(THETA1_M);
         acc_yM(i+1,1)=accM(i+1,1).*sind(THETA1_M);
   
         F_v(i+1,1)= max(0,abs(v_M+accM(i+1,1)*delT_M)); % speed of follower from IDM

         v_xM(i+1,1)=F_v(i+1,1).*cosd(THETA1_M);
         v_yM(i+1,1)=F_v(i+1,1).*sind(THETA1_M);
        
         dsM(i+1,1) = max(0, abs(v_M*delT_M+0.5*accM(i+1,1)*delT_M^2)); % distane of follower from IDM
         
            dxM(i+1,1)=dsM(i+1,1)*cosd(THETA1_M);
            dyM(i+1,1)=dsM(i+1,1)*sind(THETA1_M);

         S_M(1,1)=F_vehM(1,2);
         
         S_M(i+1,1)=F_idmM(i,7)+dsM(i+1,1); % Total distanec of subject vehicle in  IDM model
         
        
         DX(i+1,1)=min(F_idmM(i,2)+dxM(i+1,1),LONGDist); % longitudinal distance for SV in MIDM
     
         DY(i+1,1)=min(F_idmM(i,3)+dyM(i+1,1),LATDist); %lateral distance for SV in MIDM
         
         
         F_idmM(i+1,2)=DX(i+1,1);
         F_idmM(i+1,3)=DY(i+1,1);
         
          F_idmM(i+1,4)=F_v(i+1,1)*cosd(THETA1_M);
          F_idmM(i+1,5)=F_v(i+1,1)*sind(THETA1_M);
         
         F_idmM(i+1,6)=F_v(i+1,1);
        
         F_idmM(i+1,7)=S_M(i+1,1);
         F_idmM(i+1,8)=k; % Vehicle ID
         
         elseif dataM(i,CID_Pred_AltM)==4
        
         spacingM(i,1)=dataM(i,24);
      
         rspeedM(i,1)=dataM(i,25); % relative speed
     
         v_M=F_idmM(i,6).*cosd(atand(F_idmM(i,5)./F_idmM(i,4))-THETA4_M);
         second_termM(i,1)=v_M*t0_M+v_M*rspeedM(i,1)/(2*(a_max_M*b_M)^0.5);
         s_starM(i,1)=s0_M+max(0,second_termM(i,1));
        
         
         extra_term = s1_M*sqrt(v_M/v0_M);
         extra_term = 0;
         s_starM(i,1)=s_starM(i,1)+extra_term; % modification is added here
         
         accM(i+1,1)=a_max_M*(1-(v_M/v0_M).^delta_M-(s_starM(i,1)/spacingM(i,1)).^2);  % acceleration of follower from IDM 
         acc_xM(i+1,1)=accM(i+1,1).*cosd(THETA4_M);
         acc_yM(i+1,1)=accM(i+1,1).*sind(THETA4_M);
         
   
         F_v(i+1,1)= max(0,abs(v_M+accM(i+1,1)*delT_M)); % speed of follower from IDM

         v_xM(i+1,1)=F_v(i+1,1).*cosd(THETA4_M);
         v_yM(i+1,1)=F_v(i+1,1).*sind(THETA4_M);
        
         dsM(i+1,1) = max(0, abs(v_M*delT_M+0.5*accM(i+1,1)*delT_M^2)); % distane of follower from IDM
        
            dxM(i+1,1)=dsM(i+1,1)*cosd(THETA4_M);
             dyM(i+1,1)=dsM(i+1,1)*sind(THETA4_M);

         S_M(1,1)=F_vehM(1,2);
         
         S_M(i+1,1)=F_idmM(i,7)+dsM(i+1,1); % Total distanec of subject vehicle in  IDM model
         
        
         DX(i+1,1)=min(F_idmM(i,2)+dxM(i+1,1),LONGDist); % longitudinal distance for SV in MIDM
     
         DY(i+1,1)=min(F_idmM(i,3)+dyM(i+1,1),LATDist); %lateral distance for SV in MIDM
         
         
         F_idmM(i+1,2)=DX(i+1,1);
         F_idmM(i+1,3)=DY(i+1,1);
         
          F_idmM(i+1,4)=F_v(i+1,1)*cosd(THETA4_M);
          F_idmM(i+1,5)=F_v(i+1,1)*sind(THETA4_M);
         
         F_idmM(i+1,6)=F_v(i+1,1);
        
         F_idmM(i+1,7)=S_M(i+1,1);
         F_idmM(i+1,8)=k; % Vehicle ID
         
         elseif dataM(i,CID_Pred_AltM)==5
        
         spacingM(i,1)=dataM(i,27);
      
         rspeedM(i,1)=dataM(i,28); % relative speed
     
        
         v_M=F_idmM(i,6).*cosd(atand(F_idmM(i,5)./F_idmM(i,4))-THETA5_M);
         second_termM(i,1)=v_M*t0_M+v_M*rspeedM(i,1)/(2*(a_max_M*b_M)^0.5);
         s_starM(i,1)=s0_M+max(0,second_termM(i,1));
         
         
         extra_term = s1_M*sqrt(v_M/v0_M);
         
         s_starM(i,1)=s_starM(i,1)+extra_term; % modification is added here
         
         accM(i+1,1)=a_max_M*(1-(v_M/v0_M).^delta_M-(s_starM(i,1)/spacingM(i,1)).^2);  % acceleration of follower from IDM 
         acc_xM(i+1,1)=accM(i+1,1).*cosd(THETA5_M);
         acc_yM(i+1,1)=accM(i+1,1).*sind(THETA5_M);
   
         F_v(i+1,1)= max(0,abs(v_M+accM(i+1,1)*delT_M)); % speed of follower from IDM

         v_xM(i+1,1)=F_v(i+1,1).*cosd(THETA5_M);
         v_yM(i+1,1)=F_v(i+1,1).*sind(THETA5_M);
        
         dsM(i+1,1) = max(0, abs(v_M*delT_M+0.5*accM(i+1,1)*delT_M^2)); % distane of follower from IDM
         
        
             dxM(i+1,1)=dsM(i+1,1)*cosd(THETA5_M);
             dyM(i+1,1)=dsM(i+1,1)*sind(THETA5_M);

         S_M(1,1)=F_vehM(1,2);
         
         S_M(i+1,1)=F_idmM(i,7)+dsM(i+1,1); % Total distanec of subject vehicle in  IDM model
         
        
         DX(i+1,1)=min(F_idmM(i,2)+dxM(i+1,1),LONGDist); % longitudinal distance for SV in MIDM
     
         DY(i+1,1)=min(F_idmM(i,3)+dyM(i+1,1),LATDist); %lateral distance for SV in MIDM
         
         
         F_idmM(i+1,2)=DX(i+1,1);
         F_idmM(i+1,3)=DY(i+1,1);
         
          F_idmM(i+1,4)=F_v(i+1,1)*cosd(THETA5_M);
          F_idmM(i+1,5)=F_v(i+1,1)*sind(THETA5_M);
         
         F_idmM(i+1,6)=F_v(i+1,1);
        
         F_idmM(i+1,7)=S_M(i+1,1);
         F_idmM(i+1,8)=k; % Vehicle ID  
     end
      
     % Objective function for calibration of model parameters
     
    ErrorM(i+1,1)= (F_vehM(i+1,2)-F_idmM(i+1,2)).^2+(F_vehM(i+1,3)-F_idmM(i+1,3)).^2;
     
 end
 RMSE_M=sqrt(mean(ErrorM))
 
 
 plot(F_idmM(:,1),F_idmM(:,2),'Marker','.','MarkerSize',8,'Color','b','linestyle','-.')  % Simulated Trajectory
 hold on
 
%  plot3(F_idmM(:,2),F_idmM(:,3), F_idmM(:,1),'Marker','.','MarkerSize',8,'Color','b','linestyle','-.')  % Simulated Trajectory
%  hold on
%  plot3(F_vehM(:,2),F_vehM(:,3), F_vehM(:,1),'Marker','.','MarkerSize',8,'Color','r','linestyle','-.') % Observed Trajectory
 
% plot(F_vehM(:,1),F_vehM(:,2),'Marker','.','MarkerSize',8,'Color','r','linestyle','-.')  % Simulated Trajectory

 end
F_idm_datasetM =[F_idm_datasetM;F_idmM]; % Simulated dataset build-up
  end
  
  
SimulatedData=[];
SimulatedData=[SimulatedData;F_idm_dataset];
SimulatedData=[SimulatedData;F_idm_datasetM];


 hold off
  xlabel('Time (second)');
  ylabel('Simulated Longitudinal Distance (meter)');
  
%   ylabel('Lateral Distance (meter)');
%   xlabel('Longitudinal Distance (meter)');
%   ylabel('Lateral Distance (meter)');
%   zlabel('Time (second)');
%   legend('Simulated','Observed')
  
  
  %% Macroscopic properties for simulated dataset
  
 TableUpdate=[];
 To=min(SimulatedData(:,1));
 dt=60;
 Tend=max((SimulatedData(:,1)))-dt;
 d1=50;
 d2=195;
 
 for t=To:dt:Tend
 t1=t;
 t2=t+dt;
 dt=t2-t1;
 Time =t2;
 area=(t2-t1)*(d2-d1);
 data2region=SimulatedData(SimulatedData(:,1)>= t1 & SimulatedData(:,1) <= t2 & SimulatedData(:,2)>= d1 & SimulatedData(:,2)<=d2,:);
 VehID = data2region(:,8);
 UniqueVehID = unique(data2region(:,8)); % Unique vehicle ID
 N =numel(UniqueVehID); 
 TT=[];
 TD=[];

 for i=min(VehID):max(VehID)
 data2eachVeh=data2region(data2region(:,8)==i,:); % Individual vehicle data
 TimeEnd = max(data2eachVeh(:,1));
 TimeStart = min(data2eachVeh(:,1));
 TravelTime = TimeEnd -TimeStart;
 DistEnd = max(data2eachVeh(:,2));
 DistStart = min(data2eachVeh(:,2));
 TravelDist = DistEnd-DistStart;
 TT=[TT;TravelTime]; % Travel time for each vehicle
 TD=[TD;TravelDist]; % Travel distance for each vehicle
 end
 TTT = sum(TT);
 TDT= sum(TD);
 density = TTT/area; %(veh/m)
 Density=(density*1000); %(Veh/km)
 flow = TDT/area;
 Flow =flow*3600; % (Veh/h)
 Table =[];
 Table=[Time TTT TDT N Density Flow]; % time, total time taken, total distance taken, No. of vehicles, Density, Flow
 TableUpdate=[TableUpdate;Table];
 realData=real(TableUpdate);
 end
Speed =TableUpdate(:,6)./TableUpdate(:,5);

% %%Plotting macroscopic properties
% 
% figure
% plot(TableUpdate(:,5),TableUpdate(:,6),'.','MarkerSize', 20) % 2D plot
% hold on
% plot(TableUpdate_E(:,5),TableUpdate_E(:,6),'.','MarkerSize', 20) % 2D plot
% hold off
% axis([0 300 0 10000])
% xlabel('Density (veh/km)')
% ylabel('Flow (veh/h)')
% legend('Simulated','Observed')
% 
% figure
% plot(TableUpdate(:,5),Speed,'.','MarkerSize', 20)
% hold on
% plot(TableUpdate_E(:,5),Speed_E,'.','MarkerSize', 20)
% hold off
% axis([0 300 0 100])
% xlabel('Density (veh/km)')
% ylabel('Speed (km/h)')
% legend('Simulated','Observed')
% 
% figure
% plot(TableUpdate(:,6),Speed,'.','MarkerSize', 20)
% hold on
% plot(TableUpdate_E(:,6),Speed_E,'.','MarkerSize', 20)
% axis([0 10000 0 100])
% xlabel('Flow (veh/h)')
% ylabel('Speed (km/h)')
% legend('Simulated','Observed')
% 
% % figure
% % plot(TableUpdate(:,2),TableUpdate(:,3),'o')
% % xlabel('Total time taken (TTT) (s)')
% % ylabel('Total distance travel (TDT) (m)')
% % 
% % figure
% %  plot(TableUpdate(:,1),TableUpdate(:,5),'-r.','MarkerSize', 12)% time vs density graph
% %   xlabel('Time (s)')
% %   ylabel('Density (veh/km)')
% %   figure
% %   plot(TableUpdate(:,1),TableUpdate(:,6),'-r.','MarkerSize', 12)% time vs flow graph
% %   xlabel('Time (s)')
% %   ylabel('Flow (veh/h)')
% %   figure
% %   plot(TableUpdate(:,1),Speed,'-r.','MarkerSize', 12)% time vs speed graph
% %   xlabel('Time (s)')
% %   ylabel('Speed (km/h)')
  
  

  