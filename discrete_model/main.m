clear all
global m n alpha Ktheta Ks
global damp_theta damp_u
global Time load_Time delta_y
global progress gamma hori_shift vert_shift
global strain stress count ang_b m0 w0 k_contact

rc_m = 2; rc_n = 2; % a two-by-two unit cell

%% User defined parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file = 'Fig2c_iv'; % input design name
p_m = 5; p_n = 4; % structure size (number of unit cells)
% -----------------------------------
alpha = 1.73; % normalized inertia
% hinge stiffness
kl = 470;
Ktheta = 0.0075; % normalized torsional stiffness
Ks = 0.33; % normalized shearing stiffness
gamma = -0.2;
k_contact = 0.05;
damp = 0.025; % damping coefficient
damp_theta = damp; damp_u = damp;
% loading conditions
comp_per = 0.10; % maximal compression strain
FrameN  = 100; % output frame number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Derived parameters
n = rc_n*p_n+2; m = rc_m*p_m; % total size of the system (number of quadrilaterals)
load_Time = 500*n; % loading time
input = load(['input_design/',file,'.txt']);
rc_rand = reshape(input,[4,2,2]);
rc_rand(:,1:rc_n,rc_m+1) = rc_rand(:,1:rc_n,1);
rc_rand(:,rc_n+1,1:rc_m) = rc_rand(:,1,1:rc_m);

hori_shift =zeros([2,n,m+1]);
vert_shift = zeros([2,n+1,m]);

for p_i = 1:p_n
    for p_j = 1:p_m
        vert_shift(:,(p_i-1)*rc_n+2:(p_i)*rc_n+2,(p_j-1)*rc_n+1:(p_j)*rc_n+1) = rc_rand(1:2,:,:);
        hori_shift(:,(p_i-1)*rc_n+2:(p_i)*rc_n+2,(p_j-1)*rc_n+1:(p_j)*rc_n+1) = rc_rand(3:4,:,:);
    end
end
hori_shift(:,n,:) = 0;

delta_y = -comp_per*(n-1);
Time = load_Time;
y0 = 0.00001.*(rand(m*n*6,1)-0.5);%Initial value for variables


% edge vectors
for i = 1:n
    for j = 1:m
        R0_right = [+0.5,0]' + hori_shift(:,i,j+1);
        R0_left = [-0.5,0]' + hori_shift(:,i,j);
        R0_up = [0,+0.5]' + vert_shift(:,i+1,j);
        R0_down = [0,-0.5]' + vert_shift(:,i,j);
        ev(:,i,j,1) = R0_up - R0_right;
        ev(:,i,j,2) = R0_left - R0_up;
        ev(:,i,j,3) = R0_down - R0_left;
        ev(:,i,j,4) = R0_right - R0_down;
    end
end

ang_b = NaN([n,m,4,2]);
for i = 1:n
    for j = 1:m
        if j<m % there is a right neighbour
            ang_b(i,j,1,1) = ang_bet(ev(:,i,j,4),-ev(:,i,j+1,3));
            ang_b(i,j,1,2) = ang_bet(ev(:,i,j,1),-ev(:,i,j+1,2));
        end
        if i<n % there is a up neighbour
            ang_b(i,j,2,1) = ang_bet(ev(:,i,j,1),-ev(:,i+1,j,4));
            ang_b(i,j,2,2) = ang_bet(ev(:,i,j,2),-ev(:,i+1,j,3));
        end
        if j>1 % there is a left neighbour
            ang_b(i,j,3,1) = ang_bet(ev(:,i,j,2),-ev(:,i,j-1,1));
            ang_b(i,j,3,2) = ang_bet(ev(:,i,j,3),-ev(:,i,j-1,4));
        end
        if i>1 % there is a bot neighbour
            ang_b(i,j,4,1) = ang_bet(ev(:,i,j,3),-ev(:,i-1,j,2));
            ang_b(i,j,4,2) = ang_bet(ev(:,i,j,4),-ev(:,i-1,j,1));
        end
    end
end

%% ODE45 solver
progress = 0; count = 0;
tic;
[T,Yout] = ode45('ODE',[0,Time], y0);                            %Test Mark
[Timelength, Distance]=size(Yout);
timeelapse=toc

%% Data processing
clear DispX DispY thetaR theta_R
for i = 1:n
    for j = 1:m
        posi = ((i-1)*m + j-1)*6;
        DispX(i,j,:) = Yout(:,posi+1);
        DispY(i,j,:) = Yout(:,posi+2);
        theta_R(i,j,:) = Yout(:,posi+3);
    end
end

% Original positions
for i=1:n
    for j=1:m
        X0(i,j) = j;
        Y0(i,j) = i;
    end
end

% find array of time point with total number of FrameN
timepoints = round(interp1(T,1:length(T),linspace(T(1),T(end-1),FrameN)));
for t=1:FrameN
    for i=1:n
        for j=1:m
            UX(i,j,t) = DispX(i,j,timepoints(t));%x diplacement for the center
            UY(i,j,t) = DispY(i,j,timepoints(t));%y displacment for the center
            thetaR(i,j,t) = theta_R(i,j,timepoints(t));%Rotational angle
        end
    end
end
% save data
mkdir('data');
save(['data/raw_data_',file,'.mat'],'UX','UY','thetaR','strain','stress','Ktheta','Ks','gamma','damp_u','damp_theta');

%% Video: Absolute rotation plot
l=0.5;
u=[0,0]';
count1 = 0;
for t=1:FrameN
    count1 = count1 + 1;
    for i=1:n
        for j=1:m
            R0_right = [+0.5,0]' + hori_shift(:,i,j+1);
            R0_left = [-0.5,0]' + hori_shift(:,i,j);
            R0_up = [0,+0.5]' + vert_shift(:,i+1,j);
            R0_down = [0,-0.5]' + vert_shift(:,i,j);
            
            Point1 = R0_left;
            Point2 = R0_up;
            Point3 = R0_right;
            Point4 = R0_down;
            thet = thetaR(i,j,t);
            RoMatrix = [cos(thet),sin(thet); -sin(thet),cos(thet)];
            u(1)=UX(i,j,t) + X0(i,j);
            u(2)=UY(i,j,t) + Y0(i,j);
            Point(i,j,t,1,:) = u + (RoMatrix*Point1);
            Point(i,j,t,2,:) = u + (RoMatrix*Point2);
            Point(i,j,t,3,:) = u + (RoMatrix*Point3);
            Point(i,j,t,4,:) = u + (RoMatrix*Point4);
        end
    end
end

close all
clf
fig = figure(1);
set(fig,'position', [0, 0, m*100,n*90])
framerate = 20;
mkdir('animation');
filename = ['animation/',file,'_ab_ang'];
myVideo = VideoWriter(filename,'MPEG-4');
myVideo.FrameRate = framerate;
open(myVideo);

timemark = 0;
count1 = 0; kk = 0;
base_color = [0.5,0.5,0.5];

for t=1:FrameN
    hold off
    count1 = count1 + 1;  
    % plot the base grid
    for j = 1:m
        plot(2*l.*[j,j],[0,2*(n+1)*l],'--k','linewidth',0.5);hold on
    end
    for i = 1:n
        plot([0,2*(m+1)*l],2*l.*[i,i],'--k','linewidth',0.5);hold on
    end
    
    for i=1:n
        for j=1:m 
            % set the vectors to plot squares
            X=[Point(i,j,t,1,1),Point(i,j,t,2,1),Point(i,j,t,3,1),...
                Point(i,j,t,4,1)];
            Y=[Point(i,j,t,1,2),Point(i,j,t,2,2),Point(i,j,t,3,2),...
                Point(i,j,t,4,2)];
            set(gcf,'Color',[1,1,1]);
            
            scale = 1/(pi/4);
            colorvalue = (abs(thetaR(i,j,t))*scale);            
            colorvalue = min(max(colorvalue,0),1);cm = jet;
            colorID = max(1,sum(colorvalue > [0:1/length(cm(:,1)):1])); myColor = cm(colorID, :);
            fill(X,Y,myColor,'FaceAlpha',0.9); hold on          
            axis([-1,m+2,0,n+1]); axis off
        end
    end
    
    colorvalue = abs(thetaR(1,1,t)*scale);colorvalue = min(max(colorvalue,0),1);cm = jet;
    colorID = max(1,sum(colorvalue > [0:1/length(cm(:,1)):1])); myColor = cm(colorID, :);
    
    X_rec = [l*1,2*l*m+l,2*l*m+l,1*l];
    Y_rec = [2*l,2*l,l,l];
    fill(X_rec,Y_rec,myColor,'FaceAlpha',1); 
    
    Y_rec = Y_rec + (n-1/2)*2*l + UY(n,1,t);
    fill(X_rec,Y_rec,myColor,'FaceAlpha',1);
    
    frame = getframe(gcf);
    writeVideo(myVideo, frame);
    pause(0.00001)
end
close(myVideo)


%% Video: Rotation induced shrinkage
% record initial height of the quadrilaterals
t = 1;
for i=1:n
    for j=1:m
        % set the vectors to plot squares
        X=[Point(i,j,t,1,1),Point(i,j,t,2,1),Point(i,j,t,3,1),...
            Point(i,j,t,4,1)];
        Y=[Point(i,j,t,1,2),Point(i,j,t,2,2),Point(i,j,t,3,2),...
            Point(i,j,t,4,2)];
        rs0(i,j) = Y(2) - Y(4);
    end
end

clf
fig = figure(1);
set(fig,'position', [0, 0, m*100,n*90])
framerate = 20;
mkdir('animation');
filename = ['animation/',file,'_RIS'];
myVideo = VideoWriter(filename,'MPEG-4');
myVideo.FrameRate = framerate;
open(myVideo);

timemark = 0;
count1 = 0; kk = 0;
base_color = [0.5,0.5,0.5];

for t=1:FrameN
    hold off
    count1 = count1 + 1;  
    % plot the base grid
    for j = 1:m
        plot(2*l.*[j,j],[0,2*(n+1)*l],'--k','linewidth',0.5);hold on
    end
    for i = 1:n
        plot([0,2*(m+1)*l],2*l.*[i,i],'--k','linewidth',0.5);hold on
    end
    
    for i=1:n
        for j=1:m 
            % set the vectors to plot squares
            X=[Point(i,j,t,1,1),Point(i,j,t,2,1),Point(i,j,t,3,1),...
                Point(i,j,t,4,1)];
            Y=[Point(i,j,t,1,2),Point(i,j,t,2,2),Point(i,j,t,3,2),...
                Point(i,j,t,4,2)];
            set(gcf,'Color',[1,1,1]);
            
            rs = (rs0(i,j) - (Y(2)-Y(4)))/1;
            colorvalue = (rs/0.2)^1;      
            colorvalue = min(max(colorvalue,0),1);cm = jet;
            colorID = max(1,sum(colorvalue > [0:1/length(cm(:,1)):1])); myColor = cm(colorID, :);
            fill(X,Y,myColor,'FaceAlpha',0.9); hold on          
            axis([-1,m+2,0,n+1]); axis off
        end
    end
    
    colorvalue = (0);colorvalue = min(max(colorvalue,0),1);cm = jet;
    colorID = max(1,sum(colorvalue > [0:1/length(cm(:,1)):1])); myColor = cm(colorID, :);
    
    X_rec = [l*1,2*l*m+l,2*l*m+l,1*l];
    Y_rec = [2*l,2*l,l,l];
    fill(X_rec,Y_rec,myColor); 
    
    Y_rec = Y_rec + (n-1/2)*2*l + UY(n,1,t);
    fill(X_rec,Y_rec,myColor);
    
    frame = getframe(gcf);
    writeVideo(myVideo, frame);
    pause(0.00001)
end
close(myVideo)

%% Stress-strain curve
fig = figure(2);
set(fig,'position', [0, 0, 500,400]);
set(fig,'color','w');

area = (0.02*0.1); % crosssection of the sample
plot(-strain,-stress.*kl/area*1e-3,'-k','linewidth',2);  hold off
xlim([0,0.1]);
ylim([-max(stress.*kl/area*1e-3).*1.2,0]);
%     ylim([0,20]);
%     yticks([0,0.05,0.1]);
xticklabels({'0','-0.05','-0.1'})
set(gca,'Ydir','reverse')
set(gca,'fontsize',24);
ylabel('Nominal stress \sigma [kPa]');
xlabel('Applied strain \epsilon');
xticks([0,0.05,0.1]);
xticklabels({'0','-0.05','-0.1'})
mkdir('figures');
saveas(gca,['figures/stress_strain_',file,'.png']);



