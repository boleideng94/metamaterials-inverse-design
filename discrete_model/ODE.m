function dy = ODE(T,y)

global m n alpha Ktheta Ks 
global damp_theta damp_u Time
global progress gamma delta_y load_Time
global hori_shift vert_shift
global strain stress count ang_b

flag = 0;
if T/Time > progress
    progress = progress + 0.01;
    disp([num2str(progress*100),'%']);
    flag = 1;
    count = count + 1; stress(count) = 0; strain(count) = 0;
end
% if T/Time<0.5
Veloy = delta_y/load_Time;
% else
%  Veloy = -delta_y/load_Time*2;   
% end

mass = 1;
J = 1/alpha^2*mass*0.5^2;

% Initial and boundary conditions
dy = zeros(m*n*6,1);
F = zeros(m*n*6,1);

% Main part of ODE
for i=1:n
    for j=1:m
        % i-th row and j-th column
        % Each cell have 3 DOFs: U V Th
        % Define coodinates
        posi = ((i-1)*m + j-1)*6;
        up_posi = ((i-1+1)*m + j-1)*6;
        down_posi = ((i-1-1)*m + j-1)*6;
        left_posi = ((i-1)*m + j-1-1)*6;
        right_posi = ((i-1)*m + j-1+1)*6;
        %Coordinates of the center mass
        U_cen=posi+1;
        V_cen=posi+2;
        Th_cen=posi+3;
        %Coordinates of the top mass
        U_up=up_posi+1;
        V_up=up_posi+2;
        Th_up=up_posi+3;
        %Coordinates of the bottom mass
        U_down=down_posi+1;
        V_down=down_posi+2;
        Th_down=down_posi+3;
        %Coordinates of the right mass
        U_right=right_posi+1;
        V_right=right_posi+2;
        Th_right=right_posi+3;
        %Coordinates of the left mass
        U_left=left_posi+1;
        V_left=left_posi+2;
        Th_left=left_posi+3;
        %% Direction vectors
        %Vectors from center of mass
        R0_right = [+0.5,0]' + hori_shift(:,i,j+1);
        R0_left = [-0.5,0]' + hori_shift(:,i,j);
        R0_up = [0,+0.5]' + vert_shift(:,i+1,j);
        R0_down = [0,-0.5]' + vert_shift(:,i,j);
        
        R0_right_neighbor = [-0.5,0]' + hori_shift(:,i,j+1);
        R0_left_neighbor = [+0.5,0]' + hori_shift(:,i,j);
        R0_up_neighbor = [0,-0.5]' + vert_shift(:,i+1,j);
        R0_down_neighbor = [0,+0.5]' + vert_shift(:,i,j);
        
        R_right = Rot_matrix(y(Th_cen))*R0_right;
        R_left = Rot_matrix(y(Th_cen))*R0_left;
        R_up = Rot_matrix(y(Th_cen))*R0_up;
        R_down = Rot_matrix(y(Th_cen))*R0_down;
   
%         R_right = 0.5.*[cos(y(Th_cen)), -sin(y(Th_cen))];
%         R_up = 0.5.*[sin(y(Th_cen)), cos(y(Th_cen))];
%         R_left = 0.5.*[-cos(y(Th_cen)), sin(y(Th_cen))];
%         R_down = 0.5.*[-sin(y(Th_cen)), -cos(y(Th_cen))];
        
        %Elongations (forces) of springs
        D_right=[0,0];D_up=[0,0];D_left=[0,0];D_down=[0,0];
        
        if j~=m %otherwise nothing on the right
%             D_right = [(y(U_right)-y(U_cen))+0.5*(2-cos(y(Th_cen))-cos(y(Th_right))),...
%                 (y(V_right)-y(V_cen))+0.5*(sin(y(Th_right))+sin(y(Th_cen)))];           
            D_right = [1,0]' - R_right + Rot_matrix(y(Th_right))*R0_right_neighbor ...
                + [y(U_right)-y(U_cen);y(V_right)-y(V_cen)];
        end
        if i~=n %otherwise nothing on the top
%             D_up = [(+y(U_up)-y(U_cen))-0.5*(sin(y(Th_up))+sin(y(Th_cen))),...
%                 (y(V_up)-y(V_cen))+0.5*(2 -cos(y(Th_cen))-cos(y(Th_up)))];
            D_up = [0,1]' - R_up + Rot_matrix(y(Th_up))*R0_up_neighbor ...
                + [y(U_up)-y(U_cen);y(V_up)-y(V_cen)];
        end
        if j~=1 %otherwise nothing on the left
%             D_left = [(y(U_left)-y(U_cen))-0.5*(2 -cos(y(Th_cen))-cos(y(Th_left))),...
%                 (y(V_left)-y(V_cen))-0.5*(sin(y(Th_left))+sin(y(Th_cen)))];
            D_left = [-1,0]' - R_left + Rot_matrix(y(Th_left))*R0_left_neighbor ...
                + [y(U_left)-y(U_cen);y(V_left)-y(V_cen)];
        end
        if i~=1 %otherwise nothing on the bottom
%             D_down = [(y(U_down)-y(U_cen))+0.5*(sin(y(Th_down))+sin(y(Th_cen))),...
%                 (y(V_down)-y(V_cen))-0.5*(2-cos(y(Th_cen))-cos(y(Th_down)))];
            D_down = [0,-1]' - R_down + Rot_matrix(y(Th_down))*R0_down_neighbor ...
                + [y(U_down)-y(U_cen);y(V_down)-y(V_cen)];
        end
        
%         if i==4 && j == 8 
%             y(Th_cen)
%             if y(Th_cen)<-0.6
%             'test'
%             end
%         end
   
%         if i == 2 && j ==2
%             'check';
%             
%         end
        % Forces and moments
        % Moment from rotational springs
        Mm=0;
        % Moment from springs via different directions

        M_right=0;M_up=0;M_left=0;M_down=0; M_cont = 0;
        if j~=m %otherwise nothing on the right
            M_right = Ks*R_right(1)*D_right(2) - R_right(2)*D_right(1);
            d_angle = y(Th_cen)-y(Th_right);
            Mm = Mm + Ktheta/4*(d_angle + gamma*sign(d_angle)*d_angle^2);%gamma*d_angle^3
            M_cont = M_cont + contact_moment(y(Th_cen),y(Th_right),reshape(ang_b(i,j,1,:),[2,1]));
        end
        if i~=n %otherwise nothing on the top
            M_up = R_up(1)*D_up(2) - Ks*R_up(2)*D_up(1);
            d_angle = y(Th_cen)-y(Th_up);
            Mm = Mm + Ktheta/4*(d_angle +  gamma*sign(d_angle)*d_angle^2);
            M_cont = M_cont + contact_moment(y(Th_cen),y(Th_up),reshape(ang_b(i,j,2,:),[2,1]));
        end
        if j~=1 %otherwise nothing on the left
            M_left = Ks*R_left(1)*D_left(2) - R_left(2)*D_left(1);
            d_angle = y(Th_cen)-y(Th_left);
            Mm = Mm + Ktheta/4*(d_angle +  gamma*sign(d_angle)*d_angle^2);
            M_cont = M_cont +  contact_moment(y(Th_cen),y(Th_left),reshape(ang_b(i,j,3,:),[2,1]));
        end
        if i~=1 %otherwise nothing on the bottom
            M_down = R_down(1)*D_down(2) - Ks*R_down(2)*D_down(1);
            d_angle = y(Th_cen)-y(Th_down);
            Mm = Mm + Ktheta/4*(d_angle +  gamma*sign(d_angle)*d_angle^2);
            M_cont = M_cont +  contact_moment(y(Th_cen),y(Th_down),reshape(ang_b(i,j,4,:),[2,1]));
        end
       
        
        F(U_cen) = D_right(1) + Ks*D_up(1) + Ks*D_down(1) + D_left(1);
        F(V_cen) = Ks*D_right(2) + D_up(2) + D_down(2) + Ks*D_left(2);
        F(Th_cen) = +M_right+M_left+M_up+M_down+Mm + M_cont;
        
        dy(U_cen) = y(U_cen+3);
        dy(V_cen) = y(V_cen+3);
        dy(Th_cen) = y(Th_cen+3);
        dy(U_cen+3) = F(U_cen)/mass-damp_u*y(U_cen+3);
        dy(V_cen+3) = F(V_cen)/mass-damp_u*y(V_cen+3);
        dy(Th_cen+3) = -F(Th_cen)/J-damp_theta*y(Th_cen+3);
        
        % loadign boundary
        if i == n % the top row units
            %                 y(V_cen) = Dispy;
            dy(V_cen) = Veloy;
            dy(U_cen) = 0;
            dy(Th_cen) = 0;
            
            if flag == 1
                stress(count) = stress(count) + D_down(2)/m;
                strain(count) = strain(count) + y(V_cen)/m/(n-1);
            end
        end
        if i == 1 % the bottom row units
            %                 y(V_cen) = 0;
            dy(V_cen) = 0;
            dy(U_cen) = 0;
            dy(Th_cen) = 0;
        end
    end
end
end

function M_cont = contact_moment(theta1,theta2,th0)
% global m0 w0
global k_contact
% w0 = 0.03; m0 = 0.5;
if isnan(th0(1)) || isnan(th0(2))
    M_cont = 0;
else
%     M1 = (tanh((-theta1 + theta2-th0(1))/w0-w0)+1)*(-m0/2);
%     M2 = (tanh((-theta2 + theta1-th0(2))/w0-w0)+1)*(m0/2);
%     
    M1 = -k_contact*max(-theta1 + theta2-th0(1)+0.1,0);
    M2 = k_contact*max(-theta2 + theta1-th0(2)+0.1,0);
    
    M_cont = M1 + M2;
end
end

function Rot = Rot_matrix(theta)
    Rot = [cos(theta),sin(theta);
        -sin(theta),cos(theta)];
end




