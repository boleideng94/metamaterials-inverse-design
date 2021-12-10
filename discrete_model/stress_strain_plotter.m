clear all

file ='softest1';
load(['data/raw_data_',file,'.mat']);

factor = 1/(0.02*0.1)/1e3;
kl = 470;

%%
% mycolor = [102,204,0]./255;
% exp_stress = dlmread(['treated_exp_data/',file,'.txt']);
% target_stress = dlmread(['target_curves/',file,'.txt']);
% target_stress = target_stress.*1;
% plot(-strain,stress.*kl*factor,'-k','linewidth',2); hold on
% plot(-strain,exp_stress','--','linewidth',4,'color',mycolor);

%%
close all
clf
fig = figure(1);
set(fig,'position', [0, 0, 500,400]);
set(fig,'color','w');
framerate = 20;
mkdir('animation');
filename = ['animation/',file,'_stress_strain_num'];
myVideo = VideoWriter(filename,'MPEG-4');
myVideo.FrameRate = framerate;
open(myVideo);
mycolor = [102,204,0]./255;

for k = 1:100
    hold off
%     plot(-strain,-target_stress.*kl*factor,'--r','linewidth',2); hold on
%     plot(-strain(k),-exp_stress(k),'-ko','linewidth',0.5,'MarkerFaceColor',mycolor,'MarkerSize',10); hold on
%     plot(-strain(1:k),-exp_stress(1:k)','--','linewidth',3,'color',mycolor);
    plot(-strain(k),-stress(k)*kl*factor,'-ko','linewidth',0.5,'MarkerFaceColor','k','MarkerSize',10);hold on
    plot(-strain(1:k),-stress(1:k).*kl*factor,'-k','linewidth',2);  hold off
    xlim([0,0.1]);
    ylim([-max(stress.*kl*factor).*1.2,0]);
%     ylim([0,20]);
%     yticks([0,0.05,0.1]);
    xticklabels({'0','-0.05','-0.1'})
    set(gca,'Ydir','reverse')
    
    set(gca,'fontsize',24);
    ylabel('Nominal stress \sigma [kPa]');
    xlabel('Applied strain \epsilon');
    xticks([0,0.05,0.1]);
    xticklabels({'0','-0.05','-0.1'})
    frame = getframe(gcf);
    writeVideo(myVideo, frame);
    pause(0.00001)
end
close(myVideo)
% 


%%

% %%
% % plot(-strain,stress,'-k');
% clf
% hold on
% fig = figure(1);
% set(fig,'position',[0,0,360,270]);
% factor = 1/(0.02*0.1)/1e3;
% kl = 470;
% 
% plot(-strain,stress.*kl*factor ,'-k','linewidth',2); hold on
% 
% % filename0 = 'random_case_2'; 
% % filename0 = 'stiff_long_1';
% 
% 
% filename0 = '0907/super_stiff_2';
% % filename_exp = '0921/test2jump_2';
% % filename_exp = '0929/nonmono_2';
% % filename_exp = 'random_case_1';
% % filename_exp = 'super_stiff_2';
% % filename_exp ='big_jump_2';
% % filename_exp = 'perfect_squares'; 
% 
% % experiemntal results
% M = readtable(['exp_data/',filename_exp,'.csv']);
% strain_exp = smooth(M.Displacement./89);
% force_exp = smooth(M.Force.*1000/1000,20);
% strain_shift = 0.013;
% strain_exp = strain_exp - strain_shift;
% 
% [~,id1] = max(strain_exp);
% [~,id2] = min(abs(strain_exp(id1:end)));
% id2 = id2+id1;
% strain1 = strain_exp(3:id1);
% strain2 = flip(strain_exp(id1+1:id2-1));
% force1 = force_exp(3:id1);
% force2 = flip(force_exp(id1+1:id2-1));
% strain_inter = linspace(0.0001,0.1,99);
% force_inter = (interp1(strain1,force1,strain_inter') +  interp1(strain2,force2,strain_inter'))./2;
% % 
% force_exp = smooth(force_exp,20);
% exp_stress = force1.*factor*1e3*1.1;
% plot(strain1,exp_stress,'-r');
% strain_take = round(interp1(strain1,1:length(strain1),-strain));
% stress_exp = exp_stress(strain_take);
% plot(-strain,stress_exp,'-r');
% 
% xlim([0,0.10]);
% % ylim([0,20]);
% % ylim([0,2])
% 
% set(gca,'fontsize',24);
% xlabel('strain');
% 
% folder = 'treated_exp_data';
% mkdir(folder);
% dlmwrite([folder,'/',file,'.txt'],stress_exp);
% 
% %
% % plot(strain1,force1);
% % plot(strain2,force2);
% % plot(strain_exp,force_exp); hold onn/',file,'.png']);
