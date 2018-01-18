clean
load('maxEff.mat');

i_max = max(size(x));
while isnan(x(i_max))    
       i_max = i_max -1;
end


for i = 1:i_max
   
    Ehs(i)   = fc{i}{end}(1)/fc{i}{end}(2);
    Ehs_c(i) = fcopt{i}{end}(1)/fcopt{i}{end}(2);
    Exf(i)   = ff{i}{end}(1)/ff{i}{end}(2);
    
end
% 
% figure(1000);
% plot(Ehs); hold on; plot(Ehs_c,'r'); plot(Exf,'k')

for i = 2:i_max
    
    figure(1000+i-1)
    title(sprintf('%d VS %d',i-1,i));
    
    % caso i-1
    %alpha_double_check = 2*linspace(floor(x(i-1)/2),ceil(x(i-1)/2),6);
    hold on        
    plot(alpha_test,surPol.Eff{i-1},'bx--',...
         x(i-1),Ehs_c(i-1),'bo');
         %        alpha_double_check,double_check.Eff{i-1},'bs--',...
    hold on
    
    % caso i 
    %alpha_double_check = 2*linspace(floor(x(i)/2),ceil(x(i)/2),6);
    plot(alpha_test,surPol.Eff{i},'rx--',...
        x(i),Ehs_c(i),'ro')
    %        alpha_double_check,double_check.Eff{i},'rs--',...
    grid on
    
    legend(sprintf('Polare iterazione %d',i-1),...
           sprintf('Massimo trovato: [%1.2f %1.2f]',x(i-1),Ehs_c(i-1)),...
           sprintf('Polare iterazione %d',i),...
           sprintf('Massimo trovato: [%1.2f %1.2f]',x(i),Ehs_c(i)));
       
    plot(x_ms(i-1,:),-f_ms(i-1,:),'bs');
    plot(x_ms(i,:),  -f_ms(i,:)  ,'rs');
    
end



% %% pp
% clean
% load('log.mat');
% color_code = {'g','c','m','k','w'};
% % quante iterazioni?
% itc = size(x,2);

% 
% figure(10);
% plot([0:30],-14*ones(size([0:30])),'r'); hold on;
% plot([0:30],pol.DeltaP,'b','LineWidth',2); grid on;
% 
% leg_cell_10{1} = 'Valarezo Limit';
% leg_cell_10{2} = sprintf('Xfoil (intersection @ %1.2f deg)',16.85);
% 
% 
% Dmat = nan(itc,max(size(alpha_test)));
% 
% hold on
% for i = 1:itc
%     for j = 1:max(size(alpha_test))
%        Dmat(i,j) = polCor{i,j}{1};
%     end
%  
%     plot(alpha_test,Dmat(i,:),strcat(color_code{i},'x--'));
%     leg_cell_10{end+1} = sprintf('Iter %d',i);
%     plot(x(i),fcopt{i}{1},strcat(color_code{i},'s'),'LineWidth',2);
%     leg_cell_10{end+1} = sprintf('Iter %d result -> [%1.2f %1.2f]',i,x(i),fcopt{i}{1});
% end
% 
% legend(leg_cell_10,'Location','bestoutside')