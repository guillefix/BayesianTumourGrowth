x0 = ones(1,8);
x=x0;
V0=[0.031;0.015];
tspan=[0,15];
% sampled_times=1:0.5:15;
% sampled_times_individual_data=[10,15];
sampled_times_individual_data=[10,15];
sampled_times_sum_data=[3 5 7 10 13 15];
standard_deviation_noise=0.05;

odds_of_initial_values = [1 0;
        3 1;
        1 1;
        1 3;
        0 1];

% neg log likelihood
V_sum_data_average_PC3_0Gy = [0.0636679167 0.0661274167 0.0987616667 0.16010675 0.25829375 0.280966;
    0.0574364167 0.0799836667 0.1146740833 0.2577545833 0.41375175 0.4119126667;
    0.0589646667 0.0776658333 0.1201630833 0.27822025 0.5145295 0.557941;
    0.0544568333 0.08331 0.1326801667 0.35862575 0.5359255 0.523286;
    0.0503110833 0.0921101667 0.1695063333 0.4100345833	0.5936558333 0.6179516667];

V_individual_control_data_PC3_0Gy = [0.99 0.98;
    0.56 0.45;
    0.32 0.18;
    0.12 0.06;
    0.00 0.00];

V_individual_resistant_data_PC3_0Gy = [0.00 0.00;
    0.43 0.53;
    0.67 0.80;
    0.87 0.92;
    0.99 0.97];

V_sum_data_sd_PC3_0Gy = [0.009479858 0.016532594 0.0242550332 0.0714584987 0.0963189337 0.0344588741;
    0.0073424551 0.0101382147 0.0258780892 0.0847120964 0.0766670382 0.0301724037;
    0.0056577373 0.0083716803 0.0113536203 0.0565088006 0.0929378678 0.0378267744;
    0.0077258534 0.0156706135 0.0118896548 0.0923792651 0.1044882296 0.0394104073;
    0.0080204573 0.0246143359 0.0285837478 0.1238795509 0.1385946846 0.0198155617];

V_individual_control_sd_PC3_0Gy = [0.001 0.01;
    0.01 0.05;
    0.02 0.01;
    0.01 0.01;
    0.001 0.001];

V_individual_resistant_sd_PC3_0Gy = [0.001 0.001;
    0.01 0.05;
    0.02 0.01;
    0.01 0.01;
    0.001 0.01];

ratiosNumber = 5;

NLLsig = @(x) findLogLikRealerData(x,odds_of_initial_values,V_sum_data_average_PC3_0Gy,V_sum_data_sd_PC3_0Gy,...
    V_individual_control_data_PC3_0Gy,V_individual_control_sd_PC3_0Gy,...
    V_individual_resistant_data_PC3_0Gy,V_individual_resistant_sd_PC3_0Gy,...
    ratiosNumber,x(8),sampled_times_individual_data, ...
    sampled_times_sum_data, x(7),tspan);

NumSample=10000;

unnorm_posterior = @(x) exp(-NLLsig(x))*16*prod(normpdf(x(1:6),0.5*x0(1:6),0.5))*cauchypdf(x(7),0,1)*2*normpdf(x(8),0,0.1);

xs = [];
sigma=0.05;
accepted_num=0;
accepted_num_prev = accepted_num;
filename = 'scatterhist_MCMC_animated.gif';


for i=1:NumSample
   if (mod(i,NumSample/100) == 0) disp(i);
   end
   xnew = proposalSample(x,sigma);
   if (xnew(2) < 0 || xnew(5) < 0 || xnew(1) < 0 || xnew(4) < 0 || xnew(8) < 0) continue;
   end


  %  if (accepted_num > 0)
  %    plot(xs(:,1), 'b')
  %    hold on;
  %    plot(accepted_num,xs(end,1), 'g.', 'MarkerSize',50)
  %  end
  %  plot(accepted_num+1,xnew(1), 'ro', 'MarkerSize',30, 'LineWidth', 5);
  %  xlim([1,accepted_num+2])
  %  xlabel('i');
  %  ylabel('r_c')
  %  drawnow;
  %  hold off;

  index1=3;
  index2=6;
   if (accepted_num ~= accepted_num_prev)
     scatterhist(xs(:,index1), xs(:,index2),'Color','b')
     hold on;
     plot(xs(end,index1),xs(end,index2),'g.', 'MarkerSize',30)
    accepted_num_prev = accepted_num;
   end
   plot(xnew(index1),xnew(index2), 'r.', 'MarkerSize',30);
   if (accepted_num > 0)
     x1=min([0,min(xs(:,index1)), xnew(index1)]);
     x2=max([1.5,max(xs(:,index1)),xnew(index1)]);
     y1=min([0,min(xs(:,index2)), xnew(index2)]);
     y2=max([1.5,max(xs(:,index2)), xnew(index2)]);
     xlim([x1,x2]);
     ylim([y1,y2]);
   end
   drawnow;

   frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end

   r = unnorm_posterior(xnew)*proposalDist(x,xnew,sigma)...
       /(unnorm_posterior(x)*proposalDist(xnew,x,sigma));
   if (rand < r)
       x = xnew;
       xs = [xs; x];
       accepted_num=accepted_num+1;
   end
end
xs=xs(end/2:end,:);
sizes=size(xs);
sum(xs,1)/(sizes(1))
