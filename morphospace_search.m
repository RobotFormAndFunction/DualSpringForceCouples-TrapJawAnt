%%% dual spring force couple mandible simulation
%%% R. St. Pierre - University at Buffalo

clear all; close all; clc;

%%
Lratio=[.8:.01:.99]; %lever ratio
ktval=[3000]; %value of the spring constant for muscle apodeme unit
partition=[.01:.01:.99]; % energy partition 


for index1=1:length(ktval)
    for index2=1:length(partition)
        for index3=1:length(Lratio)
            %calling the EOM function, A1 and A2 are just dummy variables
            %for the power outputs from that function
            [A1,A2]=mandiblelever_energypartion_fxn(Lratio(index3),...
                partition(index2),ktval(index1));

            poweroutputmat(index1,index2,index3)=A1;
            poweroutputmat2(index1,index2,index3)=A2;
          
        end
    end
end

%% 
% if the spring constant (kmau) is a vector, this will resort the matrices
B=poweroutputmat;
C = permute(B,[2 3 1]);

B2=poweroutputmat2;
C2 = permute(B2,[2 3 1]);


for index10=1:length(ktval)
    designheatmaps{index10,1}=C(:,:,index10);
end

for index10=1:length(ktval)
    designheatmaps2{index10,1}=C2(:,:,index10);
end


% plotting the contour maps of power within the energy partition and
% mandible level spaces

figure(1)
contourf(partition(1:end),Lratio,transpose(designheatmaps{index11,1}(1:end,1:end)))
ylabel('L_c/L_t')
xlabel('\alpha')
title(strcat('k_t=',mat2str(ktval(index11)),' N/m'))
colorbar
hold on
plot(.5, .9, 'ko') %plotting animal data
hold off



figure(2)
contourf(partition(1:end),Lratio,transpose(designheatmaps2{index11,1}(1:end,1:end)))
ylabel('L_c/L_t')
xlabel('\alpha')
title(strcat('k_t=',mat2str(ktval(index11)),' N/m'))
title('Power (W)')
colorbar
hold on
plot(.5, .9, 'ko') %plotting animal data
hold off