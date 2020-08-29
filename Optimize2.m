%For loop method is very slow
%alpha01=0.3; beta01=0.4; gamma01=0.3; beta02=0.3; gamma02 = 0.2; % Designed Parameters]
alpha01Vec = 0:0.05:2;
beta01Vec = 0:0.05:1;
gamma01Vec = 0:0.05:1;

energyData = []
counter = 1
numIterations = length(alpha01Vec)*length(beta01Vec)*length(gamma01Vec);
for alpha01=alpha01Vec
    for beta01=beta01Vec
    	for gamma01 = gamma01Vec
            energy = EConsumpFunc2([alpha01, beta01, gamma01]);
            energyData(counter, :) = [alpha01, beta01, gamma01, energy];
            counter = counter+1;
            disp(['Percent Complete: ', num2str(counter/numIterations*100)]);
        end
    end
end
[minEnergy, minEnergyIdx] = min(energyData(:,4))