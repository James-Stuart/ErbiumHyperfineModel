%Erbium_hyperfine_model 3.0 Almost a new hope.
%Update log:
%changing dimensions of energy level and oscillator strength matrices to be
%8x8's with d_m = 0 on the mian diagonal.
%Changing energy level populations to be size freq_array x 8. Instead of
%freq_array x 3 x 8.
%Now runs an order of magnitude faster.

clc, clear, close all
%% Generate initial populations
res = 1; %Resolution in MHz

%Get energy level of hyperfine structure from previous work
levels = open('Erbium_hyperfine_energy_levels.mat');
levels = levels.levels();

for i = 1:15
    %Plot position of each hyperfine level
    hold on
    plot(levels(:,i),zeros(8),'x')
    xlim([-1.5e3 1.5e3])
end

%Create oscillator strengths matrix, diag is delta m = 0, off diag are d_m = -1,1.
osc_str = toeplitz([1, zeros(1,7)]);
osc_str = osc_str + diag(linspace(0.33,0.033,7),1) + diag(linspace(0.25,0.025,7),-1);

%Matrix of the energy levels for the d_m = -1,0,1
model_levels = diag(levels(:,8)) + diag(levels(2:8,9),1) + diag(levels(2:8,7),-1);
index_model_levels = round(model_levels/res);

%matrix of difference between adjact energy levels
d_model_levels = diag(levels(2:8,8) - levels(1:7,8)) + diag(levels(3:8,9) ...
    - levels(2:7,9),1) + diag(levels(3:8,7) - levels(2:7,7),-1);
index_d_model_levels = round(d_model_levels/res);

%Creates a matrix of difference in energies between dm = -1,0,1.
d_m = zeros(8,8);
for i = 1:7
    d_m(i,i+1) = d_m(i,i+1) - model_levels(i,i) + model_levels(i,i+1);
    d_m(i+1,i) = d_m(i+1,i) - model_levels(i+1,i+1) + model_levels(i+1,i);
end
index_d_m = round(d_m/res);




%Create a frequency array and pre-allocate two matrices for hyperfine
%populations
freq_array = -1.5e3:res:2e3;
ground_abs = zeros(length(freq_array),8);
ground_level = ground_abs;

%Create voigt profiles for each hyperfine level
population = [0   0   0   0   0   0.00   0.00   1];
% population = ones(1,8); %Populations of each hyperfine state.
population = population*7/sum(population);
for i = 1:8
    ground_level(:,i) = voigt(freq_array,model_levels(i,i),45,45);
    ground_level(:,i) = ground_level(:,i)/max(ground_level(:,i))*population(i);
end

total_abs = sum(permute(ground_level,[2,1]));
neg1 = zeros(length(freq_array),1);
pos1 = zeros(length(freq_array),1);
for i = 1:7
    neg1 = neg1 + shift_array(ground_level(:,i+1),index_d_m(i+1,i))*osc_str(i+1,i);
    pos1 = pos1 + shift_array(ground_level(:,i),index_d_m(i,i+1))*osc_str(i,i+1);
end
total_abs = total_abs';
total_abs = total_abs + pos1+neg1;
absorption_impurities = voigt(freq_array,180,45,45); %Adds absorption due to isotopic impurities
absorption_impurities = absorption_impurities/max(absorption_impurities)*0.68;
total_abs = total_abs + absorption_impurities;

figure(1)
hold on
plot(freq_array,total_abs,'blue')

freq_neg = zeros(length(freq_array),7);
freq_pos = freq_neg;
for i = 1:7
    freq_neg(:,i) = freq_array-d_m(i+1,i);
    freq_pos(:,i) = freq_array-d_m(i,i+1);
end




%% Create a Guassian laser
laser_width = 1; %FHWM of laser (MHz)
laser_wavelength_vec = fliplr([-0.9447,-0.8707,-0.7797,-0.6757,-0.5677,-0.4607])*1e3;

for xx = 1:length(laser_wavelength_vec)
laser_wavelength = laser_wavelength_vec(xx);    
%Round laser wavelength to nearest freq_array value (could be done with
%while loop but what ever)
for i = 1:length(freq_array)
    if laser_wavelength < freq_array(i)
        low = laser_wavelength - freq_array(i-1);
        high = freq_array(i) - laser_wavelength;
        if high < low
            laser_wavelength = freq_array(i);
            index = i;
        else
            laser_wavelength = freq_array(i-1);
            index = i-1;
        end
        break
    end
end

%This if statement determines which transition the laser is hitting.
if laser_wavelength < 750
    if laser_wavelength < -381
        laser_state = 1; %delta m = -1
    else
        laser_state = 2; %delta m = 0
    end
else
    laser_state = 3; %delta m = 1
end

ls = laser_state;

%Creates a laser with guassian dist. But only out to 4 sigma, to save
%computation time.
laser_array = [-laser_width*4:res:laser_width*4]+laser_wavelength;
laser = voigt(laser_array,laser_wavelength,laser_width,0);

if length(laser_array) < 10;
    disp('Warning: laser width small, consider a finer resolution')
end

excited_level = zeros(size(ground_level));

%Sets 'index' to be where the laser array starts within the freq array
if mod(length(laser_array),2) == 0;
   index = index - length(laser_array)/2;
else
   index =  index - round(length(laser_array)/2);
end
in = index;
inn = in+length(laser)-1; %Last index of the laser array





%% Simulation
steps = 5; %Number of time steps
ts = 5e-3; %Time per step
lt = 1e-2; %Life time of excited level
laser = laser*1000; %Change laser power

%Probability of decaying from excited state to ground state via dm = -1,0,1
%based on oscillator strengths
osc_strn = zeros(8,8);
for i = 1:8
    s = sum(osc_str(i,:));
    osc_strn(i,:) = osc_str(i,:)/s;
end

decay = exp(-ts/lt); %Ratio of the excited state that will fall to ground
decay_level = excited_level; %Pre-allocate decay array

%Difference in energies between d_m = 0;
sh = diag(index_d_model_levels);

fprintf('Simulation resolution: %.3f MHz. \nLaser wavelength: %.1f MHz. \n',res, laser_wavelength)
fprintf('Time resolution: %.2f ms. \nSimulation time: %.4f s. \n\n',ts*1e3, steps*ts)    

for j = 1:steps
    for i = 1:8
        
        %If the laser is hitting the d_m = -1 feature, we ignore |-7/2> and
        %we 'move' the laser so that it hits the d_m = 0 feature
        if ls == 1
        if i > 1 %Ignore |-7/2>
            in1 = in - d_m(i,i-1); %"moves" the laser by changing the index of the laser
            inn1 = inn - d_m(i,i-1);
            excited_level(in1-sh(i-1):inn1-sh(i-1),i-1) = excited_level(in1-sh(i-1):inn1-sh(i-1),i-1) ...
                + ground_level(in1:inn1,i).*(1-exp(-laser*ts)); %|e> gains population due to laser burning |g>
            ground_level(in1:inn1,i) = ground_level(in1:inn1,i).*exp(-laser*ts); %|g> loses population due to laser.
            
            %|g> gains population due to |e> decaying via delta m = 1
            ground_level(in1:inn1,i) = ground_level(in1:inn1,i) + decay_level(in1-sh(i-1):inn1-sh(i-1),i-1)*osc_strn(i-1,i);
            %|g> gains population due to |e> decaying via delta m = 0
            ground_level(in1-sh(i-1):inn1-sh(i-1),i-1) = ground_level(in1-sh(i-1):inn1-sh(i-1),i-1) + decay_level(in1-sh(i-1):inn1-sh(i-1),i-1)*osc_strn(i-1,i-1);
            if i > 2 %ignore |-7/2,-5/2>
            %|g> gains population due to |e> decaying via delta m = -1
            ground_level(in1-sh(i-1)-sh(i-2):inn1-sh(i-1)-sh(i-2),i-2) = ground_level(in1-sh(i-1)-sh(i-2):inn1-sh(i-1)-sh(i-2),i-2) + decay_level(in1-sh(i-1):inn1-sh(i-1),i-1)*osc_strn(i-1,i-2);
            end
        end
        
        %if the laser is hitting the d_m = 1 feature, we ignore the |7/2>
        %and like before 'move' the laser
        elseif ls == 3
        if i < 8 %Ignore |7/2>
            in1 = in - d_m(i,i+1);
            inn1 = inn - d_m(i,i+1);
            excited_level(in1+sh(i):inn1+sh(i),i+1) = excited_level(in1+sh(i):inn1+sh(i),i+1) ...
                + ground_level(in1:inn1,i).*(1-exp(-laser*ts));            
            ground_level(in1:inn1,i) = ground_level(in1:inn1,i).*exp(-laser*ts);

            
            
            ground_level(in1:inn1,i) = ground_level(in1:inn1,i) + decay_level(in1+sh(i):inn1+sh(i),i+1)*osc_strn(i+1,i);
            ground_level(in1+sh(i):inn1+sh(i),i+1) = ground_level(in1+sh(i):inn1+sh(i),i+1) + decay_level(in1+sh(i):inn1+sh(i),i+1)*osc_strn(i+1,i+1);
            if i < 7
            ground_level(in1+sh(i)+sh(i+1):inn1+sh(i)+sh(i+1),i+2) = ground_level(in1+sh(i)+sh(i+1):inn1+sh(i)+sh(i+1),i+2) + decay_level(in1+sh(i):inn1+sh(i),i+1)*osc_strn(i+1,i+2);
            end            
            
            
            
        end        
        
        else %If ls == , in this case i dont have to 'move' the laser so there is no need for in1 and inn1
            excited_level(in:inn,i) = excited_level(in:inn,i) + ground_level(in:inn,i).*(1-exp(-laser*ts));           
            ground_level(in:inn,i) = ground_level(in:inn,i).*exp(-laser*ts);
            
            if i > 1
            ground_level(in-sh(i-1):inn-sh(i-1),i-1) = ground_level(in-sh(i-1):inn-sh(i-1),i-1) + decay_level(in:inn,i)*osc_strn(i,i-1);
            end
            
            ground_level(in:inn,i) = ground_level(in:inn,i) + decay_level(in:inn,i)*osc_strn(i,i);
            
            if i < 8
            ground_level(in+sh(i):inn+sh(i),i+1) = ground_level(in+sh(i):inn+sh(i),i+1) + decay_level(in:inn,i+1)*osc_strn(i,i+1);
            end                   
            
         end
        
        
        
    end
    decay_level = excited_level*(1-decay);
    excited_level = excited_level*decay;
end

%At the end of the simulation create an array for the dm = -1,1
ground_pos = zeros(length(freq_array),1);
ground_neg = ground_pos;
for i = 1:7
    ground_neg = ground_neg + shift_array(ground_level(:,i+1),index_d_m(i+1,i))*osc_str(i+1,i);
    ground_pos = ground_pos + shift_array(ground_level(:,i),index_d_m(i,i+1))*osc_str(i,i+1);
end

final_spectrum = sum(ground_level')' + ground_neg + ground_pos + absorption_impurities;

figure()
hold on
subplot(2,1,1)
plot(freq_array,final_spectrum)
hold on
for i = 1:7
    plot(freq_array,shift_array(ground_level(:,i+1),index_d_m(i+1,i))*osc_str(i+1,i),'r')
    plot(freq_array,shift_array(ground_level(:,i),index_d_m(i,i+1))*osc_str(i,i+1),'r')
    plot(freq_array,ground_level(:,i),'r')
end
plot(freq_array,ground_level(:,8),'r')

subplot(2,1,2)
hold on
for i = 1:8
    plot(freq_array,excited_level(:,i))
end

t = trapz(total_abs);
tt = trapz(final_spectrum+sum(excited_level')');
fprintf('Initial atoms: %.1f. \nFinal atoms: %.1f. \n=======Finished=======\n',t,tt)
end























    