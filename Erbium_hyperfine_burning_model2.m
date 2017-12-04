%Erbium Model 2.0
clc, clear, close all
%% Generate initial populations
%At 1kHz res this section will take ~2 minutes
res = 1; %Resolution in MHz
laser_width = 1; %FWHM of laser (MHz)
laser_wavelength = -964+1; %Where the laser is burning the spectrum

%Get energy level of hyperfine structure from previous work
levels = open('Erbium_hyperfine_energy_levels.mat');
levels = levels.levels();

for i = 1:15
    %Plot position of each hyperfine level
    hold on
    plot(levels(:,i),zeros(8),'x')
    xlim([-1.5e3 1.5e3])
end

%Oscillator strnegths of delta m = -1,0,1
osc_str(1,:) = [0,linspace(0.25,0.025,7)]; %-1
osc_str(3,:) = [linspace(0.33,0.033,7),0]; %+1
osc_str(2,:) = ones(1,8);   


model_levels = levels(:,7:9);%takes only the levels of the delta m = -1,0,1
model_levels(1:7,3) = model_levels(2:8,3); %Shifts the d_m = 1 column to fit osc_str
model_levels(8,3) = 0;
model_levels = model_levels';

d_model_levels=  zeros(3,8);
for i = 1:3
    for j = 1:7
    
        d_model_levels(i,j) = model_levels(i,j+1)-model_levels(i,j);

    end
end

d_model_levels(1,1) = 0;
d_model_levels(3,2:7) = d_model_levels(3,1:6);
d_model_levels(2,2:8) = d_model_levels(2,1:7);
d_model_levels(3,1) = 0;
d_model_levels(2,1) = 0;


%Create a frequency array and pre-allocate two matrices for hyperfine
%populations
freq_array = -1.5e3:res:2e3;
each_level_abs = zeros(length(freq_array),3,8);
each_level = each_level_abs;

%Create voigt profiles for each hyperfine level taking into account of
%oscillator strength and eventually population.
%population = [0   0   7   0   0   0   0   0];
population = ones(1,8);
for j = 1:3
    for i = 1:8
        temp_abs = voigt(freq_array,model_levels(j,i),45,45);
        temp_abs = temp_abs./max(temp_abs);
        %Keep note of population of each state at each frequency
        each_level(:,j,i) = temp_abs*population(i); %Important for later
        each_level_abs(:,j,i) = osc_str(j,i)*temp_abs*population(i);
    end
end

%Create the total absorption spectrum including isotopic impurities.
total_abs = sum(sum(permute(each_level_abs,[3,2,1])));
total_abs = permute(total_abs,[3,2,1]);
absorption_impurities = voigt(freq_array,180,45,45); %Adds absorption due to isotopic impurities
absorption_impurities = absorption_impurities/max(absorption_impurities)*0.68;
total_abs = total_abs + absorption_impurities;

plot(freq_array,total_abs)
%Plot population of each state.
for j = 1:3
    for i = 1:8
        plot(freq_array,each_level_abs(:,j,i),'r:')
    end
end

%% Create a Guassian Laser
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

%This if statement determines which hyperfine level the laser is hitting.
if laser_wavelength < 750;
    if laser_wavelength < -381
        laser_state = 1; %delta m = -1
    else
        laser_state = 2; %delta m = 0
    end
else
    laser_state = 3; %delta m = 1
end
ls = laser_state;
%Creates a laser with guassian dist.
laser_array = [-laser_width*4:res:laser_width*4]+laser_wavelength;
laser = voigt(laser_array,laser_wavelength,laser_width,0);

if length(laser_array) < 10;
    disp('Warning: laser width small, consider a finer resolution')
end

excited_state = zeros(size(each_level_abs));
excited = zeros(length(freq_array),8);

%Sets index to be where the laser array starts in the freq array
if mod(length(laser_array),2) == 0;
   index = index - length(laser_array)/2;
else
   index =  index - round(length(laser_array)/2);
end
in = index;

%% Time to simulate
steps = 1000;
ts = 1e-4; %Time for each simulation step in s
laser = laser*100; %Allows me to change "laser power".

%Create 'normalised' oscillator strengths, this is just the probability
%an excited state will decay by d_m = +1,0,-1
for i = 1:8
    norm = sum(osc_str(:,i));
    osc_strn(:,i) = osc_str(:,i)/norm;
end

%finds the index for each energy level
model_levels_res = round(model_levels/res);
%Finds the index change for the change in each energy level
d = round(d_model_levels/res);
d(1,2:8) = d(1,1:7);

decay = exp(-ts/1e-2);% '%' of how much will decay in a given time step (10ms life time)
decay_spectrum = excited_state; %Pre-allocates a decay spectrum 

inn = in+length(laser)-1;

for k = 1:steps
for j = 1:3
for i = 1:8      
    %energy difference between the hyperfine levels of a given state
    ml = model_levels_res(j,i) - model_levels_res(ls,i) + 1;
    
    %If the laser is burning on the delta m = 0
    if ls == 2
        excited_state(in+ml:inn+ml,j,i) = excited_state(in+ml:inn+ml,j,i) ...
            + each_level(in:inn,ls,i).*(1-exp(-laser*ts));                
    %Is the laser is burning on the delta m = -1
    elseif ls == 1
        excited_state(in+ml-d(j,i):inn+ml-d(j,i),j,i) = excited_state(in+ml-d(j,i):inn+ml-d(j,i),j,i) ...
            + each_level(in:inn,ls,i).*(1-exp(-laser*ts));
    %If the laser is burning on the delta m = 1
    else
    if i < 8
        excited_state(in+ml+d(j,i+1):inn+ml+d(j,i+1),j,i) = excited_state(in+ml+d(j,i+1):inn+ml+d(j,i+1),j,i) ...
            + each_level(in:inn,ls,i).*(1-exp(-laser*ts));
    end
    end
    %Adds population back due to decay from excited state

    each_level(1:end-1,j,i) = each_level(1:end-1,j,i) + decay_spectrum(2:end,j,i)*osc_strn(2,i);
    
    
    if ls == 1
    if i > 1
    if (j == 3 && i == 3) || (j == 2 && i == 3) 
        each_level(1-d(j,i-1):end-1,j,i) = each_level(1-d(j,i-1):end-1,j,i) + decay_spectrum(2:end+d(j,i-1),j,i)*osc_strn(1,i-1);
    else
        each_level(1:end-1-d(j,i-1),j,i) = each_level(1:end-1-d(j,i-1),j,i) + decay_spectrum(2+d(j,i-1):end,j,i)*osc_strn(1,i-1);
    end
    end
    if (j == 3 && i == 2) || (j == 2 && i == 2)
        each_level(2:end+d(j,i),j,i) = each_level(2:end+d(j,i),j,i) + decay_spectrum(1-d(j,i):end-1,j,i)*osc_strn(3,i);
    else
        each_level(d(j,i)+1:end-1,j,i) = each_level(d(j,i)+1:end-1,j,i) + decay_spectrum(2:end-d(j,i),j,i)*osc_strn(3,i);
    end
    
    elseif ls == 2
    if (j == 3 && i == 2) || (j == 2 && i == 2)
        each_level(1-d(j,i):end-1,j,i) = each_level(1-d(j,i):end-1,j,i) + decay_spectrum(2:end+d(j,i),j,i)*osc_strn(1,i);
    else
        each_level(1:end-1-d(j,i),j,i) = each_level(1:end-1-d(j,i),j,i) + decay_spectrum(2+d(j,i):end,j,i)*osc_strn(1,i);
    end
    if i < 8
    if (j == 3 && i == 1) || (j == 2 && i == 1)
        each_level(2:end+d(j,i+1),j,i) = each_level(2:end+d(j,i+1),j,i) + decay_spectrum(1-d(j,i+1):end-1,j,i)*osc_strn(3,i);
    else
        each_level(d(j,i+1)+1:end-1,j,i) = each_level(d(j,i+1)+1:end-1,j,i) + decay_spectrum(2:end-d(j,i+1),j,i)*osc_strn(3,i);
    end        
    end
    
    else
    if (j == 3 && i == 1) || (j == 2 && i == 1)
    if i < 8
        each_level(1-d(j,i+1):end-1,j,i) = each_level(1-d(j,i+1):end-1,j,i) + decay_spectrum(2:end+d(j,i+1),j,i)*osc_strn(1,i);
    else
        each_level(1:end-1-d(j,i+1),j,i) = each_level(1:end-1-d(j,i+1),j,i) + decay_spectrum(2+d(j,i+1):end,j,i)*osc_strn(1,i);
    end
    end
    if i < 7
        each_level(d(j,i+2)+1:end-1,j,i) = each_level(d(j,i+2)+1:end-1,j,i) + decay_spectrum(2:end-d(j,i+2),j,i)*osc_strn(3,i);
    end        
    end     
    
    

    %Population burnt out of ground state
    each_level(in+ml-1:inn+ml-1,j,i) = each_level(in+ml-1:inn+ml-1,j,i).*exp(-laser*ts);
    
    
    decay_spectrum(:,j,i) = excited_state(:,j,i)*(1-decay);
end
end
    excited_state = excited_state - decay_spectrum;
end

for j = 1:3
    for i = 1:8 
        excited_state(:,j,i) = excited_state(:,j,i)*osc_str(j,i);
        each_level_abs(:,j,i) = each_level(:,j,i)*osc_str(j,i);  
        decay_spectrum(:,j,i) = decay_spectrum(:,j,i)*osc_str(j,i);
    end
end


final_abs = squeeze(sum(sum(permute(each_level_abs,[3,2,1]))));
final_abs = final_abs + absorption_impurities;

final_excite = squeeze(sum(sum(permute(excited_state,[3,2,1]))));

figure()

subplot(2,1,1)
hold on
plot(freq_array,final_abs)
plot(freq_array,total_abs)
for j = 1:3
    for i = 1:8
        plot(freq_array,each_level_abs(:,j,i),'g')
    end
end
%  xlim([27 35])
% ylim([0 0.9])

subplot(2,1,2)
% plot(freq_array,final_excite)
hold on
for j = 1:3
    for i = 1:8
        plot(freq_array,excited_state(:,j,i))
        %xlim([27 35])
    end
end

% figure()
% plot(freq_array,squeeze(sum(sum(permute(decay_spectrum,[3,2,1])))))

t = trapz(total_abs);
tt = trapz(final_abs+final_excite);
fprintf('Initial atoms: %.1f. \n  Final atoms: %.1f. \n',t,tt)

        










