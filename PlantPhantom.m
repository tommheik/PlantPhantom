function [ph, ph_NoRot] = PlantPhantom(N, n, delta, highContrast)
% PLANTPHANTOM creates a (2+1)D phantom of a plant stem with resolution N and n
% time frames to simulate the flow of a fluid such as iodine in the plant.
% Finally a small rotational error of at most delta is added to each layer.
%
% Tommi Heikkilä 2018
% Last edited 8.1.2020
%
% INPUTS
% N             resolution of the phantom
% n             number of time frames, has to be at least 10!
% delta         rotational error
% highContrast  High contrast 'iodine' on/off (1/0)
%
% OUTPUTS
% ph       Desired phantom with the rotational error added
% ph_NoRot Same as ph but without the rotational error for comparison

% Helper function to create the fluid spots
function m = fluidspot(m,x,y,rx,ry,angles)
    % m is a layer we want to add the fluid spots into
    % (x,y) is the center
    % rx and ry are the radii
    % angles is the array of angles used elsewhere in this code

    % M is the edited m
    for t = angles
        Dx = round(rx*cos(t));
        Dy = round(ry*sin(t));
        m(y-Dy:y+Dy, x-Dx:x+Dx) = 1;
    end
end

if n < 10
    error('At least 10 time frames are required!')
end

if N < 32
    warning('Consider increasing the resolution!')
end

ph = zeros([N N n]);
ph_NoRot = ph;
layer = zeros(N);
stagestep = floor((n-5)/5);
stages = 1:stagestep:n-4;
% Note that this vector does not scale with n!
radii = [0.1,0.12,0.15,0.2,0.22,0.25,0.27,0.3];

% Outer border
coeff = 0.9;
x0 = N/2;
y0 = N/2;
Rx = 0.91*N/2;
Ry = 0.87*N/2;

% Choose a sufficient resolution for the angular sampling.
N_theta = 4*N; 
theta = linspace(0,pi/2,N_theta);

for iii = 1:N_theta
    dx = round(Rx*cos(theta(iii)));
    dy = round(Ry*sin(theta(iii)));
    layer(y0-dy:y0+dy, x0-dx:x0+dx) = coeff;
end

coeff = 0.6;
Rx = 0.89*N/2;
Ry = 0.85*N/2;

for iii = 1:N_theta
    dx = round(Rx*cos(theta(iii)));
    dy = round(Ry*sin(theta(iii)));
    layer(y0-dy:y0+dy, x0-dx:x0+dx) = coeff;
end

% Interior
coeff = 0.3;
Rx = 0.85*N/2;
Ry = 0.8*N/2;

for iii = 1:N_theta
    dx = round(Rx*cos(theta(iii)));
    dy = round(Ry*sin(theta(iii)));
    layer(y0-dy:y0+dy, x0-dx:x0+dx) = coeff;
end

% Air gap
coeff = 0.15;
Rx = 0.52*N/2;
Ry = 0.54*N/2;

% 13 sine waves are used to make the interior have a wavy outline.
for iii = 1:N_theta
    dx = round((0.05*sin(13*theta(iii))+1)*Rx*cos(theta(iii)));
    dy = round((0.05*sin(13*theta(iii))+1)*Ry*sin(theta(iii)));
    layer(y0-dy:y0+dy, x0-dx:x0+dx) = coeff;
end

% Core
coeff = 0.31;
Rx = 0.34*N/2;
Ry = 0.32*N/2;

% 11 sine waves are used to make the interior have a wavy outline.
for iii = 1:N_theta
    dx = round((0.05*sin(11*(theta(iii)+1))+1)*Rx*cos(theta(iii)));
    dy = round((0.05*sin(11*(theta(iii)+1))+1)*Ry*sin(theta(iii)));
    layer(y0-dy:y0+dy, x0-dx:x0+dx) = coeff;
end

% Thin walls dividing the air gap
coeff = 0.22;
x0 = N/2;
y0 = N/2;
for dir = [0.15,0.7,1.9,2.45,3.35,4.4,5,5.5] % These are the directions for the dividers
    for iii = 0.32*N/2:0.05:0.54*N/2 % This is the radius
        dx = round(iii*cos(dir));
        dy = round(iii*sin(dir));
        if layer(x0-dx,y0-dy) == 0.15
            layer(x0-dx,y0-dy) = coeff;
        end
    end
end

% Fluid flow and random rotational error
delta = delta*randn(1,n);
coeff = 0.03;
s = 1;
if ~highContrast
    coeff = coeff/6;
end

nospots = layer;    % Save the phantom before adding any fluid spots.
% Convolution kernel for smoothing
p = [1,1,1,1;1,1,1,1;1,1,1,1;1,1,1,1];
    p = conv2(p,p);
    p = conv2(p,p);
    p = p / sum(p(:));

for iii = 1:n
    sublayer = zeros(N);
    % Fluid is added in 5 stages and with growing radii
    % Spot 1
    if iii > stages(1) && iii <= stages(1)+length(radii)
        x0 = round(0.88*N);
        y0 = round(0.5*N);
        sublayer = fluidspot(sublayer,x0,y0,N/10*radii(iii-stages(1)),N/12*radii(iii-stages(1)),theta);
    end
    
    % Spot 2
    if iii > stages(2) && iii <= stages(2)+length(radii)
        x0 = round(0.85*N);
        y0 = round(0.63*N);
        sublayer = fluidspot(sublayer,x0,y0,N/10*radii(iii-stages(2)),N/12*radii(iii-stages(2)),theta);
    end
    
    % Spot 3
    if iii > stages(2)+1 && iii <= stages(2)+length(radii)+1
        x0 = round(0.84*N);
        y0 = round(0.41*N);
        sublayer = fluidspot(sublayer,x0,y0,N/10*radii(iii-stages(2)-1),N/12*radii(iii-stages(2)-1),theta);
    end
    
    % Spot 4
    if iii > stages(3)-1 && iii <= stages(3)+length(radii)-1
        x0 = round(0.8*N);
        y0 = round(0.57*N);
        sublayer = fluidspot(sublayer,x0,y0,N/10*radii(iii-stages(3)+1),N/12*radii(iii-stages(3)+1),theta);
    end
    
    % Spot 5
    if iii > stages(3) && iii <= stages(3)+length(radii)
        x0 = round(0.77*N);
        y0 = round(0.68*N);
        sublayer = fluidspot(sublayer,x0,y0,N/10*radii(iii-stages(3)),N/12*radii(iii-stages(3)),theta);
    end
    
    % Finally spread all of the spots
    if iii > stages(4) && iii <= stages(4)+3
%         sublayer = layer - nospots;
%         sublayer = s*conv2(sublayer,p,'same');
%         sublayer(sublayer < 0.5) = 3*sublayer(sublayer < 0.5);
        radii = [0.6,0.95,1.1];
        s = 5; % Secret multiplier
        
        x0 = round(0.88*N); % Spot 1
        y0 = round(0.5*N);
        sublayer = fluidspot(sublayer,x0,y0,N/11*radii(iii-stages(4)),N/12*radii(iii-stages(4)),theta);
        x0 = round(0.85*N); % Spot 2
        y0 = round(0.63*N);
        sublayer = fluidspot(sublayer,x0,y0,N/11*radii(iii-stages(4)),N/12*radii(iii-stages(4)),theta);
        x0 = round(0.84*N); % Spot 3
        y0 = round(0.41*N);
        sublayer = fluidspot(sublayer,x0,y0,N/12*radii(iii-stages(4)),N/11*radii(iii-stages(4)),theta);
        x0 = round(0.8*N); % Spot 4
        y0 = round(0.57*N);
        sublayer = fluidspot(sublayer,x0,y0,N/11*radii(iii-stages(4)),N/12*radii(iii-stages(4)),theta);
        x0 = round(0.77*N); % Spot 5
        y0 = round(0.68*N);
        sublayer = fluidspot(sublayer,x0,y0,N/10*radii(iii-stages(4)),N/10*radii(iii-stages(4)),theta);
        sublayer = conv2(sublayer,p,'same');
    end
    sublayer(nospots < 0.3) = 2*sublayer(nospots < 0.3);    % Fluid spreads faster in the air gap
    sublayer(nospots == 0) = 0;
    sublayer(nospots >= 0.6) = 0;  % Fluid can't pass the outer border
    sublayer(sublayer > 0.7) = 0.7;
  
    layer = layer + s*coeff*sublayer;
    layer(layer > 1) = 1;
    
    ph(:,:,iii) = imrotate(layer,delta(iii),'crop');
    ph_NoRot(:,:,iii) = layer;
end
end

