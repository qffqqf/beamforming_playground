clc
clear all
freq = 150;
kWave = 2*pi*freq/340;

%% Get mics
X_mic = 1;
Y_mic = 1;
nMic_x = 3;
nMic_y = 3;
nMic = nMic_x*nMic_y;
rect_mic = [-X_mic,X_mic;...
            -Y_mic,Y_mic];
z_mic = 1;
[pos, ~] = get_grid(rect_mic, z_mic, nMic_x, nMic_y);
mic_pos = pos';

%% data
soc_pos = [1.2,-0.7,0.01];
dist_x = soc_pos(1)*ones(1,nMic) - mic_pos(1,:);
dist_y = soc_pos(2)*ones(1,nMic) - mic_pos(2,:);
dist_z = soc_pos(3)*ones(1,nMic) - mic_pos(3,:);
dist = sqrt(dist_x.^2 + dist_y.^2 + dist_z.^2);
soc_pos_noise1 = [1, 0.7, 0];
dist_x = soc_pos_noise1(1)*ones(1,nMic) - mic_pos(1,:);
dist_y = soc_pos_noise1(2)*ones(1,nMic) - mic_pos(2,:);
dist_z = soc_pos_noise1(3)*ones(1,nMic) - mic_pos(3,:);
dist_noise1 = sqrt(dist_x.^2 + dist_y.^2 + dist_z.^2);
soc_pos_noise2 = [-2, 0.7, 0];
dist_x = soc_pos_noise2(1)*ones(1,nMic) - mic_pos(1,:);
dist_y = soc_pos_noise2(2)*ones(1,nMic) - mic_pos(2,:);
dist_z = soc_pos_noise2(3)*ones(1,nMic) - mic_pos(3,:);
dist_noise2 = sqrt(dist_x.^2 + dist_y.^2 + dist_z.^2);
pres = - exp(1i*kWave*dist)./(4*pi*dist) - 1*exp(1i*kWave*dist_noise1)./(4*pi*dist_noise1)...
       - 1*exp(1i*kWave*dist_noise2)./(4*pi*dist_noise2);

%% Get sources
nSoc_x = 500;
nSoc_y = 500;
lSoc = 20;
rect_soc = [-lSoc,lSoc;...
            -lSoc,lSoc];
z_soc = 0;
[soc_pos_all, dA] = get_grid(rect_soc, z_soc, nSoc_x, nSoc_y);
E_all = compute_energy(soc_pos_all, mic_pos, kWave, dA);


%% Monitoring
nPix_x = 30;
nPix = nPix_x^2;
LX = 5;
LY = 5;
rect_moni = [-LX,LX;...
            -LY,LY];
[soc_moni, ~] = get_grid(rect_moni, 0, nPix_x, nPix_x);
focus = 0.1;
nInt = 50;
energy = zeros(nPix,1);
parfor iPix = 1:nPix
    iPix
    center = soc_moni(iPix,:);
    [soc_pos_focus, dA] = get_grid([center(1)-focus/2, center(1)+focus/2; center(2)-focus/2, center(2)+focus/2], ...
                                            0, nInt, nInt);
%     [soc_pos_focus, dA] = get_grid_circle([center(1)-5*focus, center(1)+5*focus; center(2)-5*focus, center(2)+5*focus], ...
%                                             0, nInt, nInt, [center(1), center(2), 0], focus/2);
    E_focus = compute_energy(soc_pos_focus, mic_pos, kWave, dA);
    [beamformer, enhance_rate] = eigs(E_focus, E_all, 1, 'largestabs');
    energy(iPix) = abs(beamformer'*pres.');
end
% energy = reshape(energy,[nPix_x,nPix_x]);

%% plot
figure
pcolor(reshape(soc_moni(:,1),[nPix_x,nPix_x]), reshape(soc_moni(:,2),[nPix_x,nPix_x]), reshape(energy,[nPix_x,nPix_x]))
shading interp;
hold on
scatter(mic_pos(1,:),mic_pos(2,:),'MarkerFaceColor',[0,0,0])
scatter(soc_pos(1),soc_pos(2), 60, 'd', 'MarkerEdgeColor',[0 .5 .5],...
                                          'MarkerFaceColor',[0 .7 .7],...
                                          'LineWidth',1.5)

scatter(soc_pos_noise1(1),soc_pos_noise1(2), 60, 's', 'MarkerEdgeColor',[0 .5 .5],...
                                          'MarkerFaceColor',[0 .7 .7],...
                                          'LineWidth',1.5)

scatter(soc_pos_noise2(1),soc_pos_noise2(2), 60, 's', 'MarkerEdgeColor',[0 .5 .5],...
                                          'MarkerFaceColor',[0 .7 .7],...
                                          'LineWidth',1.5)




function [pos, dA] = get_grid(rect, z, nX, nY)
    x = linspace(rect(1,1), rect(1,2), nX);
    y = linspace(rect(2,1), rect(2,2), nY);
    [X,Y] = meshgrid(x,y);
    pos = [X(:), Y(:), z*ones(nX*nY,1)];
    dA = (rect(1,2) - rect(1,1))*(rect(2,2) - rect(2,1))/nX/nY;
end

function [pos, dA] = get_grid_circle(rect, z, nX, nY, center, radius)
    x = linspace(rect(1,1), rect(1,2), nX);
    y = linspace(rect(2,1), rect(2,2), nY);
    [X,Y] = meshgrid(x,y);
    pos = [X(:), Y(:), z*ones(nX*nY,1)];
    nPoint = numel(X);
    dist = pos - ones(nPoint,1)*center;
    R = sqrt(dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2);
    index_delete = find(R > radius);
    pos(index_delete,:) = [];
    dA = (rect(1,2) - rect(1,1))*(rect(2,2) - rect(2,1))/nX/nY;
end

function E = compute_energy(soc_pos, mic_pos, kWave, dA)
    nMic = size(mic_pos,2);
    nSoc = size(soc_pos,1);
    %% Get distance
    dist_x = soc_pos(:,1)*ones(1,nMic) - ones(nSoc,1)*mic_pos(1,:);
    dist_y = soc_pos(:,2)*ones(1,nMic) - ones(nSoc,1)*mic_pos(2,:);
    dist_z = soc_pos(:,3)*ones(1,nMic) - ones(nSoc,1)*mic_pos(3,:);
    dist = sqrt(dist_x.^2 + dist_y.^2 + dist_z.^2);
    %% Get pressure
    pres = -exp(1i*kWave*dist)./(4*pi*dist);
    %% Integration
    E = pres'*pres*dA;
end