clc
clear all
freq = 150;
kWave = 2*pi*freq/340;

%% Get mics
X_mic = 2;
Y_mic = 2;
nMic_x = 10;
nMic_y = 10;
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

%% Monitoring
nPix_x = 10;
nPix = nPix_x^2;
LX = 5;
LY = 5;
rect_moni = [-LX,LX;...
            -LY,LY];
[soc_moni, ~] = get_grid(rect_moni, 0, nPix_x, nPix_x);
G = compute_transfer_mx(soc_moni, mic_pos, kWave);
eta = 1e-15;
energy = abs(G'*((G*G'+ eta*eye(nMic))\pres.'));
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

function G = compute_transfer_mx(soc_pos, mic_pos, kWave)
    nMic = size(mic_pos,2);
    nSoc = size(soc_pos,1);
    %% Get distance
    dist_x = soc_pos(:,1)*ones(1,nMic) - ones(nSoc,1)*mic_pos(1,:);
    dist_y = soc_pos(:,2)*ones(1,nMic) - ones(nSoc,1)*mic_pos(2,:);
    dist_z = soc_pos(:,3)*ones(1,nMic) - ones(nSoc,1)*mic_pos(3,:);
    dist = sqrt(dist_x.^2 + dist_y.^2 + dist_z.^2).'; % mxn
    %% Get pressure
    G = -(exp(1i*kWave*dist)./(4*pi*dist));
end