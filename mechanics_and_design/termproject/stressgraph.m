function createfigure(Vertices1, Faces1, FaceVertexCData1)
%CREATEFIGURE(Vertices1, Faces1, FaceVertexCData1)
%  VERTICES1:  patch vertices
%  FACES1:  patch faces
%  FACEVERTEXCDATA1:  patch facevertexcdata

%  MATLAB에서 12-Dec-2017 14:00:35에 자동 생성됨

% figure 생성
figure1 = figure;
colormap(jet);

% axes 생성
axes1 = axes('Parent',figure1);

% patch 생성
patch('Parent',axes1,'Vertices',Vertices1,'Faces',Faces1,'FaceColor','flat',...
    'FaceVertexCData',FaceVertexCData1);

% ylabel 생성
ylabel('y (m)');

% xlabel 생성
xlabel('x (m)');

% title 생성
title('Tresca Stress - 20 triangular elements');

% 다음 라인의 주석 처리를 제거하여 좌표축의 X 제한을 유지
xlim(axes1,[0 20]);
% 다음 라인의 주석 처리를 제거하여 좌표축의 Y 제한을 유지
ylim(axes1,[-5 6]);
% 다음 라인의 주석 처리를 제거하여 좌표축의 Z 제한을 유지
zlim(axes1,[0 50000]);
view(axes1,[0.5 90]);
% 나머지 axes 속성 설정
set(axes1,'OuterPosition',[0 0 0.921171171171171 1],'XGrid','on','YGrid',...
    'on');
% colorbar 생성
colorbar('peer',axes1,'Units','centimeters');

