function createfigure(Vertices1, Faces1, FaceVertexCData1)
%CREATEFIGURE(Vertices1, Faces1, FaceVertexCData1)
%  VERTICES1:  patch vertices
%  FACES1:  patch faces
%  FACEVERTEXCDATA1:  patch facevertexcdata

%  MATLAB���� 12-Dec-2017 14:00:35�� �ڵ� ������

% figure ����
figure1 = figure;
colormap(jet);

% axes ����
axes1 = axes('Parent',figure1);

% patch ����
patch('Parent',axes1,'Vertices',Vertices1,'Faces',Faces1,'FaceColor','flat',...
    'FaceVertexCData',FaceVertexCData1);

% ylabel ����
ylabel('y (m)');

% xlabel ����
xlabel('x (m)');

% title ����
title('Tresca Stress - 20 triangular elements');

% ���� ������ �ּ� ó���� �����Ͽ� ��ǥ���� X ������ ����
xlim(axes1,[0 20]);
% ���� ������ �ּ� ó���� �����Ͽ� ��ǥ���� Y ������ ����
ylim(axes1,[-5 6]);
% ���� ������ �ּ� ó���� �����Ͽ� ��ǥ���� Z ������ ����
zlim(axes1,[0 50000]);
view(axes1,[0.5 90]);
% ������ axes �Ӽ� ����
set(axes1,'OuterPosition',[0 0 0.921171171171171 1],'XGrid','on','YGrid',...
    'on');
% colorbar ����
colorbar('peer',axes1,'Units','centimeters');

