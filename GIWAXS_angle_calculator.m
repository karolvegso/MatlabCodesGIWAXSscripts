%**************************************************************************
% Author: Karol Vegso
% Affiliation: Institute of Physics, Sloavak Academy of Sciences
% 
% You define unit cell of crystal by basis vectors a1, a2, a3.
% Then, define reciprocal space vectors h1 and h2.
% Program calculates angle between reciprocal space vectors h1 and h2 in
% variable angle.
% Then, program calculates projections of reciprocal space vectors h1 and h2
% on the plane represented by reciprocal space vector h perpendicular to
% that plane. Finally, it calculates angle between those projected vectors 
% in variable angle_proj.
%**************************************************************************
clear all
close all
clc
%**************************************************************************
% length of basis vector a1
a1=8.8392; % [A]
% length of basis vector a2
a2=8.8392; % [A]
% length of basis vector a3
a3=12.6948; % [A]
% angle between a2 and a3 basis vectors
alpha=90; % [deg]
% angle between a1 and a3 basis vectors
beta=90; % [deg]
% angle between a1 and a2 basis vectors
gama=90; % [deg]
%**************************************************************************
% calculate metric tensor g in real space
% first row
g11=a1*a1;
g12=a1*a2*cosd(gama);
g13=a1*a3*cosd(beta);
% second row
g21=a2*a1*cosd(gama);
g22=a2*a2;
g23=a2*a3*cosd(alpha);
% third row
g31=a3*a1*cosd(beta);
g32=a3*a2*cosd(alpha);
g33=a3*a3;
% define metric tensor g in real space
g=[g11 g12 g13; g21 g22 g23; g31 g32 g33];
% define unity matrix
I=eye(3,3);
% calculate metric tenosr g* in reciprocal space
g_star=I(:,:)*inv(g(:,:));
%**************************************************************************
% first vector in reciprocal space, Miller indices
h1=[0; 0; 2];
% second vector in reciprocal space, Miller indices
h2=[2; 0; 2];
%**************************************************************************
% length of h1
h1_length=transpose(h1(:,1))*g_star(:,:)*h1(:,1);
h1_length=sqrt(h1_length);
% length of h2
h2_length=transpose(h2(:,1))*g_star(:,:)*h2(:,1);
h2_length=sqrt(h2_length);
%**************************************************************************
% calculate scalar multiplication between h1 and h2
scalar_multiplication=transpose(h1(:,1))*g_star(:,:)*h2(:,1);
% calculate angle between h1 and h2
angle=acosd(scalar_multiplication/(h1_length*h2_length));
%**************************************************************************
% calculate contravariant components of h1
h1_contravariant=g_star(:,:)*h1(:,1);
% calculate contravariant components of h2
h2_contravariant=g_star(:,:)*h2(:,1);
%**************************************************************************
% project vectors to (001) plane
h=[0; 0; 1];
% length of h vector
h_length=transpose(h(:,1))*g_star(:,:)*h(:,1);
h_length=sqrt(h_length);
%**************************************************************************
% project vector h1 to (001) plane
proj_vector_h1=zeros(3,1);
for index_0=1:3
    for index_1=1:3
        proj_vector_h1(index_0,1)=proj_vector_h1(index_0,1)-h(index_1,1)*h1_contravariant(index_1,1)*h(index_0,1)/h_length^2;
    end
end
proj_vector_h1(:,1)=h1(:,1)+proj_vector_h1(:,1);
%**************************************************************************
% project vector h2 to (001) plane
proj_vector_h2=zeros(3,1);
for index_0=1:3
    for index_1=1:3
        proj_vector_h2(index_0,1)=proj_vector_h2(index_0,1)-h(index_1,1)*h2_contravariant(index_1,1)*h(index_0,1)/h_length^2;
    end
end
proj_vector_h2(:,1)=h2(:,1)+proj_vector_h2(:,1);
%**************************************************************************
% length of projection of h1 to (001) plane
proj_vector_h1_length=transpose(proj_vector_h1(:,1))*g_star(:,:)*proj_vector_h1(:,1);
proj_vector_h1_length=sqrt(proj_vector_h1_length);
% length of projection of h2 to (001) plane
proj_vector_h2_length=transpose(proj_vector_h2(:,1))*g_star(:,:)*proj_vector_h2(:,1);
proj_vector_h2_length=sqrt(proj_vector_h2_length);
%**************************************************************************
% scalar multiplication between projection vectors to (001) plane
scalar_multiplication_projections=transpose(proj_vector_h1(:,1))*g_star(:,:)*proj_vector_h2(:,1);
% angle between projected vectors
angle_proj=acosd(scalar_multiplication_projections/(proj_vector_h1_length*proj_vector_h2_length));
%**************************************************************************
% end of programe
%**************************************************************************