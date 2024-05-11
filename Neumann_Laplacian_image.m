format long
I=imread('cameraman256.jpg'); %the initial image
O=imread('house256.jpg'); %the final image
%I=I(1:259,1:194,1:3);
%imshow(I) %to view the RGB image
I=im2gray(I); %transform the color image in a grayscale image
O=im2gray(O);
%figure
%imshow(G) % if you want to see the grayscale image
I=double(I)/255; %we make M to be the matrix with double values in [0,1]
O=double(O)/255;
I=I';
U=O';
U=0.8*U+0.2; %to have U far from black which gives singular points of p;
%figure(1)
%imshow(U')
Nx=length(I(:,1)); %number of pixels in a row
Ny=length(I(1,:)); %number of pixels in a column
D=15.6*25.4; % the diagonal of the monitor expressed in milimeters
dx=D/sqrt(4852800); %the length of the edge of a pixel for a monitor with aspect ratio 16:9, expressed in mm.
Nt=200; %number of iterations with respect to time that we want to produce
dt=0.01; %we choose a small time step in order to have a good accuracy
x=1:Nx;
y=1:Ny;
[X,Y]=ndgrid(x,y);


%Diffusion Coefficient
d=0.1; %it is important to choose a good diffusion coefficient in order to obtain relevant results
alpha=5; %calibration constant
Delta_U=zeros(Nx,Ny); % computing the laplacean of U in order to define p
for i=3:(Nx-2)
    for j=3:(Ny-2)
        Delta_U(i,j)=(1/(12*(dx)^2))*(-U(i+2,j)-U(i-2,j)-U(i,j+2)-U(i,j-2)+16*U(i+1,j)+16*U(i-1,j)+16*U(i,j+1)+16*U(i,j-1)-60*U(i,j));
    end
end
for i=[2,(Nx-1)]
    for j=2:(Ny-1)
        Delta_U(i,j)=(1/(dx)^2)*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1)-4*U(i,j));
    end
end

for j=[2,(Ny-1)]
     for i=3:(Nx-2)
         Delta_U(i,j)=(1/(dx)^2)*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1)-4*U(i,j));
     end
end

for i=2:(Nx-1) %setting the values of Delta_r on Gamma1
    Delta_U(i,1)=(1/(dx)^2)*(U(i+1,1)+U(i-1,1)+2*U(i,2)-4*U(i,1));
end

for i=2:(Nx-1) %setting the values of Delta_r on Gamma3
    Delta_U(i,Ny)=(1/(dx)^2)*(U(i+1,Ny)+U(i-1,Ny)+2*U(i,Ny-1)-4*U(i,1));
end

for j=2:(Ny-1) %setting the values of Delta_r on Gamma4
    Delta_U(1,j)=(1/(dx)^2)*(U(1,j+1)+U(1,j-1)+2*U(2,j)-4*U(1,j));
end

for j=2:(Ny-1) %setting the values of Delta_r on Gamma2
    Delta_U(Nx,j)=(1/(dx)^2)*(U(Nx,j+1)+U(Nx,j-1)+2*U(Nx-1,j)-4*U(Nx,j));
end

%Setting the values of Delta_r on the 4 corners of the rectangle

Delta_U(1,1)=(1/(dx)^2)*(2*U(1,2)+2*U(2,1)-4*U(1,1));
Delta_U(Nx,1)=(1/(dx)^2)*(2*U(Nx-1,1)+2*U(Nx,2)-4*U(Nx,1));
Delta_U(1,Ny)=(1/(dx)^2)*(2*U(1,Ny-1)+2*U(2,Ny)-4*U(1,Ny));
Delta_U(Nx,Ny)=(1/(dx)^2)*(2*U(Nx-1,Ny)+2*U(Nx,Ny-1)-4*U(Nx,Ny));

%Defining the diffusion effect r that we want

%r=ones(Nx,Ny);

%r=2.5-sin(X./10)-(cos((Y./10))); %BEST OF ALL sinusoidal effect
%r=sqrt(4-(3-sqrt((X./100).^2+(Y./100).^2)));
r=1.5-sin(X.^2+Y.^2); % fancy
%r=1.5-sin(X./5).*cos(Y./5);%dotted effect
%r=0.5+exp(-(X./30-1).^2-(Y./20-1).^2); % STABLE
%r=0.1+exp(-0.1.*(X./10-10).^2-0.1*(Y./10-5).^2);
%r=exp(-0.001.*(X-1).^2-0.001*(Y-1).^2)+0.5 %STABLE
%r=1./(sin(X./10)+cos(Y./20)+3.01); %STABLE
%r=0.5+heavi(2000-(X-120).^2-(Y-120).^2);
%r=1.75-sin(X./5); %Vertical effect
%r=1.75-sin(Y./5); %Horizontal effect
%r=0.2+abs(sin(X-Y./30));
%r=sin(X/5+Y/5)+1.75; %Diagonal effect


%Defining the corespondent p that will give the wanted image U;

%p=(r./U)+(d/alpha)*((Delta_U)./(U.^2)); %the theoretical p 
%p=(r./U)+(d/alpha)*((Delta_U)./(U.^2+10)); %adding something to U.^2 does the thing - very stable
p=r./U; %the practical p gives better results

%Initialization
u=zeros(Nx,Ny,Nt);
u(:,:,1)=I; %this is the initial data of our PDE (the grayscale image matrix)

F=zeros(Nx,Ny,Nt); 
%F(:,:,1)=f(X,Y,u(x,y,1)); %this is the initialization of the source
%function when we work with mathematical function f
F(:,:,1)=alpha.*u(:,:,1).*(r-p.*u(:,:,1));

%Starting the iterations
for k=1:(Nt-1)
    if k==1
        for i=3:(Nx-2)
        for j=3:(Ny-2)
                %here we set the core of our grid with a centered 5-point
                %formula for the second derivative (or 9-point formula for the
                %laplacian)
                u(i,j,k+1)=u(i,j,k)+dt.*((d./(12.*((dx)^2))).*(-u(i+2,j,k)+16.*u(i+1,j,k)-30.*u(i,j,k)+16*u(i-1,j,k)-u(i-2,j,k))+(d./(12.*((dx)^2))).*(-u(i,j+2,k)+16.*u(i,j+1,k)-30.*u(i,j,k)+16*u(i,j-1,k)-u(i,j-2,k))+F(i,j,k));
        end
        end

        for i=[2,(Nx-1)]
        for j=2:(Ny-1)
        %here we set the vertical shell of our core using a centered
        %3-point formula for the second derivative (or 5-point formula for
        %the laplacian)
        u(i,j,k+1)=u(i,j,k)+dt.*((d./((dx).^2)).*(u(i-1,j,k)-2.*u(i,j,k)+u(i+1,j,k))+(d./((dx).^2)).*(u(i,j-1,k)-2.*u(i,j,k)+u(i,j+1,k))+F(i,j,k));
        end
        end

        for j=[2,(Ny-1)]
        for i=3:(Nx-2) %for i=2 and i=Nx1-1 we previously set their values
        %here we set the horizontal shell of our core using the same
        %centered 3-point formula for the second derivative
        u(i,j,k+1)=u(i,j,k)+dt.*((d./((dx).^2)).*(u(i-1,j,k)-2.*u(i,j,k)+u(i+1,j,k))+(d./((dx).^2)).*(u(i,j-1,k)-2.*u(i,j,k)+u(i,j+1,k))+F(i,j,k));
        end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %Now we start to assign values on the boundary, but without the 4
       %corners of the rectangle 
       for i=2:(Nx-1)
        %on Gamma1 - using Forward Five-point Endpoint Formula
        u(i,1,k+1)=(0.04).*(48.*u(i,2,k+1)-36.*u(i,3,k+1)+16.*u(i,4,k+1)-3.*u(i,5,k+1));
        %on Gamma3 - using Backward Five-point Formula
        u(i,Ny,k+1)=(0.04).*(48.*u(i,Ny-1,k+1)-36.*u(i,Ny-2,k+1)+16.*u(i,Ny-3,k+1)-3.*u(i,Ny-4,k+1));
       end

       for j=2:(Ny-1)
        %on Gamma4 - using Forward Five-point Endpoint Formula
        u(1,j,k+1)=(0.04).*(48.*u(2,j,k+1)-36.*u(3,j,k+1)+16.*u(4,j,k+1)-3.*u(5,j,k+1));
        %on Gamma 2 - using Backward Five-point Formula
        u(Nx,j,k+1)=(0.04).*(48.*u(Nx-1,j,k+1)-36.*u(Nx-2,j,k+1)+16.*u(Nx-3,j,k+1)-3.*u(Nx-4,j,k+1));
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %Setting the values on the 4 corners of the rectangle (i,j)=(1,1),(1,Nx2),(Nx1,1),(Nx1,Nx2)
       
       %(i,j)=(1,1)-we use the arithmetic mean of the values that were
       %obtained on the sides Gamma1 and Gamma4 in (1,1)
        u(1,1,k+1)=(0.02).*((48.*u(1,2,k+1)-36.*u(1,3,k+1)+16.*u(1,4,k+1)-3.*u(1,5,k+1))+(48.*u(2,1,k+1)-36.*u(3,1,k+1)+16.*u(4,1,k+1)-3.*u(5,1,k+1)));
       
       %(i,j)=(1,Ny)- we make the mean between the values obtained on Gamma3 and Gamma4
        u(1,Ny,k+1)=(0.02).*((48.*u(1,Ny-1,k+1)-36.*u(1,Ny-2,k+1)+16.*u(1,Ny-3,k+1)-3.*u(1,Ny-4,k+1))+(48.*u(2,Ny,k+1)-36.*u(3,Ny,k+1)+16.*u(4,Ny,k+1)-3.*u(5,Ny,k+1)));
        
       %(i,j)=(Nx,1)- we make the mean between the values obtained on
       %Gamma1 and Gamma2
        u(Nx,1,k+1)=(0.02).*((48.*u(Nx,2,k+1)-36.*u(Nx,3,k+1)+16.*u(Nx,4,k+1)-3.*u(Nx,5,k+1))+(48.*u(Nx-1,1,k+1)-36.*u(Nx-2,1,k+1)+16.*u(Nx-3,1,k+1)-3.*u(Nx-4,1,k+1)));
       
       %(i,j)=(Nx,Ny)- we make the mean between the values obtained on
       %Gamma2 and Gamma3
        u(Nx,Ny,k+1)=(0.02).*((48.*u(Nx-1,Ny,k+1)-36.*u(Nx-2,Ny,k+1)+16.*u(Nx-3,Ny,k+1)-3.*u(Nx-4,Ny,k+1))+(48.*u(Nx,Ny-1,k+1)-36.*u(Nx,Ny-2,k+1)+16.*u(Nx,Ny-3,k+1)-3.*u(Nx,Ny-4,k+1)));
    
        %NOW WE SET THE VALUES OF THE SOURCE FOR k+1=2
        %F(:,:,k+1)=f(X,Y,u(x,y,k+1));
        F(:,:,k+1)=alpha.*u(:,:,k+1).*(r-p.*u(:,:,k+1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    elseif k==2
        for i=3:(Nx-2)
        for j=3:(Ny-2)
        %here we set the core of our grid
        u(i,j,k+1)=(4./3).*u(i,j,k)-(1./3).*u(i,j,k-1)+(2./3).*dt.*((d./(12.*((dx)^2))).*(-u(i+2,j,k)+16.*u(i+1,j,k)-30.*u(i,j,k)+16*u(i-1,j,k)-u(i-2,j,k))+(d./(12.*((dx)^2))).*(-u(i,j+2,k)+16.*u(i,j+1,k)-30.*u(i,j,k)+16*u(i,j-1,k)-u(i,j-2,k))+F(i,j,k));
        end
        end

        for i=[2,(Nx-1)]
        for j=2:(Ny-1)
        %here we set the shell of our core using a quadratic formula for the laplacian
        u(i,j,k+1)=(4./3).*u(i,j,k)-(1./3).*u(i,j,k-1)+(2./3).*dt.*((d./((dx).^2)).*(u(i-1,j,k)-2.*u(i,j,k)+u(i+1,j,k))+(d./((dx).^2)).*(u(i,j-1,k)-2.*u(i,j,k)+u(i,j+1,k))+F(i,j,k));
        end
        end

        for j=[2,(Ny-1)]
        for i=3:(Nx-2) %for i=2 and i=Nx-1 we previously set their values
        %here we set the shell of our core using a quadratic formula for the laplacian
        u(i,j,k+1)=(4./3).*u(i,j,k)-(1./3).*u(i,j,k-1)+(2./3).*dt.*((d./((dx).^2)).*(u(i-1,j,k)-2.*u(i,j,k)+u(i+1,j,k))+(d./((dx).^2)).*(u(i,j-1,k)-2.*u(i,j,k)+u(i,j+1,k))+F(i,j,k));
        end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Now we start to assign values on the boundary, but without the 4
       %corners of the rectangle 
       for i=2:(Nx-1)
        %on Gamma1 - using Forward Five-point Endpoint Formula
        u(i,1,k+1)=(0.04).*(48.*u(i,2,k+1)-36.*u(i,3,k+1)+16.*u(i,4,k+1)-3.*u(i,5,k+1));
        %on Gamma3 - using Backward Five-point Formula
        u(i,Ny,k+1)=(0.04).*(48.*u(i,Ny-1,k+1)-36.*u(i,Ny-2,k+1)+16.*u(i,Ny-3,k+1)-3.*u(i,Ny-4,k+1));
       end

       for j=2:(Ny-1)
        %on Gamma4 - using Forward Five-point Endpoint Formula
        u(1,j,k+1)=(0.04).*(48.*u(2,j,k+1)-36.*u(3,j,k+1)+16.*u(4,j,k+1)-3.*u(5,j,k+1));
        %on Gamma 2 - using Backward Five-point Formula
        u(Nx,j,k+1)=(0.04).*(48.*u(Nx-1,j,k+1)-36.*u(Nx-2,j,k+1)+16.*u(Nx-3,j,k+1)-3.*u(Nx-4,j,k+1));
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %Setting the values on the 4 corners of the rectangle (i,j)=(1,1),(1,Nx2),(Nx1,1),(Nx1,Nx2)
       
       %(i,j)=(1,1)-we use the arithmetic mean of the values that were
       %obtained on the sides Gamma1 and Gamma4 in (1,1)
        u(1,1,k+1)=(0.02).*((48.*u(1,2,k+1)-36.*u(1,3,k+1)+16.*u(1,4,k+1)-3.*u(1,5,k+1))+(48.*u(2,1,k+1)-36.*u(3,1,k+1)+16.*u(4,1,k+1)-3.*u(5,1,k+1)));
       
       %(i,j)=(1,Ny)- we make the mean between the values obtained on Gamma3 and Gamma4
        u(1,Ny,k+1)=(0.02).*((48.*u(1,Ny-1,k+1)-36.*u(1,Ny-2,k+1)+16.*u(1,Ny-3,k+1)-3.*u(1,Ny-4,k+1))+(48.*u(2,Ny,k+1)-36.*u(3,Ny,k+1)+16.*u(4,Ny,k+1)-3.*u(5,Ny,k+1)));
        
       %(i,j)=(Nx,1)- we make the mean between the values obtained on
       %Gamma1 and Gamma2
        u(Nx,1,k+1)=(0.02).*((48.*u(Nx,2,k+1)-36.*u(Nx,3,k+1)+16.*u(Nx,4,k+1)-3.*u(Nx,5,k+1))+(48.*u(Nx-1,1,k+1)-36.*u(Nx-2,1,k+1)+16.*u(Nx-3,1,k+1)-3.*u(Nx-4,1,k+1)));
       
       %(i,j)=(Nx,Ny)- we make the mean between the values obtained on
       %Gamma2 and Gamma3
        u(Nx,Ny,k+1)=(0.02).*((48.*u(Nx-1,Ny,k+1)-36.*u(Nx-2,Ny,k+1)+16.*u(Nx-3,Ny,k+1)-3.*u(Nx-4,Ny,k+1))+(48.*u(Nx,Ny-1,k+1)-36.*u(Nx,Ny-2,k+1)+16.*u(Nx,Ny-3,k+1)-3.*u(Nx,Ny-4,k+1)));
    

        %NOW WE SET THE VALUES OF THE SOURCE FOR k+1=3
        %F(:,:,k+1)=f(X,Y,u(x,y,k+1));
        F(:,:,k+1)=alpha.*u(:,:,k+1).*(r-p.*u(:,:,k+1));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    elseif k==3 %we apply here the backward formula with 4 nodes for the time derivative f'(x3)=1/6h (-2f(x0)+9f(x1)-18f(x2)+11f(x3))
        
       for i=3:(Nx-2)
       for j=3:(Ny-2)
       %here we set the core of our grid
        u(i,j,k+1)=(18./11).*u(i,j,k)-(9./11).*u(i,j,k-1)+(2./11).*u(i,j,k-2)+(6./11).*dt.*((d./(12.*((dx)^2))).*(-u(i+2,j,k)+16.*u(i+1,j,k)-30.*u(i,j,k)+16*u(i-1,j,k)-u(i-2,j,k))+(d./(12.*((dx)^2))).*(-u(i,j+2,k)+16.*u(i,j+1,k)-30.*u(i,j,k)+16*u(i,j-1,k)-u(i,j-2,k))+F(i,j,k));
       end
       end

       for i=[2,(Nx-1)]
       for j=2:(Ny-1)
       %here we set the shell of our core using a quadratic formula for the laplacian
        u(i,j,k+1)=(18./11).*u(i,j,k)-(9./11).*u(i,j,k-1)+(2./11).*u(i,j,k-2)+(6./11).*dt.*((d./((dx).^2)).*(u(i-1,j,k)-2.*u(i,j,k)+u(i+1,j,k))+(d./((dx).^2)).*(u(i,j-1,k)-2.*u(i,j,k)+u(i,j+1,k))+F(i,j,k));
       end
       end

       for j=[2,(Ny-1)]
       for i=3:(Nx-2) %for i=2 and i=Nx-1 we previously set their values
        %here we set the shell of our core using a quadratic formula for the laplacian
        u(i,j,k+1)=(18./11).*u(i,j,k)-(9./11).*u(i,j,k-1)+(2./11).*u(i,j,k-2)+(6./11).*dt.*((d./((dx).^2)).*(u(i-1,j,k)-2.*u(i,j,k)+u(i+1,j,k))+(d./((dx).^2)).*(u(i,j-1,k)-2.*u(i,j,k)+u(i,j+1,k))+F(i,j,k));
       end
       end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now we start to assign values on the boundary, but without the 4
       %corners of the rectangle 
       for i=2:(Nx-1)
        %on Gamma1 - using Forward Five-point Endpoint Formula
        u(i,1,k+1)=(0.04).*(48.*u(i,2,k+1)-36.*u(i,3,k+1)+16.*u(i,4,k+1)-3.*u(i,5,k+1));
        %on Gamma3 - using Backward Five-point Formula
        u(i,Ny,k+1)=(0.04).*(48.*u(i,Ny-1,k+1)-36.*u(i,Ny-2,k+1)+16.*u(i,Ny-3,k+1)-3.*u(i,Ny-4,k+1));
       end

       for j=2:(Ny-1)
        %on Gamma4 - using Forward Five-point Endpoint Formula
        u(1,j,k+1)=(0.04).*(48.*u(2,j,k+1)-36.*u(3,j,k+1)+16.*u(4,j,k+1)-3.*u(5,j,k+1));
        %on Gamma 2 - using Backward Five-point Formula
        u(Nx,j,k+1)=(0.04).*(48.*u(Nx-1,j,k+1)-36.*u(Nx-2,j,k+1)+16.*u(Nx-3,j,k+1)-3.*u(Nx-4,j,k+1));
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %Setting the values on the 4 corners of the rectangle (i,j)=(1,1),(1,Nx2),(Nx1,1),(Nx1,Nx2)
       
       %(i,j)=(1,1)-we use the arithmetic mean of the values that were
       %obtained on the sides Gamma1 and Gamma4 in (1,1)
        u(1,1,k+1)=(0.02).*((48.*u(1,2,k+1)-36.*u(1,3,k+1)+16.*u(1,4,k+1)-3.*u(1,5,k+1))+(48.*u(2,1,k+1)-36.*u(3,1,k+1)+16.*u(4,1,k+1)-3.*u(5,1,k+1)));
       
       %(i,j)=(1,Ny)- we make the mean between the values obtained on Gamma3 and Gamma4
        u(1,Ny,k+1)=(0.02).*((48.*u(1,Ny-1,k+1)-36.*u(1,Ny-2,k+1)+16.*u(1,Ny-3,k+1)-3.*u(1,Ny-4,k+1))+(48.*u(2,Ny,k+1)-36.*u(3,Ny,k+1)+16.*u(4,Ny,k+1)-3.*u(5,Ny,k+1)));
        
       %(i,j)=(Nx,1)- we make the mean between the values obtained on
       %Gamma1 and Gamma2
        u(Nx,1,k+1)=(0.02).*((48.*u(Nx,2,k+1)-36.*u(Nx,3,k+1)+16.*u(Nx,4,k+1)-3.*u(Nx,5,k+1))+(48.*u(Nx-1,1,k+1)-36.*u(Nx-2,1,k+1)+16.*u(Nx-3,1,k+1)-3.*u(Nx-4,1,k+1)));
       
       %(i,j)=(Nx,Ny)- we make the mean between the values obtained on
       %Gamma2 and Gamma3
        u(Nx,Ny,k+1)=(0.02).*((48.*u(Nx-1,Ny,k+1)-36.*u(Nx-2,Ny,k+1)+16.*u(Nx-3,Ny,k+1)-3.*u(Nx-4,Ny,k+1))+(48.*u(Nx,Ny-1,k+1)-36.*u(Nx,Ny-2,k+1)+16.*u(Nx,Ny-3,k+1)-3.*u(Nx,Ny-4,k+1)));
    
        %NOW WE SET THE VALUES OF THE SOURCE FOR k+1=4
        %F(:,:,k+1)=f(X,Y,u(x,y,k+1));
        F(:,:,k+1)=alpha.*u(:,:,k+1).*(r-p.*u(:,:,k+1));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    else %meaning that k>=4. Here we use a 5 nodes approximation formula for the time derivative
       
        for i=3:(Nx-2)
        for j=3:(Ny-2)
        %here we set the core of our grid using a fourth order
        %approximation for the laplacian
        u(i,j,k+1)=(48./25).*u(i,j,k)-(36./25).*u(i,j,k-1)+(16./25).*u(i,j,k-2)-(3./25).*u(i,j,k-3)+(12./25).*dt.*((d./(12.*((dx)^2))).*(-u(i+2,j,k)+16.*u(i+1,j,k)-30.*u(i,j,k)+16*u(i-1,j,k)-u(i-2,j,k))+(d./(12.*((dx)^2))).*(-u(i,j+2,k)+16.*u(i,j+1,k)-30.*u(i,j,k)+16*u(i,j-1,k)-u(i,j-2,k))+F(i,j,k));
        end
        end

        for i=[2,(Nx-1)]
        for j=2:(Ny-1)
        %here we set the shell of our core using a quadratic formula for the laplacian
        u(i,j,k+1)=(48./25).*u(i,j,k)-(36./25).*u(i,j,k-1)+(16./25).*u(i,j,k-2)-(3./25).*u(i,j,k-3)+(12./25).*dt.*((d./((dx).^2)).*(u(i-1,j,k)-2.*u(i,j,k)+u(i+1,j,k))+(d./((dx).^2)).*(u(i,j-1,k)-2.*u(i,j,k)+u(i,j+1,k))+F(i,j,k));
        end
        end

        for j=[2,(Ny-1)]
        for i=3:(Nx-2) %for i=2 and i=Nx1-1 we previously set their values
        %here we set the shell of our core using a quadratic formula for the laplacian
        u(i,j,k+1)=(48./25).*u(i,j,k)-(36./25).*u(i,j,k-1)+(16./25).*u(i,j,k-2)-(3./25).*u(i,j,k-3)+(12./25).*dt.*((d./((dx).^2)).*(u(i-1,j,k)-2.*u(i,j,k)+u(i+1,j,k))+(d./((dx).^2)).*(u(i,j-1,k)-2.*u(i,j,k)+u(i,j+1,k))+F(i,j,k));
        end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now we start to assign values on the boundary, but without the 4
       %corners of the rectangle 
       for i=2:(Nx-1)
        %on Gamma1 - using Forward Five-point Endpoint Formula
        u(i,1,k+1)=(0.04).*(48.*u(i,2,k+1)-36.*u(i,3,k+1)+16.*u(i,4,k+1)-3.*u(i,5,k+1));
        %on Gamma3 - using Backward Five-point Formula
        u(i,Ny,k+1)=(0.04).*(48.*u(i,Ny-1,k+1)-36.*u(i,Ny-2,k+1)+16.*u(i,Ny-3,k+1)-3.*u(i,Ny-4,k+1));
       end

       for j=2:(Ny-1)
        %on Gamma4 - using Forward Five-point Endpoint Formula
        u(1,j,k+1)=(0.04).*(48.*u(2,j,k+1)-36.*u(3,j,k+1)+16.*u(4,j,k+1)-3.*u(5,j,k+1));
        %on Gamma 2 - using Backward Five-point Formula
        u(Nx,j,k+1)=(0.04).*(48.*u(Nx-1,j,k+1)-36.*u(Nx-2,j,k+1)+16.*u(Nx-3,j,k+1)-3.*u(Nx-4,j,k+1));
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %Setting the values on the 4 corners of the rectangle (i,j)=(1,1),(1,Nx2),(Nx1,1),(Nx1,Nx2)
       
       %(i,j)=(1,1)-we use the arithmetic mean of the values that were
       %obtained on the sides Gamma1 and Gamma4 in (1,1)
        u(1,1,k+1)=(0.02).*((48.*u(1,2,k+1)-36.*u(1,3,k+1)+16.*u(1,4,k+1)-3.*u(1,5,k+1))+(48.*u(2,1,k+1)-36.*u(3,1,k+1)+16.*u(4,1,k+1)-3.*u(5,1,k+1)));
       
       %(i,j)=(1,Ny)- we make the mean between the values obtained on Gamma3 and Gamma4
        u(1,Ny,k+1)=(0.02).*((48.*u(1,Ny-1,k+1)-36.*u(1,Ny-2,k+1)+16.*u(1,Ny-3,k+1)-3.*u(1,Ny-4,k+1))+(48.*u(2,Ny,k+1)-36.*u(3,Ny,k+1)+16.*u(4,Ny,k+1)-3.*u(5,Ny,k+1)));
        
       %(i,j)=(Nx,1)- we make the mean between the values obtained on
       %Gamma1 and Gamma2
        u(Nx,1,k+1)=(0.02).*((48.*u(Nx,2,k+1)-36.*u(Nx,3,k+1)+16.*u(Nx,4,k+1)-3.*u(Nx,5,k+1))+(48.*u(Nx-1,1,k+1)-36.*u(Nx-2,1,k+1)+16.*u(Nx-3,1,k+1)-3.*u(Nx-4,1,k+1)));
       
       %(i,j)=(Nx,Ny)- we make the mean between the values obtained on
       %Gamma2 and Gamma3
        u(Nx,Ny,k+1)=(0.02).*((48.*u(Nx-1,Ny,k+1)-36.*u(Nx-2,Ny,k+1)+16.*u(Nx-3,Ny,k+1)-3.*u(Nx-4,Ny,k+1))+(48.*u(Nx,Ny-1,k+1)-36.*u(Nx,Ny-2,k+1)+16.*u(Nx,Ny-3,k+1)-3.*u(Nx,Ny-4,k+1)));
    
        %NOW WE SET THE VALUES OF THE SOURCE FOR k+1
        %F(:,:,k+1)=f(X,Y,u(x,y,k+1));
        F(:,:,k+1)=alpha.*u(:,:,k+1).*(r-p.*u(:,:,k+1));
    end

end

%for k=1:Nt
%for k=Nt
for k=[1 10 30 70 100 200]
    figure
    imagesc(u(:,:,k)',[0 1]), axis equal; axis off; colormap(gray)
    %title(['The image that will disappear seen after t='...
        % num2str(k) ' timesteps']);
    %title(['Initial image diffused into the final image after t='...
     %   num2str(k) ' timesteps']);
    pause(0.1)
end

% 
% P=0;
% No=0;
% for k=1:Nt
%     if psnr(U,u(:,:,k))>P
%     P=psnr(U,u(:,:,k));
%     No=k;
%     end
% end
% 
% G=0;
% [Gx_U,Gy_U]=imgradientxy(U);
% for k=1:Nt
%     [Gx,Gy]=imgradientxy(u(:,:,k));
%     if 0.5.*(psnr(Gx_U,Gx)+psnr(Gy_U,Gy))>G
%         G=0.5.*(psnr(Gx_U,Gx)+psnr(Gy_U,Gy));
%         Nr=k;
%     end
% end
% 
% [Gx_ufinal,Gy_ufinal]=imgradientxy(u(:,:,Nr));
% 
% disp(['PSNR_maximum value at ',num2str(No),' timesteps'])
% psnr(U,u(:,:,Nr))
% disp(['PSNRgrad maximum value at ',num2str(Nr),' timesteps'])
% 0.5.*(psnr(Gx_U,Gx_ufinal)+psnr(Gy_U,Gy_ufinal))
% NOISE_u_psnr=1.4826.*median(median((abs(u(:,:,No)-median(median((u(:,:,No))))))))./sqrt(2)
% NOISE_u_psnrgrad=1.4826.*median(median((abs(u(:,:,Nr)-median(median((u(:,:,Nr))))))))./sqrt(2)
% figure
% imagesc(u(:,:,Nr)',[0 1]), axis equal; axis off; colormap(gray);

% v = VideoWriter("image_to_image",'MPEG-4');
% v.Quality=100;
% v.FrameRate=10;
% open(v)
% set(gcf, 'renderer', 'zbuffer');
% 
% for k=1:Nt
%     figure(4)
%     imagesc(u(:,:,k)',[0 1]), axis equal; axis off; colormap(gray)
%    title(['Initial image diffused into the final image after t='...
%      num2str(k) ' timesteps']);
%     pause(0.1)
%     frame = getframe(gcf);
%     writeVideo(v,frame);
% end
% close(v)






