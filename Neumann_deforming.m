format long
%I=imread('house256.jpg');
I=imread('cameraman256.jpg');
I=im2gray(I); %transform the color image in a grayscale image
%figure
%imshow(G) % if you want to see the grayscale image
I=double(I)/255; %we make M to be the matrix with double values in [0,1]
I=I';
Nx=length(I(:,1)); %number of pixels in a row
Ny=length(I(1,:)); %number of pixels in a column
D=15.6*25.4; % the diagonal of the monitor expressed in milimeters
dx=D/sqrt(4852800); %the length of the edge of a pixel for a monitor with aspect ratio 16:9, expressed in mm.
Nt=5000; %number of iterations with respect to time that we want to produce
dt=0.01; %we choose a small time step in order to have a good accuracy
x=1:Nx;
y=1:Ny;
[X,Y]=ndgrid(x,y);


%Diffusion Coefficient
d=0.3; %it is important to choose a good diffusion coefficient in order to obtain relevant results



%The source function
f=@(x,y,u) x.*0.*y+0.5.*u.*(1-(0.4).*sin(x./20)-(0.4).*cos(y./20)-2*u);
%f=@(x,y,u) x.*0+y.*0+sin(u).*cos(pi.*u./2);

%Initialization
u=zeros(Nx,Ny,Nt);
u(:,:,1)=I; %this is the initial data of our PDE (the grayscale image matrix)

F=zeros(Nx,Ny,Nt); 
F(:,:,1)=f(X,Y,u(x,y,1)); %this is the initialization of the source
%function when we work with mathematical function f
%F(:,:,1)=alpha.*u(:,:,1).*(r-p.*u(:,:,1));

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
        F(:,:,k+1)=f(X,Y,u(x,y,k+1));
        %F(:,:,k+1)=alpha.*u(:,:,k+1).*(r-p.*u(:,:,k+1));


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
        F(:,:,k+1)=f(X,Y,u(x,y,k+1));
        %F(:,:,k+1)=alpha.*u(:,:,k+1).*(r-p.*u(:,:,k+1));




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
        F(:,:,k+1)=f(X,Y,u(x,y,k+1));
        %F(:,:,k+1)=alpha.*u(:,:,k+1).*(r-p.*u(:,:,k+1));




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
        F(:,:,k+1)=f(X,Y,u(x,y,k+1));
        %F(:,:,k+1)=alpha.*u(:,:,k+1).*(r-p.*u(:,:,k+1));
    end

end

% for k=1:Nt
%     figure;
%     imagesc(u(:,:,k)',[0 1]), axis equal; axis off; colormap(gray)
%      title(['The diffused cameraman image after t='...
%         num2str(k) ' timesteps']);
%     pause(0.1)
% end

v = VideoWriter("neumann_diffusion_cameraman",'MPEG-4');
v.Quality=100;
v.FrameRate=10;
open(v)
set(gcf, 'renderer', 'zbuffer');

for k=1:Nt
    figure(4)
    imagesc(u(:,:,k)',[0 1]), axis equal; axis off; colormap(gray)
    title(['The diffused cameraman image after t='...
        num2str(k) ' timesteps']);
    %title(['Initial image diffused into the final image after t='...
        %num2str(k) ' timesteps']);
    pause(0.1)
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)


