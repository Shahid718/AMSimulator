% =====
% data
Nx     = 100;
Ny     = 100;
nsteps = 3000;
dx     = 20e-9;
dy     = 20e-9;
delta  = 4.0*dx;
sigma  = 0.37;
alpha  = 0.05;
aniso  = 4.0;
theta0 = 0.0;
Tm     = 1728.0;
K      = 84.01;
Cp     = 5.42e+06;
H      = 2.35e+09;
mu     = 2.0;
T0     = Tm-H/Cp*0.3;
lambda = 0.1;
b      = 2.0*atanh(1.0-2.0*lambda);
a0     = sqrt(3.0*delta*sigma/b);
w      = 6.0*sigma*b/delta;
pmob   = b*Tm*mu/(3.0*delta*H);
dt     = dx^2 /(5.0*K/Cp);

% get the start time
time0 = clock();

% ====
% Initial microstructure
T = zeros(Nx,Ny);
p = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        if ((i-Nx/2.0)^2 + (j-Ny/2.0)^2 - 3 <= 0)
            p(i,j) = 1.0;
        end
        T(i,j) = T0;
    end
end

% ====
% Evolution
for istep =1:nsteps
    for i=1:Nx
        for j=1:Ny

            % boundary conditions
            jp = j + 1;
            jm = j - 1;
            ip = i + 1;
            im = i - 1;

            if (im == 0)
                im = Nx;
            end
            if(ip == (Nx+1))
                ip = 1;
            end
            if(jm == 0)
                jm = Ny;
            end
            if(jp == (Ny+1))
                jp = 1;
            end

            % laplacian
            p_lap(i,j)=(p(ip,j)+p(im,j)+p(i,jm)+p(i,jp)-4.0*p(i,j))/(dx*dy);
            T_lap(i,j)=(T(ip,j)+T(im,j)+T(i,jm)+T(i,jp)-4.0*T(i,j))/(dx*dy);

            % gradients and angle
            pdx(i,j) = (p(ip,j)-p(im,j)) /(2.0*dx);
            pdy(i,j) = (p(i,jp)-p(i,jm)) /(2.0*dy);
            theta = atan2(pdy(i,j),pdx(i,j));

            % anisotropy and its derivative
            a(i,j) = a0*(1.0+alpha*cos(aniso*(theta-theta0)));
            da(i,j)= -a0*aniso*alpha*sin(aniso*(theta-theta0));
            term1 = (a(i,jp)*da(i,jp)*pdx(i,jp)-a(i,jm)*da(i,jm)*pdx(i,jm))/dy;
            term2 = -(a(ip,j)*da(ip,j)*pdy(ip,j)-a(im,j)*da(im,j)*pdy(im,j))/dx;

            % driving force
            f = -H*(T(i,j)-Tm)/Tm;

            % phi and temperature
            pnew(i,j) = p(i,j) + (a(i,j)^2*p_lap(i,j) + term1 + term2 + ...
                4.0*w*p(i,j)*(1.0-p(i,j))*(p(i,j)-0.5+15.0/(2.0*w)*f* ...
                p(i,j)*(1.0-p(i,j))))*dt*pmob;
            Tnew(i,j) = T(i,j) + (K/Cp)*T_lap(i,j)*dt + ...
                30.0*p(i,j)^2*(1.0-p(i,j))^2*(H/Cp)*(pnew(i,j)-p(i,j));

        end
    end
    % update
    T = Tnew;
    p = pnew;
end

compute_time = etime(clock(),time0);
fprintf('time: %f\n',compute_time);
