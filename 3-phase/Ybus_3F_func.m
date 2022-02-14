function [Yss, Ysd, Yds, Ydd, Zr] = Ybus_3F_func(System_Data_Lines)
%     """Builds the admittance matrix Ybus in p.u. and the submatrices
%     needed @author: Juan S. Giraldo (UTwente) jnse@ieee.org Obs: Slack
%     bus needs to be numbered as 1 !! """
    
fn = System_Data_Lines(:,1);
tn = System_Data_Lines(:,2);
Raa = System_Data_Lines(:,3);
Xaa = System_Data_Lines(:,4);
Rbb = System_Data_Lines(:,5);
Xbb = System_Data_Lines(:,6);
Rcc = System_Data_Lines(:,7);
Xcc = System_Data_Lines(:,8);
Rab = System_Data_Lines(:,9);
Xab = System_Data_Lines(:,10);
Rac = System_Data_Lines(:,11);
Xac = System_Data_Lines(:,12);
Rbc = System_Data_Lines(:,13);
Xbc = System_Data_Lines(:,14);


nr = length(fn(:,1));
nb = nr + 1;
Y = zeros(3*nb,3*nb);


ADM = zeros(nr,length(System_Data_Lines(1,:)));
Yl = zeros(3,3);


% fileID = fopen('Y_barra_123_bus_3fPU_HR.dat','w');
% nbytes = fprintf(fileID,'param : Ybus : gaa             baa           gbb            bbb              gcc           bcc             gab              bab            gac             bac          gbc             bbc:=');

Zr = zeros(3*nr,3*nr);

% Zr = load('Zr_8_Dan').Zr;
% Zl2 = Zr(3*1-2:3*1,3*1-2:3*1)
%%

for i = 1:nr
    
    Zl = [Raa(i) + 1j*Xaa(i),Rab(i) + 1j*Xab(i),Rac(i) + 1j*Xac(i)
        Rab(i) + 1j*Xab(i),Rbb(i) + 1j*Xbb(i),Rbc(i) + 1j*Xbc(i)
        Rac(i) + 1j*Xac(i),Rbc(i) + 1j*Xbc(i),Rcc(i) + 1j*Xcc(i)];
    
    
    
%     Zr(3*i-2:3*i,3*i-2:3*i) = Zl; % Uncomment to obtain matrix for
%     Danilos' algorithms

    if Raa(i) == 0 && Xaa(i) == 0
        Zl(1,1) = 1000;
    end
    
    if Rbb(i) == 0 && Xbb(i) == 0
        Zl(2,2) = 1000;
    end
    
    if Rcc(i) == 0 && Xcc(i) == 0
        Zl(3,3) = 1000;
    end
      
    Yl = inv(Zl);
      
      
    if Raa(i) == 0 && Xaa(i) == 0
        Yl(1,:) = zeros;
        Yl(:,1) = zeros;
    end
    
    if Rbb(i) == 0 && Xbb(i) == 0
        Yl(2,:) = zeros;
        Yl(:,2) = zeros;
    end
    
    if Rcc(i) == 0 && Xcc(i) == 0
        Yl(3,:) = zeros;
        Yl(:,3) = zeros;
    end
    
    %               e      r       gaa             baa           gbb            bbb              gcc           bcc             gab              bab            gac             bac          gbc             bbc           Imax      
    ADM(i,1:end) = [fn(i), tn(i), real(Yl(1,1)), imag(Yl(1,1)), real(Yl(2,2)), imag(Yl(2,2)), real(Yl(3,3)), imag(Yl(3,3)),  real(Yl(1,2)), imag(Yl(1,2)), real(Yl(1,3)), imag(Yl(1,3)), real(Yl(2,3)), imag(Yl(2,3))];
%     nbytes = fprintf(fileID,'\n%d %d %2.10f %2.10f %2.10f %2.10f %2.10f %2.10f %2.10f %2.10f %2.10f %2.10f %2.10f %2.10f',fn(i), tn(i), real(Yl(1,1)), imag(Yl(1,1)), real(Yl(2,2)), imag(Yl(2,2)), real(Yl(3,3)), imag(Yl(3,3)),  real(Yl(1,2)), imag(Yl(1,2)), real(Yl(1,3)), imag(Yl(1,3)), real(Yl(2,3)), imag(Yl(2,3)));
end
% fprintf(fileID,'\n;');
% fclose(fileID);

gaa = ADM(:,3);
baa = ADM(:,4);
gbb = ADM(:,5);
bbb = ADM(:,6);
gcc = ADM(:,7);
bcc = ADM(:,8);
gab = ADM(:,9);
bab = ADM(:,10);
gac = ADM(:,11);
bac = ADM(:,12);
gbc = ADM(:,13);
bbc = ADM(:,14);
e = fn;
r=tn;

% Matrices Fuera Diagonal

for i = 1:1:nr
    l = e(i);
    c = r(i);
    if Raa(i) ~= 0 || Xaa(i) ~= 0
        Y(3*l-2,3*c-2) = -(gaa(i) + 1j*baa(i));
        Y(3*c-2,3*l-2) = -(gaa(i) + 1j*baa(i));
%         Y(3*c-2,3*c-2) = Y(3*c-2,3*c-2) - Y(3*l-2,3*c-2);
%         Y(3*l-2,3*l-2) = Y(3*l-2,3*l-2) - Y(3*c-2,3*l-2);
    end

    if Rbb(i)~= 0 || Xbb(i)~= 0
        Y(3*l-1,3*c-1) = -(gbb(i) + 1j*bbb(i));
        Y(3*c-1,3*l-1) =  -(gbb(i) + 1j*bbb(i));
    end

    if Rcc(i)~= 0 || Xcc(i)~= 0
        Y(3*l,3*c) =  -(gcc(i) + 1j*bcc(i));
        Y(3*c,3*l) =  -(gcc(i) + 1j*bcc(i));
    end

    if Rab(i)~= 0 || Xab(i)~= 0
        Y(3*l-2,3*c-1) = -(gab(i) + 1j*bab(i));
        Y(3*l-1,3*c-2) = -(gab(i) + 1j*bab(i));
        Y(3*c-1,3*l-2) = -(gab(i) + 1j*bab(i));
        Y(3*c-2,3*l-1) = -(gab(i) + 1j*bab(i));
    end

    if Rac(i)~= 0 || Xac(i)~= 0
        Y(3*l-2,3*c) =  -(gac(i) + 1j*bac(i));
        Y(3*l,3*c-2) =  -(gac(i) + 1j*bac(i));
        Y(3*c,3*l-2) =  -(gac(i) + 1j*bac(i));
        Y(3*c-2,3*l) =  -(gac(i) + 1j*bac(i));
    end

    if Rbc(i)~= 0 || Xbc(i)~= 0
        Y(3*l-1,3*c) =  -(gbc(i) + 1j*bbc(i));
        Y(3*l,3*c-1) =  -(gbc(i) + 1j*bbc(i));
        Y(3*c,3*l-1) =  -(gbc(i) + 1j*bbc(i));
        Y(3*c-1,3*l) =  -(gbc(i) + 1j*bbc(i));
    end
end

% Matrices Diagonal

for i = 1:1:nb
    for k = 1:1:nb
        if k ~= i
            Y(3*i-2,3*i-2) = Y(3*i-2,3*i-2) - Y(3*i-2,3*k-2);
            Y(3*i-2,3*i-1) = Y(3*i-2,3*i-1) - Y(3*i-2,3*k-1);
            Y(3*i-2,3*i) = Y(3*i-2,3*i) - Y(3*i-2,3*k);

            Y(3*i-1,3*i-2) = Y(3*i-1,3*i-2) - Y(3*i-1,3*k-2);
            Y(3*i-1,3*i-1) = Y(3*i-1,3*i-1) - Y(3*i-1,3*k-1);
            Y(3*i-1,3*i) = Y(3*i-1,3*i) - Y(3*i-1,3*k);

            Y(3*i,3*i-2) = Y(3*i,3*i-2) - Y(3*i,3*k-2);
            Y(3*i,3*i-1) = Y(3*i,3*i-1) - Y(3*i,3*k-1);
            Y(3*i,3*i) = Y(3*i,3*i) - Y(3*i,3*k);
        end
    end
end

Y = Y + 0.0001*diag(diag(Y)==0);

Yss = sparse(Y(1:3,1:3));
Ysd = sparse(Y(1:3,4:end));
Yds = sparse(Y(4:end,1:3));
Ydd = sparse(Y(4:end,4:end));
end

