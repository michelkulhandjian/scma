%% MPA MIMO scma decoder is for 2 transmitter (Tx) and 1 receiver (Rx)system.
%  The algorithm produces LLRs given the received vector of dimension K.
%  For two Tx system we create joint factor graph (size KxMx(2J)) by 
%  concatinating of two factor graphs F (size KxMxJ).

%   Author:
%   Michel Kulhandjian
%   Research Scientist
%   School of Electrical Engineering and Computer Science
%   University of Ottawa, Ottawa, Ontario, K1N 6N5, Canada
%   Email: mkk6@buffalo.edu or cdamours@uottawa.ca
%   Web: http://mkulhandjian.x10host.com/main/
%
%   Date: December 2017
%
%   References: 
%   1) M. Kulhandjian and C. D'Amours, "Design of Permutation-Based Sparse
%   Code Multiple Access System," accepted for publication in Proc. IEEE
%   Pers., Indoor, Mobile Radio Conf. (PIMRC) 2017, Montreal, Canada, Oct. 2017. 
%   
%   **Any feedback of comments regarding the source code are welcomed. 
%   This source code is publically shared in the hope that it will be useful
%   for researchers. If you may implement the source even partially please
%   cite the reference article above.** 

function [LLR] = mpaMIMO(Y,C, K, M, N, J, df, VN, FN, VNdf, FNn, IndRowOne, IndRowZer, N0, Nit)
%  MPA MIMO SCMA decoder (Log-MPA)
%
%  Input arguments:
%
%  Y  - received SCMA signal at Rx1 after fading channel (size Kx1)
%  C  - SCMA codebooks (size KxMx(2J)) undergone with fast fading channel 
%  K  - number of orthogonal resources (4)
%  M  - number of codewords in each codebook (4)
%  N  - number of non-zero element in the codebook (2)
%  J  - number of users (layers) (12)
%  df - is the number branches arriving to a resource node (6)
%  VN - variable node matrix (size (2J)x2), each row r of column in F
%       connected to resource nodes
%  FN - function node matrix (size Kx(2df)) each row r of row in F
%       connected to variable nodes
%  VNdf  - matrix (size (2J)x2) each row r shows r at which columns appear in FN
%  FNn  - matrix (size Kx(2df)) each row r indicates the number of connection
%         at FN(r,k) variable node 
%  IndRowOne  - position of +1 bits in all B (size Mx1)
%  IndRowZer  - position of -1 bits in all B (size Mx1)
%  N0  - variance of noise (in AWGN channel)
%  Nit - number of MPA iterations
%
%  Output arguments:
%
%  LLR - Log-Likelihood Ratio (size (log2(M)*J)x1)

EF = zeros(df, M, K);
EF2 = zeros(df, M, K);
EV = zeros(N, M, J);
FF = zeros(M, M, M, M, M, M, K);
dMf = M^(df-1);
mM = zeros(1,df);

% Step 1: Initial calculations
for k = 1:K % resourses
    for m1 = 1:M
        for m2 = 1:M
            for m3 = 1:M
                for m4 = 1:M
                    for m5 = 1:M
                        for m6 = 1:M
                            FF(m1,m2,m3,m4,m5,m6,k) = -(1/(2*N0))*abs(Y(k)-(C(k,m1,FN(k,1))+C(k,m2,FN(k,2))+C(k,m3,FN(k,3))+C(k,m4,FN(k,4))+C(k,m5,FN(k,5))+C(k,m6,FN(k,6))))^2;
                        end
                    end
                end
            end
        end
    end
end

% Step 2: Iterative procedure
for nmMPA = 1: Nit
    %% EF computaiton
    for k=1:K
        % df = 1
        for m1 = 1:M
            sIgv = zeros(1,M^5);
            for m2 = 1:M
                for m3 = 1:M
                    for m4 = 1:M
                        for m5 = 1:M
                            for m6 = 1:M
                                sIgv((m2-1)*M^4+(m3-1)*M^3+(m4-1)*M^2+(m5-1)*M+m6) = EV(FNn(k,2),m2,FN(k,2))+...
                                    EV(FNn(k,3),m3,FN(k,3))+EV(FNn(k,4),m4,FN(k,4))+EV(FNn(k,5),m5,FN(k,5))+...
                                    EV(FNn(k,6),m6,FN(k,6))+FF(m1,m2,m3,m4,m5,m6,k);
                            end
                        end
                    end
                end
            end
            EF(1,m1,k) = max(sIgv);
        end
        
        % df = 2
        for m2 = 1:M
            sIgv = zeros(1,M^5);
            for m1 = 1:M
                for m3 = 1:M
                    for m4 = 1:M
                        for m5 = 1:M
                            for m6 = 1:M
                                sIgv((m1-1)*M^4+(m3-1)*M^3+(m4-1)*M^2+(m5-1)*M+m6) = EV(FNn(k,1),m1,FN(k,1))+...
                                    EV(FNn(k,3),m3,FN(k,3))+EV(FNn(k,4),m4,FN(k,4))+EV(FNn(k,5),m5,FN(k,5))+...
                                    EV(FNn(k,6),m6,FN(k,6))+FF(m1,m2,m3,m4,m5,m6,k);
                            end
                        end
                    end
                end
            end
            EF(2,m2,k) = max(sIgv);
        end
        
        % df = 3
        for m3 = 1:M
            sIgv = zeros(1,M^5);
            for m1 = 1:M
                for m2 = 1:M
                    for m4 = 1:M
                        for m5 = 1:M
                            for m6 = 1:M
                                sIgv((m1-1)*M^4+(m2-1)*M^3+(m4-1)*M^2+(m5-1)*M+m6) = EV(FNn(k,1),m1,FN(k,1))+...
                                    EV(FNn(k,2),m2,FN(k,2))+EV(FNn(k,4),m4,FN(k,4))+EV(FNn(k,5),m5,FN(k,5))+...
                                    EV(FNn(k,6),m6,FN(k,6))+FF(m1,m2,m3,m4,m5,m6,k);
                            end
                        end
                    end
                end
            end
            EF(3,m3,k) = max(sIgv);
        end
        
        % df = 4
        for m4 = 1:M
            sIgv = zeros(1,M^5);
            for m1 = 1:M
                for m2 = 1:M
                    for m3 = 1:M
                        for m5 = 1:M
                            for m6 = 1:M
                                sIgv((m1-1)*M^4+(m2-1)*M^3+(m3-1)*M^2+(m5-1)*M+m6) = EV(FNn(k,1),m1,FN(k,1))+...
                                    EV(FNn(k,2),m2,FN(k,2))+EV(FNn(k,3),m3,FN(k,3))+EV(FNn(k,5),m5,FN(k,5))+...
                                    EV(FNn(k,6),m6,FN(k,6))+FF(m1,m2,m3,m4,m5,m6,k);
                            end
                        end
                    end
                end
            end
            EF(4,m4,k) = max(sIgv);
        end
        
        % df = 5
        for m5 = 1:M
            sIgv = zeros(1,M^5);
            for m1 = 1:M
                for m2 = 1:M
                    for m3 = 1:M
                        for m4 = 1:M
                            for m6 = 1:M
                                sIgv((m1-1)*M^4+(m2-1)*M^3+(m3-1)*M^2+(m4-1)*M+m6) = EV(FNn(k,1),m1,FN(k,1))+...
                                    EV(FNn(k,2),m2,FN(k,2))+EV(FNn(k,3),m3,FN(k,3))+EV(FNn(k,4),m4,FN(k,4))+...
                                    EV(FNn(k,6),m6,FN(k,6))+FF(m1,m2,m3,m4,m5,m6,k);
                            end
                        end
                    end
                end
            end
            EF(5,m5,k) = max(sIgv);
        end
        
        % df = 6
        for m6 = 1:M
            sIgv = zeros(1,M^5);
            for m1 = 1:M
                for m2 = 1:M
                    for m3 = 1:M
                        for m4 = 1:M
                            for m5 = 1:M
                                sIgv((m1-1)*M^4+(m2-1)*M^3+(m3-1)*M^2+(m4-1)*M+m5) = EV(FNn(k,1),m1,FN(k,1))+...
                                    EV(FNn(k,2),m2,FN(k,2))+EV(FNn(k,3),m3,FN(k,3))+EV(FNn(k,4),m4,FN(k,4))+...
                                    EV(FNn(k,5),m5,FN(k,5))+FF(m1,m2,m3,m4,m5,m6,k);
                            end
                        end
                    end
                end
            end
            EF(6,m6,k) = max(sIgv);
        end     
    end
    
    %% EV computaiton
    for j= 1: J
        for nk = 1 : N
            for m = 1 : M
                for ni = 1 : N
                    if (ni ~= nk)
                        EV(nk, m, j) =  EF(VNdf(j,ni), m, VN(j,ni));
                    end
                end
            end
        end
    end
end
IV = zeros(J,M);
Nbits = log2(M);
B = zeros(Nbits,J);
ym = zeros(1,J);
llr =[];
for j = 1 : J
    IV(j,:) = sum( EV(:, :, j), 1);
    for bi = Nbits:-1:1
        llr = [llr ; max(IV(j,IndRowOne((bi-1)*M/2+1:bi*M/2)))-max(IV(j,IndRowZer((bi-1)*M/2+1:bi*M/2))) ];
    end
end

LLR =[];
 for j = 1 : J/2
     LLL =[LLR; [llr((j-1)*2+1); llr(j*2); llr((j-1)*2+1+J); llr(j*2+J)]; ];
 end