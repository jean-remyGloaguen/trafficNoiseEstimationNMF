function NMF = algo_nmfSemiSupervisedEXP(H,W,Y,Z,V,iteration,setting)

beta = setting.beta;
sparsity = setting.sparsity;
smoothness = setting.smoothness;
J = setting.SS_sizeWrand;
K = size(W,2);

if ~isequal(beta,0) || ~isequal(beta,1) || ~isequal(beta,2)
    if beta < 1
        gamma = 1/(2-beta);
    elseif beta >= 1 && beta <= 2
        gamma = 1;
    else
        gamma = 1/(beta-1);
    end
end

Vap = W*H;
Vtot = Vap+Y*Z;

NMF = [];

switch beta
    case 2
        for iter = 1:iteration
            Y = Y.*(((V*Z')./(Vtot*Z')).^gamma);
            Y = Y./repmat(sum(Y),size(Y,1),1);
            Vtot = Vap+Y*Z;
            
            if smoothness == 0 || strcmp(smoothForm,'traffic')
                Z = Z.*(((Y'*V)./(Y'*Vtot)).^gamma);
            elseif smoothness == 1 && strcmp(smoothForm,'all')
                Z_1 = [zeros(J,1) Z(:,1:end-1)];
                Z_2 = [Z(:,2:end) zeros(J,1)];
                Z_12 = [zeros(J,1) Z(:,2:end-1) zeros(J,1)];
                Z = Z.*(Y'*V-sparsity+2*smoothness.*(Z_1+Z_2))./(Y'*Vtot+2*smoothness.*(Z+Z_12));
            end
            Vtot = Vap+Y*Z;
            
            if smoothness == 0
                H = H.*(((W'*V-sparsity)./(W'*Vtot)).^gamma);
                H(H<0)=0;
            else
                H_1 = [zeros(K,1) H(:,1:end-1)];
                H_2 = [H(:,2:end) zeros(K,1)];
                H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                H = H.*(W'*V-sparsity+2*smoothness.*(H_1+H_2))./(W'*Vap+2*smoothness.*(H+H_12));
            end
            
            Vap = W*H;
            Vtot = Vap+Y*Z;
        end
        
    case 1
        F = size(V,1);
        N = size(H,2);
        
        for iter = 1:iteration
            Y = Y.*((((V.*(Vtot).^(-1))*Z')./repmat(sum(Z,2)',F,1)).^gamma);          
            Y = Y./repmat(sum(Y), size(Y, 1), 1);
            Vtot = Vap+Y*Z;
            
            if smoothness == 0 || strcmp(smoothForm,'traffic')
                Z = Z.*(((Y'*(V.*(Vtot).^(-1)))./repmat(sum(Y,1)',1,N)).^gamma);
            elseif smoothness == 1 && strcmp(smoothForm,'all')
                Z_1 = [zeros(J,1) Z(:,1:end-1)];
                Z_2 = [Z(:,2:end) zeros(J,1)];
                Z_12 = [zeros(J,1) Z(:,2:end-1) zeros(J,1)];
                Z = Z.*(((Y'*(Vtot.^(beta-2).*V)-sparsity+2*smoothness.*(Z_1+Z_2))./(Y'*Vtot.^(0))+2*smoothness.*(Z+Z_12)).^gamma);
            end
            Vtot = Vap+Y*Z;
            
            if smoothness == 0
                H = H.*(((W.'*(V.*(Vtot).^(-1)))./(sparsity+repmat(sum(W,1)',1,N))).^gamma);
            else
                H_1 = [zeros(K,1) H(:,1:end-1)];
                H_2 = [H(:,2:end) zeros(K,1)];
                H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                H = H.*((W'*(V.*Vap.^(-1))+2*smoothness.*(H_1+H_2))./(sparsity+W'*Vap.^(0)+2*smoothness.*(H+H_12)));
            end
            
            Vap = W*H;
            Vtot = Vap+Y*Z;
        end
    case 0
        for iter = 1:iteration
            Y = Y.*((((V.*Vtot.^(-2))*Z')./(Vtot.^(-1)*Z')).^gamma);            
            Y = Y./repmat(sum(Y), size(Y, 1), 1);
            Vtot = Vap+Y*Z;
            
            if smoothness == 0 || strcmp(smoothForm,'traffic')
                Z = Z.*((Y'*(V.*(Vtot).^(-2)))./(sparsity+Y'*Vtot.^(-1))).^gamma;
            elseif smoothness == 1 && strcmp(smoothForm,'all')
                Z_1 = [zeros(J,1) Z(:,1:end-1)];
                Z_2 = [Z(:,2:end) zeros(J,1)];
                Z_12 = [zeros(J,1) Z(:,2:end-1) zeros(J,1)];
                Z = Z.*(((Y'*(Vtot.^(-2).*V)-sparsity+2*smoothness.*(Z_1+Z_2))./(Y'*Vtot.^(-1))+2*smoothness.*(Z+Z_12)).^gamma);
            end
            Vtot = Vap+Y*Z;
            
            if smoothness == 0
                H = H.*(((W'*(Vtot.^(-2).*V))./(sparsity+W'*Vtot.^(-1))).^gamma);
            else
                H_1 = [zeros(K,1) H(:,1:end-1)];
                H_2 = [H(:,2:end) zeros(K,1)];
                H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                H = H.*((W'*(V.*Vap.^(-2))+2*smoothness.*(H_1+H_2))./(sparsity+W'*Vap.^(-1)+2*smoothness.*(H+H_12))).^(1/2);
            end
            
            Vap = W*H;
            Vtot = Vap+Y*Z;
        end
    otherwise
        for iter = 1:iteration
            Y = Y.*((((V.*Vtot.^(beta-2))*Z')./(Vtot.^(beta-1)*Z')).^gamma);
            Y = Y./repmat(sum(Y), size(Y, 1), 1);
            Vtot = Vap+Y*Z;
            
            Z = Z.*(((Y'*(Vtot.^(beta-2).*V))./(Y'*Vtot.^(beta-1))).^gamma);
            Vtot = Vap+Y*Z;
            
            if smoothness == 0
                if beta < 2
                    H = H.*(((W'*(Vtot.^(beta-2).*V))./(sparsity+W'*Vtot.^(beta-1))).^gamma);
                else
                    H = H.*(((W'*(Vtot.^(beta-2).*V)-sparsity)./(W'*Vtot.^(beta-1))).^gamma);
                end
            else
                if beta < 2
                    H_1 = [zeros(K,1) H(:,1:end-1)];
                    H_2 = [H(:,2:end) zeros(K,1)];
                    H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                    H = H.*((W'*(V.*Vap.^(beta-2))+2*smoothness.*(H_1+H_2))./(sparsity+W'*Vap.^(beta-1)+2*smoothness.*(H+H_12))).^(gamma);
                else
                    H_1 = [zeros(K,1) H(:,1:end-1)];
                    H_2 = [H(:,2:end) zeros(K,1)];
                    H_12 = [zeros(K,1) H(:,2:end-1) zeros(K,1)];
                    H = H.*((W'*(V.*Vap.^(beta-2))+2*smoothness.*(H_1+H_2)-sparsity)./(W'*Vap.^(beta-1)+2*smoothness.*(H+H_12))).^(1/2);
                end
            end
            Vap = W*H;
            Vtot = Vap+Y*Z;
            
        end
end

NMF.H = H;
NMF.W = W;
NMF.Z = Z;
NMF.Y = Y;
NMF.Vtot = Vtot;