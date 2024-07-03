%%
classdef fns_get_recurPara
    methods (Static)
        %%
        function [Rmat,AuBu,AuBubar]=computeRplusAuBuMat(omega, k_vect, o_I,...
                alpha_vect, beta_vect, h_vect, mu_vect,...
                N, ha_inv_2by2,J,F_S,F_R)
            alpha_N2 = alpha_vect(N+2);
            beta_N2 = beta_vect(N+2);
            mu_N2 = mu_vect(N+2);
            [~, Pmat1] = fns_get_recurPara.construct_EPmat(alpha_N2,...
                beta_N2, mu_N2,omega, o_I, k_vect,Inf);
            Rmat1=zeros(2,2,length(k_vect));
            Emat1=zeros(2,2,length(k_vect));
            AuBu1=zeros(2,1,length(k_vect));
            AuBubar1=zeros(2,1,length(k_vect));
            FSFRVect=zeros(2,1,length(k_vect));
            FSFRVect(1,1,:)=k_vect.*F_S;
            FSFRVect(2,1,:)=k_vect.*F_R;
            % Rtil_mat1=zeros(2,2,length(k_vect));

            % for loop runs from (second to last: N+1) to (first layer)
            for i_L=N+1:-1:1
                alpha=alpha_vect(i_L);
                beta=beta_vect(i_L);
                mu=mu_vect(i_L);
                h=h_vect(i_L);
                [Emat, Pmat] = fns_get_recurPara.construct_EPmat(alpha,...
                    beta, mu, omega, o_I, k_vect,h);
                Lmat=fns_RmatGen.SliceMultiSolver(Pmat,Pmat1);
                L11mat = Lmat(1:2, 1:2, :);
                L12mat = Lmat(1:2, 3:4, :);
                L21mat = Lmat(3:4, 1:2, :);
                L22mat = Lmat(3:4, 3:4, :);

                ER_mat1=fns_RmatGen.SliceMultiProd(Emat1,Rmat1);
                Rtil_Part1=fns_RmatGen.SliceMultiProd(L11mat,ER_mat1)+L12mat;
                Rtil_Part2=fns_RmatGen.SliceMultiProd(L21mat,ER_mat1)+L22mat;
                Rtil_Part2INV=ha_inv_2by2(Rtil_Part2,length(k_vect));
                Rtil_mat=fns_RmatGen.SliceMultiProd(Rtil_Part1,Rtil_Part2INV);
                Rmat=fns_RmatGen.SliceMultiProd(Rtil_mat,Emat);
                %-----------------------------------%
                AuBu_Part11=L11mat-fns_RmatGen.SliceMultiProd(Rtil_mat,L21mat);
                AuBu_Part12=fns_RmatGen.SliceMultiProd(Emat1,AuBu1);
                AuBu_Part1=fns_RmatGen.SliceMultiProd(AuBu_Part11,AuBu_Part12);
                %-----------------------------------%
                AuBubar_Part12=fns_RmatGen.SliceMultiProd(Emat1,AuBubar1);
                AuBubar_Part1=fns_RmatGen.SliceMultiProd(AuBu_Part11,AuBubar_Part12);
                %-----------------------------------%
                if i_L==J
                    [Q12mat,Q22mat,Q22barmat,Q12barmat] =...
                        fns_get_recurPara.construct_Qmat(alpha,beta, mu,...
                        omega, o_I, k_vect);
                    AuBu_Part2=fns_RmatGen.SliceMultiProd(Rtil_mat,Q22mat)-Q12mat;
                    AuBu_Part2=fns_RmatGen.SliceMultiProd(AuBu_Part2,FSFRVect);
                    AuBubar_Part2=fns_RmatGen.SliceMultiProd(Rtil_mat,Q22barmat)-Q12barmat;
                    AuBubar_Part2=fns_RmatGen.SliceMultiProd(AuBubar_Part2,FSFRVect);
                else
                    AuBu_Part2=zeros(2,1,length(k_vect));
                    AuBubar_Part2=zeros(2,1,length(k_vect));
                end
                AuBu=AuBu_Part1+AuBu_Part2;
                AuBubar=AuBubar_Part1+AuBubar_Part2;
                Pmat1=Pmat;
                Emat1=Emat;
                Rmat1=Rmat;
                Rmat_nanidx = isnan(Rmat);
                AuBu1=AuBu;
                AuBubar1=AuBubar;
                if any(Rmat_nanidx)
                    [row_idx, col_idx] = find(Rmat_nanidx);
                    for i = 1:length(row_idx)
                        warning('NaN value detected at index (%d, %d)',...
                            row_idx(i), col_idx(i));
                    end
                end
            end
        end
        %%
        function [Q12mat,Q22mat,Q22barmat,Q12barmat] = construct_Qmat(alpha, beta, mu,...
                omega, o_I, k_vect)
            kBetaVect = (omega + 1i*o_I) ./ beta;
            kAlphaVect = (omega + 1i*o_I) ./ alpha;
            gammaVect = (kBetaVect.^2 - k_vect.^2).^(1/2);
            nuVect = (kAlphaVect.^2 - k_vect.^2).^(1/2);

            %%
            Q22mat=zeros(2,2,length(k_vect));
            Q22mat(1,1,:)=-1i.*k_vect./nuVect;
            Q22mat(1,2,:)=ones(1,1,length(k_vect));
            Q22mat(2,1,:)=ones(1,1,length(k_vect));
            Q22mat(2,2,:)=-1i.*k_vect./gammaVect;
            Q22mat=(-1/(2*mu*kBetaVect.^2)).*Q22mat;
            %-----------%
            Q12mat=zeros(2,2,length(k_vect));
            Q12mat(1,1,:)=1i.*k_vect./nuVect;
            Q12mat(1,2,:)=ones(1,1,length(k_vect));
            Q12mat(2,1,:)=ones(1,1,length(k_vect));
            Q12mat(2,2,:)=1i.*k_vect./gammaVect;
            Q12mat=(-1/(2*mu*kBetaVect.^2)).*Q12mat;
            %-----------%
            Q22barmat=zeros(2,2,length(k_vect));
            Q22barmat(1,1,:)=-k_vect;
            Q22barmat(1,2,:)=-nuVect*1i;
            Q22barmat(2,1,:)=-gammaVect*1i;
            Q22barmat(2,2,:)=-k_vect;
            Q22barmat=(-1/(2*mu*kBetaVect.^2)).*Q22barmat;
            %-----------%
            Q12barmat=zeros(2,2,length(k_vect));
            Q12barmat(1,1,:)=-k_vect;
            Q12barmat(1,2,:)=nuVect*1i;
            Q12barmat(2,1,:)=gammaVect*1i;
            Q12barmat(2,2,:)=-k_vect;
            Q12barmat=(-1/(2*mu*kBetaVect.^2)).*Q12barmat;
        end
        %%
        function [Emat, Pmat] = construct_EPmat(alpha, beta, mu,...
                omega, o_I, k_vect, h)
            kBetaVect = (omega + 1i*o_I) ./ beta;
            kAlphaVect = (omega + 1i*o_I) ./ alpha;
            gammaVect = (kBetaVect.^2 - k_vect.^2).^(1/2);
            nuVect = (kAlphaVect.^2 - k_vect.^2).^(1/2);
            chiVect = k_vect.^2 - gammaVect.^2;

            Pmat = zeros(4, 4, length(k_vect));
            Pmat(1, 1, :) = k_vect;
            Pmat(1, 3, :) = k_vect;
            Pmat(2, 2, :) = k_vect;
            Pmat(2, 4, :) = k_vect;

            Pmat(1, 2, :) = 1i.*gammaVect;
            Pmat(1, 4, :) = -1i.*gammaVect;
            Pmat(2, 1, :) = 1i.*nuVect;
            Pmat(2, 3, :) = -1i.*nuVect;

            Pmat(3, 1, :) = 2*1i.*mu.*k_vect.*nuVect;
            Pmat(3, 3, :) = -2*1i.*mu.*k_vect.*nuVect;
            Pmat(3, 2, :) = mu*chiVect;
            Pmat(3, 4, :) = mu*chiVect;

            Pmat(4, 1, :) = mu*chiVect;
            Pmat(4, 3, :) = mu*chiVect;
            Pmat(4, 2, :) = 2*1i.*mu.*k_vect.*gammaVect;
            Pmat(4, 4, :) = -2*1i.*mu.*k_vect.*gammaVect;

            %%  
            Emat=zeros(2,2,length(k_vect));
            Emat(1,1,:)=exp(1i.*nuVect*h);
            Emat(1,2,:)=zeros(1,1,length(k_vect));
            Emat(2,1,:)=zeros(1,1,length(k_vect));
            Emat(2,2,:)=exp(1i.*gammaVect*h);
        end
        %%
        function [r,Cu,Cubar]=compute_rCu(omega, k_vect, o_I,...
                beta_vect, h_vect, mu_vect,N,J,F_T)
            beta1 = beta_vect(N+2);
            mu1 = mu_vect(N+2);
            kBeta1 = (omega + 1i*o_I) ./ beta1;
            gamma1 = (kBeta1.^2 - k_vect.^2).^(1/2);
            r1=zeros(1,length(k_vect));
            e1=zeros(1,length(k_vect));
            Cu1=zeros(1,length(k_vect));
            Cubar1=zeros(1,length(k_vect));
            % FTvect=k_vect.*F_T;
            % for loop runs from (second to last: N+1) to (first layer)
            for i_L=N+1:-1:1
                beta=beta_vect(i_L);
                mu=mu_vect(i_L);
                kBeta = (omega + 1i*o_I) ./ beta;
                gamma = (kBeta.^2 - k_vect.^2).^(1/2);
                h=h_vect(i_L);
                e=exp(1i.*gamma*h);

                r=((mu.*gamma.*(1+e1.*r1)-mu1.*gamma1.*(1-e1.*r1))./...
                    (mu.*gamma.*(1+e1.*r1)+mu1.*gamma1.*(1-e1.*r1))).*e;
                if i_L==J
                    Cu=(2*mu1.*gamma1.*e1.*Cu1)+(1i*F_T*(1+e1.*r1))./...
                        (mu.*gamma.*(1+e1.*r1)+mu1.*gamma1.*(1-e1.*r1));
                    Cubar=(2*mu1.*gamma1.*e1.*Cu1)+(gamma.*F_T.*(e1.*r1-1))./...
                        (mu.*gamma.*(1+e1.*r1)+mu1.*gamma1.*(1-e1.*r1));
                else
                    Cu=(2*mu1.*gamma1.*e1.*Cu1)./...
                        (mu.*gamma.*(1+e1.*r1)+mu1.*gamma1.*(1-e1.*r1));
                    Cubar=zeros(1,length(k_vect));
                end
                mu1=mu;
                gamma1=gamma;
                e1=e;
                r1=r;
                Cu1=Cu;
                Cubar1=Cubar;
                r_nanidx = isnan(r);
                if any(r_nanidx)
                    [row_idx, col_idx] = find(r_nanidx);
                    for i = 1:length(row_idx)
                        warning('NaN value detected at index (%d, %d)',...
                            row_idx(i), col_idx(i));
                    end
                end
            end
        end
    end
end