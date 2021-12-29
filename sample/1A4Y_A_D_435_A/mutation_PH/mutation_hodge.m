clc;
clear;

formatSpec = '%d %f %f %f';
sizeA = [4,Inf];
nloop = 48;
idim = 2;
file0 = fopen(strcat('./wild_mutation_eig.txt'),'w');
file1 = fopen(strcat('./wild_mutation_simp.txt'),'w');

for j=1:3
    for k=1:3
        Elements = {'C','N','O'};
        e1 = Elements{j}; e2 = Elements{k};
        fprintf('%s %s \n', e1, e2);
        WildName = strcat('wild_mutation_',e1,'_',e2,'.pts');
        fileID = fopen(strcat('./', WildName), 'r');
        A = fscanf(fileID, formatSpec, sizeA);
        fclose(fileID);
        for iloop=1:nloop
            VRdiameter=0.25*iloop;
            distances = zeros(size(A,2),size(A,2));
            for ii=1:size(A,2)
                for jj=(ii+1):size(A,2)
                    if A(1,ii)+A(1,jj) == 1
                        dis = sqrt((A(2,ii) - A(2,jj))^2 + (A(3,ii) - A(3,jj))^2 + (A(4,ii) - A(4,jj))^2);
                        if (dis < VRdiameter)
                            distances(ii,jj) = 1;
                            distances(jj,ii) = 1;
                            fprintf('%.2f, %d, %d, %.2f\n', VRdiameter, ii, jj, dis);
                        end
                    end
                end
            end

            [kSkeletonOfVRComplex,simplexDimension]=computeVRComplex(distances,idim);

           %%  computer boundaries
            N = size(A,2);
            bdyMatricies{1}=zeros(N,1);
            for it=2:idim     
                simplexDimension=it;

                KSimplicies=kSkeletonOfVRComplex{simplexDimension};

            % if (simplexDimension==1)
            %     boundaryMatrix=zeros(size(KSimplicies,1),1)';
            %   %  fprintf('simplexdimension: %i %i \n', simplexDimension,it); 
            % end

                KMinusOneSimplicies=kSkeletonOfVRComplex{simplexDimension-1};
                nKSimplicies=size(KSimplicies,1);
                nKMinusOneSimplicies=size(KMinusOneSimplicies,1);
                boundaryMatrix=zeros(nKMinusOneSimplicies,nKSimplicies);

            % fprintf('simplexdimension: %i %i %i\n', simplexDimension,size(boundaryMatrix));
            % fprintf('bdmatrix dimension: %i %i %i\n', nKSimplicies,nKMinusOneSimplicies);

                for i=1:nKMinusOneSimplicies
                    for m=1:nKSimplicies

                        if KSimplicies(m,1)>KMinusOneSimplicies(i,1)
                            break
                        end

                        if KMinusOneSimplicies(i,simplexDimension-1)>KSimplicies(m,simplexDimension)
                            continue 
                        end

                        if ismembc(KMinusOneSimplicies(i,:),KSimplicies(m,:))
                            [a, indiciesOfMissingElements] = find(ismembc(KSimplicies(m,:), KMinusOneSimplicies(i,:))==0);
                            boundaryMatrix(i,m)=(-1)^mod(indiciesOfMissingElements+1,2);  %%%%%%
                            %  boundaryMatrix(i,m)=1;  %%%% remove sign
                        end

                    end
                end
                bdyMatricies{it}=boundaryMatrix(:,:); 
            end
           %% calculate eigenvectors
            for id=1:idim-1
                nmt=size(kSkeletonOfVRComplex{id},1);
                fprintf(file1,'%s %s %5d %6.2f %7d\n',Elements{j}, Elements{k}, id, VRdiameter, nmt);
                hodgematrix{id}=zeros(nmt,nmt);
                matI=bdyMatricies{id};
                matIadd1=bdyMatricies{id+1};
                if(id==1)
                    hodgematrix{id}=matIadd1*matIadd1';
                else
                    hodgematrix{id}=matI'*matI+matIadd1*matIadd1';
                end
            end
            for i =1:idim-1
                mathdg=hodgematrix{i};
            %     if numel(mathdg)>0
                [Vec,Dmat] = eig(mathdg);
                ndg=size(Dmat,1);
                ivct=[];
                for indg=1:ndg
                    ivct(indg)=Dmat(indg,indg);
                    fprintf(file0,'%s %s %5d %6.2f %18.3e\n',Elements{j}, Elements{k}, i, VRdiameter, ivct(indg));
                end
            end
        end
    end
end

clearvars -except indexlist fileID C idim data nloop file0 file1;
fclose(file0);
fclose(file1);

formatSpec = '%d %f %f %f';
sizeA = [4,Inf];
nloop = 48;
idim = 2;
file0 = fopen(strcat('./mut_mutation_eig.txt'),'w');
file1 = fopen(strcat('./mut_mutation_simp.txt'),'w');

for j=1:3
    for k=1:3
        Elements = {'C','N','O'};
        e1 = Elements{j}; e2 = Elements{k};
        fprintf('%s %s \n', e1, e2);
        WildName = strcat('mut_mutation_',e1,'_',e2,'.pts');
        fileID = fopen(strcat('./', WildName), 'r');
        A = fscanf(fileID, formatSpec, sizeA);
        fclose(fileID);
        for iloop=1:nloop
            VRdiameter=0.25*iloop;
            distances = zeros(size(A,2),size(A,2));
            for ii=1:size(A,2)
                for jj=(ii+1):size(A,2)
                    if A(1,ii)+A(1,jj) == 1
                        dis = sqrt((A(2,ii) - A(2,jj))^2 + (A(3,ii) - A(3,jj))^2 + (A(4,ii) - A(4,jj))^2);
                        if (dis < VRdiameter)
                            distances(ii,jj) = 1;
                            distances(jj,ii) = 1;
                            fprintf('%.2f, %d, %d, %.2f\n', VRdiameter, ii, jj, dis);
                        end
                    end
                end
            end

            [kSkeletonOfVRComplex,simplexDimension]=computeVRComplex(distances,idim);

           %%  computer boundaries
            N = size(A,2);
            bdyMatricies{1}=zeros(N,1);
            for it=2:idim     
                simplexDimension=it;

                KSimplicies=kSkeletonOfVRComplex{simplexDimension};

            % if (simplexDimension==1)
            %     boundaryMatrix=zeros(size(KSimplicies,1),1)';
            %   %  fprintf('simplexdimension: %i %i \n', simplexDimension,it); 
            % end

                KMinusOneSimplicies=kSkeletonOfVRComplex{simplexDimension-1};
                nKSimplicies=size(KSimplicies,1);
                nKMinusOneSimplicies=size(KMinusOneSimplicies,1);
                boundaryMatrix=zeros(nKMinusOneSimplicies,nKSimplicies);

            % fprintf('simplexdimension: %i %i %i\n', simplexDimension,size(boundaryMatrix));
            % fprintf('bdmatrix dimension: %i %i %i\n', nKSimplicies,nKMinusOneSimplicies);

                for i=1:nKMinusOneSimplicies
                    for m=1:nKSimplicies

                        if KSimplicies(m,1)>KMinusOneSimplicies(i,1)
                            break
                        end

                        if KMinusOneSimplicies(i,simplexDimension-1)>KSimplicies(m,simplexDimension)
                            continue 
                        end

                        if ismembc(KMinusOneSimplicies(i,:),KSimplicies(m,:))
                            [a, indiciesOfMissingElements] = find(ismembc(KSimplicies(m,:), KMinusOneSimplicies(i,:))==0);
                            boundaryMatrix(i,m)=(-1)^mod(indiciesOfMissingElements+1,2);  %%%%%%
                            %  boundaryMatrix(i,m)=1;  %%%% remove sign
                        end

                    end
                end
                bdyMatricies{it}=boundaryMatrix(:,:); 
            end
           %% calculate eigenvectors
            for id=1:idim-1
                nmt=size(kSkeletonOfVRComplex{id},1);
                fprintf(file1,'%s %s %5d %6.2f %7d\n',Elements{j}, Elements{k}, id, VRdiameter, nmt);
                hodgematrix{id}=zeros(nmt,nmt);
                matI=bdyMatricies{id};
                matIadd1=bdyMatricies{id+1};
                if(id==1)
                    hodgematrix{id}=matIadd1*matIadd1';
                else
                    hodgematrix{id}=matI'*matI+matIadd1*matIadd1';
                end
            end
            for i =1:idim-1
                mathdg=hodgematrix{i};
            %     if numel(mathdg)>0
                [Vec,Dmat] = eig(mathdg);
                ndg=size(Dmat,1);
                ivct=[];
                for indg=1:ndg
                    ivct(indg)=Dmat(indg,indg);
                    fprintf(file0,'%s %s %5d %6.2f %18.3e\n',Elements{j}, Elements{k}, i, VRdiameter, ivct(indg));
                end
            end
        end
    end
end

clearvars -except indexlist fileID C idim data nloop file0 file1;
fclose(file0);
fclose(file1);

