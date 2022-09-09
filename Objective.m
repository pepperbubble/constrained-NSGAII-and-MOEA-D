
function [Output,Boundary,Infeasible] = Objective(Operation,Problem,M,Input)
%input:N ��Ⱥ��ģ
    persistent K; %����־ñ���
    %�־ñ������������ǵĺ����ľֲ�������
    %����ֵ�����ڶԸú����ĸ��ε�����ʹ�õ��ڴ��С�
    
    Boundary = NaN;
    switch Operation  %ѡ�����ģʽ0/1
        
        %0����ʼ����Ⱥ
        case 0
            k = find(~isstrprop(Problem,'digit'),1,'last'); % �ж��м���Ӣ����ĸ��k =4
            %isstrprop:����ַ�ÿһ���ַ��Ƿ�����ָ���ķ�Χ,'digit'�ж��ǲ������֣�����0/1����
            %��~�� ��0/1����ȡ��
            switch Problem
                case 'func1'
                    MaxValue = [99 99 200 200];
                    MinValue = [1 1 10 10];
                    Boundary = [MaxValue;MinValue];
                    D = 4;
                    Population = rand(Input,D);   %ȷ����Ⱥ��С��(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;
                case 'func2'
                    MaxValue = [0.5 0.5 0.6 0.5 6];
                    MinValue = [0.05 0.2 0.2 0.35 3];
                    Boundary = [MaxValue;MinValue];
                    D = 5;
                    Population = rand(Input,D);   %ȷ����Ⱥ��С��(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;
                case 'func3'
                    MaxValue = [100 100 3];
                    MinValue = [1e-5 1e-5 1];
                    Boundary = [MaxValue;MinValue];
                    D = 3;
                    Population = rand(Input,D);   %ȷ����Ⱥ��С��(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;
                case 'func4'
                    MaxValue = [5 10 10 5];
                    MinValue = [0.125 0.1 0.1 0.125];
                    Boundary = [MaxValue;MinValue];
                    D = 4;
                    Population = rand(Input,D);   %ȷ����Ⱥ��С��(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;
                case 'func5'
                    MaxValue = [80 110 3000 20];
                    MinValue = [55 75 1000 11];
                    Boundary = [MaxValue;MinValue];
                    D = 4;
                    Population = rand(Input,D);   %ȷ����Ⱥ��С��(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;  
                case 'func8'
                    MaxValue = [1.5 1.35 1.5 1.5 2.625 1.2 1.2];
                    MinValue = [0.5 0.45 0.5 0.5 0.875 0.4 0.4];
                    Boundary = [MaxValue;MinValue];
                    D = 7;
                    Population = rand(Input,D);   %ȷ����Ⱥ��С��(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;  
                case 'func9'
                    F=10;
                    sig=10;
                    MaxValue = [(3*F/sig) (3*F/sig) (3*F/sig) (3*F/sig)];
                    MinValue = [(F/sig) (sqrt(2)*F/sig) (sqrt(2)*F/sig) (F/sig)];
                    Boundary = [MaxValue;MinValue];
                    D = 4;
                    Population = rand(Input,D);   %ȷ����Ⱥ��С��(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population; 
                case 'func16'
                    MaxValue = [0.05 1];
                    MinValue = [0.01 0.2];
                    Boundary = [MaxValue;MinValue];
                    D = 2;
                    Population = rand(Input,D);   %ȷ����Ⱥ��С��(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population; 
                case 'func21'
                    MaxValue = [1.7 3.5 1.7 1.7 1.7 1.7];
                    MinValue = [1.3 2.5 1.3 1.3 1.3 1.3];
                    Boundary = [MaxValue;MinValue];
                    D = 6;
                    Population = rand(Input,D);   %ȷ����Ⱥ��С��(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population; 
                case 'func12'
                    MaxValue = [80 50 5 5];
                    MinValue = [10 10 0.9 0.9];
                    Boundary = [MaxValue;MinValue];
                    D = 4;
                    Population = rand(Input,D);   %ȷ����Ⱥ��С��(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;       
                case 'func23'
                    MaxValue = [1 1 1 1 16 16];
                    MinValue = [0 0 0 0 1e-5 1e-5];
                    Boundary = [MaxValue;MinValue];
                    D = 6;
                    Population = rand(Input,D);   %ȷ����Ⱥ��С��(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;     
                    
            otherwise
                K = [5 10 10 10 10];
                K = K(str2double(Problem(k+1:end))); %'DTZL2'K=10
                %str2double��һ�ֺ������书���ǰ��ַ���ת����ֵ��

                D = M+K-1; %D����������
                MaxValue   = ones(1,D);
                MinValue   = zeros(1,D);
                Population = rand(Input,D);   %ȷ����Ⱥ��С��(N*D)
                Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                % %�����µĳ�ʼ��Ⱥ
                Output   = Population;
                Boundary = [MaxValue;MinValue];
            end
        %1������Ŀ�꺯��ֵ������ֻ����DTLZ1~DTZL4��������
        case 1
            Population    = Input;  %�Ѿ���ʼ����ɵ���Ⱥ
            [popSize,D] = size(Population);
            FunctionValue = zeros(size(Population,1),M);
            switch Problem
                case 'DTLZ1'
                    g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(20.*pi.*(Population(:,M:end)-0.5)),2));
                    for i = 1 : M  %�����iάĿ�꺯��ֵ
                        FunctionValue(:,i) = 0.5.*prod(Population(:,1:M-i),2).*(1+g);
                        if i > 1
                            FunctionValue(:,i) = FunctionValue(:,i).*(1-Population(:,M-i+1));
                        end
                    end
                case 'DTLZ2'
                    g = sum((Population(:,M:end)-0.5).^2,2);
                    for i = 1 : M
                        FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                        if i > 1
                            FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                        end
                    end
                case 'DTLZ3'
                    g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(20.*pi.*(Population(:,M:end)-0.5)),2));
                    for i = 1 : M
                        FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                        if i > 1
                            FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                        end
                    end
                case 'DTLZ4'
                    Population(:,1:M-1) = Population(:,1:M-1).^100;
                    g = sum((Population(:,M:end)-0.5).^2,2);
                    for i = 1 : M
                        FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                        if i > 1
                            FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                        end
                    end
                case 'DTLZ5'
                    g = sum((Population(:,M:end)-0.5).^2,2);
                    
                    for i = 1 : M
                        theta =(pi./(4*(1+g))).*(1+2.*g.*Population(:,1:M-i));
                        FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*theta),2);
                        if i > 1
                            FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*(pi./(4*(1+g))).*(1+2.*g.*Population(:,M-i+1)));
                        end
                    end
                otherwise
                    k = find(~isstrprop(Problem,'digit'),1,'last'); % �ж��м���Ӣ����ĸ��k = 4
                    a = Problem(:,k+1:end);
                    num = str2num(a);
                    [f,g,h] = CEC2021_func(Population,num);
                    [r,c] = size(g);
                    if c ~= 1
                        g2 = sum(g');
                        g2 = g2';
                    end
                    if c==1
                        g2=g;
                    end
                    Infeasible = zeros(popSize, 1);
                    for i=1:r
                        Infeasible(i) = Infeasible(i)+g2(i);
                    end
                    FunctionValue = f; 
            end
            Output = FunctionValue;
    end
end