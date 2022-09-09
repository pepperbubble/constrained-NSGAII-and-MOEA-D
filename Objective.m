
function [Output,Boundary,Infeasible] = Objective(Operation,Problem,M,Input)
%input:N 种群规模
    persistent K; %定义持久变量
    %持久变量是声明它们的函数的局部变量；
    %但其值保留在对该函数的各次调用所使用的内存中。
    
    Boundary = NaN;
    switch Operation  %选择操作模式0/1
        
        %0：初始化种群
        case 0
            k = find(~isstrprop(Problem,'digit'),1,'last'); % 判断有几个英文字母，k =4
            %isstrprop:检测字符每一个字符是否属于指定的范围,'digit'判断是不是数字，返回0/1数组
            %“~” 对0/1数组取反
            switch Problem
                case 'func1'
                    MaxValue = [99 99 200 200];
                    MinValue = [1 1 10 10];
                    Boundary = [MaxValue;MinValue];
                    D = 4;
                    Population = rand(Input,D);   %确定种群大小：(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;
                case 'func2'
                    MaxValue = [0.5 0.5 0.6 0.5 6];
                    MinValue = [0.05 0.2 0.2 0.35 3];
                    Boundary = [MaxValue;MinValue];
                    D = 5;
                    Population = rand(Input,D);   %确定种群大小：(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;
                case 'func3'
                    MaxValue = [100 100 3];
                    MinValue = [1e-5 1e-5 1];
                    Boundary = [MaxValue;MinValue];
                    D = 3;
                    Population = rand(Input,D);   %确定种群大小：(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;
                case 'func4'
                    MaxValue = [5 10 10 5];
                    MinValue = [0.125 0.1 0.1 0.125];
                    Boundary = [MaxValue;MinValue];
                    D = 4;
                    Population = rand(Input,D);   %确定种群大小：(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;
                case 'func5'
                    MaxValue = [80 110 3000 20];
                    MinValue = [55 75 1000 11];
                    Boundary = [MaxValue;MinValue];
                    D = 4;
                    Population = rand(Input,D);   %确定种群大小：(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;  
                case 'func8'
                    MaxValue = [1.5 1.35 1.5 1.5 2.625 1.2 1.2];
                    MinValue = [0.5 0.45 0.5 0.5 0.875 0.4 0.4];
                    Boundary = [MaxValue;MinValue];
                    D = 7;
                    Population = rand(Input,D);   %确定种群大小：(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;  
                case 'func9'
                    F=10;
                    sig=10;
                    MaxValue = [(3*F/sig) (3*F/sig) (3*F/sig) (3*F/sig)];
                    MinValue = [(F/sig) (sqrt(2)*F/sig) (sqrt(2)*F/sig) (F/sig)];
                    Boundary = [MaxValue;MinValue];
                    D = 4;
                    Population = rand(Input,D);   %确定种群大小：(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population; 
                case 'func16'
                    MaxValue = [0.05 1];
                    MinValue = [0.01 0.2];
                    Boundary = [MaxValue;MinValue];
                    D = 2;
                    Population = rand(Input,D);   %确定种群大小：(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population; 
                case 'func21'
                    MaxValue = [1.7 3.5 1.7 1.7 1.7 1.7];
                    MinValue = [1.3 2.5 1.3 1.3 1.3 1.3];
                    Boundary = [MaxValue;MinValue];
                    D = 6;
                    Population = rand(Input,D);   %确定种群大小：(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population; 
                case 'func12'
                    MaxValue = [80 50 5 5];
                    MinValue = [10 10 0.9 0.9];
                    Boundary = [MaxValue;MinValue];
                    D = 4;
                    Population = rand(Input,D);   %确定种群大小：(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;       
                case 'func23'
                    MaxValue = [1 1 1 1 16 16];
                    MinValue = [0 0 0 0 1e-5 1e-5];
                    Boundary = [MaxValue;MinValue];
                    D = 6;
                    Population = rand(Input,D);   %确定种群大小：(N*D)
                    Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                    Output = Population;     
                    
            otherwise
                K = [5 10 10 10 10];
                K = K(str2double(Problem(k+1:end))); %'DTZL2'K=10
                %str2double是一种函数，其功能是把字符串转换数值，

                D = M+K-1; %D：变量数量
                MaxValue   = ones(1,D);
                MinValue   = zeros(1,D);
                Population = rand(Input,D);   %确定种群大小：(N*D)
                Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
                % %产生新的初始种群
                Output   = Population;
                Boundary = [MaxValue;MinValue];
            end
        %1：计算目标函数值；这里只包含DTLZ1~DTZL4函数问题
        case 1
            Population    = Input;  %已经初始化完成的种群
            [popSize,D] = size(Population);
            FunctionValue = zeros(size(Population,1),M);
            switch Problem
                case 'DTLZ1'
                    g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(20.*pi.*(Population(:,M:end)-0.5)),2));
                    for i = 1 : M  %计算第i维目标函数值
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
                    k = find(~isstrprop(Problem,'digit'),1,'last'); % 判断有几个英文字母，k = 4
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