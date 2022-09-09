%NSGA-II
function NSGAII(Problem,M)
clc;
format compact;%�ո����
tic;%��¼����ʱ�� 


%�����趨
Generations = 1000;   %��������
if M == 2  %Ŀ������
    N = 100; %��Ⱥ��С
else M == 3
    N = 105;
end

    %��ʼ����Ⱥ�����س�ʼ����Ⱥ�;��߱���������
    [Population,Boundary] = Objective(0,Problem,M,N);
    [FunctionValue,null,Infeasible] = Objective(1,Problem,M,Population);%����Ŀ�꺯��ֵ
    
    % ���з�֧������
    FrontValue = NonDominateSort(FunctionValue,0); 
    CrowdDistance = CrowdDistances(FunctionValue,FrontValue);%����ۼ�����
    
    %��ʼ����
    for Gene = 1 : Generations    
        %�����Ӵ���
        MatingPool = Mating(Population,FrontValue,CrowdDistance); %�����ѡ��2�Ľ�����ѡ��ʽ
        Offspring = NSGA_Gen(MatingPool,Boundary,N); %����,���죬Խ�紦�������µ���Ⱥ

        Population = [Population;Offspring];  %��Ⱥ�ϲ�
        
        [FunctionValue,null,Infeasible] = Objective(1,Problem,M,Population);%����Ŀ�꺯��ֵ
        [FrontValue,MaxFront] = NonDominateSort(FunctionValue,1); % ���з�֧������
        CrowdDistance = CrowdDistances(FunctionValue,FrontValue);%����ۼ�����

        
        %ѡ����֧��ĸ���        
        Next = zeros(1,N);
        NoN = numel(FrontValue,FrontValue<MaxFront);
        Next(1:NoN) = find(FrontValue<MaxFront);
        
        %ѡ�����һ����ĸ���
        Last = find(FrontValue==MaxFront);
        [~,Rank] = sort(CrowdDistance(Last),'descend');
        Next(NoN+1:N) = Last(Rank(1:N-NoN));
        
        %��һ����Ⱥ
        Population = Population(Next,:);
        FrontValue = FrontValue(Next);
        CrowdDistance = CrowdDistance(Next);
        
		[FunctionValue,null,Infeasible] = Objective(1,Problem,M,Population);
        
        cla;%cla �ӵ�ǰ������ɾ�������ɼ����������ͼ�ζ���
        DrawGraph(FunctionValue);%�Ѻ���ֵ�ĵ㻭����ά����ϵ�
        hold on;%������֪Pareto������,�������ŵ�
        switch Problem
            case 'DTLZ1'
                if M == 2
                    pareto_x = linspace(0,0.5);
                    pareto_y = 0.5 - pareto_x;
                    plot(pareto_x, pareto_y, 'r');
                elseif M == 3
                    [pareto_x,pareto_y]  = meshgrid(linspace(0,0.5));
                    pareto_z = 0.5 - pareto_x - pareto_y;
                    axis([0,1,0,1,0,1]);
                    mesh(pareto_x, pareto_y, pareto_z);
                end
            otherwise
                if M == 2
                    k = find(~isstrprop(Problem,'digit'),1,'last'); % �ж��м���Ӣ����ĸ��k = 4
                    a = Problem(:,k+1:end);
                    num = str2num(a);
                    point = load(['RefPoints\nadir_' num2str(num) '.txt']);
                    scatter(point(1),point(2));
                    title(Problem);
                elseif M == 3
                    k = find(~isstrprop(Problem,'digit'),1,'last'); % �ж��м���Ӣ����ĸ��k = 4
                    a = Problem(:,k+1:end);
                    num = str2num(a);
                    point = load(['RefPoints\nadir_' num2str(num) '.txt']);
                    scatter(point(1),point(2),point(3));
                    title(Problem);
                end
        end
        pause(0.01);
        %clc;
    end
    HV = 0;
    k = find(~isstrprop(Problem,'digit'),1,'last'); % �ж��м���Ӣ����ĸ��k = 4
    a = Problem(:,k+1:end);
    num = str2num(a);
    HV = HV_score(FunctionValue, num);
    disp(HV);
end
