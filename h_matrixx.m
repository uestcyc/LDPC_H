%����һ�����ҽڵ�����жԣ����볤 ����һ��h��
%1.0�汾ֻ����һ������ľ����û������Ȧ�Ĵ���
%����ÿһ���ظ��Ĳ��ܳ���1֮��û���κ�����


%************************����������********************
clear all;

   xx_1=input('\n','s');
   if isempty(xx_1)
       lamda(:,1)=[4;3];
       lamda(:,2)=[0.8884;0.1116];
   else 
   lamda1=sscanf(xx_1,'%f');
   lamda(:,1)=lamda1(1:2:length(lamda1));
   lamda(:,2)=lamda1(2:2:length(lamda1));
end;
   

  xx_2=input('\n','s');
  if isempty(xx_2)
      ro(:,1)=[8;7];
      ro(:,2)=[0.7395;0.2605];
  else
   ro1=sscanf(xx_2,'%f');
   ro(:,1)=ro1(1:2:length(ro1));
   ro(:,2)=ro1(2:2:length(ro1));
end;

  code_length=input('');
  if isempty(code_length)
      code_length=100;
  end;
  


lamda_c=0;
for i1=1:length(lamda(:,1))
lamda_c=lamda_c+lamda(i1,2)/lamda(i1,1); 
end;
e_num=code_length/lamda_c%�ߵ���Ŀ
ro_c=0;
for i1=1:length(ro(:,1))
   ro_c=ro_c+ro(i1,2)/ro(i1,1); 
end;
r=1-ro_c/lamda_c %����

for i1=1:length(lamda(:,1))
l_i_node(i1,:)=[lamda(i1,1),e_num*lamda(i1,2)/lamda(i1,1)]%��ڵ�Ķ��Լ��˶ȵĽڵ���
end;

for i1=1:length(ro(:,1))
r_i_node(i1,:)=[ro(i1,1),e_num*ro(i1,2)/ro(i1,1)]%��ڵ�Ķ��Լ��˶ȵĽڵ���
end;



   l_array1=ones(l_i_node(1,1),1)*[1:round(l_i_node(1,2))];
   length1=l_i_node(1,1)*round(l_i_node(1,2));
   l_array(1:length1)=l_array1;
   clear l_array1;%��ս�ſռ�Ϊ�´�ʹ��׼��
   
   i1=2;
  while i1<length(l_i_node(:,1))+1
    index1=l_array(length(l_array));
    l_array1=ones(l_i_node(i1,1),1)*[index1+1:index1+round(l_i_node(i1,2))];
    length1=l_i_node(i1,1)*round(l_i_node(i1,2));
    l_array2(1:length1)=l_array1;
    l_array=[l_array,l_array2];
    clear l_array2;
    clear l_array1;%Ϊ�´�ʹ�����Ҫ��Ȼ������ Ҳ��ʡ���ڴ�
    i1=i1+1;
end;
 
 

   r_array1=ones(r_i_node(1,1),1)*[1:round(r_i_node(1,2))];
   length1=r_i_node(1,1)*round(r_i_node(1,2));
   r_array(1:length1)=r_array1;
   clear r_array1;%��ս�ſռ�Ϊ�´�ʹ��׼��
   
   i1=2;
  while i1<length(r_i_node(:,1))+1
    index1=r_array(length(r_array));
    r_array1=ones(r_i_node(i1,1),1)*[index1+1:index1+round(r_i_node(i1,2))];
    length1=r_i_node(i1,1)*round(r_i_node(i1,2));
    r_array2(1:length1)=r_array1;
    r_array=[r_array,r_array2];
    clear r_array2;
    clear r_array1;%Ϊ�´�ʹ�����Ҫ��Ȼ������ Ҳ��ʡ���ڴ�
    i1=i1+1;
end;




len_l_array=length(l_array);

h_coordinate=[0,0];%������ʼ�� ������� ������һ��ȥ��
i1=1;
while i1<len_l_array
    
    index2=find(l_array==l_array(i1));%�ҵ�r_array(i1) ��l_array�ķֲ� Ϊÿһ��Ԫ����r_array���������һ����Ӧֵ����Щ��Ӧ��ֵ������ͬ
    
    i1=i1+length(index2);%�����Ĳ������ǹ̶�����Ϊ l_array =*1     1     1     1     1     1     1     1     *2     2     2     2     2     2     2     2     *3     3     3     
    
    randnum=randint(1,length(index2),[1,length(r_array)]);%����һ���������
    coordinate=[r_array(randnum)',l_array(index2)'];
    
    
    hh=sparse(ones(length(coordinate(:,1)),1),coordinate(:,1),ones(length(coordinate(:,1)),1));%Ҫ������ͬ��ֵ��hh�ض��д���1��ֵ��ֹ���رߵ�
    i2=0;
    while max(hh)>1&(i2<100)%���������׵ľͲ����ر����Ҫ�趨һ��ѭ�������ֵ������ѭ��
        i2=i2+1
        randnum=randint(1,length(index2),[1,length(r_array)]);%����һ���������
        coordinate=[r_array(randnum)',l_array(index2)'];
        hh=sparse(ones(length(coordinate(:,1)),1),coordinate(:,1),ones(length(coordinate(:,1)),1));
    end;
  %���Ѿ������Ĵ�r_array��ȥ��
  r_array(randnum)=0;
  r_array=sort(r_array);
  length_r_array=length(r_array);
  r_array=r_array(length(randnum)+1:length_r_array);
  %������
    h_coordinate=[h_coordinate; coordinate];
    clear index2;%��Ϊ�´μ����index2���ȿ��ܻ�仯����´εĳ���С�Ļ� �� �����±��μ����ֵ ���Ҫ�����
    clear randnum;%��Ϊ�´μ����index2���ȿ��ܻ�仯����´εĳ���С�Ļ� �� �����±��μ����ֵ ���Ҫ�����
    clear hh;%ͬ��
end;
h_coordinate=h_coordinate(2:length(h_coordinate),:);
h=sparse(h_coordinate(:,1),h_coordinate(:,2),ones(length(h_coordinate(:,1)),1));%��û�������ص�
%hh=pailie(h);

%dlmwrite('G:\xingang\ldpc matlab\H_martrix generater\random generater\first\r_h.txt',h);
%dlmwrite('G:\xingang\ldpc matlab\H_martrix generater\random generater\first\r_hh.txt',hh);

spy(h);





