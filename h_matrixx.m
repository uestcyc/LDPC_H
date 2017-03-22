%给定一个左右节点的序列对，和码长 生成一个h阵
%1.0版本只生成一个随机的矩阵而没有消除圈的处理
%除了每一列重复的不能超过1之外没有任何限制


%************************参数的输入********************
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
e_num=code_length/lamda_c%边的数目
ro_c=0;
for i1=1:length(ro(:,1))
   ro_c=ro_c+ro(i1,2)/ro(i1,1); 
end;
r=1-ro_c/lamda_c %码率

for i1=1:length(lamda(:,1))
l_i_node(i1,:)=[lamda(i1,1),e_num*lamda(i1,2)/lamda(i1,1)]%左节点的度以及此度的节点数
end;

for i1=1:length(ro(:,1))
r_i_node(i1,:)=[ro(i1,1),e_num*ro(i1,2)/ro(i1,1)]%左节点的度以及此度的节点数
end;



   l_array1=ones(l_i_node(1,1),1)*[1:round(l_i_node(1,2))];
   length1=l_i_node(1,1)*round(l_i_node(1,2));
   l_array(1:length1)=l_array1;
   clear l_array1;%清空解放空间为下次使用准备
   
   i1=2;
  while i1<length(l_i_node(:,1))+1
    index1=l_array(length(l_array));
    l_array1=ones(l_i_node(i1,1),1)*[index1+1:index1+round(l_i_node(i1,2))];
    length1=l_i_node(i1,1)*round(l_i_node(i1,2));
    l_array2(1:length1)=l_array1;
    l_array=[l_array,l_array2];
    clear l_array2;
    clear l_array1;%为下次使用清空要不然会出错的 也节省了内存
    i1=i1+1;
end;
 
 

   r_array1=ones(r_i_node(1,1),1)*[1:round(r_i_node(1,2))];
   length1=r_i_node(1,1)*round(r_i_node(1,2));
   r_array(1:length1)=r_array1;
   clear r_array1;%清空解放空间为下次使用准备
   
   i1=2;
  while i1<length(r_i_node(:,1))+1
    index1=r_array(length(r_array));
    r_array1=ones(r_i_node(i1,1),1)*[index1+1:index1+round(r_i_node(i1,2))];
    length1=r_i_node(i1,1)*round(r_i_node(i1,2));
    r_array2(1:length1)=r_array1;
    r_array=[r_array,r_array2];
    clear r_array2;
    clear r_array1;%为下次使用清空要不然会出错的 也节省了内存
    i1=i1+1;
end;




len_l_array=length(l_array);

h_coordinate=[0,0];%迭代初始化 否则出错 最后把这一行去掉
i1=1;
while i1<len_l_array
    
    index2=find(l_array==l_array(i1));%找到r_array(i1) 在l_array的分布 为每一个元素在r_array中随机的找一个对应值但这些对应的值不能相同
    
    i1=i1+length(index2);%递增的步长不是固定的因为 l_array =*1     1     1     1     1     1     1     1     *2     2     2     2     2     2     2     2     *3     3     3     
    
    randnum=randint(1,length(index2),[1,length(r_array)]);%生成一个随机数列
    coordinate=[r_array(randnum)',l_array(index2)'];
    
    
    hh=sparse(ones(length(coordinate(:,1)),1),coordinate(:,1),ones(length(coordinate(:,1)),1));%要是有相同的值则hh必定有大于1得值防止有重边的
    i2=0;
    while max(hh)>1&(i2<100)%到最后很容易的就产生重边因此要设定一个循环的最大值以免死循环
        i2=i2+1
        randnum=randint(1,length(index2),[1,length(r_array)]);%生成一个随机数列
        coordinate=[r_array(randnum)',l_array(index2)'];
        hh=sparse(ones(length(coordinate(:,1)),1),coordinate(:,1),ones(length(coordinate(:,1)),1));
    end;
  %把已经连过的从r_array中去掉
  r_array(randnum)=0;
  r_array=sort(r_array);
  length_r_array=length(r_array);
  r_array=r_array(length(randnum)+1:length_r_array);
  %清除完毕
    h_coordinate=[h_coordinate; coordinate];
    clear index2;%因为下次计算的index2长度可能会变化如果下次的长度小的话 则 会留下本次计算的值 因此要清除掉
    clear randnum;%因为下次计算的index2长度可能会变化如果下次的长度小的话 则 会留下本次计算的值 因此要清除掉
    clear hh;%同上
end;
h_coordinate=h_coordinate(2:length(h_coordinate),:);
h=sparse(h_coordinate(:,1),h_coordinate(:,2),ones(length(h_coordinate(:,1)),1));%还没有消除重叠
%hh=pailie(h);

%dlmwrite('G:\xingang\ldpc matlab\H_martrix generater\random generater\first\r_h.txt',h);
%dlmwrite('G:\xingang\ldpc matlab\H_martrix generater\random generater\first\r_hh.txt',hh);

spy(h);





