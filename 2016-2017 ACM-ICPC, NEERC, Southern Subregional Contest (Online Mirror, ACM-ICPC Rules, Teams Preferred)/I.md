# 省赛训练
recording problems </br>

## [contest730 I][贪心] Olympiad in Programming and Sports
## 题目描述
一共有n个人，每个人都有编程能力a(i)和运动能力b(i)</br>
要求选出p人参加编程队，s人参加运动队</br>
使得p人的a(i)值总和加上s人的b(i)值总和最大</br>
并且输出两队分别选择的人员编号</br>
$2\le n\le 3e3$ </br>
$0\le p+s\le n$ </br>
$1\le a(i)\le 3e3$ </br>
$1\le b(i)\le 3e3$ </br>

## 思路分析
主要思想是贪心</br>
假设现在1在编程队，2在运动队，当a_1+b_2 < b_1+a_2 时，即b_1-a_1>b_2-b_1时需要交换1、2从属的队伍</br>
故可以根据b_i-a_i的值来排序，以确定是否交换人员。

具体实现</br>
使用优先队列p1，p2分别存放a，b数组的值和ID。</br>
现假设编程队中所有人都取尽可能大的值，并保存其(b_i-a_i)至优先队列排序
其中，这部分人要么确实在编程队中，要么在运动队中更优。</br>
其被放入运动队的条件在于：当p2.top() (即b值最强且未被选过的人)放入运动队后(贡献为a_i+p2.top)，没有将(b_i-a_i)值最大的人换入运动队时(此时需在编程队新添加一人,即贡献为p1.top+b_i)更优。</br>
即当关系p1.top+b_i>a_i+p2.top时需要将p1.top放入编程队，将i放入运动队。</br>
由此循环直至将运动队填满位置，可达到贡献最大值。</br>


## 核心代码
``` c++
/* http://codeforces.com/contest/730 */
#include"bits/stdc++.h"
#include<iostream>
#include<string.h>
#include<algorithm>
#include<cstdio>
using namespace std;
const int N=3030;
#define pii pair<int,int>
int a[N],b[N],c[N];
int n,p,s;
int ans;
priority_queue<pii>p1,p2,tmp;
int main()
{
	scanf("%d%d%d",&n,&p,&s);
	for(int i=1;i<=n;i++){
		scanf("%d",&a[i]);
		p1.push(make_pair(a[i],i));
	}
	for(int i=1;i<=n;i++){
		scanf("%d",&b[i]);
		p2.push(make_pair(b[i],i));
	}
	int pp=p,ss=s;
	while(pp--){//先按照a[]中从大到小排满编程队(队1)
		int u=p1.top().second;
		p1.pop();
		c[u]=1;
		ans+=a[u];
		tmp.push(make_pair(b[u]-a[u],u));
        //b[u]-a[u]差值越大说明放在运动队(队2)中会更好
	}
	while(ss--){
        //去除已经排入对的点
		while(!p1.empty()&&c[p1.top().second])p1.pop();
        while(!p2.empty()&&c[p2.top().second])p2.pop();
        //表示p1.top+b[v]>a[v]+p2.top 说明将v点放入2中，并插入一个新的1，比在2中直接插入p2.top()更优
        if(p1.top().first+tmp.top().first>p2.top().first){
        	int u=p1.top().second;
        	int v=tmp.top().second;
        	p1.pop();	tmp.pop();
        	c[u]=1;c[v]=2;
        	tmp.push(make_pair(b[u]-a[u],u));
        	ans=ans-a[v]+a[u]+b[v];
		}
		else{//此时将p2.top()放在2中为最优情况
			int v=p2.top().second;
			p2.top();
			c[v]=2;
			ans+=b[v];
		}
	}
	printf("%d\n",ans);
	for(int i=1;i<=n;i++)
		if(c[i]==1)printf("%d ",i);
	printf("\n");
	for(int i=1;i<=n;i++)
		if(c[i]==2)printf("%d ",i);
}

```
