# 省赛训练
recording problems </br>

## [contest730 D][贪心] Running Over The Bridges
## 题目描述
小明要跨过n座桥，每座桥长度为l(i),停留时间不能超过t(i)，不然桥会塌</br>
小明每秒钟可以移动0.5个单位,他也可以嗑药，在r秒内获得1单位/秒的速度</br>
求小明最少嗑多少药才能安全跨国所有桥</br>
$1\le n\le 2e5$ </br>
$1\le r\le 1e12$ </br>
$1\le l(i)\le 5e6$ </br>
$1\le t(i)\le 1e7$ </br>

## 思路分析
补题的时候明确是贪心了，直接从贪心开始想题
考虑能安全跑过每座桥就行，计算出每座桥要跑步多久，走路多久
但是不能先跑步再走路，跑步的时间放在最后，零头可以用在下一座桥上
所以
先计算上座桥的零头能跑过的距离
再计算走路的距离
再计算跑步的距离，记录嗑药点
计算零头

## 核心代码
``` c++
#define _CRT_SECURE_NO_WARNINGS
#include<cstring>
#include<vector>
#include<cstdio>
#include<iostream>
#include<algorithm>
#include<map>
#include<queue>
#include<string>
#include<stack>
using namespace std;
typedef long long LL;
const int maxn = 200010;

LL n, r;
LL L[maxn], T[maxn];
vector<LL> ans;

int main() {
	scanf("%lld %lld", &n, &r);
	for (int i = 0; i<n; i++) scanf("%d", &L[i]);
	for (int i = 0; i<n; i++) scanf("%d", &T[i]);
	int fg = 1;
	for (int i = 0; i<n; i++) {
		if (T[i]<L[i]) {
			fg = 0;
			break;
		}
	}
	if (!fg) {
		puts("-1");
		return 0;
	}
	LL rm = 0; //剩余加速时间
	LL nans = 0; //加速点数目
	LL poi = 0; //地址
	LL tim = 0; //时间
	for (int i = 0; i<n; i++) {
		if (rm >= L[i]) { //全程冲刺
			rm -= L[i];
			poi += L[i];
			continue;
		}
		// 处理掉剩余的冲刺
		L[i] -= rm;
		T[i] -= rm;
		poi += rm;
		rm = 0;

		//慢慢走可以走完
		if (L[i] * 2 <= T[i]) {
			poi += L[i] * 2;
			continue;
		}

		LL t = T[i] - L[i]; //需要跑步的时间
		T[i] -= t*2;
		L[i] -= t;
		poi += 2 * t;
		LL tmp = poi;
		while (ans.size()<=100000 && tmp < poi + T[i]) {
			ans.push_back(tmp);
			tmp += r;
		}
		nans += L[i] / r;
		rm = L[i] % r;
		poi += L[i];
		if (rm) {
			nans++;
			rm = r - rm;
		}
	}
	printf("%lld\n", nans);
	if (nans <= 100000) {
		sort(ans.begin(), ans.end());
		for (int i = 0; i<nans; i++) {
			printf("%lld ", ans[i]);
		}puts("");
	}
}

```

