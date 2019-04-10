### 3. 拉格朗日插值法求幂

$O(N)$时间计算 [横坐标连续的点集]的拟合多项式 在$X$点的取值
由 $N+1$ 个点可以唯一确定一个 $N$次多项式
注意 需要初始化 $y[1..n]​$.
$$
\begin{equation}
P(x) = \sum_{i=1}^ny_i\prod_{j=1,j\ne i}^n\frac{x-x_j}{x_i-x_j}
\end {equation}
$$

用于解决类似求解 $f(k) = \sum_{i=a}^{b}i^k$的值

```C++
struct PolyInter {
	enum {MOD = 1e9+7,N = 1015}; 	/*类内常量定义MOD,N*/
	int y[N], k, finv[N];
	int l[N], r[N];
	void init(int n, int Y[]) {
		memcpy(y,Y,sizeof(int)*(n+2));
		k=n;
		finv[0]=finv[1]=1;
		for(int i=2; i<=k; i++) {
			finv[i]=1LL*(MOD-MOD/i)*finv[MOD%i]%MOD;
		}
		int tmp = 1;
		for(int i=2; i<=k; i++) {
			tmp = 1LL*tmp*i%MOD;
			finv[i]=1LL*finv[i]*finv[i-1]%MOD;
		}
	}
	int calc(long long x) {
		int ret=0,cur=0;
		l[0]=r[k+1]=1;
		for(int i=1; i<=k; i++) {
			l[i]=1LL*l[i-1]*(x-i+MOD)%MOD;
		}
		for(int i=k; i>0; i--) {
			r[i]=1LL*r[i+1]*(x-i+MOD)%MOD;
		}
		for(int i=1; i<=k; i++) {
			cur=1LL*l[i-1]*finv[i-1]%MOD;
			cur=1LL*cur*r[i+1]%MOD*finv[k-i]%MOD;
			if((k-i)%2) cur=1LL*cur*(MOD-1)%MOD;
			ret=(ret+1LL*cur*y[i]%MOD)%MOD;
		}
		return (ret%MOD+MOD)%MOD;
	}
} poly;
```

###二进制最后一个1的位置

```c++
int ffs(unsigned long long word) {
  int num = 0;
  if ((word & 0xffffffff) == 0) {
    num += 32;
    word >>= 32;
  }
  if ((word & 0xffff) == 0) {
    num += 16;
    word >>= 16;
  }
  if ((word & 0xff) == 0) {
    num += 8;
    word >>= 8;
  }
  if ((word & 0xf) == 0) {
    num += 4;
    word >>= 4;
  }
  if ((word & 0x3) == 0) {
    num += 2;
    word >>= 2;
  }
  if ((word & 0x1) == 0)
    num += 1;
  return num;
}
```

### 快速沃尔什变换

用于快速计算二进制卷积，复杂度 $O(N\log N)$

注意，数组的长度应为二的幂次$2^n$

```c++
void FWT(int a[], int n)
{
  for (int d = 1; d < n; d <<= 1)
    for (int m = d << 1, i = 0; i < n; i += m)
      for (int j = 0; j < d; j++)
      {
        int x = a[i + j], y = a[i + j + d];
        a[i + j] = (x + y) % MOD;
        a[i + j + d] = (x - y + MOD) % MOD;
        //xor:a[i+j]=x+y,a[i+j+d]=(x-y+MOD)%MOD;
        //and:a[i+j]=x+y;
        //or:a[i+j+d]=x+y;
      }
}
void UFWT(int a[], int n)
{
  for (int d = 1; d < n; d <<= 1)
    for (int m = d << 1, i = 0; i < n; i += m)
      for (int j = 0; j < d; j++)
      {
        int x = a[i + j], y = a[i + j + d];
        a[i + j] = 1LL * (x + y) * inv2 % MOD;
        a[i + j + d] = (1LL * (x - y) * inv2 % MOD + MOD) % MOD;
        //xor:a[i+j]=(x+y)/2,a[i+j+d]=(x-y)/2;
        //and:a[i+j]=x-y;
        //or:a[i+j+d]=y-x;
      }
}
//n取2的整数次幂
void solve(int a[],int b[],int n)  
{  
    FWT(a,n);  
    FWT(b,n);  
    for(int i=0;i<n;i++) a[i]=1LL*a[i]*b[i]%mod;  
    UFWT(a,n);  
} 
```

### 简易大整数

支持加、减、比较，注意只适用于正整数，大数减小数

```c++
struct BigInt { // 仅支持非负整数
  typedef long long LL;
  static const int N = 200;
  static const int base = 1e6; // 修改它时记得修改输入输出格式

  int a[N];
  int length;

  BigInt(): length(0) {memset(a, 0, sizeof(a));}

  BigInt(LL p) {
    memset(a, 0, sizeof(a));
    length = 0;
    if (!p) return;
    for (LL x = std::abs(p); x; x /= base) {
      a[length ++] = x % base;
    }
  }

  int &operator [](int sit) {return a[sit];}

  bool operator < (const BigInt &q)const {
    if (length != q.length) return length < q.length;
    for (int i = length - 1; ~i; -- i) {
      if (a[i] != q.a[i]) return a[i] < q.a[i];
    }
    return false;
  }

  BigInt operator + (const BigInt &p)const {
    BigInt ret;
    ret.length = std::max(length, p.length) + 1;
    for (int i = 0; i < ret.length; ++ i) {
      ret.a[i] += a[i] + p.a[i];
      if (ret.a[i] >= base) ret.a[i] -= base, ++ ret.a[i + 1];
    }
    for ( ; ret.length && !ret.a[ret.length - 1]; -- ret.length)
      ;
    return ret;
  }

  BigInt operator - (const BigInt &p)const {
    BigInt ret;
    ret.length = length;
    for (int i = 0; i < ret.length; ++ i) {
      ret.a[i] += a[i] - p.a[i];
      if (ret.a[i] < 0) ret.a[i] += base, -- ret.a[i + 1];
    }
    for ( ; ret.length && !ret.a[ret.length - 1]; -- ret.length)
      ;
    return ret;
  }

  BigInt operator * (const BigInt &p)const {
    static LL aux[N << 1];
    memset(aux, 0, sizeof(LL) * (length + p.length));
    for (int i = 0; i < length; ++ i) {
      for (int j = 0; j < p.length; ++ j) {
        aux[i + j] += (LL) a[i] * p.a[j];
      }
    }
    BigInt ret;
    ret.length = p.length + length;
    for (int i = 0; i < ret.length; ++ i) {
      aux[i + 1] += aux[i] / base;
      aux[i] %= base;
    }
    for ( ; ret.length && !aux[ret.length - 1]; -- ret.length)
      ;
    for (int i = 0; i < ret.length; ++ i) ret.a[i] = aux[i];
    return ret;
  }

  void write() {
    if (!length) return (void) printf("0");
    printf("%d", a[length - 1]);
    for (int i = length - 2; ~i; -- i) {
      printf("%06d", a[i]);
    }
  }
};
```
###生成函数

#### 正五边形定理优化整数拆分

##### 朴素整数拆分将n划分成不大于m的正整数之和

1)划分多个整数可以相同	 $dp[n][m] = dp[n][m-1] + dp[n-m][m]$

2)划分多个整数不可相同 	 $dp[n][m] = dp[n][m-1] + dp[n-m][m-1]$  

3)将n划分为k个数			 $dp[n][k] = dp[n-k][k] + dp[n-1][k-1]$

#### 生成函数表示整数拆分数

设整数拆分数生成函数P(x)
$$
\begin {equation}
P(x) = \sum_{n=0}^{\infin}p(n)x^n \\
P(x) = (1+x+x^2+\cdots ) (1+x^2+x^4+\cdots)\cdots \\
P(x)= \prod_{i=1}^{\infty} \frac{1}{(1-x^i)}
\end {equation}
$$

其中$x^m$项的系数即为$m$的拆分方案数，可使用五边形定理优化递推  $O(N\sqrt{N})$
$$
\begin{equation}
\phi(x)= \prod_{i=1}^{\infin}(1-x^i)= \sum_{-\infin}^{\infin}  (-1)^kx^{\frac {k(3k-1)}{2}} = 1+\sum_{k=1}^{\infin} (-1)^kx^{\frac{k(3k\pm1)}{2}}\\
P(x) = 1/\phi(x) \\

\end{equation}
$$
手算一下系数即可。若每一整数用的次数要小于k次，那么修改生成函数即可。

```c++
#include <bits/stdc++.h>
using namespace std;
typedef long long LL;
const int N = 100005;
const int MOD = 1e9 + 7;
inline calc(int x) {
  return (3 * x * x - x) / 2;
}
int exf[N], tot, n, k, p[N], T;
int main() {
  tot = 0;
  for (int i = 1; i < N; i++) {
    exf[tot++] = calc(i);
    exf[tot++] = calc(-i);
    if (exf[tot - 1] > N) break;
  }
  p[0] = 1;
  for (int i = 1; i < N; i++) {
    for (int j = 0, mu = 1; j < tot; j++) {
      if (exf[j] > i) break;
      p[i] = (p[i] + mu * p[i - exf[j]]) % MOD;
      p[i] = p[i] % MOD;
      if (j % 2 == 1) mu = -mu;
    }
  }
  scanf("%d", &T);
  while (T--) {
    scanf("%d%d", &n, &k);
    int ans = p[n];
    for (int i = 0, mu = -1; i < tot; i++) {
      if (exf[i] * k > n) break;
      ans += mu * p[n - exf[i] * k];
      ans %= MOD;
      if (i % 2 == 1) mu = -mu;
    }
    printf("%d\n", (ans % MOD + MOD) % MOD);
  }
  return 0;
}
```



### 一些特殊计数数列

#### Large Schröder numbers [A006318]  

> a(n) is the **number** of subdiagonal paths from (0, 0) to (n, n) consisting of steps East (1, 0), North (0, 1) and Northeast (1, 1) (sometimes called royal paths) 
>
> a(n) is the **number** of separable permutations, i.e.,	 permutations avoiding 2413 and 3142  
>
> ```c++
> 1, 4, 10, 33, 91, 298, 910, 3017, 9945, 34207, 119369, 429250, 1574224, 5916148, 22699830, 89003059, 356058540, 1453080087, 6044132794, 25612598436, 110503627621, 485161348047, 2166488899642, 9835209912767, 45370059225318
> ```


$$
\begin{equation}
a_n =  \sum_{k=0}^n {{n+k \choose n}{n \choose k}\over (k+1)} \\
a_n = \frac 1n\sum_{k=0}^n 2^k{n\choose k}{n \choose k-1}
\end{equation}
$$

#### little Schröder numbers [A001003] 

>in number theory, the Schröder–Hipparchus numbers form an integer sequence that can be used to count the number of plane trees with a given set of leaves, the number of ways of inserting parentheses into a sequence, and the number of ways of dissecting a convex polygon into smaller polygons by inserting diagonals. 
>
>```c++
>	1, 1, 3, 11, 45, 197, 903, 4279, 20793, 103049, 518859, 2646723, 13648869, 71039373, 372693519, 1968801519, 10463578353, 55909013009, 300159426963, 1618362158587, 8759309660445, 47574827600981, 259215937709463, 1416461675464871 
>```

$$
\begin{equation}
S(n)=\frac{1}{n}((6n-9)S(n-1)-(n-3)S(n-2))
\end{equation}
$$

### 树上倍增LCA，维护树链最值

```c++
/*
[倍增LCA]
维护树上树链u->v最大值
维护点权最大时注意还要加上lca的贡献*
*/
const int N = 100005;
const int LOG = 19;
typedef pair<int, int> Pii;
struct Edge {
  int u, v, w;
  Edge(int u = 0, int v = 0, int w = 0):
    u(u), v(v), w(w) {}
  bool operator < (const Edge& r) const {
    return w < r.w;
  }
};
vector<Edge>E, G[N];
int parent[LOG][N];
int mx[LOG][N];
int root;
int depth[N];

void dfs(int u, int fa, int d, int c) {
  parent[0][u] = fa;
  mx[0][u] = c;
  depth[u] = d;
  for (int i = 0; i < G[u].size(); i++) {
    int v = G[u][i].v;
    if (v == fa) continue;
    dfs(v, u, d + 1, G[u][i].w);
  }
}
void init(int n) { //G[1..n] 下标从1开始
  dfs(root, -1, 0, 0);
  for (int k = 0; k + 1 < LOG; k++) {
    for (int v = 1; v <= n; v++) {
      if (parent[k][v] < 0) {
        parent[k + 1][v] = -1;
        mx[k + 1][v] = mx[k][v];
      }
      else {
        mx[k + 1][v] = max(mx[k][v], mx[k][parent[k][v]]);
        parent[k + 1][v] = parent[k][parent[k][v]];
      }
    }
  }
}
//返回u与v的lca
int lca(int u, int v) {
  if (depth[u] > depth[v]) swap(u, v);
  for (int k = 0; k < LOG; k++) {
    if ((depth[v] - depth[u]) >> k & 1) {
      v = parent[k][v];
    }
  }
  if (u == v) return u;
  for (int k = LOG - 1; k >= 0; k--) {
    if (parent[k][u] != parent[k][v]) {
      u = parent[k][u];
      v = parent[k][v];
    }
  }
  return parent[0][u];
}
//返回u->v树链最值,要将边权下放到点权
//处理点权最值时要注意额外计算lca的贡献
int rmq(int u, int v) {
  int ret = 0;
  //根据情况设定初始值
  if (depth[u] > depth[v]) swap(u, v);
  for (int k = 0; k < LOG; k++) {
    if ((depth[v] - depth[u]) >> k & 1) {
      ret = max(ret, mx[k][v]);
      v = parent[k][v];
    }
  }
  if (u == v) return ret;
  for (int k = LOG - 1; k >= 0; k--) {
    if (parent[k][u] != parent[k][v]) {
      ret = max({ret, mx[k][u], mx[k][v]});
      u = parent[k][u];
      v = parent[k][v];
    }
  }
  return max({ret, mx[0][u], mx[0][v]});
}
```

### 树上欧拉序列(BFS序)+最近公共祖先

```c++
/*
HDU6393[欧拉序列+LCA+BIT]
环套树，任意去除环上一边
分类讨论图上任意两点最短路
单点到树根的距离转化为欧拉序列前缀和
利用树状数组的差分区间修改，单点求值
*/
#include <bits/stdc++.h>
using namespace std;
typedef long long LL;
const int N = 100005;
int head[N], to[N << 1], weight[N << 1], nxt[N << 1], tot;
/*链式前向星部分*/
int depth[N << 1], seq[N << 1], in[N], out[N], val[N], be[N], pos;
/*深度序列,欧拉序列,入栈时间,出栈时间,点权,边对应点编号*/
int st[N << 1][20], lg[N << 1];
/*ST表,log2表*/
LL bit[N << 1];
int T, n, q, X, Y, Z, W, x, y, z;
void add(int x, int k) {
  while (x <= pos) {
    bit[x] += k;
    x += x & -x;
  }
}
LL sum(int x) {
  LL ret = 0;
  while (x) {
    ret += bit[x];
    x -= x & -x;
  }
  return ret;
}
void clear() {
  for (int i = 0; i <= n; i++) {
    head[i] = -1;
    bit[i << 1] = bit[i << 1 | 1] = 0;
    in[i] = 0;
  } pos = tot = 0;
}
void add_edge(int u, int v, int w) {
  to[tot] = v;
  weight[tot] = w;
  nxt[tot] = head[u];
  head[u] = tot++;
}
void dfs(int u, int fa, int deep) {
  seq[++pos] = u;
  depth[pos] = deep;
  in[u] = pos;
  for (int i = head[u]; ~i; i = nxt[i]) {
    int v = to[i];
    if (v == fa) continue;
    if (in[v]) {
      X = u;
      Y = v;
      Z = i / 2 + 1;
      W = weight[i];
      continue;
    }
    be[i / 2 + 1] = v; /*绑定边权和点权，下放到叶子方向*/
    val[v] = weight[i];
    dfs(v, u, deep + 1);
    seq[++pos] = u;
    depth[pos] = deep;
  }
  out[u] = pos;
}
void init_rmq(int size) {
  lg[0] = -1;
  for (int i = 1; i <= size; i++) {
    lg[i] = (i & (i - 1)) == 0 ? lg[i - 1] + 1 : lg[i - 1];
  }
  for (int i = 1; i <= size; i++) {
    st[i][0] = i;
  }
  for (int j = 1; (1 << j) <= size; j++) {
    for (int i = 1; i + (1 << j) - 1 <= size; i++) {
      if (depth[st[i][j - 1]] < depth[st[i + (1 << (j - 1))][j - 1]]) {
        st[i][j] = st[i][j - 1];
      } else {
        st[i][j] = st[i + (1 << (j - 1))][j - 1];
      }
    }
  }
}
int lca(int l, int r) {
  l = in[l], r = in[r];
  if (l > r) swap(l, r);
  int k = lg[r - l + 1];
  if (depth[st[l][k]] < depth[st[r - (1 << k) + 1][k]]) {
    return seq[st[l][k]];
  } else {
    return seq[st[r - (1 << k) + 1][k]];
  }
}
LL dist(int a, int b) {
  return sum(in[a]) + sum(in[b]) - 2 * sum(in[lca(a, b)]);
}
int main() {
  scanf("%d", &T);
  while (T--) {
    scanf("%d%d", &n, &q);
    clear();
    for (int i = 0, u, v, w; i < n; i++) {
      scanf("%d%d%d", &u, &v, &w);
      add_edge(u, v, w);
      add_edge(v, u, w);
    }
    dfs(1, -1, 0);
    init_rmq(pos);
    for (int i = 1; i <= n; i++) {
      add(in[i], val[i]);
      add(out[i] + 1, -val[i]);
    }
    while (q--) {
      scanf("%d%d%d", &z, &x, &y);
      if (z) {
        LL ans = dist(x, y);
        ans = min(ans, dist(x, X) + dist(Y, y) + W);
        ans = min(ans, dist(x, Y) + dist(X, y) + W);
        printf("%lld\n", ans);
      } else {
        if (x == Z) {
          W = y;
        } else {
          int cur = be[x];
          add(in[cur], y - val[cur]);
          add(out[cur] + 1, val[cur] - y);
          val[cur] = y;
        }
      }
    }
  }
  return 0;
}
```

### 三维空间向量仿射变换

```c++
double mat[3][3];
// 绕过原点和（x，y，z）的轴，转 a 角度的 旋转矩阵，注意传的点向量为单位向量，a为弧度
// 方向为从(x,y,z) 看向 原点 逆时针
void rotate_mat(double x,double y,double z, double a) {
	double c=cos(a);
	double s=sin(a);
	double xy = x*y, xz = x*z, yz = y*z, xx = x*x, yy = y*y , zz = z*z ;
	mat[0][0] = xx + (1 - xx) * c; 
    mat[0][1] = xy * (1-c) - z * s; 
    mat[0][2] = xz *  (1-c) + y * s;
	mat[1][0] = xy * (1 - c) + z * s; 
    mat[1][1] = yy + (1 - yy) * c; 
    mat[1][2] = yz * (1-c) -x * s ;
	mat[2][0] = xz * (1 - c) - y * s; 
    mat[2][1] = yz * (1-c) + x * s; 
    mat[2][2] = zz +(1 - zz) * c;
}
```



#### 基本空间仿射变换

##### 平移变换

$$
\begin {pmatrix}
&1 &0 &0 &\Delta x\\
&0&1&0 &\Delta y \\
&0&0&1 &\Delta z\\
&0&0&0&1
\end{pmatrix} 
\begin{pmatrix}
x \\ y \\z\\1
\end{pmatrix}
$$

##### 缩放变换

$$
\begin {pmatrix}
&kx&0&0&0\\
&0&ky&0&0 \\
&0&0&kz&0\\
&0&0&0&1
\end{pmatrix} 
\begin{pmatrix}
x\\y\\z\\1
\end{pmatrix}
$$

#### 旋转变换

$$
rotate(x,y,z,d) \\
\begin {pmatrix}
&(1-\cos(d))x^2+\cos(d)&(1-\cos(d))xy-\sin(d)z& (1-\cos(d))xz+\sin(d)y&0\\
&(1-\cos(d))yx+\sin(d)z&(1-\cos(d))y^2+\cos(d)& (1-\cos(d))yz-\sin(d)x&0 \\
&(1-\cos(d))zx-\sin(d)y  & (1-\cos(d))zy+\sin(d)x & (1-cos(d))z^2+\cos(d)&0\\
&0&0&0&1
\end{pmatrix} 
\begin{pmatrix}
x\\y\\z\\1
\end{pmatrix}
$$





### 计算几何整理

#### 单位圆覆盖最多的点

##### 用大小R为圆去覆盖一个点集使得被包含的点最多。

```c++
/*用大小R为圆去覆盖一个点集使得被包含的点最多。*/
int rac_cirtocir(Circle c1, Circle c2, pair<db, db>& rac) {
  double d = len(c1.o - c2.o);
  if (!sgn(d)) {
    if (!sgn(c1.r - c2.r)) return -1;
    return 0;
  }
  if (sgn(c1.r + c2.r - d) < 0) return 0;
  if (sgn(fabs(c1.r - c2.r) - d) > 0) return 0;
  db a = angle(c2.o - c1.o);
  db da = acos((c1.r * c1.r + d * d - c2.r * c2.r) / (2 * c1.r * d));
  rac.first = a - da;
  rac.second = a + da;
  return 1;
}
int MaxCircleCover(double R, Point *p, int n) {
  int res = 1;
  for (int i = 0; i < n; i++) {
    int sum = 1;
    Circle now(p[i], R);
    vector<pair<db, bool> > mp;
    pair<db, db> rac;
    for (int j = 0; j < n; j++) {
      if (i == j) continue;
      rac = make_pair(0, 2e18);
      if (rac_cirtocir(now, Circle(p[j], R), rac)) {
        mp.push_back(make_pair(rac.first, false));
        mp.push_back(make_pair(rac.second, true));
      }
    }
    sort(mp.begin(), mp.end());
    for (int i = 0; i < mp.size(); i++) {
      if (!mp[i].second) sum++;
      else sum--;
      res = max(sum, res);
    }
  }
  return res;
}
```



#### 线段求交[可支持整点判断]

```c++
using LL = long long;
bool is_Seginterseg(Line A, Line B) {
  return
    max(A.s.x, A.e.x) >= min(B.s.x, B.e.x) &&
    max(A.s.y, A.e.y) >= min(B.s.y, B.e.y) &&
    min(A.s.x, A.e.x) <= max(B.s.x, B.e.x) &&
    min(A.s.y, A.e.y) <= max(B.s.y, B.e.y) &&
    sgn(cross(B.s - B.e, A.s - A.e)) * sgn(cross(B.e-A.e, A.s - A.e)) <= 0 &&
    sgn(cross(A.s - B.e, B.s - B.e)) * sgn(cross(A.e-B.e, B.s - B.e)) <= 0;
}
bool is_Pointonseg(Point p, Line l) {
  return sgn(cross(l.s - p, l.e - p)) == 0 && sgn(dot(p - l.e, p - l.s)) <= 0;
}
bool cross_Linetoline(Line u, Line v, Point& P) {
  if (!is_Seginterseg(u, v)) return false;
  LL p = cross(v.e - v.s, u.s - v.s);
  LL q = cross(v.e - v.s, u.e - v.s);
  P = Point(
        (u.s.x * q - u.e.x * p) / (q - p),
        (u.s.y * q - u.e.y * p) / (q - p)
      );
// printf(" cross at (%lld, %lld) ck: %s\n ", P.x, P.y, is_Pointonseg(P, u) && is_Pointonseg(P, v)?"Yes":"No");
  return is_Pointonseg(P, u) && is_Pointonseg(P, v);
}
```

#### 凸包缩小 shrink_poly

```c++
#include <bits/stdc++.h>
using namespace std;
const int N = 1005;
const double eps = 1e-7;
int T, n, r;
int sgn(double x) {
  return (x > eps) - (x < -eps);
}
struct Point {
  double x, y;
  Point(double x = 0, double y = 0):
    x(x), y(y) {}
  Point operator - (const Point& r) const {
    return Point(x - r.x, y - r.y);
  }
  Point operator + (const Point& r) const {
    return Point(x + r.x, y + r.y);
  }
  Point operator * (const double& r) const {
    return Point(x * r, y * r);
  }
  Point operator / (const double& r) const {
    return Point(x / r, y / r);
  }
  bool operator < (const Point& r) const {
    return sgn(x - r.x) == 0 ? sgn(y - r.y) < 0 : x < r.x;
  }
  double len() {
    return sqrt(x * x + y * y);
  }
  Point norm() {
    return Point(-y, x) / len();
  }
};
Point p[N], ch[N];
typedef Point Vector;
double cross(Vector u, Vector v) {
  return u.x * v.y - u.y * v.x;
}
double dot(Vector u, Vector v) {
  return u.x * v.x + u.y * v.y;
}
struct Line {
  Point s, e;
  double ang;
  Line() {}
  Line(Point s, Point e): s(s), e(e) {
    ang = atan2((s - e).y, (s - e).x);
  }
  bool operator < (const Line& l) const {
    if (sgn(ang - l.ang) != 0) return ang < l.ang;
    return cross(s - l.s, l.e-l.s) < 0;
  }
};
bool on_left(Line l, Point p) {
  return cross(l.e-l.s, p - l.s) >= 0;
}
Point calc_intersection(Line A, Line B) {
  Point P = A.s;
  double t = cross(B.e-B.s, A.s - B.s) / cross(A.e-A.s, B.e-B.s);
  return  P + (A.e-A.s) * t;
}
int half_plane_intersection(Line l[], int n, Point poly[]) {
  sort(l, l + n);
  int first, last;
  Point *p = new Point[n];
  Line *q = new Line[n];
  q[first = last = 0] = l[0];
  for (int i = 1; i < n; i++) {
    while (first < last && !on_left(l[i], p[last - 1])) last--;
    while (first < last && !on_left(l[i], p[first])) first++;
    q[++last] = l[i];
    if (sgn(cross(q[last].e-q[last].s, q[last - 1].e-q[last - 1].s)) == 0) {
      last--;
      if (on_left(q[last], l[i].s)) q[last] = l[i];
    }
    if (first < last) {
      p[last - 1] = calc_intersection(q[last - 1], q[last]);
    }
  }
  while (first < last && !on_left(q[first], p[last - 1])) last--;
  if (last - first <= 1) return 0;
  p[last] = calc_intersection(q[last], q[first]);
  int m = 0;
  for (int i = first; i <= last; i++) poly[m++] = p[i];
  delete []p;
  delete []q;
  return m;
}
int shrink_poly(Point* p, int n, double r, Point* np) {
  Line *l = new Line[n];
  for (int i = 0; i < n; i++) {
    Vector offset = (p[(i + 1) % n] - p[i]).norm() * r;
    l[i] = Line(p[i] + offset, p[(i + 1) % n] + offset);
  }
  int ret = half_plane_intersection(l, n, np);
  delete[] l;
  return ret;
}
/*注意方向统一为逆时针*/
```

### 大整数[1e18]开根

```c++
const LL INF = 1e18 + 300;
const LL TT = (LL)1 << 31;
LL multi(LL a, LL b) {
  LL ans = 1;
  while (b) {
    if (b & 1) {
      double judge = 1.0 * INF / ans;
      if (a > judge) return -1;
      ans *= a;
    }
    b >>= 1;
    if (a > TT && b > 0) return -1;
    a = a * a;
  }
  return ans;
}
LL find(LL x, LL k) {
  LL r = (LL)pow(x, 1.0 / k);
  LL t, p;
  p = multi(r, k);
  if (p == x) return r;
  else if (p > x || p == -1) {
    do {
      r--;
      p = multi(r, k);
    } while (p > x || p == -1);
  }
  else {
    t = multi(r + 1, k);
    while (t != -1 && t <= x) {
      r++;
      t = multi(r + 1, k);
    }
  }
  return r;
}
```



###  第k短路(Dijkstra+A*近似算法)

```c++
#include <bits/stdc++.h>
using namespace std;
const int N = 1005;
const int inf = 0x3f3f3f3f;
typedef long long LL;
struct Edge {
  int v, w;
  Edge(int v = 0, int w = 0):
    v(v), w(w) {}
};
struct Path {
  int u;
  LL w;
  Path(int u = 0, LL w = 0):
    u(u), w(w) {}
  bool operator < (const Path& p) const {
    return w > p.w;
  }
};
struct Node {
  LL g, h;
  int v;
  Node(LL g = 0, LL h = 0, int v = 0): g(g), h(h), v(v) {}
  bool operator < (const Node& r) const {
    return g + h > r.g + r.h;
  }
};
vector<Edge> G[N], rG[N];
int h[N];
int n, m, s, e, k, t;
void dijkstra(int src) {
  memset(h, 0x3f, sizeof(int) * (n + 1));
  h[src] = 0;
  priority_queue<Path> Q;
  Q.push(Path(src, 0LL));
  while (!Q.empty()) {
    Path p = Q.top();
    Q.pop();
    int u = p.u;
    LL w = p.w;
    if (w > h[u]) continue;
    for (int i = 0; i < (int)rG[u].size(); i++) {
      int v = rG[u][i].v;
      if (h[v] > h[u] + rG[u][i].w) {
        h[v] = h[u] + rG[u][i].w;
        Q.push(Path(v, h[v]));
      }
    }
  }
}
LL kth_shortest() {
  if (s == t) k++;
  if (h[s] == inf) return inf;
  priority_queue<Node> Q;
  Q.push(Node(0, h[s], s));
  while (!Q.empty()) {
    Node now = Q.top();
    Q.pop();
    LL g = now.g;
    int u = now.v;
    if (u == e) {
      if (k > 1) k--;
      else return g;
    }
    for (int i = 0; i < (int)G[u].size(); i++) {
      int v = G[u][i].v;
      int w = G[u][i].w;
      Q.push(Node(g + w, h[v], v));
    }
  }
  return inf;
}
int main() {
  while (~scanf("%d%d", &n, &m)) {
    for (int i = 0; i <= n; i++) {
      G[i].clear();
      rG[i].clear();
    }
    scanf("%d%d%d%d", &s, &e, &k, &t);
    for (int i = 0, u, v, w; i < m; i++) {
      scanf("%d%d%d", &u, &v, &w);
      G[u].push_back(Edge(v, w));
      rG[v].push_back(Edge(u, w));
    }
    dijkstra(e); //反向图跑dij预处理估价函数
    LL kth = kth_shortest();
    if (kth - t <= 0) {
      puts("yareyaredawa");
    } else puts("Whitesnake!");
  }
  return 0;
}
```

### MillerRabin黑科技

```c++
#include <bits/stdc++.h>
using namespace std;
typedef long long LL;
const int maxp = 1e6 + 1, maxv = 25, maxc = (int)1e4 + 1;
int ptot, pr[maxp], d[maxp], cnt;
LL n, p[maxc];
LL mod_add(LL x, LL y, LL p) {
	return (x += y) < p ? x : x - p;
}
LL mod_mul(LL x, LL y, LL p) {
	LL ret = x * y - (LL)((long double)x * y / p + 0.5) * p;
	return ret < 0 ? ret + p : ret;
}
LL mod_pow(LL x, LL k, LL p) {
	LL ret = 1 % p;
	for( ; k > 0; k >>= 1, x = mod_mul(x, x, p))
		(k & 1) && (ret = mod_mul(ret, x, p));
	return ret;
}
bool miller_rabin(LL n) {
	if(n == 2) return 1;
	if(n < 2 || !(n & 1))
		return 0;
	LL s = 0, r = n - 1;
	for( ; !(r & 1); r >>= 1, ++s);
	for(int i = 0; pr[i] < n && pr[i] < maxv; ++i) {
		LL cur = mod_pow(pr[i], r, n), nxt;
		for(int j = 0; j < s; ++j) {
			nxt = mod_mul(cur, cur, n);
			if(nxt == 1 && cur != 1 && cur != n - 1) return 0;
			cur = nxt;
		}
		if(cur != 1) return 0;
	}
	return 1;
}
LL gcd(LL a, LL b) {
	int ret = 0;
	while(a) {
		for( ; !(a & 1) && !(b & 1); ++ret, a >>= 1, b >>= 1);
		for( ; !(a & 1); a >>= 1);
		for( ; !(b & 1); b >>= 1);
		if(a < b)
			swap(a, b);
		a -= b;
	}
	return b << ret;
}
LL pollard_rho(LL n) {
	static LL seq[maxp];
	while(1) {
		LL x = rand() % n, y = x, c = rand() % n;
		LL *px = seq, *py = seq, tim = 0, prd = 1;
		while(1) {
			*py++ = y = mod_add(mod_mul(y, y, n), c, n);
			*py++ = y = mod_add(mod_mul(y, y, n), c, n);
			if((x = *px++) == y) break;
			LL tmp = prd;
			prd = mod_mul(prd, abs(y - x), n);
			if(!prd) return gcd(tmp, n);
			if((++tim) == maxv) {
				if((prd = gcd(prd, n)) > 1 && prd < n) return prd;
				tim = 0;
			}
		}
		if(tim && (prd = gcd(prd, n)) > 1 && prd < n) return prd;
	}
}
void decompose(LL n) {
	for(int i = 0; i < cnt; ++i)
		if(n % p[i] == 0) {
			p[cnt++] = p[i];
			n /= p[i];
		}
	if(n < maxp) {
		for( ; n > 1; p[cnt++] = d[n], n /= d[n]);
	} else if(miller_rabin(n)) {
		p[cnt++] = n;
	} else {
		LL fact = pollard_rho(n);
		decompose(fact), decompose(n / fact);
	}
} // prepare pr(prime) and d(minimal factor)
```

###Zeller 公式

```c++
int zeller(int y, int m, int d) {
	if (m<=2) y--,m+=12; int c=y/100; y%=100;
	int w=((c>>2)-(c<<1)+y+(y>>2)+(13*(m+1)/5)+d-1)%7;
	if (w<0) w+=7; return(w);
}
```

### 子集反演

```c++

for (int i = 0; i < n; i++)
  for (int s = 0; s < (1 << n); i++)
    if (s >> i & 1) f[s] += f[s ^ 1 << i];
for (int i = 0; i < n; i++)
  for (int s = 0; s < (1 << n); s++)
    if (s >> i & 1) f[s] -= f[s ^ 1 << i];
```

### 左偏堆，可合并堆

插入、合并、删除$\log(N)$

```c++
const int N = 500005;
int l[N],r[N],d[N],v[N],heap[N],tot,
int merge(int x,int y) {
  if (!x) return y;
  if (!y) return x;
  if (v[x]>v[y]) swap(x,y);
  r[x]=merge(r[x],y);
  if (d[l[x]]<d[r[x]]) {
    swap(l[x],r[x]);
  }
  d[x]=d[r[x]]+1;
  return x;
}
int init(int key) {
  v[++tot]=key;
  l[tot]=r[tot]=d[tot]=0;
  return tot;
}
int insert(int x,int y) {
  return merge(x,init(y));
}
int pop(int x) {
  return merge(l[x],r[x]);
}
int top(int x) {
  return v[x];
}
bool empty(int x) {
  return x==0;
}
```

