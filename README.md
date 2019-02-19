# wintercamp
recording problems
记录使用标准markdown语法


# 标准记录格式
[题目oj+题号][知识点]题目标题(*可省略标题)

## 题目描述
blablablablabla
## 思路分析(套路分析)
blablablablabla
## 核心代码
``` c++
#include <iostream>
using namespace std;
int main() {
  return 0;
}
```
e.g.

[codeforcesXXXXC][模拟]a_plus_b_problem
## 题目描述
输出a+b的值。\\
$0\le a,b\le 1e18$
(简单描述即可，起码要自己能看懂，也可以直接把题面粘过来)
## 思路分析
直接输出a+b, 小心溢出。
## 核心代码
``` c++
#include <iostream>
using namespace std;
int main() {
  int a, b;
  cin >> a >> b;
  cout << a + b << endl;
  return 0;
}
```

