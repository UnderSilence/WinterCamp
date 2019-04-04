# 标准记录格式
## e.g.[codeforcesXXXXC][模拟]a_plus_b_problem
## 题目描述
输出a+b的值。</br>
$0\le a,b\le 1e18$ </br>
(简单描述即可，起码要自己能看懂，也可以直接把题面粘过来，题目链接也可)
## 思路分析
直接输出a+b, 小心溢出。
## 核心代码
``` c++
#include <iostream>
using namespace std;
int main() {
  long long a, b;
  cin >> a >> b;
  cout << a + b << endl;
  return 0;
}
```

