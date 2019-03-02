/*
gym101982F 线段树，离散化，扫描线模板题，注意前提端点，线段树维护区间
*/
#include <bits/stdc++.h>
using namespace std;
typedef long long LL;
const int N = 200005;
const int inf = 0x3f3f3f3f;
int n, x[N], tot;
int len[N << 2], sxor[N << 2], tag[N << 2];
void push_down(int rt) {
    if (tag[rt]) {
        tag[rt << 1] ^= 1;
        tag[rt << 1 | 1] ^= 1;
        sxor[rt << 1] = len[rt << 1] - sxor[rt << 1];
        sxor[rt << 1 | 1] = len[rt << 1 | 1] - sxor[rt << 1 | 1];
        tag[rt] = 0;
    }
}
void push_up(int rt) { sxor[rt] = sxor[rt << 1] + sxor[rt << 1 | 1]; }
void build(int l, int r, int rt) {
    len[rt] = (x[r + 1] - x[l]);
    //cout << l << ' ' << r << ' ' << len[rt] << endl;
    sxor[rt] = 0;
    tag[rt] = 0;
    if (l == r) return;
    int mid = (l + r) / 2;
    build(l, mid, rt << 1);
    build(mid + 1, r, rt << 1 | 1);
}
void flip(int ql, int qr, int l, int r, int rt) {
    if (ql <= l && r <= qr) {
        tag[rt] ^= 1;
        sxor[rt] = len[rt] - sxor[rt];
        return;
    }
    push_down(rt);
    int mid = (l + r) / 2;
    if (ql <= mid) flip(ql, qr, l, mid, rt << 1);
    if (qr > mid) flip(ql, qr, mid + 1, r, rt << 1 | 1);
    push_up(rt);
}
struct Line {
    int x1, x2, y;
    Line(int x1 = 0, int x2 = 0, int y = 0) : x1(x1), x2(x2), y(y) {}
};
vector<Line> l;
int main() {
    scanf("%d", &n);
    for (int i = 0, x1, y1, x2, y2; i < n; i++) {
        scanf("%d%d%d%d", &x1, &y1, &x2, &y2);
        l.emplace_back(x1, x2, y1);
        l.emplace_back(x1, x2, y2);
        x[tot++] = x1;
        x[tot++] = x2;
    }
    x[tot++] = -inf;
    x[tot++] = inf;
    sort(x, x + tot);
    tot = unique(x, x + tot) - x;
    build(1, tot - 2, 1);
    sort(l.begin(), l.end(), [&](auto a, auto b) { return a.y < b.y; });
    LL ans = 0, last;
    for (int i = 0; i < l.size(); i++) {
        //cout << sxor[1] << "*(" << l[i].y << "-" << last << ")\n";
        ans += i ? 1LL * (l[i].y - last) * sxor[1] : 0;
        if (l[i].x1 > l[i].x2) swap(l[i].x1, l[i].x2);
        int L = lower_bound(x, x + tot, l[i].x1) - x;
        int R = lower_bound(x, x + tot, l[i].x2) - x;
        flip(L, R - 1, 1, tot - 2, 1);
        last = l[i].y;
    }
    printf("%lld\n", ans);
    return 0;
}
/*
2
0 0 3 3
1 1 3 3
4
0 0 10 10
1 1 11 11
2 2 12 12
3 3 13 13
*/