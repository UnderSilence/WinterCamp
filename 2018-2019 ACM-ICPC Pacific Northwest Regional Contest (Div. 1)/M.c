/*
gym101982m, 数学,优化，函数最值
*/
#include <stdio.h>
#include <stdlib.h>
#define max(a, b) (a) > (b) ? (a) : (b);
#define min(a, b) (a) < (b) ? (a) : (b);
#define N 30005
#define eps 1e-7
typedef long long LL;
struct point {
    double x, y;
};
inline struct point sub(struct point a, struct point b) {
    return (struct point){.x = a.x - b.x, .y = a.y - b.y};
}
inline double cross(struct point a, struct point b) {
    return a.x * b.y - a.y * b.x;
}
inline int sgn(double x) { return (x > eps) - (x < eps); }
int cmp(const void* a, const void* b) {
    struct point* p = (struct point*)a;
    struct point* q = (struct point*)b;
    if (!sgn(p->x - q->x)) {
        return sgn(p->y - q->y);
    } else
        return sgn(p->x - q->x);
}
double calc(struct point u, struct point v) {
    double A = (u.x - v.x) * (u.y - v.y);
    double B = (u.x - v.x) * v.y + (u.y - v.y) * v.x;
    double C = v.x * v.y;
    double x = -0.5 * B / A;
    if (sgn(x - 0.0) <= 0)
        return 0;
    else
        x = min(x, 1.0);
    return (A * x + B) * x + C;
}
int main() {
    struct point v[N], ch[N];
    double c[N], h[N], p[N], ans = 0.0;
    int n, C, m = 0;

    scanf("%d%d", &n, &C);
    for (int i = 0; i < n; i++) {
        scanf("%lf%lf%lf", &c[i], &h[i], &p[i]);
        v[i] = (struct point){.x = h[i] / c[i] * C, .y = p[i] / c[i] * C};
        ans = max(ans, v[i].x * v[i].y);
    }
    qsort(v, n, sizeof(v[0]), cmp);
    /*  for (int i = 0; i < n; i++) {
          printf("point %.3f %.3f\n", v[i].x, v[i].y);
      }*/
    for (int i = 0; i < n; i++) {
        while (m > 1 &&
               sgn(cross(sub(ch[m - 1], ch[m - 2]), sub(v[i], ch[m - 2]))) > 0)
            m--;
        ch[m++] = v[i];
    }
    /*  for (int i = 0; i < m; i++) {
          printf("convex hull %.3f %.3f\n", ch[i].x, ch[i].y);
      }*/
    for (int i = 0; i < m - 1; i++) {
        ans = max(ans, calc(ch[i], ch[i + 1]));
    }
    printf("%.2f\n", ans);
    return 0;
}
/*
4 100000
300 1 0.02
500 0.2 1
250 0.3 0.1
1000 1 0.1
*/