# [20210912 2021 BUAA Autumn Training 1](https://codeforces.com/gym/344373)

| 排名 | 当场过题数 | 至今过题数 | 总题数 |
| ---- | ---------- | ---------- | ------ |
| 2/3  | 8          | 9          | 10     |

## **A**

**solved by JLK**

### 题意

签到，略

### 题解



## **B**

**solved by JLK**

### 题意

签到，略

### 题解



## **C**

**upsolved by TYB**

### 题意

给一个字符串$S$，从包含其所有回文子串的集合中选出两个不同子串，使得一个是另外一个的子串，求这样的回文子串对数。

$|S|\le10^5$

### 题解

建出回文自动机，并求出回文自动机上每个节点表示的回文串第一次出现的区间，对每个这样的区间询问该区间内本质不同的回文串有多少个，那么所有询问的和即为本题答案。区间本质不同回文串个数是经典问题，~~但是我还不会~~，直接套模板即可。

## **D**

**solved by JLK**

### 题意

$n$个物品要放到至多$K$个容积相同的容器里。策略是：每次找一个最大能放进容器的物品放进去，直到放不进之后放下一个容器。求最小的容器容积。

$1 \le n,K,v_i \le 1000$

### 题解

（乱搞）二分答案。对于一个容积，可以直接$O(n\log n)$模拟判断。显然答案在大体上满足单调性，但中间一部分会有01交错。再观察可以发现这部分交错不会很多。先二分到一个最小1的位置，然后暴力向前找比$n$​多一点的位置，就可以找到了。

正经做法是，可以推出答案的下界为$\lceil \frac {sum}k\rceil$，上界为$\lceil \frac {sum}k\rceil+maxV$，直接枚举判断即可。

$O(nv\log n)$​

## **E**

**solved by YZW**

### 题意



### 题解



## **F**

**upsolved by **

### 题意



### 题解



## **G**

**solved by YZW**

### 题意



### 题解



## **H**

**solved by TYB**

### 题意

给出$n$个点$m$条边的无向连通图，边权全为$1$，给出点集$A,B$，从$A$中等概率随机选出一个点$a$，从$B$中等概率随机选出一个点$b$，从$n$个点中等概率随机选出一个点$c$，再根据这三个点选出一个点$d$，使得$dis(a,d)+dis(b,d)+dis(c,d)$最小，其中$dis(x,y)$为$x$到$y$的最短路长度，求$dis(a,d)+dis(b,d)+dis(c,d)$​的期望。

$|A|,|B|\le20,n,m\le10^5$

### 题解

所有情况概率相等，只需要考虑总和。$|A|,|B|$较小，考虑直接枚举其中的点$a,b$。分别以$a,b$为起点跑最短路，然后把每个点的初始权值设为该点到$a$的最短路和到$b$的最短路之和，再跑一次多源最短路，这次得到每个点$c$的最短路即为$dis(a,d)+dis(b,d)+dis(c,d)$。正确性显然。由于边权为$1$，只需要bfs即可，注意最后一次跑多源最短路时需要用vector或桶排序之类的方法对初始权值排序，复杂度才正确。

复杂度$\mathcal{O}(|A||B|n)$。

## **I**

**solved by YZW**

### 题意



### 题解



## **J**

**solved by TYB**

### 题意

签到。

### 题解

略。

## **记录**

开局JLK签A(8)B(35)。

看CD，感觉D可以二分，然后WA了好多次。

YZW写E，WA一次后AC(73)。

JLK讲了H做法，TYB开写，WA一次后AC(144)。

然后YZW写G，WA三次后AC(189)。

发现J是签到，TYB写，WA一次后AC(206)。

然后YZW开始写I，JLK和TYB想CD。

TYB觉得C可以抄板子，决定最后乱搞一下D。

YZW WA一次后AC I(277)。

TYB开始抄板子。

时间快到了，JLK试了几次D，然后AC(291)。

TYB抄完了，WA5次最终也没AC。

## **总结**

JLK：缺乏在榜歪的时候发现可做题的能力。

TYB：被多case坑了，要注意输出格式和初始化。

## **Dirt**

D(-3)

E(-1)

G(-3)

H(-1)

I(-1)

J(-1)
