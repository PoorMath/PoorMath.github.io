## 用途

筛积性函数$f(x)$的前缀和（某些不是积性函数，但可以按质因子拆贡献的$f(x)$也可筛）。

要求

①$f(p^k)$​​能够快速计算。

②能够构造出一个完全积性函数$h(x)$​​​，使$f(p)=h(p)$​​​。

③质因子对函数值的贡献能够快速计算。

（若无特殊说明，$p$​都代表质数）。

## 怎么做

分两步：

1、对所有$x=\lfloor\frac{n}{i}\rfloor$筛出$\le x$的所有质数的$f$值之和。

2、对所有$x=\lfloor\frac{n}{i}\rfloor$求出$\le x$的所有数的$f$值之和。

## 第一步

设$prime_i$​​表示从小到大第$i$​​​个质数（特别地，规定$prime_0=1$），$minp(i)$表示$i$​​的最小质因子。

在这一步我们需要用到性质②：构造一个完全积性函数$h(x)$，使$h(p)=f(p)$。具体原因后面会讲。

设$g(x,j)=\sum_{i=2}^x[i\in Prime\ or\ minp(i)>prime_j]h(i)$​​​​。

即统计质数和最小质因子大于$prime_j$​​​​​​所有的数的$h$​​​​​​​值之和（注意$prime_0=1$，所以$g$是不包括$1$的）。感性理解一下就是在做一个埃氏筛法，第$j$轮就用$prime_j$筛去剩下的数中最小质因子为$prime_j$的数，$g(x,j)$就是$1$到$x$这些数在筛了$j$轮后，剩下的数的$h$​值之和。

最后我们要求的是所有质数的$f$​值，注意$h(p)=f(p)$​，所以为$g(x,tot)$​，$tot$​为$[1,x]$​范围内质数个数（对每个$x=\lfloor\frac{n}{i}\rfloor$​范围内的质数个数一般会在min_25筛的过程中顺便求出）。在这一步中，我们只是为了得到质数的$f$​值之和，合数的$f$​无关紧要，而且在最后会被筛掉，因此为了有快速的计算方法，即使$f($合数$)\neq h($合数$)$，我们也用$h$​计算。

因为所有合数的最小质因子都$\le\sqrt n$​，所以我们要先筛出$\sqrt n$​范围内的质数。

初始值$g(x,0)=\sum_{i=2}^xh(i)$​。

从小到大枚举$j$，考虑如何转移。

若$prime_j^2>x$，则我们在第$j$轮不会筛掉任何额外的数，所以$g(x,j)=g(x,j-1)$。

若$prime_j^2\le x$​​，最小质因子为$prime_j$​的数会被筛掉，考虑减去$h(prime_j)\times g(x/prime_j,j-1)$​，注意$h$是完全积性函数，所以不需要除掉所有$prime_j$这个质因子，但因为$g$中还包含了所有质数，所以这样多减去了一些最小质因子小于$prime_j$的值，所以还要加上这一部分。这一部分的方程为$g(x,j)=g(x,j-1)-h(prime_j)\times g(x/prime_j,j-1)+h(prime_j)\sum_{i=1}^{j-1}h(prime_i)$​。

##### 第一步具体实现

以求区间质数个数为例，直接上代码：

```C++
int gid(LL x)
{
	if(x<=sn)return id1[x];
	return id2[n/x];
}
void init()
{
	tot=lw=0;sn=(LL)sqrt(n);
	for(int i=0;i<=sn;i++)mark[i]=false;
	for(int i=2;i<=sn;i++)
	{
		if(!mark[i])prime[++tot]=i;
		for(int j=1;prime[j]*i<=sn&&j<=tot;j++)
		{
			mark[prime[j]*i]=true;
			if(i%prime[j]==0)break;
		}
	}
	LL p;
	for(LL i=1;i<=n;i=p+1)
	{
		w[++lw]=n/i;p=n/(n/i);
		if(n/i<=sn)id1[n/i]=lw;else id2[n/(n/i)]=lw;
		g[lw]=w[lw]-1;
	}
	for(int j=1;j<=tot;j++)
	for(int i=1;i<=lw;i++)
    {
        if((LL)prime[j]*prime[j]>w[i])break;
        int id=gid(w[i]/prime[j]);g[i]=g[i]-g[id]+j-1;
	}
}
```

几个要点：

1、注意我们只需要求所有$x=\lfloor\frac{n}{i}\rfloor$的$g(x)$​的取值，只有根号级别个值（数组要开到$2\sqrt n$），我们的$x$不用具体的值而用一个编号代替。储存编号的时候用到一个常用的小技巧，大于和小于根号的分别放到两个数组存，这样空间是$O(\sqrt n)$的。

2、我们只需要用到$g(x,tot)$，所以只用开一维数组，类似背包DP，每一轮从大到小转移即可。

3、最后求$[1,x]$​内质数个数，就是$g[gid(x)]$​。

这一部分的复杂度是$O(\frac{n^{\frac{3}{4}}}{\log n})$，注意如果筛的函数就是完全积性函数（如质数个数、质数和）那么只要这一步就足够了。

## 第二步

设$S(x,j)=\sum_{i=2}^x[minp(i)\ge prime_j]f(i)$。

即所有最小质因子$\ge prime_j$的数的$f$值之和。

考虑将这个拆成质数和合数两部分来算。

质数的部分比较简单，即$g(x,j)-\sum_{i=1}^{j-1}f(prime_i)$​。

合数部分，考虑枚举最小质因子及其次数，可以得到$\sum_{k=j}^{prime_k^2\le x}\sum_{e=1}^{prime_k^{e+1}\le x}[S(x/prime_k^e,j+1)\times f(prime_k^e)+f(prime_k^{e+1})]$​​，正确性显然。

合起来就是$S(x,j)=g(x,j)-\sum_{i=1}^{j-1}f(prime_i)+\sum_{k=j}^{prime_k^2\le x}\sum_{e=1}^{prime_k^{e+1}\le x}[S(x/prime_k^e,j+1)\times f(prime_k^e)+f(prime_k^{e+1})]$​。

##### 第二步具体实现

```C++
LL S(LL n,int m)
{
	if(prime[m]>n)return 0;
	LL re;//re初始值设为质数部分贡献
	for(int j=m;j<=tot&&prime[j]*prime[j]<=n;j++)
	{
		LL t=prime[j];
		for(int k=1;t*prime[j]<=n;k++)
		{
            //计算合数部分贡献
			t*=prime[j];
		}
	}
	return re;
}
```

程序中的变量名和上面公式中的略有不同。

为什么非积性函数也能筛呢？我们用到积性函数的性质其实是$S(x/prime_k^e,j+1)\times f(prime_k^e)$这个地方，实际上，只要能快速由$S(x/prime_k^e,j+1)$转移到$S(x,j)$的函数就能筛，对着模板改一下就好了。

注意我们始终没有考虑$1$，所以$\sum_{i=1}^nf(i)=S(n,1)+f(1)$。

这一部分复杂度也是$O(\frac{n^{\frac{3}{4}}}{\log n})$​。

但如果计算出每个$x=\lfloor\frac{n}{i}\rfloor$处的值，直接用递推写，复杂度貌似是$O(n^{1-\epsilon})$​的，亚线性复杂度（这个不太了解，可能有误）。

##### update：补充递推的写法（算出在所有$\lfloor\frac{n}{i}\rfloor$处函数前缀和）

设$S(x,j)=\sum_{i=2}^x[minp(i)\ge prime_j\ or\ i\in Prime]f(i)$​（注意定义略有不同）。

显然$S(x,tot+1)$的值就是我们第一步筛出来的$g(x)$（只包含了质数）。

考虑$S(x,j)$​比$S(x,j+1)$​多了哪些数的贡献，类似地，我们将其分为两类，一类是至少有两个不同的质因子，且最小质因子为$prime_j$​，一类是$prime_j^k$​，其中$k\ge2$​。所以式子为$S(x,j)=S(x,j+1)+\sum_{e=1}^{prime_j^{e+1}\le x}\{[S(x/{prime_j^e},j+1)-\sum_{k=1}^{j-1}f(prime_k)]\times f(prime_j^e)+f(prime_j^{e+1})\}$。

##### 具体实现（以筛$\phi$为例）

```C++
for(int j=tot;j;j--)
{
    for(int i=1;i<=lw;i++)
    {
        if((LL)prime[j]*prime[j]>w[i])break;
        int phi=prime[j]-1;LL t=prime[j];
        while(t*prime[j]<=w[i])
        {
            up(sphi[i],(LL)phi*(sphi[gid(w[i]/t)]-sprime[j]+j)%P);//至少两个不同质因子合数
            phi=(ll)phi*prime[j]%P;
            up(sphi[i],phi);//质数的幂
            t*=(ll)prime[j];
        }
    }
}
```

## 小结

做了几个题后，感觉min_25筛是个很灵活的东西，很多跟质因数有关的题都可以用这个做。

~~完结撒花~~

