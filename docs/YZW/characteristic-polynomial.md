### 特征多项式

对于 $n$ 阶矩阵 $A$，其特征多项式定义为 $\varphi_A(\lambda)=\det(\lambda I-A)$，易见 $\varphi_A(\lambda)$ 为关于 $\lambda$ 的 $n$ 次首一多项式。

给定 $A$，计算 $\varphi_A(\lambda)$。



### Lagrange 插值

由于 $\varphi_A(\lambda)$ 为 $n$ 次多项式，我们只需要 $n+1$ 个点的点值即可 Lagrange 插值求得 $\varphi_A(\lambda)$。

> 实际上，我们知道 $[\lambda^n]\varphi_A(\lambda)=1$，可以少计算一项，但这对减少复杂度并没有什么作用。

分别代入 $\lambda=0, 1,\dots, n$ 计算 $\det(\lambda I-A)$ 再 Lagrange 插值即可。

时间复杂度 $O(n^4)$。



### Upper Hessenberg 矩阵

> ——写成英文就可以规避分词错误了。

我们希望能求得与 $A$ 相似的 $P^{-1}AP$，使之具有可以简单计算特征多项式的形式。

> 相似矩阵具有相同的特征多项式。
>
> 更本质地，相似矩阵只是同一线性变换 $\mathscr{A}$ 在不同基下的矩阵，「特征」只是 $\mathscr{A}$ 的特征。

当我们对 $A$ 左乘 $P^{-1}$ 进行初等行变换时，同时右乘 $P$ 进行初等列变换即可保证变换前后矩阵相似。

容易发现，单次变换只会影响第 $i,j$ 行与第 $i,j$ 列，于是可以以第 $i+1$ 行消去第 $i+2$ 行至第 $n$ 行的第 $i$ 列，得到与 $A$ 相似的 Upper Hessenberg 矩阵 $B=P^{-1}AP$：

$$B=\begin{bmatrix}a_{11}&a_{12}&a_{13}&\cdots&a_{1n}\\a_{21}&a_{22}&a_{23}&\cdots&a_{2n}\\&a_{32}&a_{33}&\cdots&a_{3n}\\&&\ddots&\ddots&\vdots\\&&&a_{n,n-1}&a_{nn}\end{bmatrix}$$

此时计算 $\det(\lambda I-B)$ 将简化不少。设 $(\lambda I-B)_i$ 为 $\lambda I-B$ 第 $i$ 至 $n$ 行与第 $i$ 至 $n$ 列组成的 $n-i+1$ 阶主子式，对其按第 $i$ 列展开：

$$D_i=\det (\lambda I-B)_i=(\lambda-a_{ii})D_{i+1}-\sum_{j=i+1}^n a_{i,j}D_{j+1}\prod_{k=i+1}^j a_{k,k-1}$$

注意到 $D_i$ 为关于 $\lambda$ 的 $n-i+1$ 次首一多项式，递推计算 $D_i$ 时仅出现单个一次项 $\lambda-a_{ii}$，因此可以依次递推计算 $D_i$，$D_1=\det(\lambda I-B)=\varphi_B(\lambda)$ 即为 $B$ 的特征多项式。

$A$ 与 $B$ 相似，所以 $\varphi_A(\lambda)=\varphi_B(\lambda)$ 即为所求。

将 $A$ 变换为 Upper Hessenberg 矩阵复杂度为 $O(n^3)$，递推计算 $\varphi_B(\lambda)$ 复杂度也为 $O(n^3)$，总时间复杂度 $O(n^3)$。
