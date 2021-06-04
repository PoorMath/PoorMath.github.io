### 常系数线性齐次递推

我们关注这样一类问题：已知 $f_1,\dots,f_k$，并且有 $f_i=\sum_{j=1}^k c_i f_{j-k}$，求 $f_n$。

将递推式写为矩阵形式：

$$\begin{bmatrix}f_{i+k}\\f_{i+k-1}\\\vdots\\\vdots\\f_{i+1}\end{bmatrix}=\begin{bmatrix}c_1&c_2&\cdots&\cdots&c_k\\1&&&&\\&1&&&\\&&\ddots&&\\&&&1&\end{bmatrix}\begin{bmatrix}f_{i+k-1}\\f_{i+k-2}\\\vdots\\\vdots\\f_i\end{bmatrix}$$

设转移矩阵为 $A$，那么所求的 $f_n$ 可以通过计算下式得到：

$$\begin{bmatrix}f_{n}\\f_{n-1}\\\vdots\\\vdots\\f_{n-k+1}\end{bmatrix}=A^{n-k}\begin{bmatrix}f_{k}\\f_{k-1}\\\vdots\\\vdots\\f_1\end{bmatrix}$$

问题转化为快速求得 $A$ 的幂次（或许需要在 $\bmod m$ 意义下）。

由 Cayley-Hamilton 定理，$A$ 的特征多项式 $\varphi(\lambda)=\vert \lambda I-A\vert$ 满足 $\varphi(A)=0$，这说明 $\varphi(A)$ 是 $A$ 的一个零化多项式，因此可以求 $A^{n-k}\bmod \varphi(A)$ 以避免多次矩阵乘法。

将 $\vert\lambda I-A\vert$ 按行按列展开：

$$\begin{aligned}\begin{vmatrix}\lambda I-A\end{vmatrix}&=\begin{vmatrix}\lambda-c_1&-c_2&\cdots&\cdots&-c_k\\-1&\lambda&&&\\&-1&\lambda&&\\&&\ddots&\ddots&\\&&&-1&\lambda\end{vmatrix}\\&=(\lambda-c_1)\lambda^{k-1}-c_2\lambda^{k-2}-\cdots-c_k\\&=\lambda^k-\sum_{i=1}^kc_i\lambda^{k-i}\end{aligned}$$

因此 $\varphi(A)=A^k-\sum_{i=1}^k c_i A^{k-i}=0$，移项即有 $A^k=\sum_{i=1}^k c_i A^{k-i}$。

对于两个 $\deg<k$ 的 $A$ 的多项式，其相乘结果对 $\varphi(A)$ 取模就可以暴力地在 $O(k^2)$ 内完成，那么快速幂计算 $A^n \bmod \varphi(A)$ 时间复杂度为 $O(k^2\log n)$。

若 $A\in(\mathrm{Z}_p)^{k\times k}$，$p$ 为 NTT 质数，那么取模这一步可以利用 NTT 在 $O(k\log k)$ 内完成，时间复杂度为 $O(k\log k\log n)$。

若 $A^{n-k}\bmod \varphi(A)=b_0I+b_1A+b_2A^2+\cdots+b_{k-1}A^{k-1}$，可得下式：

$$\begin{aligned}\begin{bmatrix}f_{n}\\f_{n-1}\\\vdots\\\vdots\\f_{n-k+1}\end{bmatrix}&=(A^{n-k}\bmod \varphi(A))\begin{bmatrix}f_{k}\\f_{k-1}\\\vdots\\\vdots\\f_1\end{bmatrix}\\&=\sum_{i=0}^{k-1}b_iA^i\begin{bmatrix}f_{k}\\f_{k-1}\\\vdots\\\vdots\\f_1\end{bmatrix}\\&=\sum_{i=0}^{k-1}b_i\begin{bmatrix}f_{k+i}\\f_{k+i-1}\\\vdots\\\vdots\\f_{i+1}\end{bmatrix}\end{aligned}$$

因此：

$$f_n=\sum_{i=0}^{k-1}b_if_{k+i}$$

也可以左乘 $A^{n-1}$ 取最后一行，这是类似的。
