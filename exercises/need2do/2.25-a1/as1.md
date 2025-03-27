# Assignment 1  

## Question 1  

Determine if each set below is convex.  

1. $\{(x,y)\in R_{++}^{2}|x/y\le1\}$   
2. $\{(x,y)\in R_{++}^{2}|x/y\geq1\}$   
3. $\{(x,y)\in R_{++}^{2}|x y\geq1\}$   
4. $\{(x,y)\in R_{++}^{2}|x y\leq1\}$   
5. $S=\{x\in R^{n}|x^{T}y\leq1$ for all $x\in C\}$ , where $C$ is a set (may not be convex).   
6. The ellipsoid $\{x|x^{T}P^{-1}x\leq1\}$ where $P\in S_{++}^{n}$ .  

## Question 2  

Give an example of two closed convex sets that are disjoint but cannot be strictly separated.  

## Question 3  

Supporting hyperplanes.  

1. Express the closed convex set $\{x\in\mathbb{R}_{+}^{2}|x_{1}x_{2}\geq1\}$ as an intersection of halfspaces.  

2. Let $C=\{x\in\mathbb{R}^{n}|\|x\|_{\infty}\leq1\}$ , the $l_{\infty}$ -norm unit ball in $\mathbb{R}^{n}$ , and let $\hat{x}$ be a point in the boundary of $C$ . Identify the supporting hyperplanes of $C$ at $\hat{x}$ explicitly.  

3. Express the ball $B=\{x\in R^{n}|\|x\|_{2}\leq1\}$ as an intersection of halfspaces.  

## Question 4  

Perspective function $P:R^{n+1}\to R^{n}$ is defined as:  

$$
P(x,t)=x/t,\quad\mathrm{dom}P=\{(x,t)\mid t>0\}
$$  

where $x\in R^{n},t>0$ (i.e. $d o m(P)=\{(x,t)|x\in R^{n},t>0\}$ .) Prove:  

• If $C\subset d o m(P)$ is convex, then $P(C)$ is convex.   
• If $D\subset R^{n}$ is convex, then $P^{-1}(D)$ is convex.  