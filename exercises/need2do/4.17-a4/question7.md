# Solutions to additional exercises

# Standard form LP barrier method

In the following three exercises, you will implement a barrier method for solving the standard form LP  

$$
{\begin{array}{l l}{{\mathrm{minimize}}}&{c^{T}x}\\ {{\mathrm{subject~to}}}&{A x=b,\quad x\succeq0,}\end{array}}
$$  

with variable $x\in\mathbf{R}^{n}$ , where $A\in\mathbf{R}^{m\times n}$ , with $m\:<\:n$ . Throughout this exercise we will assume that $A$ is full rank, and the sublevel sets $\{x\mid A x=b,x\succeq0,c^{T}x\leq\gamma\}$ are all bounded. (If this is not the case, the centering problem is unbounded below.)  

1. Centering step. Implement Newton’s method for solving the centering problem  

$$

\begin{array}{l l}{\mathrm{minimize~}}&{c^{T}x-\sum_{i=1}^{n}\log x_{i}}\\ {\mathrm{subject~to}}&{A x=b,}\end{array}
$$  

with variable $x$ , given a strictly feasible starting point $x_{0}$ .  

Your code should accept $A$ , $b$ , $c$ , and $x_{0}$ , and return $x^{\star}$ , the primal optimal point, $\nu^{\star}$ , a dual optimal point, and the number of Newton steps executed.  

Use the block elimination method to compute the Newton step. (You can also compute the Newton step via the KKT system, and compare the result to the Newton step computed via block elimination. The two steps should be close, but if any $x_{i}$ is very small, you might get a warning about the condition number of the KKT matrix.)  

Plot $\lambda^{2}/2$ versus iteration $k$ , for various problem data and initial points, to verify that your implementation gives asymptotic quadratic convergence. As stopping criterion, you can use $\lambda^{2}/2\le10^{-6}$ . Experiment with varying the algorithm parameters $\alpha$ and $\beta$ , observing the effect on the total number of Newton steps required, for a fixed problem instance. Check that your computed $x^{\star}$ and $\nu^{\star}$ (nearly) satisfy the KKT conditions.  

To generate some random problem data $(i.e.,A,b,c,x_{0})$ ), we recommend the following approach. First, generate $A$ randomly. (You might want to check that it has full rank.) Then generate a random positive vector $x_{0}$ , and take $b=A x_{0}$ . (This ensures that $x_{0}$ is strictly feasible.) The parameter $c$ can be chosen randomly. To be sure the sublevel sets are bounded, you can add a row to $A$ with all positive elements. If you want to be able to repeat a run with the same problem data, be sure to set the state for the uniform and normal random number generators.  

Here are some hints that may be useful.  

• We recommend computing $\lambda^{2}$ using the formula $\lambda^{2}=-\Delta x_{\mathrm{nt}}^{T}\nabla f(x)$ . You don’t really need $\lambda$ for anything; you can work with $\lambda^{2}$ instead. (This is important for reasons described below.)  

There can be small numerical errors in the Newton step $\Delta x_{\mathrm{nt}}$ that you compute. When $x$ is nearly optimal, the computed value of $\lambda^{2}$ , i.e., $\lambda^{2}=-\Delta x_{\mathrm{nt}}^{T}\nabla f(x)$ , can actually be (slightly) negative. If you take the squareroot to get $\lambda$ , you’ll get a complex number, and you’ll never recover. Moreover, your line search will never exit. However, this only happens when $x$ is nearly optimal. So if you exit on the condition $\lambda^{2}/2\le10^{-6}$ , everything will be fine, even when the computed value of $\lambda^{2}$ is negative. • For the line search, you must first multiply the step size $t$ by $\beta$ until $x+t\Delta x_{\mathrm{nt}}$ is feasible (i.e., strictly positive). If you don’t, when you evaluate $f$ you’ll be taking the logarithm of negative numbers, and you’ll never recover.  


# Solution

1. The Newton step $\Delta x_{\mathrm{nt}}$ is defined by the KKT system:  

$$
\left[\begin{array}{c c}{H}&{A^{T}}\\ {A}&{0}\end{array}\right]\left[\begin{array}{c}{\Delta x_{\mathrm{nt}}}\\ {w}\end{array}\right]=\left[\begin{array}{c}{-g}\\ {0}\end{array}\right],
$$  

where $H=\mathbf{diag}(1/x_{1}^{2},\dots,1/x_{n}^{2})$ , and $g=c-(1/x_{1},\ldots,1/x_{n})$ . The KKT system can be efficiently solved by block elimination, i.e., by solving  

$$

A H^{-1}A^{T}w=-A H^{-1}g,
$$  

and setting $\Delta x_{\mathrm{nt}}=-H^{-1}(A^{T}w+g)$ . The KKT optimality condition is  

$$
A^{T}\nu^{\star}+c-\left(1/x_{1}^{\star},\ldots,1/x_{n}^{\star}\right)=0.
$$  

When the Newton method converges, i.e., $\Delta x_{\mathrm{nt}}\approx0$ , $w$ is the dual optimal point $\nu^{\star}$ .   
The following function computes the analytic center using Newton’s method.   

```
function [x_star, nu_star, lambda_hist] = lp_acent(A,b,c,x_0)
% solves problem
% minimize c’*x - sum(log(x))
% subject to A*x = b
% using Newton’s method, given strictly feasible starting point x0
% input (A, b, c, x_0)
% returns primal and dual optimal points
% lambda_hist is a vector showing lambda^2/2 for each newton step
% returns [], [] if MAXITERS reached, or x_0 not feasible
% algorithm parameters
ALPHA = 0.01;
BETA = 0.5;
EPSILON = 1e-6;
MAXITERS = 100;

if (min(x_0) <= 0) || (norm(A*x_0 - b) > 1e-3) % x0 not feasible
    fprintf('FAILED');
    nu_star = [];
    x_star = [];
    lambda_hist=[];
    return;
end
m = length(b);
n = length(x_0);
x = x_0;
lambda_hist = [];
for iter = 1:MAXITERS
    H = diag(x.^(-2));
    g = c - x.^(-1);
    % newton step by elimination method
    w = (A*diag(x.^2)*A’)\(-A*diag(x.^2)*g);
    dx = -diag(x.^2)*(A’*w + g);
    lambdasqr = -g’*dx; 
    lambda_hist = [lambda_hist lambdasqr/2];
    if lambdasqr/2 <= EPSILON
        break;
    end
    % backtracking line search
    % first bring the point inside the domain
    t = 1;
    while min(x+t*dx) <= 0
        t = BETA*t;
    end
    % now do backtracking line search
    while c’*(t*dx)-sum(log(x+t*dx))+sum(log(x))-ALPHA*t*g’*dx> 0
        t = BETA*t;
    end
    x = x + t*dx;
end
if iter == MAXITERS % MAXITERS reached
    fprintf('ERROR: MAXITERS reached.\n');
    x_star = [];
    nu_star = [];
else
    x_star = x;
    nu_star = w;
end
```

The random data is generated as given in the problem statement, with $A\in\mathbf{R}^{100\times500}$ . The Newton decrement versus number of Newton steps is plotted below. Quadratic convergence is clear. The Newton direction computed by the two methods are very close. The KKT optimality condtions are verified for the points returned by the function.  
 
