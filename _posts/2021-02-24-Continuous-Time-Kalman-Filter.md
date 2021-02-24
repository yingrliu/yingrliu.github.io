---
layout: post
title: Continuous-Time Kalman Filter
date: 2021-02-24
Author: Yingru 
tags: [Study Notes, Statistical Signal Processing]
comments: true
toc: true
---

> Up till now, I dont see any detailed deviation for the Continuous-Time Kalman Filter, the content of this note is all from the Chapter 8 of the Book: **Optimal State Estimation: Kalman, H Infinity, and Nonlinear Approaches.**

# Discretization

Recalled that for a simplified state-space model
$\frac{d}{dt} x_t = F_t x_t + B_tu_t,\quad z_t = H_tx_t + v_t, \quad t\geq 0$,
we have

$$
x_t = \Phi(t, 0) x_0 + \int^t_0 \Phi(t, s)B_tu_s ds,
$$


where `$\Phi(t, s)$` is the state transition function.
Then by setting `$0$` to `$t_{k-1}$`, `$t$` to `$t_k$` and assuming `$B_t, u_t$` are approximately constant, we have

$$
x_{t_k} = \Phi(t_k, t_{k-1}) x_{t_{k-1}} + \int^{t_k}_{t_{k-1}} \Phi(t_{k}, s)ds\cdot B_{t_{k-1}}u_{t_{k-1}}.
$$

The previous equation is exactly a discrete time state-space model:

$$
x_k = A_{k-1}x_{k-1} + G_{k-1} u_{k-1},
$$

where `$A_{k-1}=\Phi(t_{k-1}+\Delta t, t_{k-1})$` and `$G_{k-1}=\int^{t_{k-1}+\Delta t}_{t_{k-1}} \Phi(t_{k}, s)ds\cdot B_{t_{k-1}}$`.

## Discretization for a full state-space model
Consider the following state-space model:

$$
\frac{d}{dt}x = F_t x_t + G_tu_t + w_t,\,\, w_t \sim \mathcal{N}(0, Q_t)\\
y_t = H_t x_t + v_t,\,\, v_t \sim \mathcal{N}(0, R_t),
$$

its discretization is given as

$$
x_{k}= A_{k-1}x_{k-1} + B_{k-1}u_{k-1} + w_{k-1},\\
w_{k-1} = \Lambda_{k-1}w_{(k-1)\Delta t}\\
y_{k}=C_{k}x_{k} + v_{k},
$$

where

$$
A_{k-1} = \Phi(t_{k-1}+\Delta t, t_{k-1}),\\
B_{k-1} = \int^{t_{k-1}+\Delta t}_{t_{k-1}} \Phi(t_{k}, s)ds\cdot G_{t_{k-1}}\approx \Delta t I G_{t_{k-1}},\\
\Lambda_{k-1}=\int^{t_{k-1}+\Delta t}_{t_{k-1}} \Phi(t_{k}, s)ds\approx \Delta t I,\\
C_k = H_{(k-1)\Delta t},\\
w_{k-1} \sim \mathcal{N}(0, Q_{(k-1)\Delta t}\Delta t),\\
v_{k-1} \sim \mathcal{N}(0, \frac{R_{(k-1)\Delta t}}{\Delta t}).
$$

For simplicity, we define the notation `$Q_{k-1}=Q_{(k-1)\Delta t}\Delta t$`,  `$R_{k-1}=\frac{R_{(k-1)\Delta t}}{\Delta t}$` and `$\Phi_{k-1}=\Phi(t_{k-1}+\Delta t, t_{k-1})$`.

## Evaluating the Covariance after Discretization
By the previous step, we can discretize the differential equation into an discrete difference equation. However, the covariances of white noise `$w_t$` and `$w_t$` should be carefully evaluated, so that the stochasititicy is identical.

**(i).** The discretization of `$w_t$` and `$w_t$` is a little bit tricky (at leat in my opinion). In the book, If we discretize the state-space model, we can define a discrete-time Kalman Filter and the error of this book should be independent with `$\Delta t$`. To my understanding, we are discretizing the rectangle volume `$\int_0^{\Delta t}R_{\tau}d\tau$` in the curve of `$Q_t$` into a point value, to assure that the `$\Delta t \times R$` is equal to the integral, we set `$v_{t_{k-1}}\sim \mathcal{N}(0, R_{t_{k-1}}/\Delta t)$`. Similarly, we set
`$w_{t_{k-1}}\sim \mathcal{N}(0, Q_{t_{k-1}}/\Delta t)$`.

**(i).** We replace `$w_t$` by `$\Lambda_{k-1}w_{t_{k-1}}$`. As `$\Lambda_{k-1}\approx \Delta t I$` and `$w_{t_{k-1}}\sim \mathcal{N}(0, Q_{t_{k-1}}/\Delta t)$`, `$w_{k-1}\sim \mathcal{N}(0, \Delta t Q_{t_{k-1}})$`. 



# Deviation of Continuous-Time Kalman Filter
After the discretization, we can define a discrete-time Kalman Filter whose `$K_n$` is denoted as

$$
K_k = P_k^fC_k^TD_k^{-1}=P_k^fC_k^T(C_k P_k^f C_k^T + R_k)^{-1}\\
= P_k^f H_{k-1}^T(H_{k-1} P_k^fH_{k-1}^T + R_{k-1}/\Delta t)^{-1},
$$

and therefore,

$$
\lim_{\Delta t\rightarrow 0} \frac{K_k}{\Delta t} = \lim_{\Delta t\rightarrow 0} P_k^f H_{k-1}^T(\Delta t H_{k-1} P_k^fH_{k-1}^T + R_{k-1})^{-1} \\
= P_k^f H_{k-1}^T(R_{k-1})^{-1}.
$$

We also have
$$
P_k = (I - K_k C_k)P_k^f,\\
P_k^f = A_{k-1}P_{k-1}A_{k-1}^T + Q_{k-1}.
$$

Substituting `$A_{k-1}$` into it, we have

$$
P_{k}^f = \Phi_{k-1}P_{k-1}P_{k-1}^f\Phi_{k-1}^T + Q_{k-1}\Delta t\\
$$

By the property of STM, we further have

$$
P_{k}^f = (I + \Delta tF_{k-1}\Phi_{k-1})P_{k-1}(I + \Delta tF_{k-1}\Phi_{k-1})^T + Q_{k-1}\Delta t\\
= P_{k-1} + \Delta tF_{k-1}\Phi_{k-1}P_{k-1} + \Delta t P_{k-1}\Phi_{k-1}^TF_{k-1}^T\\ + (\Delta t)^2F_{k-1}\Phi_{k-1}P_{k-1}\Phi_{k-1}^TF_{k-1}^T + Q_{k-1}\Delta t,
$$

Substituting `$P_{k-1}^f$` into it, we have

$$
P_{k}^f = (I - K_{k-1} C_{k-1})P_{k-1}^f + \Delta tF_{k-1}\Phi_{k-1}(I - K_{k-1} C_{k-1})P_{k-1}^f \\
+ \Delta t (I - K_{k-1} C_{k-1})P_{k-1}^f\Phi_{k-1}^TF_{k-1}^T \\
+ (\Delta t)^2F_{k-1}\Phi_{k-1}(I - K_{k-1} C_{k-1})P_{k-1}^f\Phi_{k-1}^TF_{k-1}^T \\
+ Q_{(k-1)\Delta t}\Delta t,
$$

and

$$
\frac{P_{k}^f-P_{k-1}^f}{\Delta t} = -\frac{- K_{k-1} C_{k-1}P_{k-1}^f}{\Delta t} + F_{k-1}\Phi_{k-1}(I - K_{k-1} C_{k-1})P_{k-1}^f\\
+ (I - K_{k-1} C_{k-1})P_{k-1}^f\Phi_{k-1}^TF_{k-1}^T \\
+ \Delta tF_{k-1}\Phi_{k-1}(I - K_{k-1} C_{k-1})P_{k-1}^f\Phi_{k-1}^TF_{k-1}^T.
$$

Taking the limitation yields:

$$
\lim_{\Delta t\rightarrow 0} \frac{P_{k}^f-P_{k-1}^f}{\Delta t} =-P_k^f H_{k-1}^TR_{k-1}^{-1}H_{k-1}^TP_{k-1}^f\\
+ F_{k-1}(I - K_{k-1} C_{k-1})P_{k-1}^f + (I - K_{k-1} C_{k-1})P_{k-1}^fF_{k-1}^T + Q_{k-1}\\
=\lim_{\Delta t\rightarrow 0} -P_k^f H_{k-1}^TR_{k-1}^{-1}H_{k-1}^TP_{k-1}^f + F_{k-1}P_{k-1}^f + P_{k-1}^fF_{k-1}^T\\
- F_{k-1}K_{k-1} C_{k-1}P_{k-1}^f - K_{k-1} C_{k-1}P_{k-1}^fF_{k-1}^T + Q_{k-1}\\
$$

Recall that there is a `$\Delta t$` in computing `$K_{k-1}$`, which leads to `$\lim_{\Delta t\rightarrow 0}F_{k-1}K_{k-1} C_{k-1}P_{k-1}^f=0$` and `$\lim_{\Delta t\rightarrow 0}K_{k-1} C_{k-1}P_{k-1}^fF_{k-1}^T=0$`, we have

$$
\frac{d}{dt}P_t=-P_t H_t^TR_t^{-1}H_t^TP_t + F_tP_t + P_t F_t^T + Q_t,
$$

where `$P_t$` is the continuous counterpart of prior covariance `$P_k^f$` (not the posterior one). This equation is also an instance of the Differential Riccati Equation.

For the estimate `$x_k^a$`, when `$\Delta t\rightarrow 0$`, we also have

$$
x_k^a = A_{k-1}x_{k-1}^a + B_{k-1} u_{k-1} + K_k(z_k - C_kA_{k-1}x_{k-1}^a - C_kB_{k-1} u_{k-1})\\
=(I + \Delta tF_{k-1}\Phi_{k-1})x_{k-1}^a + B_{k-1} u_{k-1} + K_k(z_k - C_k(I + \Delta tF_{k-1}\Phi_{k-1})x_{k-1}^a - C_kB_{k-1} u_{k-1})\\
= x_{k-1}^a + \Delta tF_{k-1}\Phi_{k-1}x_{k-1}^a + B_{k-1} u_{k-1} + K_kz_k-K_kC_k(I + \Delta tF_{k-1}\Phi_{k-1})x_{k-1}^a\\
-K_k C_kB_{k-1} u_{k-1},\\
= x_{k-1}^a + \Delta tF_{k-1}\Phi_{k-1}x_{k-1}^a + \Delta t G_{k-1} u_{k-1} + \Delta tP_k^f H_{k-1}^T(R_{k-1})^{-1}z_k\\
-\Delta tP_k^f H_{k-1}^T(R_{k-1})^{-1}H_k(I + \Delta tF_{k-1}\Phi_{k-1})x_{k-1}^a\\
-\Delta t P_k^f H_{k-1}^T(R_{k-1})^{-1} H_kB_{k-1} u_{k-1}.
$$

Hence, we have

$$
\lim_{\Delta t\rightarrow 0}\frac{x_k^a - x_{k-1}^a}{\Delta t} = F_{k-1}x_{k-1}^a+ G_{k-1}u_{k-1} + P_k^f H_{k-1}^T(R_{k-1})^{-1}(z_k-H_kx_{k-1}^a).
$$

That is

$$
\frac{d}{dt}\hat{x}_t=F_t \hat{x}_t + G_t u_t + P_k H_t^TR_t^{-1}(Z-H_t\hat{x}_t).
$$

# Summary of Continuous-Time KF

$$
\frac{d}{dt}\hat{x}_t=F_t \hat{x}_t + G_t u_t + K_t(Z-H_t\hat{x}_t),\\
\frac{d}{dt}P_t=-P_t H_t^TR_t^{-1}H_t^TP_t + F_tP_t + P_t F_t^T + Q_t,\\
K_t = P_k H_t^TR_t^{-1}.
$$

# Accelerated Method for Riccati Equation
> The ODE of estimate covariance `$P$` is a Differential Riccati Equation, which is computationally costly. To acclerate the algorithm, several algorithms have been proposed, including *Transition Matrix Approach* and *Chandrasekhar Algorithm*. Details can be found in the Chapter 8.3 in the Book: **Optimal State Estimation: Kalman, H Infinity, and Nonlinear Approaches.**

# Continuous-Time KF for correlated Noise
Consider the following continuous-time system:

$$
\frac{d}{dt}x = F_t x_t + G_tu_t + w_t,\,\, w_t \sim \mathcal{N}(0, Q_t)\\
y_t = H_t x_t + v_t,\,\, v_t \sim \mathcal{N}(0, R_t),\\
\mathbb{E}(w_t v_\tau^T) = M_t\delta(t-\tau).
$$

As `$w_t$` and `$v_\tau$` are correlated in time `$t$`, we need to decorrelate them. We can rewrite the system dynamics as

$$
\frac{d}{dt}x = F_t x_t + G_tu_t + w_t + M_tR_t^{-1}(y_t - H_t x_t - v_t)\\
= (F_t-M_tR_t^{-1}H_t)x_t + M_t R_{t}^{-1}y_t + (w_t - M_tR_t^{-1}v_t)\\
= \tilde{F}_t x_t + \tilde{u}_t + \tilde{w}_t.
$$

where

$$
 \tilde{F}_t = F_t-M_tR_t^{-1}H_t,\,\, \tilde{u}_t=M_t R_{t}^{-1}y_t,\,\,\tilde{w}_t=w - M_tR_t^{-1}v_t.
$$

We have

$$
\mathbb{E}(\tilde{w}_t v_\tau^T) = \mathbb{E}(\tilde{w}_t v_\tau^T)= \mathbb{E}((w_t - M_tR_t^{-1}v_t) v_\tau^T) \\
= \mathbb{E}(w_t v_\tau^T) - M_tR_t^{-1}\mathbb{E}(v_t v_t^T) = M - M = 0
$$

and 

$$
\tilde{Q}_t = \mathbb{E}(\tilde{w}_t \tilde{w}_t^T) = Q_t - M_t R_t^{-1} M_t^T.
$$

We can define the continuous-time Kalman Filter for this equivalent dynamics system. In the book, it also introduces how to tackle the non-white noise.

# The Steady-State Continuous-Time Kalman Filter

\[To Be Continue\]