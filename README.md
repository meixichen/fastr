<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

\newcommand{\bm}[1]{\boldsymbol{#1}}
\newcommand{\tx}[1]{\mathrm{#1}}
\newcommand{\xx}{{\bm{x}}}
\newcommand{\yy}{{\bm{y}}}
\newcommand{\XX}{{\bm{X}}}
\newcommand{\YY}{{\bm{Y}}}
\newcommand{\ZZ}{{\bm{Z}}}
\newcommand{\tth}{{\bm{\theta}}}
\newcommand{\pps}{{\bm{\psi}}}
\newcommand{\uu}{{\bm{u}}}
\newcommand{\SSi}{{\bm{\Sigma}}}
\newcommand{\LLam}{{\bm{\Lambda}}}
\newcommand{\PPsi}{{\bm{\Psi}}}
\newcommand{\VV}{{\bm{V}}}
\newcommand{\iid}{{\overset{iid}{\sim}}}

# MNFA

*Meixi Chen, Martin Lysy*

---

### Description 

***MNFA*** (**M**ulti-**N**euron **F**actor **A**nalysis) is an R package for efficient factor analysis of spike trains simultaneously recorded from multiple neurons. Details about the model is given below.

### Model

The experimental window is discretized into $n$ time bins each of length $dt$. For a dataset of $q$ simultaneously recorded neurons, we define the following notations:

1. Let $y_{i,t}$ denote whether a spike occurs at time $t$ for neuron $i$. 

2. Let $x_{i,t}$ denote the value of the latent process at time $t$ for neuron $i$. Let $\xx_t \in \mathbb{R}^{q\times 1}$ represent the value of all neurons' latent paths at time $t$.

3. Let $dx_{i,t}=x_{i,t}-x_{i,t-1}$ denote the value of increment at time $t$ for neuron $i$. Similar to 2, we define $d\xx_t \in \mathbb{R}^{q\times 1}$ as its vector representation for all $q$ neurons at time $t$.

We assume that the values of $\xx_t$ only depends on the values of $\xx_{t-1}$ and some neuron-specific drift parameters $\bm{\alpha}=(\alpha_1,\ldots, \alpha_q)$. Suppose we believe that the latent path increments of the $q$ neurons can be represented using $d \ (d<q)$ independent factors: $\bm{f}_1, \bm{f}_2, \ldots, \bm{f}_d \sim \mathcal{N}(0, \bm{I})$. That is
$$d\xx_t = \bm{\Lambda} \bm{f} + \bm{\epsilon},$$
where $\LLam_{q\times d}$ is the loading matrix and $\bm{\epsilon}\sim \mathcal{N}(0,\PPsi)$ is the vector of neurons' uniqueness with a diagonal covariance matrix $\bm{\Psi} = \tx{diag}(\sigma_1^2,\ldots, \sigma_q^2)$.

Then the model is written as follows.

- **Data layer**: for $i=1,2,\ldots,q$ and $t=1,2,\ldots, n$,
\begin{equation}
y_{i,t} \sim \mathrm{Bernoulli}\left(\mathrm{logit}\left(x_{i,t} - k_i\sum_{j=1}^{t-1}y_{i,t}\right)\right)
\end{equation}
where $k_i$ is the threshold parameter for neuron $i$.

- **Latent layer**:
\begin{equation}
d\xx_t \sim \tx{MVN}(\bm{\alpha} dt, \SSi dt)
\end{equation}
where the diagonal elements of $\SSi=\LLam \LLam' + \PPsi$ are constrained to 1.

Thus, the parameters to estimate are 

- Drift parameters: $\bm{\alpha} \in \mathbb{R}^{q\times 1}$.

- Threshold parameters: $\bm{k} \in \mathbb{R}^{q\times 1}$.

- Loading matrix: $\LLam \in \mathbb{R}^{q\times d}$.




