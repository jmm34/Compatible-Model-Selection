# compatible-model-selection`


We consider a random variable $X$ defined on a probability space $(\Omega, {\cal B}, \mu)$.
The probability distribution of $X$ is assumed to be absolutely continuous with respect to the measure $\mu$.
Its density with respect to $\mu$ is denoted 
$$
g(\bx)\mu(\text{d}\bx)\,.
$$
A finite family of $M$ models $\mathcal{M}=\{1,\ldots,M\}$ is considered to approximate $g$. 
Each model $m$ is characterized by a parameterized density $f_m(\cdot;\btheta_m)$ with respect to the measure $\mu$,
with $\btheta_m\in\Theta_m$. For any $m$, the dimension of $\Theta_m$ is assumed to be finite. 

%La loi $\P^i_{\btheta_i}$ admet la densité de probabilité 
%par rapport à la mesure de référence $\mu$ :
%$$
%\P^i_{\btheta_i}(\text{d}\bx)=f_i(\bx;\btheta_i)\mu(\text{d}\bx)\,.
%$$

% Par construction les distributions de proba considérée ne sont pas étrangères.
We define 
$$
\text{KL}\left(g,f_m(\btheta_m)\right)=\int \log\left(\frac{g(\bx)}{f_m(\bx;\btheta_m)}\right)g(\bx)\mu(\text{d}\bx)
$$
%$$
%\text{KL}\left(\P,\P^i_{\btheta_i}\right)=+\infty\,.
%$$
and $\arg\min_{\btheta_m\in\Theta_m}\text{KL}\left(g,f_m(\btheta_m)\right)$ is assumed to be non empty for any $m \in {\cal M}$. 
The dimension of model $m$ is  $d_m=\dim(\Theta_m)$ and we define 
$$
\btheta_m^*\in\arg\min_{\btheta_m}\text{KL}\left(g,f_m(\btheta_m)\right)\,,
$$
and $K_m=\text{KL}\left(g,f_m(\btheta_m^*)\right)$.

%\newpage

\begin{definition}
The best model in the family ${\cal M}$ is defined as the most parsimonious model $m^*$ minimizing the 
Kullback-Leibler distance $K_m $:
$$
m^*\in\arg\min_{m\in\mathcal{M}}K_m\,.
$$
and 
$$
d_{m^*}\in\arg\min_{j\in\left\{\arg\min_{m\in\mathcal{M}}K_m\right\}}d_j\,.
$$
\end{definition}

\section{The notion of compatible models}

The aim is to select $m^*$ from a $n$-sample $X_1,\ldots,X_n$ of independent data arising from the density $g$.
Thus, $m^*$ can be estimated with the BIC criterion (\cite{Schwarz}) which under regularity conditions is expected to be consistent
(see for instance Raftery, 1995). But in practice, especially when $n$ is so large, it can happen that the BIC values of several models
are close to the highest value. In such cases, the choice of $m^*$ is questionable and we could wonder if it is not too depending of the sample at hand.
Thus, it is sensible to select a short list of models rather than choosing a unique model. This short list of models
would be referred to as the {\em compatible} models with the data. 
In this purpose, there is the need of assessing the variability on the BIC differences provided by the family of models ${\cal M}$.

Let $q_\alpha$ be the smallest value such that
\begin{equation} \label{obj}
\P_{g^{\otimes n}}\left[\left\{\mbox{BIC}(\hat m)-\mbox{BIC}(m^*)\right\}\leq q_\alpha\right]\geq 1-\alpha\,,
\end{equation}
where
$\hat m$ is the model selected with BIC, namely
$$
\hat m\in\arg\max_{m\in\{1,\ldots,M\}} 2\sum_{i=1}^n\log\left(f_m(\bx;\hat\btheta_m)\right)-\log(n)d_m\,,
$$
with
$$
\hat\btheta_m\in\arg\max_{\btheta_m}\sum_{i=1}^n\log\left(f_m(\bx;\btheta_m)\right)\,,
$$
it is then natural not to reject a model $m$ such that
\begin{equation} \label{compa}
\mbox{BIC}(\hat m)-\mbox{BIC}(m)\leq q_\alpha\,.
\end{equation}
The subset of models $m \in {\cal M}$ verifying the inequality (\ref{compa}) includes $m^*$
(the optimal model for the Kullback-Leibler distance) with a probability greater than $1-\alpha$.
It is called the subset of $\alpha$-compatible models and is denoted $C_\alpha$.
The size of $C_\alpha$ depends of $\alpha$ and $n$. For small values of  $\alpha$ and $n$, the cardinal of $C_\alpha$ will have a tendency to be large. For $n$ very large,
$C_\alpha$ will be reduced to the singleton $m^*$ since BIC is consistent.
Moreover, the larger $C_\alpha$ is, the more difficult the model selection problem is. 


