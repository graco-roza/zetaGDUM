# Phase 3 mathematics

## Target and notation

Let $S$ be a uniform unordered subset of $k$ distinct sites from $n$ sites. Write $A\subseteq S$ for an inclusion condition and $a=|A|$. Let $I_{is}=1(x_{is}>0)$, $o_s=\sum_iI_{is}$, $r_i=\sum_sI_{is}$, and $q_{is}=x_{is}/\sum_sx_{is}$. All exact expectations in this phase condition on the observed community and environmental values. They quantify the finite-community target, while the regression uncertainty experiments separately draw independent ecological sites.

## Exact conditional shared composition and dissimilarity

The probability that all remaining sites come from a species' occupied sites, given that all forced sites are occupied, is
$$
p_a(o_s)=\frac{\binom{o_s-a}{k-a}}{\binom{n-a}{k-a}}.
$$
Set this probability to zero when $o_s<k$. Multiplication by $\prod_{i\in A}I_{is}$ and summation over species gives $\zeta_I(k\mid A)$. In particular, the off-diagonal incidence pair matrix is the weighted cross-product $I\,\mathrm{diag}\{p_2(o_s)\}\,I^\top$.

For abundance, let $t_{s\ell}$ be the ordered distinct positive values of $q_{is}$, with $t_{s0}=0$, and $O_{s\ell}=\{i:q_{is}\ge t_{s\ell}\}$ with $m_{s\ell}=|O_{s\ell}|$. The minimum within a set admits the pointwise expansion
$$
\min_{i\in S}q_{is}=\sum_\ell(t_{s\ell}-t_{s,\ell-1})1\{S\subseteq O_{s\ell}\}.
$$
Taking conditional expectations replaces the indicator by $1\{A\subseteq O_{s\ell}\}p_a(m_{s\ell})$. Accumulated increments through the lowest abundance rank in $A$ give the conditional expectation for that species. An absent forced site has rank zero. For $A=\varnothing$, use all thresholds. Tied values share one threshold and its occupancy includes every tied site. The C++ implementation evaluates binomial ratios on the log scale and skips species occurring on fewer than $k$ sites.

For a random $d$-subset $A$ of sites, $d\le k$, every set $S$ contains exactly $\binom{k}{d}$ possible $A$. Consequently averaging the site, distinct pair or higher-degree conditional expectations over all $A$ recovers the corresponding marginal expectation.

For raw incidence dissimilarity, define $D_I(S)=\bar r_S-Z_I(S)$. For $a<n$,
$$
C_A=E[\bar r_S\mid A\subseteq S]=
\frac{\sum_{i\in A}r_i+\frac{k-a}{n-a}\sum_{i\notin A}r_i}{k}.
$$
For $a=n$, $C_A=\bar r$. The conditional expectation of $D_I$ is $C_A-\zeta_I(k\mid A)$. For relative abundance it is $1-\zeta_A(k\mid A)$. Richness-normalised incidence $D_N(S)=1-Z_I(S)/\bar r_S$ uses a weighted sample mean of the set-level response given site inclusion and positive set richness. An expectation of a ratio generally differs from a ratio of expectations.

### Computational cost

The implementation scans the dense matrix once per order and sorts positive abundance values within each species. The abundance sort costs $O(\sum_s o_s\log o_s)$ after the $O(np)$ scan. For each requested conditioning size, threshold accumulation uses work linear in the number of distinct positive values. Arbitrary explicit subsets add at most $O(\sum_A|A|)$ work per eligible species. The optional pair accumulation uses $O(\sum_{s:o_s\ge k}o_s^2)$ work and $O(n^2)$ output. It avoids an $n\times n\times p$ array. The incidence case needs no sorting. Site-only routines require $O(n)$ ordinary temporary storage after the input and do no pair accumulation. A compressed-column input could eliminate scanning zeros in the incidence path and reduce input memory for highly sparse data. The current R API materialises dense input, even for sparse R objects. At high occupancy or when a pair matrix is requested, compressed input alone would leave the pair accumulation and output unchanged.

## LCZD and regression descriptions

For a chosen response, $U_i(k)=E[D(S)\mid i\in S]$, and $\mathrm{LCZD}_i(k)=U_i(k)/\sum_j U_j(k)$ when the total is positive. The relative index is $n\,\mathrm{LCZD}_i(k)$. The centred version subtracts $E[D]$, which equals the mean of the site expectations for exact dissimilarities. The sampled ratio estimates each site's conditional mean with its own denominator, so their empirical mean can differ from the pooled sample mean. Undefined sampled conditional values propagate to normalised LCZD. At $k=2$, $U_i$ is mean dissimilarity over pairs containing $i$. Variance-based LCBD in the validation uses squared Euclidean distance from the centroid after Hellinger transformation, divided by the sum of those squared distances. Its scale and dissimilarity differ from those of LCZD.

The model response-scale site quantities use $U_{\text{pred},i}=E[\mu(S)\mid i\in S]$, evaluated on a uniform set design or by site-specific importance-weighted ratio on a proposal. The residual is $U_{\text{obs},i}-U_{\text{pred},i}$, with both quantities on the original dissimilarity scale. Working weights govern estimation of the regression, whereas proposal weights govern averaging of its predictions. Define $f(e,p)=\operatorname{logit}^{-1}(\theta_0+e+p)$, with $e$ the directional component including geographic distance when used, and $p$ the position component. For each set, the directional Shapley value is
$$
\tfrac12\{f(e,0)-f(0,0)+f(e,p)-f(0,p)\}.
$$
The position value is obtained by exchanging the roles. Their sum is $f(e,p)-f(0,0)$, so the identity also holds after averaging conditionally on each site. These are model contrasts against the intercept-only baseline. Residuals and Shapley values describe this model and its predictors.

## GDUM response cross-products and safeguards

Let $c(S)$ denote the working weight, $X(S)$ the existing design vector including intercept, and $Y(S)$ the fractional dissimilarity. Then
$$
L(\theta)=E[c\,\operatorname{softplus}(X^\top\theta)]-b^\top\theta,\qquad b=E[cYX].
$$
For unit-weight relative abundance, $E[cY\mid A]=1-\zeta_A(k\mid A)$. For mean-richness-weighted total incidence, $cY=\bar r_S-Z_I(S)$, including the zero contribution to the loss from a zero-richness set. In that case $E[cY\mid A]=C_A-\zeta_I(k\mid A)$. This cancellation applies to the default richness-weighted model. Constant-weight normalised incidence remains a distinct sampled target.

With $u_A=E[cY\mid A]$, $z_i$ a site predictor standardised on the reference basis and $B_b(z_i)$ its stored I-spline value, the intercept, position and spread entries respectively are
$$
u_\varnothing,\qquad \frac{k}{n}\sum_i z_i u_{\{i\}},\qquad
\frac{1}{\binom{n}{2}}\sum_{i<j}|B_b(z_i)-B_b(z_j)|u_{\{i,j\}}.
$$
The position identity uses $P(i\in S)=k/n$. The spread identity uses $P(i,j\in S)=\binom{k}{2}/\binom n2$ and the existing division of the set-level sum by $\binom{k}{2}$. Replacing $u_A$ with $E[c\mid A]$ gives companion moments $d=E[cX]$. The existing geographic term applies splines to the mean pairwise raw distance of each set. Its response and predictor moments remain sampled.

For a uniform design or a proposal with raw likelihood ratio $w=P/Q$, the hybrid estimator uses $\widehat{E}[f]=M^{-1}\sum_m w_m f_m$. Its gradient and Hessian are $\widehat E[c\mu X]-b$ and $\widehat E[c\mu(1-\mu)XX^\top]$, respectively, retaining sampled entries of $b$ for geographic predictors. The control-variate objective subtracts $\mu_0(X^\top\theta-\eta_0)$ from the sampled softplus for exactly integrable columns and adds the corresponding population linear term back. Equivalently, use $b_j+\mu_0(\widehat d_j-d_j)$ in the original hybrid for exact columns only. Here $\mu_0=E[cY]/E[c]$ is the intercept-only population response and $\eta_0$ its logit when interior. The implementation also centres the full softplus at $\eta_0$ for numerical stability. The population estimating equations agree for all three variants. For a fixed $\theta$, the raw-ratio objectives agree in expectation up to additive constants before any fallback, whereas minimisation, certification and fallback can affect finite-sample bias.

A sufficient certificate seeks adjusted fractional values $v_m\in(0,1)$ with $\sum_m w_mc_mX_mv_m/M=b_\mathrm{used}$ and full column rank among positive-weight rows. Starting at clipped conventional fitted probabilities $u_m$, solve the weighted Hessian system for a shift and set $v_m=u_m+u_m(1-u_m)X_m^\top\text{shift}$. Check strict interior and the moment residual. An exact certificate makes the objective equivalent up to a constant to fractional logistic loss with interior responses, which grows along every unbounded sequence of coefficients for a full-rank finite design. In floating point the implementation checks numerical tolerances. The certificate is sufficient and can reject compatible designs. The constrained optimiser retains nonnegative spread and free position/intercept, verifies stationarity, and returns the conventional fit with a reason if the certificate or optimisation fails. No penalty was added.

## Randomised incomplete designs

Independent uniform sets use unit weights. The existing species mixture uses the original base-to-proposal ratios. For design-based balance, site blocks form disjoint groups by cutting a uniform site permutation into $k$-sets. An affine pair cycle at prime $n$ uses all maps $t\mapsto at+b\pmod n$ with $a\ne0$, applied to $\{0,\ldots,k-1\}$, with a fresh random site relabelling per cycle. Affine maps send any ordered distinct pair to each ordered distinct pair exactly once. Hence a full labelled cycle includes each unordered site pair $k(k-1)$ times. At $k=2$, a shuffled exhaustive pair cycle is supported at any $n$. For environmental stratification, each row draws stratum counts with the uniform-set multivariate hypergeometric law. Cyclic allocation in each randomly permuted stratum yields uniform subsets conditional on these counts. Individual rows in each design are therefore marginally uniform. Dependence affects the joint distribution and the precision of nonlinear estimators, especially for truncated cycles. These constructions are inspired by design-based incomplete U-statistics. They do not invoke the variance guarantees of the original orthogonal array construction. [Kong & Zheng (2021)](https://www3.stat.sinica.edu.tw/statistica/j31n3/j31n322/j31n322.html)

Incomplete empirical risk and optimal logistic subsampling motivate evaluating variance and cost on the regression target itself. Here multisite sets reuse ecological sites, so pointwise calculations over uniform sets and independent ecological replication answer different questions. [Clémençon, Colin & Bellet (2016)](https://www.jmlr.org/papers/v17/15-012.html), [Wang, Zhu & Ma (2018)](https://doi.org/10.1080/01621459.2017.1292914)

## Overlap-aware covariance

For a fixed symmetric kernel $g$ of iid sites, let $C_j=\operatorname{Cov}[g(S),g(T)\mid |S\cap T|=j]$ over ecological samples. Disjoint sets have $C_0=0$. The complete U-statistic has covariance $V=\sum_{j=0}^k\pi_j C_j$, where $\pi_j=\binom{k}{j}\binom{n-k}{k-j}/\binom nk$. Matched sampling uses $B$ independent groups conditional on the sites, each containing $m\le\lfloor n/k\rfloor$ mutually disjoint sets drawn by random permutation. Let $N=Bm$, $W$ the mean covariance within groups with denominator $m-1$, and $T$ the pooled sample covariance with denominator $N-1$. With $\delta=(m-1)/(N-1)$, $E[W]=C_k$ and $E[T]=(1-\delta)(C_k-V)+\delta C_k$. Thus
$$
\widehat V=\frac{W-T}{1-\delta},\qquad
\widehat V_{\mathrm{matched}}=W-\frac{N-1}{N}T.
$$
The second expression targets $(1-1/B)V+C_k/N$, the unconditional covariance of the incomplete matched mean. The pooled-versus-within correction accounts for estimating the kernel mean and for finite matched sampling. The package uses the matrix extension of the authors' scalar result. [Xu, Zhu & Shao (2024), §§3.2–3.4](https://doi.org/10.1214/24-EJS2247)

The separate overlap-resolved experiment draws an ordered union of $2k-j$ distinct sites, using the same first $j$ labels in both sets. Define $D_j=E[(g(S)-g(T))(g(S)-g(T))^\top/2\mid j]$. Exchangeability gives $D_j=C_k-C_j$ for iid sites and a fixed kernel. Using empirical $\widehat D_0-\widehat D_j$ avoids inserting a sample estimate of the population score mean. $D_k=0$ exactly. The complete covariance estimate is $\sum_j\pi_j(\widehat D_0-\widehat D_j)$. The sampled contrasts can be noisy, and the common disjoint term enters every contrast. For the matched point estimate, adding $\widehat C_k/N$ and multiplying the complete term by $1-1/B$ accounts for its matched Monte Carlo variance in the fixed-kernel calculation.

For regression, the package allows sandwich propagation $\widehat H^{-1}\widehat V\widehat H^{-1}$ using a positive-definite Hessian. The calibration uses $g(S;\widehat\theta)=cX(Y-\mu_{\widehat\theta})$. Fitting makes the kernel depend on the same observed data, so finite-sample unbiasedness of the fixed-kernel estimate does not transfer. Both variants require $2k\le n$ for disjoint comparisons. They preserve signed covariance estimates and report negative eigenvalues. Regression intervals additionally require a regular interior model, appropriate linearisation, and ecological independence. Spatially dependent sites, boundary parameters, near singularity and high-order asymptotics require separate justification.
