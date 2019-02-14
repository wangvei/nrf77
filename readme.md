# Numerical Recipes in Fortran 77

## Computer Programs by Chapter and Section

| sect  | program  | description                          |
| ----- | -------- | ------------------------------------ |
|  1.0  | [`flmoon`](flmoon.f) | calculate phases of the moon by date |
|  1.1  | [`julday`](julday.f) | Julian Day number from calendar date |
|  1.1  | [`badluk`](badluk.f) | Friday the 13th when the moon is full |
|  1.1  | [`caldat`](caldat.f) | calendar date from Julian day number |
|  2.1  | [`gaussj`](gaussj.f) | Gauss-Jordan matrix inversion and linear equation solution |
|  2.3  | [`ludcmp`](ludcmp.f) | linear equation solution, LU decomposition |
|  2.3  | [`lubksb`](lubksb.f) | linear equation solution, backsubstitution |
|  2.4  | [`tridag`](tridag.f) | solution of tridiagonal systems |
|  2.4  | [`banmul`](banmul.f) | multiply vector by band diagonal matrix |
|  2.4  | [`bandec`](bandec.f) | band diagonal systems, decomposition |
|  2.4  | [`banbks`](banbks.f) | band diagonal systems, backsubstitution |
|  2.5  | [`mprove`](mprove.f) | linear equation solution, iterative improvement |
|  2.6  | [`svbksb`](svbksb.f) | singular value backsubstitution |
|  2.6  | [`svdcmp`](svdcmp.f) | singular value decomposition of a matrix |
|  2.6  | [`pythag`](pythag.f) | calculate (a2 + b2)1=2 without overflow |
|  2.7  | [`cyclic`](cyclic.f) | solution of cyclic tridiagonal systems |
|  2.7  | [`sprsin`](sprsin.f) | convert matrix to sparse format |
|  2.7  | [`sprsax`](sprsax.f) | product of sparse matrix and vector |
|  2.7  | [`sprstx`](sprstx.f) | product of transpose sparse matrix and vector |
|  2.7  | [`sprstp`](sprstp.f) | transpose of sparse matrix |
|  2.7  | [`sprspm`](sprspm.f) | pattern multiply two sparse matrices |
|  2.7  | [`sprstm`](sprstm.f) | threshold multiply two sparse matrices |
|  2.7  | [`linbcg`](linbcg.f) | biconjugate gradient solution of sparse systems |
|  2.7  | [`snrm  `](snrm.f)   | used by linbcg for vector norm |
|  2.7  | [`atimes`](atimes.f) | used by linbcg for sparse multiplication |
|  2.7  | [`asolve`](asolve.f) | used by linbcg for preconditioner |
|  2.8  | [`vander`](vander.f) | solve Vandermonde systems |
|  2.8  | [`toeplz`](toeplz.f) | solve Toeplitz systems |
|  2.9  | [`choldc`](choldc.f) | Cholesky decomposition |
|  2.9  | [`cholsl`](cholsl.f) | Cholesky backsubstitution |
|  2.10 | [`qrdcmp`](qrdcmp.f) | QR decomposition |
|  2.10 | [`qrsolv`](qrsolv.f) | QR backsubstitution |
|  2.10 | [`rsolv `](rsolv.f)  | right triangular backsubstitution |
|  2.10 | [`qrupdt`](qrupdt.f) | update a QR decomposition |
|  2.10 | [`rotate`](rotate.f) | Jacobi rotation used by qrupdt |
|  3.1  | [`polint`](polint.f) | polynomial interpolation |
|  3.2  | [`ratint`](ratint.f) | rational function interpolation |
|  3.3  | [`spline`](spline.f) | construct a cubic spline |
|  3.3  | [`splint`](splint.f) | cubic spline interpolation |
|  3.4  | [`locate`](locate.f) | search an ordered table by bisection |
|  3.4  | [`hunt  `](hunt.f)   | search a table when calls are correlated |
|  3.5  | [`polcoe`](polcoe.f) | polynomial coefficients from table of values |
|  3.5  | [`polcof`](polcof.f) | polynomial coefficients from table of values |
|  3.6  | [`polin2`](polin2.f) | two-dimensional polynomial interpolation |
|  3.6  | [`bcucof`](bcucof.f) | construct two-dimensional bicubic |
|  3.6  | [`bcuint`](bcuint.f) | two-dimensional bicubic interpolation |
|  3.6  | [`splie2`](splie2.f) | construct two-dimensional spline |
|  3.6  | [`splin2`](splin2.f) | two-dimensional spline interpolation |
|  4.2  | [`trapzd`](trapzd.f) | trapezoidal rule |
|  4.2  | [`qtrap `](qtrap.f)  | integrate using trapezoidal rule |
|  4.2  | [`qsimp `](qsimp.f)  | integrate using Simpson’s rule |
|  4.3  | [`qromb `](qromb.f)  | integrate using Romberg adaptive method |
|  4.4  | [`midpnt`](midpnt.f) | extended midpoint rule |
|  4.4  | [`qromo `](qromo.f)  | integrate using open Romberg adaptive method |
|  4.4  | [`midinf`](midinf.f) | integrate a function on a semi-infinite interval |
|  4.4  | [`midsql`](midsql.f) | integrate a function with lower square-root singularity |
|  4.4  | [`midsqu`](midsqu.f) | integrate a function with upper square-root singularity |
|  4.4  | [`midexp`](midexp.f) | integrate a function that decreases exponentially |
|  4.5  | [`qgaus `](qgaus.f)  | integrate a function by Gaussian quadratures |
|  4.5  | [`gauleg`](gauleg.f) | Gauss-Legendre weights and abscissas |
|  4.5  | [`gaulag`](gaulag.f) | Gauss-Laguerre weights and abscissas |
|  4.5  | [`gauher`](gauher.f) | Gauss-Hermite weights and abscissas |
|  4.5  | [`gaujac`](gaujac.f) | Gauss-Jacobi weights and abscissas |
|  4.5  | [`gaucof`](gaucof.f) | quadrature weights from orthogonal polynomials |
|  4.5  | [`orthog`](orthog.f) | construct nonclassical orthogonal polynomials |
|  4.6  | [`quad3d`](quad3d.f) | integrate a function over a three-dimensional space |
|  5.1  | [`eulsum`](eulsum.f) | sum a series by Eulervan Wijngaarden algorithm |
|  5.3  | [`ddpoly`](ddpoly.f) | evaluate a polynomial and its derivatives |
|  5.3  | [`poldiv`](poldiv.f) | divide one polynomial by another |
|  5.3  | [`ratval`](ratval.f) | evaluate a rational function |
|  5.7  | [`dfridr`](dfridr.f) | numerical derivative by Ridders’ method |
|  5.8  | [`chebft`](chebft.f) | fit a Chebyshev polynomial to a function |
|  5.8  | [`chebev`](chebev.f) | Chebyshev polynomial evaluation |
|  5.9  | [`chder `](chder.f)  | derivative of a function already Chebyshev fitted |
|  5.9  | [`chint `](chint.f)  | integrate a function already Chebyshev fitted |
|  5.10 | [`chebpc`](chebpc.f) | polynomial coefficients from a Chebyshev fit |
|  5.10 | [`pcshft`](pcshft.f) | polynomial coefficients of a shifted polynomial |
|  5.11 | [`pccheb`](pccheb.f) | inverse of chebpc; use to economize power series |
|  5.12 | [`pade  `](pade.f)   | Pad´e approximant from power series coefficients |
|  5.13 | [`ratlsq`](ratlsq.f) | rational fit by least-squares method |
|  6.1  | [`gammln`](gammln.f) | logarithm of gamma function |
|  6.1  | [`factrl`](factrl.f) | factorial function |
|  6.1  | [`bico  `](bico.f)   | binomial coefficients function |
|  6.1  | [`factln`](factln.f) | logarithm of factorial function |
|  6.1  | [`beta  `](beta.f)   | beta function |
|  6.2  | [`gammp `](gammp.f)  | incomplete gamma function |
|  6.2  | [`gammq `](gammq.f)  | complement of incomplete gamma function |
|  6.2  | [`gser  `](gser.f)   | series used by gammp and gammq |
|  6.2  | [`gcf   `](gcf.f)    | continued fraction used by gammp and gammq |
|  6.2  | [`erf   `](erf.f)    | error function |
|  6.2  | [`erfc  `](erfc.f)   | complementary error function |
|  6.2  | [`erfcc `](erfcc.f)  | complementary error function, concise routine |
|  6.3  | [`expint`](expint.f) | exponential integral En |
|  6.3  | [`ei    `](ei.f)     | exponential integral Ei |
|  6.4  | [`betai `](betai.f)  | incomplete beta function |
|  6.4  | [`betacf`](betacf.f) | continued fraction used by betai |
|  6.5  | [`bessj0`](bessj0.f) | Bessel function J0 |
|  6.5  | [`bessy0`](bessy0.f) | Bessel function Y0 |
|  6.5  | [`bessj1`](bessj1.f) | Bessel function J1 |
|  6.5  | [`bessy1`](bessy1.f) | Bessel function Y1 |
|  6.5  | [`bessy `](bessy.f)  | Bessel function Y of general integer order |
|  6.5  | [`bessj `](bessj.f)  | Bessel function J of general integer order |
|  6.6  | [`bessi0`](bessi0.f) | modified Bessel function I0 |
|  6.6  | [`bessk0`](bessk0.f) | modified Bessel function K0 |
|  6.6  | [`bessi1`](bessi1.f) | modified Bessel function I1 |
|  6.6  | [`bessk1`](bessk1.f) | modified Bessel function K1 |
|  6.6  | [`bessk `](bessk.f)  | modified Bessel function K of integer order |
|  6.6  | [`bessi `](bessi.f)  | modified Bessel function I of integer order |
|  6.7  | [`bessjy`](bessjy.f) | Bessel functions of fractional order |
|  6.7  | [`beschb`](beschb.f) | Chebyshev expansion used by bessjy |
|  6.7  | [`bessik`](bessik.f) | modified Bessel functions of fractional order |
|  6.7  | [`airy  `](airy.f)   | Airy functions |
|  6.7  | [`sphbes`](sphbes.f) | spherical Bessel functions jn and yn |
|  6.8  | [`plgndr`](plgndr.f) | Legendre polynomials, associated (spherical harmonics) |
|  6.9  | [`frenel`](frenel.f) | Fresnel integrals S(x) and C(x) |
|  6.9  | [`cisi  `](cisi.f)   | cosine and sine integrals Ci and Si |
|  6.10 | [`dawson`](dawson.f) | Dawson’s integral |
|  6.11 | [`rf    `](rf.f)     | Carlson’s elliptic integral of the first kind |
|  6.11 | [`rd    `](rd.f)     | Carlson’s elliptic integral of the second kind |
|  6.11 | [`rj    `](rj.f)     | Carlson’s elliptic integral of the third kind |
|  6.11 | [`rc    `](rc.f)     | Carlson’s degenerate elliptic integral |
|  6.11 | [`ellf  `](ellf.f)   | Legendre elliptic integral of the first kind |
|  6.11 | [`elle  `](elle.f)   | Legendre elliptic integral of the second kind |
|  6.11 | [`ellpi `](ellpi.f)  | Legendre elliptic integral of the third kind |
|  6.11 | [`sncndn`](sncndn.f) | Jacobian elliptic functions |
|  6.12 | [`hypgeo`](hypgeo.f) | complex hypergeometric function |
|  6.12 | [`hypser`](hypser.f) | complex hypergeometric function, series evaluation |
|  6.12 | [`hypdrv`](hypdrv.f) | complex hypergeometric function, derivative of |
|  7.1  | [`ran0  `](ran0.f)   | random deviate by Park and Miller minimal standard |
|  7.1  | [`ran1  `](ran1.f)   | random deviate, minimal standard plus shuffle |
|  7.1  | [`ran2  `](ran2.f)   | random deviate by L’Ecuyer long period plus shuffle |
|  7.1  | [`ran3  `](ran3.f)   | random deviate by Knuth subtractive method |
|  7.2  | [`expdev`](expdev.f) | exponential random deviates |
|  7.2  | [`gasdev`](gasdev.f) | normally distributed random deviates |
|  7.3  | [`gamdev`](gamdev.f) | gamma-law distribution random deviates |
|  7.3  | [`poidev`](poidev.f) | Poisson distributed random deviates |
|  7.3  | [`bnldev`](bnldev.f) | binomial distributed random deviates |
|  7.4  | [`irbit1`](irbit1.f) | random bit sequence |
|  7.4  | [`irbit2`](irbit2.f) | random bit sequence |
|  7.5  | [`psdes `](psdes.f)  | “pseudo-DES” hashing of 64 bits |
|  7.5  | [`ran4  `](ran4.f)   | random deviates from DES-like hashing |
|  7.7  | [`sobseq`](sobseq.f) | Sobol’s quasi-random sequence |
|  7.8  | [`vegas `](vegas.f)  | adaptive multidimensionalMonte Carlo integration |
|  7.8  | [`rebin `](rebin.f)  | sample rebinning used by vegas |
|  7.8  | [`miser `](miser.f)  | recursive multidimensional Monte Carlo integration |
|  7.8  | [`ranpt `](ranpt.f)  | get random point, used by miser |
|  8.1  | [`piksrt`](piksrt.f) | sort an array by straight insertion |
|  8.1  | [`piksr2`](piksr2.f) | sort two arrays by straight insertion |
|  8.1  | [`shell `](shell.f)  | sort an array by Shell’s method |
|  8.2  | [`sort  `](sort.f)   | sort an array by quicksort method |
|  8.2  | [`sort2 `](sort2.f)  | sort two arrays by quicksort method |
|  8.3  | [`hpsort`](hpsort.f) | sort an array by heapsort method |
|  8.4  | [`indexx`](indexx.f) | construct an index for an array |
|  8.4  | [`sort3 `](sort3.f)  | sort, use an index to sort 3 or more arrays |
|  8.4  | [`rank  `](rank.f)   | construct a rank table for an array |
|  8.5  | [`select`](select.f) | find the Nth largest in an array |
|  8.5  | [`selip `](selip.f)  | find the Nth largest, without altering an array |
|  8.5  | [`hpsel `](hpsel.f)  | find M largest values, without altering an array |
|  8.6  | [`eclass`](eclass.f) | determine equivalence classes from list |
|  8.6  | [`eclazz`](eclazz.f) | determine equivalence classes from procedure |
|  9.0  | [`scrsho`](scrsho.f) | graph a function to search for roots |
|  9.1  | [`zbrac `](zbrac.f)  | outward search for brackets on roots |
|  9.1  | [`zbrak `](zbrak.f)  | inward search for brackets on roots |
|  9.1  | [`rtbis `](rtbis.f)  | find root of a function by bisection |
|  9.2  | [`rtflsp`](rtflsp.f) | find root of a function by false-position |
|  9.2  | [`rtsec `](rtsec.f)  | find root of a function by secant method |
|  9.2  | [`zriddr`](zriddr.f) | find root of a function by Ridders’ method |
|  9.3  | [`zbrent`](zbrent.f) | find root of a function by Brent’s method |
|  9.4  | [`rtnewt`](rtnewt.f) | find root of a function by Newton-Raphson |
|  9.4  | [`rtsafe`](rtsafe.f) | find root of a function by Newton-Raphson and bisection |
|  9.5  | [`laguer`](laguer.f) | find a root of a polynomial by Laguerre’s method |
|  9.5  | [`zroots`](zroots.f) | roots of a polynomial by Laguerre’s method with deflation |
|  9.5  | [`zrhqr `](zrhqr.f)  | roots of a polynomial by eigenvalue methods |
|  9.5  | [`qroot `](qroot.f)  | complex or double root of a polynomial, Bairstow |
|  9.6  | [`mnewt `](mnewt.f)  | Newton’s method for systems of equations |
|  9.7  | [`lnsrch`](lnsrch.f) | search along a line, used by newt |
|  9.7  | [`newt  `](newt.f)   | globally convergent multi-dimensionalNewton’s method |
|  9.7  | [`fdjac `](fdjac.f)  | finite-difference Jacobian, used by newt |
|  9.7  | [`fmin  `](fmin.f)   | norm of a vector function, used by newt |
|  9.7  | [`broydn`](broydn.f) | secant method for systems of equations |
| 10.1  | [`mnbrak`](mnbrak.f) | bracket the minimum of a function |
| 10.1  | [`golden`](golden.f) | find minimum of a function by golden section search |
| 10.2  | [`brent `](brent.f)  | find minimum of a function by Brent’s method |
| 10.3  | [`dbrent`](dbrent.f) | find minimum of a function using derivative information |
| 10.4  | [`amoeba`](amoeba.f) | minimize in N-dimensions by downhill simplex method |
| 10.4  | [`amotry`](amotry.f) | evaluate a trial point, used by amoeba |
| 10.5  | [`powell`](powell.f) | minimize in N-dimensions by Powell’s method |
| 10.5  | [`linmin`](linmin.f) | minimum of a function along a ray in N-dimensions |
| 10.5  | [`f1dim `](f1dim.f)  | function used by linmin |
| 10.6  | [`frprmn`](frprmn.f) | minimize in N-dimensions by conjugate gradient |
| 10.6  | [`df1dim`](df1dim.f) | alternative function used by linmin |
| 10.7  | [`dfpmin`](dfpmin.f) | minimize in N-dimensions by variable metric method |
| 10.8  | [`simplx`](simplx.f) | linear programming maximization of a linear function |
| 10.8  | [`simp1 `](simp1.f)  | linear programming, used by simplx |
| 10.8  | [`simp2 `](simp2.f)  | linear programming, used by simplx |
| 10.8  | [`simp3 `](simp3.f)  | linear programming, used by simplx |
| 10.9  | [`anneal`](anneal.f) | traveling salesman problem by simulated annealing |
| 10.9  | [`revcst`](revcst.f) | cost of a reversal, used by anneal |
| 10.9  | [`revers`](revers.f) | do a reversal, used by anneal |
| 10.9  | [`trncst`](trncst.f) | cost of a transposition, used by anneal |
| 10.9  | [`trnspt`](trnspt.f) | do a transposition, used by anneal |
| 10.9  | [`metrop`](metrop.f) | Metropolis algorithm, used by anneal |
| 10.9  | [`amebsa`](amebsa.f) | simulated annealing in continuous spaces |
| 10.9  | [`amotsa`](amotsa.f) | evaluate a trial point, used by amebsa |
| 11.1  | [`jacobi`](jacobi.f) | eigenvalues and eigenvectors of a symmetric matrix |
| 11.1  | [`eigsrt`](eigsrt.f) | eigenvectors, sorts into order by eigenvalue |
| 11.2  | [`tred2 `](tred2.f)  | Householder reduction of a real, symmetric matrix |
| 11.3  | [`tqli  `](tqli.f)   | eigensolution of a symmetric tridiagonal matrix |
| 11.5  | [`balanc`](balanc.f) | balance a nonsymmetric matrix |
| 11.5  | [`elmhes`](elmhes.f) | reduce a general matrix to Hessenberg form |
| 11.6  | [`hqr   `](hqr.f)    | eigenvalues of a Hessenberg matrix |
| 12.2  | [`four1 `](four1.f)  | fast Fourier transform (FFT) in one dimension |
| 12.3  | [`twofft`](twofft.f) | fast Fourier transform of two real functions |
| 12.3  | [`realft`](realft.f) | fast Fourier transform of a single real function |
| 12.3  | [`sinft `](sinft.f)  | fast sine transform |
| 12.3  | [`cosft1`](cosft1.f) | fast cosine transform with endpoints |
| 12.3  | [`cosft2`](cosft2.f) | “staggered” fast cosine transform |
| 12.4  | [`fourn `](fourn.f)  | fast Fourier transform in multidimensions |
| 12.5  | [`rlft3 `](rlft3.f)  | FFT of real data in two or three dimensions |
| 12.6  | [`fourfs`](fourfs.f) | FFT for huge data sets on external media |
| 12.6  | [`fourew`](fourew.f) | rewind and permute files, used by fourfs |
| 13.1  | [`convlv`](convlv.f) | convolution or deconvolution of data using FFT |
| 13.2  | [`correl`](correl.f) | correlation or autocorrelation of data using FFT |
| 13.4  | [`spctrm`](spctrm.f) | power spectrum estimation using FFT |
| 13.6  | [`memcof`](memcof.f) | evaluate maximum entropy (MEM) coefficients |
| 13.6  | [`fixrts`](fixrts.f) | reflect roots of a polynomial into unit circle |
| 13.6  | [`predic`](predic.f) | linear prediction using MEM coefficients |
| 13.7  | [`evlmem`](evlmem.f) | power spectral estimation from MEM coefficients |
| 13.8  | [`period`](period.f) | power spectrum of unevenly sampled data |
| 13.8  | [`fasper`](fasper.f) | power spectrum of unevenly sampled larger data sets |
| 13.8  | [`spread`](spread.f) | extirpolate value into array, used by fasper |
| 13.9  | [`dftcor`](dftcor.f) | compute endpoint corrections for Fourier integrals |
| 13.9  | [`dftint`](dftint.f) | high-accuracy Fourier integrals |
| 13.10 | [`wt1   `](wt1.f)    | one-dimensional discrete wavelet transform |
| 13.10 | [`daub4 `](daub4.f)  | Daubechies 4-coefficient wavelet filter |
| 13.10 | [`pwtset`](pwtset.f) | initialize coefficients for pwt |
| 13.10 | [`pwt   `](pwt.f)    | partial wavelet transform |
| 13.10 | [`wtn   `](wtn.f)    | multidimensional discrete wavelet transform |
| 14.1  | [`moment`](moment.f) | calculate moments of a data set |
| 14.2  | [`ttest `](ttest.f)  | Student’s t-test for difference of means |
| 14.2  | [`avevar`](avevar.f) | calculate mean and variance of a data set |
| 14.2  | [`tutest`](tutest.f) | Student’s t-test for means, case of unequal variances |
| 14.2  | [`tptest`](tptest.f) | Student’s t-test for means, case of paired data |
| 14.2  | [`ftest `](ftest.f)  | F-test for difference of variances |
| 14.3  | [`chsone`](chsone.f) | chi-square test for difference between data and model |
| 14.3  | [`chstwo`](chstwo.f) | chi-square test for difference between two data sets |
| 14.3  | [`ksone `](ksone.f)  | Kolmogorov-Smirnov test of data against model |
| 14.3  | [`kstwo `](kstwo.f)  | Kolmogorov-Smirnov test between two data sets |
| 14.3  | [`probks`](probks.f) | Kolmogorov-Smirnov probability function |
| 14.4  | [`cntab1`](cntab1.f) | contingency table analysis using chi-square |
| 14.4  | [`cntab2`](cntab2.f) | contingency table analysis using entropy measure |
| 14.5  | [`pearsn`](pearsn.f) | Pearson’s correlation between two data sets |
| 14.6  | [`spear `](spear.f)  | Spearman’s rank correlation between two data sets |
| 14.6  | [`crank `](crank.f)  | replaces array elements by their rank |
| 14.6  | [`kendl1`](kendl1.f) | correlation between two data sets, Kendall’s tau |
| 14.6  | [`kendl2`](kendl2.f) | contingency table analysis using Kendall’s tau |
| 14.7  | [`ks2d1s`](ks2d1s.f) | KS test in two dimensions, data vs. model |
| 14.7  | [`quadct`](quadct.f) | count points by quadrants, used by ks2d1s |
| 14.7  | [`quadvl`](quadvl.f) | quadrant probabilities, used by ks2d1s |
| 14.7  | [`ks2d2s`](ks2d2s.f) | KS test in two dimensions, data vs. data |
| 14.8  | [`savgol`](savgol.f) | Savitzky-Golay smoothing coefficients |
| 15.2  | [`fit   `](fit.f)    | least-squares fit data to a straight line |
| 15.3  | [`fitexy`](fitexy.f) | fit data to a straight line, errors in both x and y |
| 15.3  | [`chixy `](chixy.f)  | used by fitexy to calculate a _2 |
| 15.4  | [`lfit  `](lfit.f)   | general linear least-squares fit by normal equations |
| 15.4  | [`covsrt`](covsrt.f) | rearrange covariance matrix, used by lfit |
| 15.4  | [`svdfit`](svdfit.f) | linear least-squares fit by singular value decomposition |
| 15.4  | [`svdvar`](svdvar.f) | variances from singular value decomposition |
| 15.4  | [`fpoly `](fpoly.f)  | fit a polynomial using lfit or svdfit |
| 15.4  | [`fleg  `](fleg.f)   | fit a Legendre polynomial using lfit or svdfit |
| 15.5  | [`mrqmin`](mrqmin.f) | nonlinear least-squares fit, Marquardt’s method |
| 15.5  | [`mrqcof`](mrqcof.f) | used by mrqmin to evaluate coefficients |
| 15.5  | [`fgauss`](fgauss.f) | fit a sum of Gaussians using mrqmin |
| 15.7  | [`medfit`](medfit.f) | fit data to a straight line robustly, least absolute deviation |
| 15.7  | [`rofunc`](rofunc.f) | fit data robustly, used by medfit |
| 16.1  | [`rk4   `](rk4.f)    | integrate one step of ODEs, fourth-order Runge-Kutta |
| 16.1  | [`rkdumb`](rkdumb.f) | integrate ODEs by fourth-order Runge-Kutta |
| 16.2  | [`rkqs  `](rkqs.f)   | integrate one step of ODEs with accuracy monitoring |
| 16.2  | [`rkck  `](rkck.f)   | Cash-Karp-Runge-Kutta step used by rkqs |
| 16.2  | [`odeint`](odeint.f) | integrate ODEs with accuracy monitoring |
| 16.3  | [`mmid  `](mmid.f)   | integrate ODEs by modified midpoint method |
| 16.4  | [`bsstep`](bsstep.f) | integrate ODEs, Bulirsch-Stoer step |
| 16.4  | [`pzextr`](pzextr.f) | polynomial extrapolation, used by bsstep |
| 16.4  | [`rzextr`](rzextr.f) | rational function extrapolation, used by bsstep |
| 16.5  | [`stoerm`](stoerm.f) | integrate conservative second-order ODEs |
| 16.6  | [`stiff `](stiff.f)  | integrate stiff ODEs by fourth-order Rosenbrock |
| 16.6  | [`jacobn`](jacobn.f) | sample Jacobian routine for stiff |
| 16.6  | [`derivs`](derivs.f) | sample derivatives routine for stiff |
| 16.6  | [`simpr `](simpr.f)  | integrate stiff ODEs by semi-implicit midpoint rule |
| 16.6  | [`stifbs`](stifbs.f) | integrate stiff ODEs, Bulirsch-Stoer step |
| 17.1  | [`shoot `](shoot.f)  | solve two point boundary value problem by shooting |
| 17.2  | [`shootf`](shootf.f) | ditto, by shooting to a fitting point |
| 17.3  | [`solvde`](solvde.f) | two point boundary value problem, solve by relaxation |
| 17.3  | [`bksub `](bksub.f)  | backsubstitution, used by solvde |
| 17.3  | [`pinvs `](pinvs.f)  | diagonalize a sub-block, used by solvde |
| 17.3  | [`red   `](red.f)    | reduce columns of a matrix, used by solvde |
| 17.4  | [`sfroid`](sfroid.f) | spheroidal functions by method of solvde |
| 17.4  | [`difeq `](difeq.f)  | spheroidal matrix coefficients, used by sfroid |
| 17.4  | [`sphoot`](sphoot.f) | spheroidal functions by method of shoot |
| 17.4  | [`sphfpt`](sphfpt.f) | spheroidal functions by method of shootf |
| 18.1  | [`fred2 `](fred2.f)  | solve linear Fredholm equations of the second kind |
| 18.1  | [`fredin`](fredin.f) | interpolate solutions obtained with fred2 |
| 18.2  | [`voltra`](voltra.f) | linear Volterra equations of the second kind |
| 18.3  | [`wwghts`](wwghts.f) | quadrature weights for an arbitrarily singular kernel |
| 18.3  | [`kermom`](kermom.f) | sample routine for moments of a singular kernel |
| 18.3  | [`quadmx`](quadmx.f) | sample routine for a quadrature matrix |
| 18.3  | [`fredex`](fredex.f) | example of solving a singular Fredholm equation |
| 19.5  | [`sor   `](sor.f)    | elliptic PDE solved by successive overrelaxation method |
| 19.6  | [`mglin `](mglin.f)  | linear elliptic PDE solved by multigrid method |
| 19.6  | [`rstrct`](rstrct.f) | half-weighting restriction, used by mglin, mgfas |
| 19.6  | [`interp`](interp.f) | bilinear prolongation, used by mglin, mgfas |
| 19.6  | [`addint`](addint.f) | interpolate and add, used by mglin |
| 19.6  | [`slvsml`](slvsml.f) | solve on coarsest grid, used by mglin |
| 19.6  | [`relax `](relax.f)  | Gauss-Seidel relaxation, used by mglin |
| 19.6  | [`resid `](resid.f)  | calculate residual, used by mglin |
| 19.6  | [`copy  `](copy.f)   | utility used by mglin, mgfas |
| 19.6  | [`fill0 `](fill0.f)  | utility used by mglin |
| 19.6  | [`maloc `](maloc.f)  | memory allocation utility used by mglin, mgfas |
| 19.6  | [`mgfas `](mgfas.f)  | nonlinear elliptic PDE solved by multigrid method |
| 19.6  | [`relax2`](relax2.f) | Gauss-Seidel relaxation, used by mgfas |
| 19.6  | [`slvsm2`](slvsm2.f) | solve on coarsest grid, used by mgfas |
| 19.6  | [`lop   `](lop.f)    | applies nonlinear operator, used by mgfas |
| 19.6  | [`matadd`](matadd.f) | utility used by mgfas |
| 19.6  | [`matsub`](matsub.f) | utility used by mgfas |
| 19.6  | [`anorm2`](anorm2.f) | utility used by mgfas |
| 20.1  | [`machar`](machar.f) | diagnose computer’s floating arithmetic |
| 20.2  | [`igray `](igray.f)  | Gray code and its inverse |
| 20.3  | [`icrc1 `](icrc1.f)  | cyclic redundancy checksum, used by icrc |
| 20.3  | [`icrc  `](icrc.f)   | cyclic redundancy checksum |
| 20.3  | [`decchk`](decchk.f) | decimal check digit calculation or verification |
| 20.4  | [`hufmak`](hufmak.f) | construct a Huffman code |
| 20.4  | [`hufapp`](hufapp.f) | append bits to a Huffman code, used by hufmak |
| 20.4  | [`hufenc`](hufenc.f) | use Huffman code to encode and compress a character |
| 20.4  | [`hufdec`](hufdec.f) | use Huffman code to decode and decompress a character |
| 20.5  | [`arcmak`](arcmak.f) | construct an arithmetic code |
| 20.5  | [`arcode`](arcode.f) | encode or decode a character using arithmetic coding |
| 20.5  | [`arcsum`](arcsum.f) | add integer to byte string, used by arcode |
| 20.6  | [`mpops `](mpops.f)  | multiple precision arithmetic, simpler operations |
| 20.6  | [`mpmul `](mpmul.f)  | multiple precision multiply, using FFT methods |
| 20.6  | [`mpinv `](mpinv.f)  | multiple precision reciprocal |
| 20.6  | [`mpdiv `](mpdiv.f)  | multiple precision divide and remainder |
| 20.6  | [`mpsqrt`](mpsqrt.f) | multiple precision square root |
| 20.6  | [`mp2dfr`](mp2dfr.f) | multiple precision conversion to decimal base |
| 20.6  | [`mppi  `](mppi.f)   | multiple precision example, compute many digits of _ |
