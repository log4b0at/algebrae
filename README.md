# Algebrae
Javascript implementations of various mathematical tools

## Functions

### Special Functions

| Name 		| Formula	| Function 		| Complexity| Average Error*	| Worst Error*	|
|-----------|-----------|---------------|-----------|-------------------|---------------|
| Gamma		| ![FGamma]	| `gamma[1p]`	| O(1)		| -16.023			| -15.352		|
| Pi		| ![FSinPI] | `[sin]pi`		| O(1)		| N/A				| N/A			|

### Discrete Functions

| Name 					| Formula	| Function 			| Complexity		| BigInt	|
|-----------------------|-----------|-------------------|-------------------|-----------|
| BinomialCoefficient	| ![FBinCo]	| `binomial`		| O(min(k, n - k))	| yes		|

<!-- Website for formulas: http://www.sciweavers.org/free-online-latex-equation-editor -->

[FGamma]: ./docs/FGamma.png "Gamma"

[FSinPi]: ./docs/FSinPi.png "SinPi"
[FBinCo]: ./docs/FBinomialCoefficient.png "Binomial"