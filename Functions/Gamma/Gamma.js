const { sinpi, pi } = require('../Pi/Pi')

function binet(x)
{
	const ix = 1 / x;
	const ix2 = ix * ix;
	const ix4 = ix2 * ix2;
	const ix6 = ix4 * ix2;
	const ix8 = ix4 * ix4;
	const ix10 = ix4 * ix6;
	const ix12 = ix6 * ix6;

	return (
		(ix / 12.) + ix * (- ix2 / 360. + (ix4 / 1260. + (
		- ix6 / 1680. + (ix8 / 1188. + (-(ix10*691) / 360360 
		+ ix12 / 156.)))))
	);
}

function pole(x)
{
	return 1 / (((((((((((((((((((
		(7.782263439905071e-12*x + 1.0434267116911005e-10) * x + -1.18127457048702e-9) * x + 5.002007644469223e-9) * x 
		+ 6.116095104481416e-9) * x + -2.056338416977607e-7) * x + 0.000001133027231981696) * x + -0.0000012504934821426706) * x 
		+ -0.00002013485478078824) * x + 0.0001280502823881162) * x + -0.00021524167411495098) * x + -0.0011651675918590652) * x 
		+ 0.0072189432466631) * x + -0.009621971527876973) * x + -0.04219773455554433) * x + 0.16653861138229148) * x 
		+ -0.04200263503409524) * x + -0.6558780715202539) * x + 0.5772156649015329) * x + 1) * x)
}

function positive(x) {
	let o = x;
	if (x < 10)
		x += 10;
	let b = binet(x), offset = x, power = x + .5;
	if (o > 142)
		power /= 2;
	let result = Math.pow(x, power);
	let last = result;
	result = result * Math.expm1(b) + result;
	result /= Math.exp(offset);
	if (o > 142)
		result *= last;
	result = result * 0.506628274631000502415765 + result * 2;
	if (o < 10)
		result /= (o + 1) * (o + 2) * (o + 3) * (o + 4) * (o + 5) * (o + 6) * (o + 7) * (o + 8) * (o + 9) * (o + 10);
	return result;
}

function reflection(x)
{
	return pi(x) / (sinpi(x) * positive(-x));
}

/**
 * Compute `°𝚪(x)` where 𝚪 is Euler gamma function
 * @param {number} x A float64 number
 */
function gamma(x)
{
	if (0 < x && x < 0.5)
		return pole(x);
	else
		return gamma1p(x - 1);
}

/**
 * Compute `°𝚪(x + 1)` where 𝚪 is Euler gamma function
 * @param {number} x A float64 number
 */
function gamma1p(x) {
	if (x === -1)
		return NaN;
	else if (x < 0)
		return reflection(x);
	else
		return positive(x);
}

module.exports = { gamma, gamma1p }