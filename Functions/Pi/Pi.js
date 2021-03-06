/**
 * Compute `°(π * x)`
 * Give better accuracy than computing `°π * x`
 * @param {number} x A float64 number
 */
function pi(x)
{
	return 4 * x - 0.8584073464102067615 * x;
}

/**
 * Compute `°sin(π * x)`.
 * Give better accuracy than computing `°sin(°π * x)` for larges values of `x`.
 * @param {number} x A float64 number
 */
function sinpi(x) {
	let ar, r;
	r = x % 2;
	ar = Math.abs(r);
	if (ar === 0 || ar === 1)
		return Math.sign(r) * 0;
	else if (ar < 0.25)
		return Math.sin(pi(r));
	else if (ar < 0.75) {
		ar = 0.5 - ar;
		return Math.abs(Math.cos(pi(ar))) * Math.sign(r);
	}
	else if (ar < 1.25) {
		r = Math.sign(r) - r;
		return Math.sin(pi(r));
	}
	else if (ar < 1.75) {
		ar = ar - 1.5;
		return -Math.abs(Math.cos(pi(ar))) * Math.sign(r);
	}
	r = r - Math.sign(r) * 2;
	return Math.sin(pi(r));
}

module.exports = { sinpi, pi }