/**
 * Compute [binomial coefficients](https://en.wikipedia.org/wiki/Binomial_coefficient) for `n` and `k`.
 * It is the coefficient of the `x^k` term in the polynomial expansion of the binomial power `(1 + x)^n`
 * @template {number|bigint} T An integer type, either bigint or number
 * @param {T} n An integer
 * @param {T} k An integer
 * @returns {T} 
 * @complexity `O(min(k, n - k))`
 */
function binomial(n, k)
{
	if (typeof n === "number")
		return int_binomial(n, k);
	else if (typeof n === "bigint")
		return big_binomial(n, k);
	else
		throw new TypeError("'n' paremeter type isn't compatible with template for this function, should be number or bigint.");
}

function int_binomial(n, k) {
	if (k * 2 > n)
		k = n - k;
	if (k <= 0)
		return k === 0 ? 1 : NaN;
	let r = n++; let d = 1;
	for (let i = 2; i <= k; i++) {
		let pr = r * (n - i)
		if (pr > Number.MAX_SAFE_INTEGER) {
			r = r / d * (n - i);
			d = i;
			continue;
		}
		r = pr;
		d *= i;
	}
	return r / d;
}

function big_binomial(n, k) {
	if ((k << 1n) > n)
		k = n - k;
	if (k <= 0n)
		return k === 0n ? 1n : NaN;
	let r = n++;
	let d = 1n;
	for (let i = 2; i <= k; i++) {
		const ii = BigInt(i)
		r *= (n - ii)
		d *= ii;
	}
	return r / d;
}

module.exports = { binomial }