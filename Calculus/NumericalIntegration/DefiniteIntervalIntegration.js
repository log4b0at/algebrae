const { pi } = require('../../Functions/Pi/Pi')
let function_call;

//Abscissas and weights for the quadrature
const aw = new Array(2000);

Initialization(2000, 1e-307, 1e-16, aw);

function Initialization(lenaw, tiny, eps) {
	/* ---- adjustable parameter ---- */
	let efs = 0.1, hoff = 8.5 * 2;
	/* ------------------------------ */
	let noff, nk, k, j;
	let pi2, tinyln, epsln, h0, ehp, ehm, h, t, ep, em, xw, wg;

	pi2 = 2 * Math.atan(1.0);
	tinyln = -Math.log(tiny);
	epsln = 1 - Math.log(efs * eps);
	h0 = hoff / epsln;
	ehp = Math.exp(h0);
	ehm = 1 / ehp;
	aw[2] = eps;
	aw[3] = Math.exp(-ehm * epsln);
	aw[4] = Math.sqrt(efs * eps);
	noff = 5;
	aw[noff] = 0.5;
	aw[noff + 1] = h0;
	aw[noff + 2] = pi(h0 / 4);
	h = 2;
	nk = 0;
	k = noff + 3;
	do {
		t = h * 0.5;
		do {
			em = Math.exp(h0 * t);
			ep = pi(em) / 2;
			em = pi(1/(2*em));
			j = k;
			do {
				const ui = Math.exp(ep)
				const uj = Math.exp(em);
				xw = uj / (uj + ui);
				wg = xw * (1 - xw) * h0;
				aw[j] = xw;
				aw[j + 1] = wg * 4;
				aw[j + 2] = wg * ep + wg * em;
				ep *= ehp;
				em *= ehm;
				j += 3;
			} while (ep < tinyln && j <= lenaw - 3);
			t += h;
			k += nk;
		} while (t < 1);
		h *= 0.5;
		if (nk == 0) {
			if (j > lenaw - 6) j -= 3;
			nk = j - noff;
			k += nk;
			aw[1] = nk;
		}
	} while (2 * k - noff - 3 <= lenaw);
	aw[0] = k - 3;
}

function DefiniteIntervalIntegration(f, a, b, out) {
	
	function call(f, x, d)
	{
		const c = Math.abs(x - a) < Math.abs(b - x);
		function delta(value)
		{
			return c ? d + (value + a) : (value + b) - d;
		}
		const v = f(x, delta)
		if (!Number.isFinite(v))
			return null;
		function_call++;
		return v;
	}

	let noff, lenawm, nk, k, j, jtmp, jm, m, klim;
	let epsh, ba, ir, xa, fa, fb, errt, errh, errd, h, iback, irback;

	function_call = 0;

	noff = 5;
	lenawm = Math.floor(aw[0] + 0.5);
	nk = Math.floor(aw[1] + 0.5);
	epsh = aw[4];
	ba = b - a;
	out.approximation = call(f, (a + b) * aw[noff], ba * (aw[noff]));
	ir = out.approximation * aw[noff + 1];
	out.approximation *= aw[noff + 2];
	out.relative_error = Math.abs(out.approximation);
	k = nk + noff;
	j = noff;
	do {
		j += 3;
		xa = ba * aw[j];
		if ((fa = call(f, a + xa, xa)) === null) break;
		if ((fb = call(f, b - xa, xa)) === null) break;
		ir += (fa + fb) * aw[j + 1];
		fa *= aw[j + 2];
		fb *= aw[j + 2];
		out.approximation += fa + fb;
		out.relative_error += Math.abs(fa) + Math.abs(fb);
	} while (aw[j] > epsh && j < k);
	errt = out.relative_error * aw[3];
	errh = out.relative_error * epsh;
	errd = 1 + 2 * errh;
	jtmp = j;
	while (Math.abs(fa) > errt && j < k) {
		j += 3;
		const d = ba * aw[j]
		if ((fa = call(f, a + d, d)) === null) break;
		ir += fa * aw[j + 1];
		fa *= aw[j + 2];
		out.approximation += fa;
	}
	jm = j;
	j = jtmp;
	while (Math.abs(fb) > errt && j < k) {
		j += 3;
		const d = ba * aw[j];
		if ((fb = call(f, b - d, d)) === null) break;
		ir += fb * aw[j + 1];
		fb *= aw[j + 2];
		out.approximation += fb;
	}
	if (j < jm) jm = j;
	jm -= noff + 3;
	h = 1;
	m = 1;
	klim = k + nk;
	while (errd > errh && klim <= lenawm) {
		iback = out.approximation;
		irback = ir;
		do {
			jtmp = k + jm;
			for (j = k + 3; j <= jtmp; j += 3) {
				xa = ba * aw[j];
				if ((fa = call(f, a + xa, xa)) === null) break;
				if ((fb = call(f, b - xa, xa)) === null) break;
				ir += (fa + fb) * aw[j + 1];
				out.approximation += (fa + fb) * aw[j + 2];
			}
			k += nk;
			j = jtmp;
			do {
				j += 3;
				const d = ba * aw[j]
				if ((fa = call(f, a + d, d)) === null) break;
				ir += fa * aw[j + 1];
				fa *= aw[j + 2];
				out.approximation += fa;
			} while (Math.abs(fa) > errt && j < k);
			j = jtmp;
			do {
				j += 3;
				const d = ba * aw[j];
				if ((fb = call(f, b - d, d)) === null) break;
				ir += fb * aw[j + 1];
				fb *= aw[j + 2];
				out.approximation += fb;
			} while (Math.abs(fb) > errt && j < k);
		} while (k < klim);
		errd = h * (Math.abs(out.approximation - 2 * iback) + Math.abs(ir - 2 * irback));
		h *= 0.5;
		m *= 2;
		klim = 2 * klim - noff;
	}
	out.call_count = function_call;
	out.approximation *= h * ba;
	if (errd > errh) {
		out.relative_error = -errd * (m * Math.abs(ba));
	} else {
		out.relative_error = out.relative_error * aw[2] * (m * Math.abs(ba));
	}
}

module.exports = { DefiniteIntervalIntegration }