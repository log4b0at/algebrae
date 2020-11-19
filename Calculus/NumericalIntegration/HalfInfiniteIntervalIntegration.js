
let function_call;

//Abscissas and weights for the quadrature
const aw = new Array(2000);

Initialization(2000, 1e-307, 1e-16, aw);

function Initialization(lenaw, tiny, eps, aw) {
	/* ---- adjustable parameter ---- */
	let efs = 0.1, hoff = 11.0;
	/* ------------------------------ */
	let noff, nk, k, j;
	let pi4, tinyln, epsln, h0, ehp, ehm, h, t, ep, em, xp, xm,
		wp, wm;

	pi4 = Math.atan(1.0);
	tinyln = -Math.log(tiny);
	epsln = 1 - Math.log(efs * eps);
	h0 = hoff / epsln;
	ehp = Math.exp(h0);
	ehm = 1 / ehp;

	aw[2] = eps;
	aw[3] = Math.exp(-ehm * epsln);
	aw[4] = Math.sqrt(efs * eps);
	noff = 5;
	aw[noff] = 1;
	aw[noff + 1] = 4 * h0;
	aw[noff + 2] = 2 * pi4 * h0;
	h = 2;
	nk = 0;
	k = noff + 6;
	do {
		t = h * 0.5;
		do {
			em = Math.exp(h0 * t);
			ep = pi4 * em;
			em = pi4 / em;
			j = k;
			do {
				xp = Math.exp(ep - em);
				xm = 1 / xp;
				wp = xp * ((ep + em) * h0);
				wm = xm * ((ep + em) * h0);
				aw[j] = xm;
				aw[j + 1] = xp;
				aw[j + 2] = xm * (4 * h0);
				aw[j + 3] = xp * (4 * h0);
				aw[j + 4] = wm;
				aw[j + 5] = wp;
				ep *= ehp;
				em *= ehm;
				j += 6;
			} while (ep < tinyln && j <= lenaw - 6);
			t += h;
			k += nk;
		} while (t < 1);
		h *= 0.5;
		if (nk == 0) {
			if (j > lenaw - 12) j -= 6;
			nk = j - noff;
			k += nk;
			aw[1] = nk;
		}
	} while (2 * k - noff - 6 <= lenaw);
	aw[0] = k - 6;
}


function HalfInfiniteIntervalIntegration(f, a, out)
{
	let noff, lenawm, nk, k, j, jtmp, jm, m, klim;
	let epsh, ir, fp, fm, errt, errh, errd, h, iback, irback;

	function call(f, x, d)
	{
		function delta(value)
		{
			return d + (value + a);
		}
		const v = f(x, delta)
		if (!Number.isFinite(v))
			return null;
		function_call++;
		return v;
	}

	noff = 5;
	lenawm = Math.floor(aw[0] + 0.5);
	nk = Math.floor(aw[1] + 0.5);
	epsh = aw[4];
	out.approximation = f(a + aw[noff]);
	ir = out.approximation * aw[noff + 1];
	out.approximation *= aw[noff + 2];
	out.relative_error = Math.abs(out.approximation);
	k = nk + noff;
	j = noff;
	do {
		j += 6;
		if ((fm = call(f, a + aw[j])) === null) break;
		if ((fp = call(f, a + aw[j + 1])) === null) break;
		ir += fm * aw[j + 2] + fp * aw[j + 3];
		fm *= aw[j + 4];
		fp *= aw[j + 5];
		out.approximation += fm + fp;
		out.relative_error += Math.abs(fm) + Math.abs(fp);
	} while (aw[j] > epsh && j < k);
	errt = out.relative_error * aw[3];
	errh = out.relative_error * epsh;
	errd = 1 + 2 * errh;
	jtmp = j;
	while (Math.abs(fm) > errt && j < k) {
		j += 6;
		if ((fm = call(f, a + aw[j])) === null) break;
		ir += fm * aw[j + 2];
		fm *= aw[j + 4];
		out.approximation += fm;
	}
	jm = j;
	j = jtmp;
	while (Math.abs(fp) > errt && j < k) {
		j += 6;
		if ((fp = call(f, a + aw[j + 1])) === null) break;
		ir += fp * aw[j + 3];
		fp *= aw[j + 5];
		out.approximation += fp;
	}
	if (j < jm) jm = j;
	jm -= noff + 6;
	h = 1;
	m = 1;
	klim = k + nk;
	while (errd > errh && klim <= lenawm) {
		iback = out.approximation;
		irback = ir;
		do {
			jtmp = k + jm;
			for (j = k + 6; j <= jtmp; j += 6) {
				if ((fm = call(f, a + aw[j])) === null) break;
				if ((fp = call(f, a + aw[j + 1])) === null) break;
				ir += fm * aw[j + 2] + fp * aw[j + 3];
				out.approximation += fm * aw[j + 4] + fp * aw[j + 5];
			}
			k += nk;
			j = jtmp;
			do {
				j += 6;
				if ((fm = call(f, a + aw[j])) === null) break;
				ir += fm * aw[j + 2];
				fm *= aw[j + 4];
				out.approximation += fm;
			} while (Math.abs(fm) > errt && j < k);
			j = jtmp;
			do {
				j += 6;
				if ((fp = call(f, a + aw[j + 1])) === null) break;
				ir += fp * aw[j + 3];
				fp *= aw[j + 5];
				out.approximation += fp;
			} while (Math.abs(fp) > errt && j < k);
		} while (k < klim);
		errd = h * (Math.abs(out.approximation - 2 * iback) + Math.abs(ir - 2 * irback));
		h *= 0.5;
		m *= 2;
		klim = 2 * klim - noff;
	}
	out.approximation *= h;
	if (errd > errh) {
		out.relative_error = -errd * m;
	} else {
		out.relative_error *= aw[2] * m;
	}
}

module.exports = { HalfInfiniteIntervalIntegration }