const { HalfInfiniteIntervalIntegration } = require("./HalfInfiniteIntervalIntegration");
const { DefiniteIntervalIntegration } = require("./DefiniteIntervalIntegration")

function NumericalIntegration(f, a, b, out = {})
{
	if (b === Infinity)
		HalfInfiniteIntervalIntegration(f, a, out);
	else if (a === -Infinity)
		HalfInfiniteIntervalIntegration(x => f(-x), b, out);
	else
		DefiniteIntervalIntegration(f, a, b, out);
	return out.approximation;
}

module.exports = { NumericalIntegration }