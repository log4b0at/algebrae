const { performance } = require('perf_hooks')
const { NumericalIntegration } = require('./NumericalIntegration');

const { sqrt, pow, log, exp, abs, sin, round, asin } = Math;

function quadradure(f, a, b, out)
{
	return NumericalIntegration(f, a, b, out)
}

let i = 0, log_average_err = 0, functions_evaluations = [], functions_types = [];

function add(f, a, b, exact) {
	const start = performance.now();
	const out = {}
	const result = quadradure(f, a, b, out);
	const end = performance.now();
	const error = abs(result/exact - 1) || 2 ** -54;
	const time = end - start;
	functions_evaluations.push({
		f:f.toString(), 
		from: a,
		to: b,
		exact, 
		result,
		error,
		time,
		i,
		call_count: out.call_count,
		measured_error: out.relative_error
	});
	log_average_err += log(error)/log(10);
	i++;
	return error;
}

function end_function_type(title)
{
	console.log('\n********  ',title, '  ********\n')
	console.table(functions_evaluations.map( ({f, from, to, exact, result, error, time, call_count, measured_error }) => { 
		let str = f.toString();
		str = f.toString().substr(str.indexOf('=>') + 3, 40) + (str.length > 40 ? "..." : "")
		return {
			"function": str , 
			from,
			to,
			"closest ieee754":exact, 
			"quadrature result":result,
			"real error": parseFloat(error.toPrecision(3)),
			"µs": round(time*1000),
			"abs error": parseFloat(measured_error.toPrecision(3)),
			"rel error": parseFloat((measured_error/result).toPrecision(3)),
			"calls": call_count,
		};}
	));
	
	const min = functions_evaluations.reduce((acc, e) => e.error < acc.error ? e : acc, functions_evaluations[0])
	const max = functions_evaluations.reduce((acc, e) => e.error > acc.error ? e : acc, functions_evaluations[0])
	const average_time = functions_evaluations.reduce((acc, e) => acc + e.time, 0) / functions_evaluations.length;
	const average_log_error = functions_evaluations.reduce((acc, e) => acc + log(e.error), 0) / functions_evaluations.length;
	const variance =  functions_evaluations.reduce((acc, e) => acc + (log(e.error) - average_log_error)**2, 0 ) / functions_evaluations.length
	functions_types.push({
		title, min, max, average_time, average_log_error, variance
	})
	functions_evaluations = [];
	i = 0;
}

function end_test()
{
	console.log('\n\n\n\n');
	console.table(functions_types.map(({variance, average_log_error, min, max, average_time, title}) => {return {
		"function type": title,
		"average log10 error": parseFloat((average_log_error / log(10)).toFixed(3)),
		"average µs": round(average_time * 1000),
		"best": min.i,
		"worst": max.i,
		"σ log10 error": parseFloat(sqrt(variance).toPrecision(3))
	}}))
}


// this quadrature is used for node.js cache
quadradure(x => log(x + 1/2), 0, 1)

add(x => sqrt(x), 0, 1, 2/3)


add((x, X) => 1 / sqrt( X(1) * -X(-1) ), -1, 1, Math.PI)
//add((x, X, P) => 1 / sqrt( P(-1, 0, 1) ), -1, 1, Math.PI)
//add((x, a_x, b_x) => 1 / sqrt( -a_x * b_x ), 0, 1, Math.PI/2)
// (x, a_x, b_x) => 1 / sqrt( a_x * -b_x )

//add((x, a_x, b_x) => sqrt(x) / sqrt(x * b_x + b_x), 0, 1, 1.19814023473559220744)
add(x => sqrt(1 - x**2), -1, 1, Math.PI/2)

add(x => log(x), 0, 1, -1)
add(x => x * log(x), 0, 1, -0.25)
add(x => pow(x, 1/20), 0, 1, 1 / (1 + 1/20) )
add((x, X) => pow(-X(-1), -1/5) + pow(x, -1/5), 0, 1, 2.5 )
add(x => 1 / sqrt(3*x*x - 5*x), -1, 0, 0.822964855789271638)
add(x => 1/x, 0.000001, 0.000002, Math.LN2 )
add(x => 1/x, 10000000, 20000000, Math.LN2 )
add((x, X) => exp(x)/sqrt(-X(-1)), 0, 1, 4.060156938557409951078179851 )
//add((x, a_x, b_x) =>  exp(-x)/sqrt(-b_x), 0, -1, -4.060156938557409951078179851 )

add(x => 1/(sqrt(1/100 + pow(x, 1/3))), 0, 1, 1.1901948695605892764)
//add(x => 1/sqrt(25-x*x) + sin(3*x)/2 + exp(-x*x), -5, 5, 4.914046504492584189227564493);
end_function_type("Singularity & Infinite derivative")

add(x => exp(x**4), -5, 5, 1.08801149989862368707e269)
end_function_type("High derivative")

add(x => (x-5)*(x+3)*(x+1)*(x+.3)*(x-.7), 0, 10, 94488.33333333333333)

end_function_type("Polynomials")

add(x => exp(-x), 0, Infinity, 1)
add(x => exp(-100000*x), 0, Infinity, 0.00001)
add(x => exp(x), -Infinity, 0, 1)
add(x => exp(-(x**2) ), -10, 10, sqrt(Math.PI))
end_function_type("Integration on real line")

end_test()