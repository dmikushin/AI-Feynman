#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

class OpFunctions
{
public :
	typedef double (*opfunc_t)(char op, double* stack, double* x);

private :
	opfunc_t opfuncs[3][CHAR_MAX];

public :

	inline const opfunc_t& operator()(char arity, char op) const
	{
		return opfuncs[arity][op];
	}

	OpFunctions()
	{
		// Nonary function.
		opfuncs[0][0] = [](char op, double* stack, double* x) -> double { return x[(int)op - 97]; };
		for (char i = 1; i < CHAR_MAX; i++)
			opfuncs[0][i] = opfuncs[0][0];
		opfuncs[0]['0'] = [](char op, double* stack, double* x) -> double { return 0.0; };
		opfuncs[0]['1'] = [](char op, double* stack, double* x) -> double { return 1.0; };
                opfuncs[0]['P'] = [](char op, double* stack, double* x) -> double { return M_PI; };
	
		// Unary function.
		opfuncs[1][0] = [](char op, double* stack, double* x) -> double
		{ 
			fprintf(stderr, "Invalid op = '%c' for arity 1\n", op);
			exit(-1);
		};
		for (char i = 1; i < CHAR_MAX; i++)
			opfuncs[1][i] = opfuncs[1][0];
		opfuncs[1]['>'] = [](char op, double* stack, double* x) -> double { return stack[0] + 1; };
                opfuncs[1]['<'] = [](char op, double* stack, double* x) -> double { return stack[0] - 1; };
                opfuncs[1]['~'] = [](char op, double* stack, double* x) -> double { return -stack[0]; };
                opfuncs[1]['\\'] = [](char op, double* stack, double* x) -> double { return 1.0 / stack[0]; };
                opfuncs[1]['L'] = [](char op, double* stack, double* x) -> double { return log(stack[0]); };
                opfuncs[1]['E'] = [](char op, double* stack, double* x) -> double { return exp(stack[0]); };
                opfuncs[1]['S'] = [](char op, double* stack, double* x) -> double { return sin(stack[0]); };
                opfuncs[1]['C'] = [](char op, double* stack, double* x) -> double { return cos(stack[0]); };
                opfuncs[1]['A'] = [](char op, double* stack, double* x) -> double { return abs(stack[0]); };
                opfuncs[1]['N'] = [](char op, double* stack, double* x) -> double { return asin(stack[0]); };
                opfuncs[1]['T'] = [](char op, double* stack, double* x) -> double { return atan(stack[0]); };
                opfuncs[1]['R'] = [](char op, double* stack, double* x) -> double { return sqrt(stack[0]); };

		// Binary function.
		opfuncs[2][0] = [](char op, double* stack, double* x) -> double
                {
                        fprintf(stderr, "Invalid op = '%c' for arity 2\n", op);
                        exit(-1);
                };
                for (char i = 1; i < CHAR_MAX; i++)
                        opfuncs[2][i] = opfuncs[2][0];
                opfuncs[2]['+'] = [](char op, double* stack, double* x) -> double { return stack[-1] + stack[0]; };
                opfuncs[2]['-'] = [](char op, double* stack, double* x) -> double { return stack[-1] - stack[0]; };
                opfuncs[2]['*'] = [](char op, double* stack, double* x) -> double { return stack[-1] * stack[0]; };
                opfuncs[2]['/'] = [](char op, double* stack, double* x) -> double { return stack[-1] / stack[0]; };
	}
};

static OpFunctions opfuncs;

extern "C" double f3_(int* nops_, int* arities, char* ops, double* x)
{
	int nops = *nops_, j = -1;
	std::vector<double> stack(nops);
	for (int i = 0; i < nops; i++)
	{
		int arity = arities[i];
		char op = ops[i];
		double y = opfuncs(arity, op)(op, &stack[j], x);
		j += 1 - arity; // Number of operands on the stack
		stack[j] = y;
	}

	if (j != 0)
	{
		fprintf(stderr, "f3: unbalanced stack!\n");
		exit(-1);
	}

	return stack[0];
}

