#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

extern "C" double f3_(int* nops_, int* arities, char* ops, double* x)
{
	int nops = *nops_, j = -1;
	std::vector<double> stack(2 * nops);
	for (int i = 0; i < nops; i++)
	{
		double y;
		int arity = arities[i];
		char op = ops[i];
		if (arity == 0) // This is a nonary function
		{
			switch (op)
			{
			case '0': y = 0.0; break;
			case '1': y = 1.0; break;
			case 'P': y = M_PI; break;
			default : y = x[(int)op - 97]; break;
			}
		}
		else if (arity == 1) // This is a unary function
		{
			switch (op)
			{
			case '>': y = stack[j] + 1; break;
			case '<': y = stack[j] - 1; break;
			case '~': y = -stack[j]; break;
			case '\\': y = 1.0 / stack[j]; break;
			case 'L': y = log(stack[j]); break;
			case 'E': y = exp(stack[j]); break;
			case 'S': y = sin(stack[j]); break;
			case 'C': y = cos(stack[j]); break;
			case 'A': y = abs(stack[j]); break;
			case 'N': y = asin(stack[j]); break;
			case 'T': y = atan(stack[j]); break;
			default : y = sqrt(stack[j]); break;
			}
		}
		else // This is a binary function
		{
			switch (op)
			{
			case '+': y = stack[j - 1] + stack[j]; break;
			case '-': y = stack[j - 1] - stack[j]; break;
			case '*': y = stack[j - 1] * stack[j]; break;
			default : y = stack[j - 1] / stack[j]; break;
			}
		}
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

