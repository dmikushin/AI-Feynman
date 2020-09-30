#include <cstdio>
#include <stdint.h>
#include <vector>

typedef void (*loop_body_t)(int* i, int64_t iformula, int64_t* iformula_out,
	double* prefactor_out, double* minloss_out, double* rmsloss_out, char* ops);

typedef void (*report_body_t)(int64_t iformula,
	double prefactor, double minloss, double rmsloss, char* ops);

using namespace std;

extern "C" void multiloop_(int* n_, int* bases, loop_body_t loop_body, report_body_t report_body,
	double* minloss, int64_t* nformulas_)
{
	int n = *n_;
        int64_t count = 1;
        for (int k = 0; k < n; k++)
                count *= bases[k];

	vector<int> i(n), i1(n);

	int64_t nformulas = *nformulas_;
	int64_t iformula_local = (int64_t)-1;
	double prefactor_local, minloss_local = *minloss, rmsloss_local;
	char ops_local[60];
	for (int64_t iformula = 0; iformula < count; iformula++)
	{
		// Decode multi-index from a planar index.
		int64_t base = iformula, divisor = count;
		for (int k = 0; k < n; k++)
		{
			divisor /= bases[k];
			i1[k] = base / divisor;
			base -= i1[k] * divisor;
		}

		loop_body(reinterpret_cast<int*>(&i1[0]), nformulas + iformula, &iformula_local,
			&prefactor_local, &minloss_local, &rmsloss_local, ops_local);
	}

	if (*minloss > minloss_local)
	{
		// Report new best fit. 
		report_body(iformula_local, prefactor_local, minloss_local, rmsloss_local, ops_local);
		*minloss = minloss_local;
	}

	*nformulas_ += count;
}

