#include <cstdio>
#include <stdint.h>
#include <vector>

typedef int (*loop_body_t)(int* i, int64_t iformula, double* minloss);

using namespace std;

extern "C" void multiloop_(int* n_, int* bases, loop_body_t loop_body, double* minloss, int64_t* nformulas_)
{
	int n = *n_;
        int64_t count = 1;
        for (int k = 0; k < n; k++)
                count *= bases[k];

	vector<int> i(n), i1(n);

	int64_t nformulas = *nformulas_;
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

		if (loop_body(reinterpret_cast<int*>(&i1[0]), nformulas + iformula, minloss) == 1)
			break;
	}

	*nformulas_ += count;
}

