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
		// 678:
		// base <- 678
		// divisor <- 100
		// i1[0] <- 678 / 100 = 6
		// base <- 78
		// divisor <- 10
		// i1[1] <- 78 / 10 = 7

		// TODO Decode multi-index from a planar index.
		int64_t base = iformula, divisor = count;
		for (int k = 0; k < n; k++)
		{
			divisor /= bases[k];
			i1[k] = base / divisor;
			base -= i1[k] * divisor;
		}

		if (loop_body(reinterpret_cast<int*>(&i1[0]), nformulas + iformula, minloss) == 1)
			break;

		// Handles <n> nested loops with loop variables i(1),...i(n).
		// Example: With n=3, bases=2, repeated calls starting with i=(000) will return
		// 001, 010, 011, 100, 101, 110, 111, 000 (and done=.true. the last time).
		// All it's doing is counting in mixed radix specified by the array <bases>.
		for (int k = 0; k < n; k++)
		{
			i[k] = i[k] + 1;
			if (i[k] < bases[k]) break;
			i[k] = 0;
		}
#if 0
		printf("%ld -> i : ", iformula);
		for (int k = 0; k < n; k++)
			printf("%d ", i[k]);
		printf("\n");

		printf("%ld -> i1: ", iformula);
		for (int k = 0; k < n; k++)
			printf("%d ", i1[k]);
		printf("\n");
#endif
	}

	*nformulas_ += count;
}

