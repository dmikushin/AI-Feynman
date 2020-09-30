#include <cstdio>
#include <memory>
#include <stdint.h>
#include <thrust/reduce.h>
#include <vector>

typedef void (*loop_body_t)(int* i,
	double* prefactor_out, double* minloss_out, double* rmsloss_out, char* ops);

typedef void (*report_body_t)(int64_t iformula,
	double prefactor, double minloss, double rmsloss, char* ops);

using namespace std;

struct LossDataCommon
{
	int n;
	int64_t count;
	int* bases;
	double minloss;

	loop_body_t loop_body;

	LossDataCommon(int n_,  int64_t count_, int* bases_, double minloss_, loop_body_t loop_body_) :
		n(n_), count(count_), bases(bases_), minloss(minloss_), loop_body(loop_body_) { }
};

unique_ptr<LossDataCommon> c;

struct LossData
{
	int64_t iformula;
	double prefactor;
	double minloss;
	double rmsloss;
	char ops[60];

	LossData(int64_t iformula_) : minloss(c->minloss), iformula(iformula_)
	{
		vector<int> i(c->n);

		// Decode multi-index from a planar index.
                int64_t base = iformula, divisor = c->count;
                for (int k = 0; k < c->n; k++)
                {
                        divisor /= c->bases[k];
                        i[k] = base / divisor;
                        base -= i[k] * divisor;
                }

                c->loop_body(reinterpret_cast<int*>(&i[0]), &prefactor, &minloss, &rmsloss, ops);
	}

	operator int64_t() { return 0; }

	LossData operator++() { return LossData(iformula + 1); }

	LossData operator+=(int n) { return LossData(iformula + n); }
};

struct LossDataReduce
{
	LossData operator()(const LossData& a, const LossData& b) const
	{
		if (a.minloss < b.minloss)
			return a;

		return b;
	}
};

extern "C" void multiloop_(int* n_, int* bases, loop_body_t loop_body, report_body_t report_body,
	double* minloss, int64_t* nformulas_)
{
	int n = *n_;
        int64_t count = 1;
        for (int k = 0; k < n; k++)
                count *= bases[k];

	int64_t nformulas = *nformulas_;
	c.reset(new LossDataCommon(n, count, bases, *minloss, loop_body));
	LossData result = thrust::reduce(
		thrust::make_counting_iterator(LossData(0)),
	       	thrust::make_counting_iterator(LossData(count)),
		LossData(0), LossDataReduce());

	if (*minloss > result.minloss)
	{
		// Report new best fit. 
		report_body(result.iformula + nformulas, result.prefactor, result.minloss, result.rmsloss, result.ops);
		*minloss = result.minloss;
	}

	*nformulas_ += count;
}

