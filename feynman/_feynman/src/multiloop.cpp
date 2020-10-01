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
	vector<int64_t> divisors;
	double minloss;

	loop_body_t loop_body;

	LossDataCommon(int n_,  int64_t count_, int* bases, double minloss_, loop_body_t loop_body_) :
		n(n_), count(count_), divisors(n_), minloss(minloss_), loop_body(loop_body_)
	{
		int64_t divisor = count;
		for (int k = 0; k < n; k++)
		{
			divisor /= bases[k];
			divisors[k] = divisor;
		}
	}
};

unique_ptr<LossDataCommon> c;

struct LossData
{
	int64_t iformula;
	double prefactor;
	double minloss;
	double rmsloss;
	char ops[60];
	bool evaluated = false;

	LossData(int64_t iformula_) : minloss(c->minloss), iformula(iformula_) { }

	LossData& eval()
	{
		if (evaluated) return *this;

		vector<int> i(c->n);

		// Decode multi-index from a planar index.
		int64_t nominator = iformula;
		for (int k = 0; k < c->n; k++)
		{
			int64_t divisor = c->divisors[k];
			i[k] = nominator / divisor;
			nominator -= i[k] * divisor;
		}

		c->loop_body(reinterpret_cast<int*>(&i[0]), &prefactor, &minloss, &rmsloss, ops);

		evaluated = true;
		return *this;
	}

	LossData& operator++()
	{
		iformula++;
		evaluated = false;
		return *this;
	}

	LossData operator+(const LossData& other) const
	{
		if (this->minloss < other.minloss)
			return *this;

		return other;
	}
};

struct MultiloopIterator;

namespace thrust {

template<>
struct iterator_system<MultiloopIterator>
{
	using type = thrust::device_system_tag;
};

template<>
struct iterator_traits<MultiloopIterator>
{
	typedef std::ptrdiff_t difference_type;
	typedef MultiloopIterator value_type;
	typedef MultiloopIterator* pointer;
	typedef MultiloopIterator& reference;
	typedef std::random_access_iterator_tag iterator_category;
};

}

struct MultiloopIterator
{
	LossData lossData;

	MultiloopIterator(int64_t iformula) : lossData(iformula) { }

	LossData& operator*()
	{
		return lossData.eval();
	}

	MultiloopIterator& operator++()
	{
		++lossData;
		return *this;
	}

	friend typename thrust::iterator_traits<MultiloopIterator>::difference_type operator-(MultiloopIterator a, MultiloopIterator b)
	{
		return (*a).iformula - (*b).iformula;
	}

	friend typename thrust::iterator_traits<MultiloopIterator>::difference_type operator+(MultiloopIterator a, int64_t b)
	{
		return (*a).iformula + b;
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
		MultiloopIterator(0), MultiloopIterator(count),
		LossData(0), thrust::plus<LossData>());

	if (*minloss > result.minloss)
	{
		// Report new best fit. 
		report_body(result.iformula + nformulas, result.prefactor, result.minloss, result.rmsloss, result.ops);
		*minloss = result.minloss;
	}

	*nformulas_ += count;
}

