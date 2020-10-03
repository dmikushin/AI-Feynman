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

template<typename T>
struct LossDataCommon
{
	int n;
	T count;
	vector<T> divisors;
	double minloss;

	loop_body_t loop_body;

	LossDataCommon(int n_, T count_, int* bases, double minloss_, loop_body_t loop_body_) :
		n(n_), count(count_), divisors(n_), minloss(minloss_), loop_body(loop_body_)
	{
		T divisor = count;
		for (int k = 0; k < n; k++)
		{
			divisor /= bases[k];
			divisors[k] = divisor;
		}
	}
};

unique_ptr<LossDataCommon<int32_t> > c32;
unique_ptr<LossDataCommon<int64_t> > c64;

template<typename T>
struct LossData
{
	T iformula;
	double prefactor;
	double minloss;
	double rmsloss;
	char ops[60];
	bool evaluated = false;

	LossData(T iformula_) : iformula(iformula_) { }

	LossData<T>& eval()
	{
		if (evaluated) return *this;

		LossDataCommon<T>* c = nullptr;
		if (sizeof(T) == sizeof(int32_t))
			c = reinterpret_cast<LossDataCommon<T>*>(c32.get());
		else
			c = reinterpret_cast<LossDataCommon<T>*>(c64.get());

		minloss = c->minloss;

		vector<int> i(c->n);

		// Decode multi-index from a planar index.
		T nominator = iformula;
		for (int k = 0; k < c->n; k++)
		{
			T divisor = c->divisors[k];
			i[k] = nominator / divisor;
			nominator -= i[k] * divisor;
		}

		c->loop_body(reinterpret_cast<int*>(&i[0]), &prefactor, &minloss, &rmsloss, ops);

		evaluated = true;
		return *this;
	}

	LossData<T>& operator++()
	{
		iformula++;
		evaluated = false;
		return *this;
	}

	const LossData<T>& operator+(const LossData<T>& other) const
	{
		if (this->minloss < other.minloss)
			return *this;

		return other;
	}
};

template<typename T>
struct MultiloopIterator;

namespace thrust {

template<typename T>
struct iterator_system<MultiloopIterator<T> >
{
	using type = thrust::device_system_tag;
};

template<typename T>
struct iterator_traits<MultiloopIterator<T> >
{
	typedef T difference_type;
	typedef MultiloopIterator<T> value_type;
	typedef MultiloopIterator<T>* pointer;
	typedef MultiloopIterator<T>& reference;
	typedef std::random_access_iterator_tag iterator_category;
};

}

template<typename T>
struct MultiloopIterator
{
	using difference_type = typename thrust::iterator_traits<MultiloopIterator<T> >::difference_type;

	LossData<T> lossData;

	MultiloopIterator(T iformula) : lossData(iformula) { }

	LossData<T>& operator*()
	{
		return lossData.eval();
	}

	MultiloopIterator<T>& operator++()
	{
		++lossData;
		return *this;
	}

	friend difference_type operator-(const MultiloopIterator<T>& a, const MultiloopIterator<T>& b)
	{
		return a.lossData.iformula - b.lossData.iformula;
	}

	friend difference_type operator+(const MultiloopIterator<T>& a, T b)
	{
		return a.lossData.iformula + b;
	}
};

template<typename T>
LossData<T> multiloop(int n, T count, int* bases, loop_body_t loop_body, report_body_t report_body,
	double* minloss, int64_t* nformulas_)
{
        LossData<T> result = thrust::reduce(
                MultiloopIterator<T>(0), MultiloopIterator<T>(count),
                LossData<T>(0).eval(), thrust::plus<LossData<T> >());

        if (*minloss > result.minloss)
        {
                // Report new best fit. 
                report_body(result.iformula + *nformulas_, result.prefactor, result.minloss, result.rmsloss, result.ops);
                *minloss = result.minloss;
        }

	*nformulas_ += count;
}

extern "C" void multiloop_(int* n_, int* bases, loop_body_t loop_body, report_body_t report_body,
	double* minloss, int64_t* nformulas_)
{
	int n = *n_;
	int64_t count = 1;
	for (int k = 0; k < n; k++)
		count *= bases[k];

	if (count < INT_MAX)
	{
		c32.reset(new LossDataCommon<int32_t>(n, count, bases, *minloss, loop_body));
		multiloop<int32_t>(*n_, (int32_t)count, bases, loop_body, report_body, minloss, nformulas_);
	}
	else
	{
		c64.reset(new LossDataCommon<int64_t>(n, count, bases, *minloss, loop_body));
		multiloop<int64_t>(*n_, count, bases, loop_body, report_body, minloss, nformulas_);
	}
}

