// Handles <n> nested loops with loop variables i(1),...i(n).
// Example: With n=3, bases=2, repeated calls starting with i=(000) will return
// 001, 010, 011, 100, 101, 110, 111, 000 (and done=.true. the last time).
// All it's doing is counting in mixed radix specified by the array <bases>.
extern "C" void multiloop_(int* n, int* bases, int* i, int* done)
{
	*done = 0;
	for (int k = 0, ke = *n; k < ke; k++)
	{
		i[k] = i[k] + 1;
		if (i[k] < bases[k]) return;
		i[k] = 0;
	}
	*done = 1;
}

