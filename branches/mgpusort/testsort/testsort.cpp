#include "mgpusort.hpp"
#include "benchmark.h"

// Need to sort fewer elements as the number of values increases, to fit
// everything in video memory.
#ifdef _DEBUG

const int ElementCounts[7] = {
	35000000,
	1<< 19,
	1<< 19,
	1<< 19,
	1<< 19,
	1<< 19,
	1<< 19
};
const int NumIterations = 1;
const int NumTests = 1;

#else


const int ElementCounts[7] = {
	40000000,
	27000000,
	16000000,
	12000000,
	10000000,
	80000000,
	70000000
};
const int NumIterations = 6;
const int NumTests = 1;
/*
const int ElementCounts[7] = {
	500000,
	500000,
	500000,
	500000,
	500000,
	500000,
	500000
};
const int NumIterations = 300;
const int NumTests = 5;
*/
#endif

bool TestSorted(CuDeviceMem** a, CuDeviceMem** b, int numElements,
	int valueCount) {
	std::vector<uint> hostA(numElements), hostB(numElements);

	for(int i(0); i < valueCount; ++i) {
		a[i]->ToHost(&hostA[0], sizeof(uint) * numElements);
		b[i]->ToHost(&hostB[0], sizeof(uint) * numElements);

		if(a != b) return false;
	}
	return true;
}


struct Throughput {
	double elementsPerSec;
	double bytesPerSec;
	double normElementsPerSec;
	double normBytesPerSec;

	void Max(const Throughput& rhs) {
		elementsPerSec = std::max(elementsPerSec, rhs.elementsPerSec);
		bytesPerSec = std::max(bytesPerSec, rhs.bytesPerSec);
		normElementsPerSec = std::max(normElementsPerSec, 
			rhs.normElementsPerSec);
		normBytesPerSec = std::max(normBytesPerSec, rhs.normBytesPerSec);
	}
};


Throughput CalcThroughput(int numBits, int numElements, int valueCount, 
	int iterations, double elapsed) {

	Throughput throughput;
	throughput.elementsPerSec = numElements * iterations / elapsed;
	throughput.bytesPerSec = sizeof(uint) * (1 + abs(valueCount)) * 
		throughput.elementsPerSec;
	throughput.normElementsPerSec = 
		(numBits / 32.0) * throughput.elementsPerSec;
	throughput.normBytesPerSec = (numBits / 32.0) * throughput.bytesPerSec;

	return throughput;
}


////////////////////////////////////////////////////////////////////////////////
// Terms for setting up a benchmark run for MGPU

struct BenchmarkTerms {
	CuContext* context;
	sortEngine_t engine;
	int count;
	int numBits;
	int bitPass;
	int numThreads;
	int valueCount;
	int numIterations;
	int numTests;
};


static uint rand_seed = 1111111;
static uint rand_carry = 2222222;

static uint
get_rand(uint &rand_seed, uint &rand_carry) {
   
	/* A multiply-with-carry generator by George Marsaglia.
	   The period is about 2^63. */

	#define RAND_MULT 2131995753

	uint64 temp;

	temp = (uint64)(rand_seed) * 
		       (uint64)RAND_MULT + 
		       (uint64)(rand_carry);
	rand_seed = (uint)temp;
	rand_carry = (uint)(temp >> 32);
	return (uint)temp;
}


bool Benchmark(BenchmarkTerms& terms, Throughput& mgpu) {

	int capacity = RoundUp(terms.count, 2048);

	std::vector<uint> keysHost(terms.count), valuesHost[6];
	DeviceMemPtr keysDevice, valuesDevice[6], 
		sortedKeysDevice, sortedValuesDevice[6];

	MgpuTerms mgpuTerms = { 0 };

	mgpuTerms.numBits = terms.numBits;
	mgpuTerms.iterations = terms.numIterations;
	mgpuTerms.reset = true;
	mgpuTerms.valueCount = terms.valueCount;
	mgpuTerms.count = terms.count;
	mgpuTerms.numThreads = terms.numThreads;
	mgpuTerms.bitPass = terms.bitPass;

	// Generate random numbers for the keys and assign to the terms structs.
	std::vector<uint> r(terms.count);
	uint mask = 0xffffffff >> (32 - terms.numBits);
	for(int i(0); i < terms.count; ++i)
		keysHost[i] = get_rand(rand_seed, rand_carry) & mask;

	// Allocate space for the random keys and sorted keys.
	terms.context->MemAlloc<uint>(capacity, &keysDevice);
	keysDevice->FromHost(keysHost);
	terms.context->MemAlloc<uint>(capacity, &sortedKeysDevice);

	mgpuTerms.randomKeys = keysDevice.get();
	mgpuTerms.sortedKeys = sortedKeysDevice.get();


	if(terms.valueCount) {
		valuesHost[0].resize(terms.count);
		for(int i(0); i < terms.count; ++i)
			valuesHost[0][i] = i;

		for(int i(0); i < abs(terms.valueCount); ++i) {
			// Allocate space for the random values and sorted values.
			terms.context->MemAlloc<uint>(capacity, &valuesDevice[i]);
			terms.context->MemAlloc<uint>(capacity, &sortedValuesDevice[i]);
			valuesDevice[i]->FromHost(valuesHost[0]);

			mgpuTerms.randomVals[i] = valuesDevice[i].get();
			mgpuTerms.sortedVals[i] = sortedValuesDevice[i].get();
		}
	}


	double elapsed;
	Throughput throughput;
	for(int test(0); test < NumTests; ++test) {
		// MGPU benchmark
		sortStatus_t status = MgpuBenchmark(mgpuTerms, terms.engine, &elapsed);
		if(SORT_STATUS_SUCCESS != status) {
			printf("Error in MGPU sort on numBits = %d: %s\n", terms.numBits,
				sortStatusString(status));
			return false;
		}
		throughput = CalcThroughput(terms.numBits, terms.count, 
			terms.valueCount, terms.numIterations, elapsed);
		mgpu.Max(throughput);
	}
#if 1
	// Read the MGPU results into host memory,
	// verify the keys are increasing
	mgpuTerms.sortedKeys->ToHost(&keysHost[0], terms.count);
	for(int i(1); i < terms.count; ++i) {
		if (keysHost[i] < keysHost[i-1]) {
			printf("sort failed at position %d of %d\n",
					i, terms.count);
			break;
		}

	}
#endif

	return true;
}



////////////////////////////////////////////////////////////////////////////////
// BenchmarkBitPass benchmarks the individual bit pass speeds. The results are
// returned in a simple format that can be parsed by tablegen to create optimal
// multi-pass algorithms for sorting keys of any size.

bool BenchmarkBitPass(CuContext* context, sortEngine_t engine, 
	const int* testSizes, int numIterations, int numTests,
	const char* tableSuffix) {

	for(int valueCount(-1); valueCount <= 6; ++valueCount) {
		for(int numThreads(128); numThreads <= 256; numThreads *= 2) {
			
			// Formulate a table name like sort_128_8_key_simple_table
			printf("sort_%d_8_", numThreads);
			switch(valueCount) {
				case -1: printf("index_"); break;
				case 0: printf("key_"); break;
				case 1: printf("single_"); break;
				default: printf("multi_%d_", valueCount); break;
			}
			// Only benchmark simple storage for now
			printf("simple_");

			printf("%s\n", tableSuffix);

			for(int bitPass(1); bitPass <= 6; ++bitPass) {
				BenchmarkTerms terms;
				terms.context = context;
				terms.engine = engine;
				terms.count = testSizes[abs(valueCount)];
				terms.numBits = (32 % bitPass) ? (32 - (32 % bitPass)) : 32;
				terms.bitPass = bitPass;
				terms.numThreads = numThreads;
				terms.valueCount = valueCount;
				terms.numIterations = numIterations;
				terms.numTests = numTests;

				Throughput mgpu = { 0 };
				bool success = Benchmark(terms, mgpu);
				if(!success) return false;

				printf("%7.3lf\n", mgpu.normElementsPerSec / 1.0e6);
			}
		}
	}
	return true;
}

void BenchmarkBitPassLarge(CuContext* context, sortEngine_t engine) {
	const int LargePass[7] = {
		35000000,
		27000000,
		16000000,
		12000000,
		10000000,
		8000000,
		7000000
	};
	BenchmarkBitPass(context, engine, LargePass, 12, 1, "large");
}

void BenchmarkBitPassSmall(CuContext* context, sortEngine_t engine) {
	const int SmallPass[7] = {
		500000,
		500000,
		500000,
		500000,
		500000,
		500000,
		500000
	};
	BenchmarkBitPass(context, engine, SmallPass, 100, 5, "small");
}




int main(int argc, char** argv) {

	cuInit(0);
	
	DevicePtr device;
	CreateCuDevice(0, &device);

	ContextPtr context;
	CreateCuContext(device, 0, &context);

	sortEngine_t engine;
	sortStatus_t status = sortCreateEngine(argv[1], &engine);
	if(SORT_STATUS_SUCCESS != status) {
		printf("Error creating MGPU sort engine: %s\n",
			sortStatusString(status));
		return 0;
	}
	BenchmarkBitPassSmall(context, engine);
	BenchmarkBitPassLarge(context, engine);

	sortReleaseEngine(engine);
}
