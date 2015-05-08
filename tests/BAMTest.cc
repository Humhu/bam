#include "bam/bam.h"
#include <cmath>
#include <ctime>
#include <cstdlib>

int main(int argc, char *argv[]) {

    using namespace bam;

    BAM32 b0, b1(0.0), b2(M_PI), b3(-M_PI), b4(2.0*M_PI), b5(-2.0*M_PI);

    std::cout << "b0 (default constructor) " << b0 << std::endl;
    std::cout << "b1 (0.0) " << b1 << std::endl;
    std::cout << "b2 (PI) " << b2 << std::endl;
    std::cout << "b3 (-PI) " << b3 << std::endl;
    std::cout << "b4 (2PI) " << b4 << std::endl;
    std::cout << "b5 (-2PI) " << b5 << std::endl;

    std::cout << "b0 (default constructor) " << b0 << std::endl;
    std::cout << "b1++ " << b1++ << std::endl;
    std::cout << "++b1 " << ++b1 << std::endl;
    std::cout << "2*b1 " << 2*b1 << std::endl;
    std::cout << "2.0*b1 " << 2.0*b1 << std::endl;
    std::cout << "-2*b1 " << -2*b1 << std::endl;
    std::cout << "-2.0*b1 " << -2.0*b1 << std::endl;
    std::cout << "b1/2 " << b1/2 << std::endl;
    std::cout << "b1/2.0 " << b1/2.0 << std::endl;
    std::cout << "b1/-2 " << b1/-2 << std::endl;
    std::cout << "b1/-2.0 " << b1/-2.0 << std::endl;

	std::cout << "abs(b3) " << bamAbs(b3) << std::endl;
	
    BAM32 b6(M_PI/3.0), b7(2*M_PI/3.0), b8(4*M_PI/3.0), b9(5*M_PI/3.0), b10(M_PI/2.0);
    std::cout << "b6 (pi/3) " << b6 << std::endl;
    std::cout << "b7 (2pi/3) " << b7 << std::endl;
    std::cout << "b8 (4pi/3) " << b8 << std::endl;
    std::cout << "b9 (5pi/3) " << b9 << std::endl;
    std::cout << "b10 (pi/2) " << b10 << std::endl;

    std::cout << "sinFast(b6) " << bamSinFast(b6) << " sin(b6.ToDouble) " << sin(b6.ToDouble()) << std::endl;
    std::cout << "sinFast(b7) " << bamSinFast(b7) << " sin(b7.ToDouble) " << sin(b7.ToDouble()) << std::endl;
    std::cout << "sinFast(b8) " << bamSinFast(b8) << " sin(b8.ToDouble) " << sin(b8.ToDouble()) << std::endl;
    std::cout << "sinFast(b9) " << bamSinFast(b9) << " sin(b9.ToDouble) " << sin(b9.ToDouble()) << std::endl;

    std::cout << "sinFast(b2) " << bamSinFast(b2) << " sin(b2.ToDouble) " << sin(b2.ToDouble()) << std::endl;
    std::cout << "sinFast(b4) " << bamSinFast(b4) << " sin(b4.ToDouble) " << sin(b4.ToDouble()) << std::endl;
    std::cout << "sinFast(b10) " << bamSinFast(b10) << " sin(b10.ToDouble) " << sin(b10.ToDouble()) << std::endl;

    std::cout << "bamAsinFast(0.3) " << bamAsinFast(0.3) << " bamAsin(0.3) " << bamAsin(0.3) << " asin(0.3) " << asin(0.3) << std::endl;
    std::cout << "bamAsinFast(1.0) " << bamAsinFast(1.0) << " bamAsin(1.0) " << bamAsin(1.0) << " asin(0.6) " << asin(1.0) << std::endl;

    std::cout << "bamAtan2(0.4, 3.0) " << bamAtan2(0.4, 3.0) << " atan2(0.4, 3.0) " << atan2(0.4, 3.0) << std::endl;
    std::cout << "bamAtan2(3.0, 0.4) " << bamAtan2(3.0, 0.4) << " atan2(3.0, 0.4) " << atan2(3.0, 0.4) << std::endl;
    
    time_t now;
    time(&now);
    std::cout << "now: " << now << std::endl;
    srand(now);

    unsigned int numRuns = 1E6;
    
    double rsum = 0;
    time_t benchStart = clock();
    for(unsigned int i = 1; i < numRuns; ++i) {
        double value = M_PI*rand()/((float)RAND_MAX);
        rsum += value;
    }
    time_t benchStop = clock();
    std::cout << "rsum = " << rsum << std::endl;
    std::cout << "creating(rand) took " << (benchStop - benchStart)/( (float) CLOCKS_PER_SEC ) << " seconds." << std::endl;
    
    double sum = 0;
    benchStart = clock();
    for(unsigned int i = 0; i < numRuns; ++i) {
        double value = M_PI*rand()/((float)RAND_MAX);
        sum += sin(value);
        
    }
    benchStop = clock();
    std::cout << "sum = " << sum << std::endl;
    std::cout << "sin(double) took " << (benchStop - benchStart)/( (float) CLOCKS_PER_SEC ) << " seconds." << std::endl;

    BAM32 bsum = 0;
    benchStart = clock();
    for(unsigned int i = 0; i < numRuns; ++i) {
        double value = M_PI*rand()/((float)RAND_MAX);
        BAM32 b(value);
        bsum += b;
    }
    benchStop = clock();
    std::cout << "bsum = " << bsum << std::endl;
    std::cout << "creating(BAM32) took " << (benchStop - benchStart)/( (float) CLOCKS_PER_SEC ) << " seconds." << std::endl;
    
    sum = 0;
    benchStart = clock();
    for(unsigned int i = 0; i < numRuns; ++i) {
        double value = M_PI*rand()/((float)RAND_MAX);
        BAM32 b(value);
        sum += bamSinFast(b);
    }
    benchStop = clock();
    std::cout << "sum = " << sum << std::endl;
    std::cout << "sinFast(BAM32) took " << (benchStop - benchStart)/( (float) CLOCKS_PER_SEC ) << " seconds." << std::endl;

    sum = 0;
    benchStart = clock();
    for(unsigned int i = 0; i < numRuns; ++i) {
        double value = M_PI*rand()/((float)RAND_MAX);
        BAM32 b(value);
        sum += bamSin(b);
    }
    benchStop = clock();
    std::cout << "sum = " << sum << std::endl;
    std::cout << "sin(BAM32) took " << (benchStop - benchStart)/( (float) CLOCKS_PER_SEC ) << " seconds." << std::endl;

    sum = 0;
    benchStart = clock();
    for(unsigned int i = 0; i < numRuns; ++i) {
        double value = rand()/((float)RAND_MAX);
        sum += asin(value);
    }
    benchStop = clock();
    std::cout << "sum = " << sum << std::endl;
    std::cout << "asin(double) took " << (benchStop - benchStart)/( (float) CLOCKS_PER_SEC ) << " seconds." << std::endl;
    
    bsum = 0;
    benchStart = clock();
    for(unsigned int i = 0; i < numRuns; ++i) {
        double value = rand()/((float)RAND_MAX);
        bsum += bamAsinFast(value);
    }
    benchStop = clock();
    std::cout << "sum = " << sum << std::endl;
    std::cout << "asinFast(BAM32) took " << (benchStop - benchStart)/( (float) CLOCKS_PER_SEC ) << " seconds." << std::endl;

    bsum = 0;
    benchStart = clock();
    for(unsigned int i = 0; i < numRuns; ++i) {
        double value = rand()/((float)RAND_MAX);
        bsum += bamAsin(value);
    }
    benchStop = clock();
    std::cout << "sum = " << sum << std::endl;
    std::cout << "asin(BAM32) took " << (benchStop - benchStart)/( (float) CLOCKS_PER_SEC ) << " seconds." << std::endl;

    sum = 0;
    benchStart = clock();
    for(unsigned int i = 0; i < numRuns; ++i) {
        double value0 = rand()/((float)RAND_MAX);
        double value1 = rand()/((float)RAND_MAX);
        sum += atan2(value0, value1);
    }
    benchStop = clock();
    std::cout << "sum = " << sum << std::endl;
    std::cout << "atan2(double) took " << (benchStop - benchStart)/( (float) CLOCKS_PER_SEC ) << " seconds." << std::endl;
    
    bsum = 0;
    benchStart = clock();
    for(unsigned int i = 0; i < numRuns; ++i) {
        double value0 = rand()/((float)RAND_MAX);
        double value1 = rand()/((float)RAND_MAX);
        bsum += bamAtan2(value0, value1);
    }
    benchStop = clock();
    std::cout << "sum = " << sum << std::endl;
    std::cout << "atan2(BAM32) took " << (benchStop - benchStart)/( (float) CLOCKS_PER_SEC ) << " seconds." << std::endl;

}
