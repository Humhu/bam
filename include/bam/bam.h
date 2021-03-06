/*
 * Binary Angle Measurement System (BAMS) Implementation
 *
 * by Humphrey Hu
 *
 * v.alpha
 *
 * Revisions:
 * Humphrey Hu 2011-10-23 Initial implementation (at UC Berkeley)
 * Humphrey Hu 2013-01-28 Ported to C++
 * Humphrey Hu 2014-02-12 Updated in C++
 */

/*
 * Two's complement (binary angles) are convenient because they
 * take apply integer math wrap-around to angle representation.
 * Also, being fixed-point, addition and subtraction of angles is
 * extremely fast and condusive to fast lookup tables.
 *
 * 16-bit binary angles map a range of [-pi, pi) to
 * [0x8000, 0x7FFF]. Note that this range does not include
 * positive pi, such that there is only one value of pi.
 * Though they are interpreted as signed integers,
 * they are equivalent to signed angles.
 *
 * 32-bit binary angles map a range of [-pi, pi) to
 * [0x80000000, 0x7FFFFFFF]. Most 32-bit functions have not been
 * implemented yet.
 */

#ifndef _BAM_H_
#define _BAM_H_

#include <stdint.h> // Needed for uintX_t typedefs
#include <iostream>

// TODO: bams16 is a holdover. Switch to BAM32 and bams64
// TODO: Add 16 <-> 32 <-> 64 conversion?
// TODO: BAMS superclass for polymorphism?
// TODO: Implement ==, >=, etc.
// TODO: Re-implement trig function lookups

namespace bam {

inline double WrapToPi(double angle);

class BAM32 {
public:

    typedef int32_t BinaryType; // The internal representation

    static const BinaryType PI = 0x80000000;
    static const BinaryType PI_2 = 0x40000000;
    static const BinaryType PI_4 = 0x80000000;
    
    BAM32();
    explicit BAM32(BinaryType b);
    BAM32(const BAM32& other);
    BAM32(double d);

    // Conversion functions
    double ToDouble() const; // Implicit conversion might be trouble...
    //BAMS16 BAMS16() const; // TODO: Implement!
    //BAMS64 BAMS64() const;
    inline BinaryType ToBinary() const { return val; }
    
    // Unary operators
    BAM32 operator+() const;
    BAM32 operator-() const;
    BAM32& operator++();
    BAM32& operator--();
    BAM32 operator++(int);
    BAM32 operator--(int);
    
    BAM32& operator=(const BAM32& rhs);
    BAM32& operator=(double rhs);
    BAM32& operator+=(const BAM32& rhs);
    BAM32& operator-=(const BAM32& rhs);
    BAM32& operator*=(int n);
    BAM32& operator*=(double d);
    BAM32& operator/=(int n);
    BAM32& operator/=(double d); // Watch out for int roundoff
    
    friend BAM32 operator+(BAM32 lhs, const BAM32& rhs);
    friend BAM32 operator-(BAM32 lhs, const BAM32& rhs);
    friend BAM32 operator*(BAM32 b, int n);
    friend BAM32 operator*(BAM32 b, double d);
    friend BAM32 operator*(int n, BAM32 b);
    friend BAM32 operator*(double d, BAM32 b);
    friend BAM32 operator/(BAM32 b, double d);
    friend BAM32 operator/(BAM32 b, int n);
	friend BAM32 abs(BAM32 b);

    BinaryType val;
    
private:
    void SetFromDouble(double d);
};

// Global operators
BAM32 operator+(BAM32 lhs, const BAM32& rhs);
BAM32 operator-(BAM32 lhs, const BAM32& rhs);
BAM32 operator*(BAM32 b, int n);
BAM32 operator*(BAM32 b, double d);
BAM32 operator*(int n, BAM32 b);
BAM32 operator*(double d, BAM32 b);
BAM32 operator/(BAM32 b, double d);
BAM32 operator/(BAM32 b, int n);
BAM32 bamAbs(BAM32 b);

std::ostream& operator<<(std::ostream& os, const BAM32& b);

// Trig functions
double bamSin(const BAM32& bam);
double bamCos(const BAM32& bam);
double bamTan(const BAM32& bam);

BAM32 bamAsin(double a);
BAM32 bamAcos(double a);
BAM32 bamAtan(double a);
BAM32 bamAtan2(double y, double x);

// Fast non-interpolating versions. Used in the normal interpolating versions.
double bamSinFast(const BAM32& bam);
double bamCosFast(const BAM32& bam);
double bamTanFast(const BAM32& bam);

BAM32 bamAsinFast(double a);
BAM32 bamAcosFast(double a);
// BAM32 atan(double a);
BAM32 bamAtan2(double y, double x);

} // end namespace BAM

#endif
