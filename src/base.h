#ifndef BASE_H
#define BASE_H

#include <utility>
#include <iterator>
#include <cmath>

const float EPS           = 1e-5;
const float ONE_EPS       = 1 - EPS;
const float MINUS_ONE_EPS = -1 + EPS;
const float EPS_SQUARE    = EPS * EPS;
const float EPS_CUBE      = EPS * EPS * EPS;

const float PI     = 3.14159265359;
const float PI_2   = 1.57079632679;
const float PI_4   = 0.78539816339;
const float PI_3_4 = 2.35619449019;

const float COS_0  = 1;
const float COS_15 = 0.965926;
const float COS_30 = 0.866025;
const float COS_45 = 0.707107;
const float COS_60 = 0.5;
const float COS_75 = 0.258819;
const float COS_90 = 0;

/**
 * Circular list function begin
 */
template<typename T>
inline static typename T::iterator next(T& collection, typename T::iterator it) {
    typename T::iterator r(it);
    ++r;
    if (r == collection.end())
        r = collection.begin();
    return r;
}

template<typename T>
inline static typename T::const_iterator next(const T& collection, typename T::const_iterator it) {
    typename T::const_iterator r(it);
    ++r;
    if (r == collection.end())
        r = collection.begin();
    return r;
}

template<typename T>
inline static typename T::iterator prev(T& collection, typename T::iterator it) {
    typename T::iterator r(it);
    if (r == collection.begin())
        r = collection.end();
    --r;
    return r;
}

template<typename T>
inline static typename T::const_iterator prev(const T& collection, typename T::const_iterator it) {
    typename T::const_iterator r(it);
    if (r == collection.begin())
        r = collection.end();
    --r;
    return r;
}

/**
 * Circular list end
 */


/**
 * Tell if Set {a, b} intersects with Set {c, d}
 */
template <typename T>
inline static bool is_set_intersected(T a, T b, T c, T d) {
    return a == c || b == d || a == d || b == c;
}

template<typename S, typename T>
inline static int associative_compare(const std::pair<S, T>& left, const std::pair<S, T>& right) {
    return left.first < right.first;
}

template<typename T, typename V = typename std::iterator_traits<typename T::const_iterator>::value_type>
inline static V average(const T& collection) {
    V val = 0;
    for (const V& v: collection)
        val += v;
    return val / collection.size();
}

template<typename T, typename V = typename std::iterator_traits<typename T::const_iterator>::value_type>
inline static V standard_deviation(const T& collection) {
    V val = 0;
    V avg = average(collection);
    for (const V& v: collection)
        val += (v - avg) * (v - avg);
    return sqrt(val / collection.size());
}

template<typename T, typename V = typename std::iterator_traits<typename T::const_iterator>::value_type>
inline static V filtered_average(const T& collection) {
    V val = 0;
    V avg = average(collection);
    V sdv = standard_deviation(collection);
    V upper_limit = avg + sdv;
    V lower_limit = avg - sdv;
    for (const V& v: collection)
        if (lower_limit < v && v < upper_limit)
            val += v;
    return val / collection.size();
}

template<typename T, typename V = typename std::iterator_traits<typename T::const_iterator>::value_type>
inline static V filtered_standard_deviation(const T& collection) {
    V val = 0;
    V avg = average(collection);
    V sdv = standard_deviation(collection);
    V upper_limit = avg + sdv;
    V lower_limit = avg - sdv;
    for (const V& v: collection)
        if (lower_limit < v && v < upper_limit)
            val += (v - avg) * (v - avg);
    return sqrt(val / collection.size());
}

#define DEBUGGABLE
#ifdef DEBUGGABLE
#define DEBUG() cerr << __FILE__ << ':' << __LINE__ << ' ' << __func__ << "()" << endl;
#endif // DEBUGGABLE

#endif // BASE_H
