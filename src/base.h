#ifndef BASE_H
#define BASE_H

const float eps = 0.00001;

#define DEBUGGABLE

#ifdef DEBUGGABLE
#define DEBUG() cerr << __FILE__ << ':' << __LINE__ << ' ' << __func__ << "()" << endl;
#endif // DEBUGGABLE

#define PI     3.14159265359
#define PI_2   1.57079632679
#define PI_4   0.78539816339
#define PI_3_4 2.35619449019

#define COS_0  1
#define COS_15 0.965926
#define COS_30 0.866025
#define COS_45 0.707107
#define COS_60 0.5
#define COS_75 0.258819
#define COS_90 0
#define SIN_0 COS_90
#define SIN_15 COS_75
#define SIN_30 COS_60
#define SIN_45 COS_45
#define SIN_60 COS_30
#define SIN_75 COS_15
#define SIN_90 COS_0

#endif // BASE_H
