#ifndef BASE_H
#define BASE_H

const float eps = 0.0001;

#define DEBUGGABLE

#ifdef DEBUGGABLE
#define DEBUG() cerr << __FILE__ << ':' << __LINE__ << ' ' << __func__ << "()" << endl;
#endif // DEBUGGABLE

#endif // BASE_H
