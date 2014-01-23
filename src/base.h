#ifndef BASE_H
#define BASE_H

const float eps = 0.0001;

#define LOGGABLE

#ifdef LOGGABLE
#define __LOG() cerr << __FILE__ << ':' << __LINE__ << ' ' << __func__ << "()" << endl;
#endif // LOGGABLE

#endif // BASE_H
