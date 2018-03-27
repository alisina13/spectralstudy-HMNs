#ifndef TYPES_HPP
#define TYPES_HPP
#include <map>
#include <string>



// This typedef makes it possible to switch between float and double accuracy
// please do not use "float" or "double" directly in your code, but use real instead
typedef double real;
typedef std::map<std::string,int>::const_iterator sii;
typedef std::map<std::string,real>::const_iterator sri;
typedef std::map<std::string,std::string>::const_iterator ssi;




#endif //TYPES_HPP
