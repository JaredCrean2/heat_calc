#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

#include <string>
#include <stdexcept>

void assertAlways(const bool condition, const std::string& msg);

#endif