#ifndef UTILS_INITIALIZATION_H
#define UTILS_INITIALIZATION_H

using ArgvType = char*[];
void initialize(int& argc, ArgvType argv);

void finalize();

#endif