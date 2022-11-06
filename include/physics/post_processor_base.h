#ifndef PHYSICS_POST_PROCESSOR_BASE_H
#define PHYSICS_POST_PROCESSOR_BASE_H

#include <string>
#include <vector>
#include "discretization/disc_vector.h"

namespace physics {

class PostProcessorBase
{
  public:
    virtual ~PostProcessorBase() {}

    // returns number of values this postprocessor returns
    virtual int numValues() const = 0;

    virtual std::vector<std::string> getNames() const = 0;

    virtual std::vector<double> getValues(DiscVectorPtr u, double t) = 0;
};

using PostProcessorPtr = std::shared_ptr<PostProcessorBase>;

}

#endif