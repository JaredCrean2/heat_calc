#ifndef PHYSICS_POST_PROCESSOR_SCHEDULER_H
#define PHYSICS_POST_PROCESSOR_SCHEDULER_H

#include <memory>

namespace physics {

class PostProcessorScheduler
{
  public:
    virtual ~PostProcessorScheduler() {}

    virtual bool shouldOutput(int timestep, double t) = 0;    
};

using PostProcessorSchedulerPtr = std::shared_ptr<PostProcessorScheduler>;


class PostProcessorScheduleFixedInterval : public PostProcessorScheduler
{
  public:
    PostProcessorScheduleFixedInterval(int interval) :
      m_interval(interval)
    {}

    bool shouldOutput(int timestep, double t) override { return (timestep % m_interval) == 0; }
  
  private:
    int m_interval;
};

} // namespace

#endif