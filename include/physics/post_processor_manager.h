#ifndef PHYSICS_POST_PROCESSOR_MANAGER_H
#define PHYSICS_POST_PROCESSOR_MANAGER_H

#include <memory>
#include <vector>
#include <fstream>
#include "post_processor_base.h"
#include "post_processor_scheduler.h"


namespace physics {

class PostProcessorManager
{
  public:
    PostProcessorManager(std::shared_ptr<PostProcessorScheduler> scheduler, const std::string& fname, MPI_Comm comm, int flush_interval=1);

    void addPostProcessor(PostProcessorPtr postproc);

    void runPostProcessors(int timestep, DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t);

    // Note: the destructor will close the file, this function is only required
    //       to close the file earlier
    void close();

  private:

    enum class State
    {
      Initializing,
      Writing,
      Closed
    };

    void writePostProcessorNames(PostProcessorPtr postproc);

    std::shared_ptr<PostProcessorScheduler> m_scheduler;
    std::vector<PostProcessorPtr> m_postprocessors;
    std::ofstream m_file;
    int m_flush_interval;
    MPI_Comm m_comm;
    bool m_am_i_root;
    State m_state = State::Initializing;
};

using PostProcessorManagerPtr = std::shared_ptr<PostProcessorManager>;

}  // namespace

#endif