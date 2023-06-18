#include "physics/post_processor_manager.h"

#include <iomanip>


namespace physics {

PostProcessorManager::PostProcessorManager(std::shared_ptr<PostProcessorScheduler> scheduler, const std::string& fname, MPI_Comm comm, int flush_interval) :
  m_scheduler(scheduler),
  m_flush_interval(flush_interval),
  m_comm(comm),
  m_am_i_root(commRank(comm) == 0)
{
  if (m_am_i_root)
  {
    m_file.open(fname);
    m_file << std::setprecision(16);
    m_file << "timestep time";
  }
}


void PostProcessorManager::addPostProcessor(PostProcessorPtr postproc)
{
  assertAlways(m_state == State::Initializing, "cannot add new postprocessors after file writing has started");

  auto it = std::find(m_postprocessors.begin(), m_postprocessors.end(), postproc);
  if (it != m_postprocessors.end())
    throw std::runtime_error("attempted to add duplicate postprocessor");

  writePostProcessorNames(postproc);
  m_postprocessors.push_back(postproc);
}


void PostProcessorManager::runPostProcessors(int timestep, DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t)
{
  assertAlways(m_state == State::Initializing || m_state == State::Writing, "cannot write post processors after file has been closed");
  
  m_state = State::Writing;

  if (m_scheduler->shouldOutput(timestep, t))
  {
    m_file << "\n" << timestep << " " << t;
    for (auto& postproc : m_postprocessors)
    {
      auto vals = postproc->getValues(u, u_aux, t);
      for (int i=0; i < postproc->numValues(); ++i)
        m_file << " " <<  vals[i];
    }
  }

  if ((timestep % m_flush_interval) == 0)
    m_file.flush();
}

void PostProcessorManager::close()
{
  m_state = State::Closed;
  m_file.close();
}



void PostProcessorManager::writePostProcessorNames(PostProcessorPtr postproc)
{
  if (m_am_i_root)
    for (auto& name : postproc->getNames())
      m_file << " " << name;
}
}