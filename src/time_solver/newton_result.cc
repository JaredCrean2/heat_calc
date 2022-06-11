#include "time_solver/newton_result.h"

namespace timesolvers {

std::ostream& operator<<(std::ostream& os, const NewtonResult& result)
{

  if (result.isConverged())
    os << "Newton solver converged in : " << result.getNIters() << " iterations: ";
  else
    os << "Newton solver failed to converge in : " << result.getNIters() << " iterations: ";


  if (result.isConverged())
  {
    if (result.isAbstolSatisfied())
      os << "Abstol satisfied " << result.getAbsNorm() << " < " << result.getAbsTol();

    if (result.isReltolSatisfied())
    {
      if (result.isAbstolSatisfied())
        os << ", ";
      
      os << "Reltol satisfied " << result.getRelNorm() << " < " << result.getRelTol();
    }
  } else
  {
    os << "Abstol not satisfied " << result.getAbsNorm() << " >= " << result.getAbsTol();
    os << ", reltol not satisfied " << result.getRelNorm() << " >= " << result.getRelTol();
  }

  os << std::endl;

  return os;
}

} // namespace