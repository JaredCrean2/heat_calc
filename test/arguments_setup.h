#ifndef TESTING_COMMAND_LINE_ARGUMENTS_H
#define TESTING_COMMAND_LINE_ARGUMENTS_H

#include <vector>
#include <string>

class CommandLineArguments
{
  public:
    explicit CommandLineArguments(const std::vector<std::string>& args) :
      m_args(args.size()),
      m_argv(args.size() + 1)
    {
      for (size_t i=0; i < m_args.size(); ++i)
      {
        // until C++17, there is no way to get a non-const char* from a string. 
        // So put the data in a vector and get a pointer from that
        int size_with_null = args[i].size() + 1;
        m_args[i].resize(size_with_null);
        for (int j=0; j < size_with_null; ++j)
          m_args[i][j] = args[i].c_str()[j];
        
        m_argv[i+1] = m_args[i].data();
      }        
    }

    int getArgc() { return m_argv.size(); }

    char** getArgv() { return m_argv.data(); }

  private:
    std::vector<std::vector<char>> m_args;
    std::vector<char*> m_argv;
};

#endif