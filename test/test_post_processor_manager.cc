#include "gtest/gtest.h"
#include "physics/post_processor_manager.h"
#include <vector>

namespace {

class PostProcessorConstant : public physics::PostProcessorBase
{
  public:
    explicit PostProcessorConstant(int id) :
      m_id(id)
    {}

    int numValues() const override { return 1; }

    std::vector<std::string> getNames() const override { return {std::string("PostProcessorConstant") + std::to_string(m_id)}; }

    std::vector<double> getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t) override { return {double(m_id)}; }

  private:
    int m_id;
};

class PostProcessorManagerTester : public ::testing::Test
{
  protected:
    PostProcessorManagerTester() :
      manager(std::make_shared<physics::PostProcessorScheduleFixedInterval>(3), "post_processor_test.txt")
    {
      manager.addPostProcessor(std::make_shared<PostProcessorConstant>(1));
      manager.addPostProcessor(std::make_shared<PostProcessorConstant>(2));      
    }

    physics::PostProcessorManager manager;
};

std::vector<std::string> split(const std::string& str, char delim)
{
  std::istringstream iss(str);
  std::string item;
  std::vector<std::string> strings;
  while (std::getline(iss, item, delim))
    strings.push_back(item);

  return strings;
}

}

TEST_F(PostProcessorManagerTester, Header)
{
  DiscVectorPtr u = nullptr;
  AuxiliaryEquationsStoragePtr u_aux = nullptr;
  manager.runPostProcessors(0, u, u_aux, 0.0);

  std::ifstream file;
  file.open("post_processor_test.txt");

  std::vector<std::string> expected_names = {"timestep", "time", "PostProcessorConstant1", "PostProcessorConstant2"};
  for (int i=0; i < 4; ++i)
  {
    std::string name;
    file >> name;
    std::cout << "name = " << name << std::endl;
    EXPECT_EQ(name, expected_names[i]);
  }

  char next_char;
  file.get(next_char);
  EXPECT_EQ(next_char, '\n');
}

TEST_F(PostProcessorManagerTester, Values)
{
  DiscVectorPtr u = nullptr;
  AuxiliaryEquationsStoragePtr u_aux = nullptr;
  double delta_t = 0.1;
  for (int i=0; i < 7; ++i)
    manager.runPostProcessors(i, u, u_aux, i*delta_t);
  manager.close();

  std::ifstream file;
  file.open("post_processor_test.txt");

  std::string line;
  std::getline(file, line); // header

  std::getline(file, line);
  std::vector<std::string> strings = split(line, ' ');
  EXPECT_EQ(strings.size(), 4u);
  EXPECT_EQ(std::stoi(strings[0]), 0);
  EXPECT_NEAR(std::stod(strings[1]), 0.0, 1e-13);
  EXPECT_NEAR(std::stod(strings[2]), 1.0, 1e-13);
  EXPECT_NEAR(std::stod(strings[3]), 2.0, 1e-13);

  
  std::getline(file, line);
  strings = split(line, ' ');
  EXPECT_EQ(strings.size(), 4u);
  EXPECT_EQ(std::stoi(strings[0]), 3);
  EXPECT_NEAR(std::stod(strings[1]), 0.3, 1e-13);
  EXPECT_NEAR(std::stod(strings[2]), 1.0, 1e-13);
  EXPECT_NEAR(std::stod(strings[3]), 2.0, 1e-13);

  std::getline(file, line);
  strings = split(line, ' ');
  EXPECT_EQ(strings.size(), 4u);
  EXPECT_EQ(std::stoi(strings[0]), 6);
  EXPECT_NEAR(std::stod(strings[1]), 0.6, 1e-13);
  EXPECT_NEAR(std::stod(strings[2]), 1.0, 1e-13);
  EXPECT_NEAR(std::stod(strings[3]), 2.0, 1e-13);

  EXPECT_FALSE(std::getline(file, line));
}