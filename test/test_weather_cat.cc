#include "gtest/gtest.h"
#include "file/WeatherCat.h"

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

TEST(WeatherCat, CommandLineArgs)
{
  CommandLineArguments args({"--file", "file1.wea"});

  WeatherCatParsedData data = parseData(args.getArgc(), args.getArgv());

  EXPECT_EQ(data.filenames.size(), 1u);
  EXPECT_EQ(data.date_ranges.size(), 1u);
  EXPECT_EQ(data.filenames[0], "file1.wea");
  EXPECT_TRUE(data.date_ranges[0].is_in_range(DateTime{Date{1, 1, 1}, Time{0, 0}}));
  
  args = CommandLineArguments({"--file", "file1.wea", "1/2/2020", "3/2/2020"});
  data = parseData(args.getArgc(), args.getArgv());

  EXPECT_EQ(data.filenames.size(), 1u);
  EXPECT_EQ(data.date_ranges.size(), 1u);
  EXPECT_EQ(data.filenames[0], "file1.wea");
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{1, 1, 1}, Time{0, 0}}));
  EXPECT_TRUE(data.date_ranges[0].is_in_range(DateTime{Date{2, 1, 2020}, Time{0, 0}}));  
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{3, 3, 2020}, Time{0, 0}}));  

  args = CommandLineArguments({"--file", "file1.wea", "1/2/2020"});
  EXPECT_ANY_THROW(parseData(args.getArgc(), args.getArgv()));

  args = CommandLineArguments({"--file", "file1.wea", "1/2/2020", "3/2/2020", "--file", "file2.wea"});
  data = parseData(args.getArgc(), args.getArgv());

  EXPECT_EQ(data.filenames.size(), 2u);
  EXPECT_EQ(data.date_ranges.size(), 2u);
  EXPECT_EQ(data.filenames[0], "file1.wea");
  EXPECT_EQ(data.filenames[1], "file2.wea");
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{1, 1, 1}, Time{0, 0}}));
  EXPECT_TRUE(data.date_ranges[0].is_in_range(DateTime{Date{2, 1, 2020}, Time{0, 0}}));  
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{3, 3, 2020}, Time{0, 0}}));  
  EXPECT_TRUE(data.date_ranges[1].is_in_range(DateTime{Date{1, 1, 1}, Time{0, 0}}));


  args = CommandLineArguments({"--file", "file1.wea", "1/2/2020", "3/2/2020", "--file", "file2.wea", "2/1/2021", "2/4/2021"});
  data = parseData(args.getArgc(), args.getArgv());

  EXPECT_EQ(data.filenames.size(), 2u);
  EXPECT_EQ(data.date_ranges.size(), 2u);
  EXPECT_EQ(data.filenames[0], "file1.wea");
  EXPECT_EQ(data.filenames[1], "file2.wea");
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{1, 1, 1}, Time{0, 0}}));
  EXPECT_TRUE(data.date_ranges[0].is_in_range(DateTime{Date{2, 1, 2020}, Time{0, 0}}));  
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{3, 3, 2020}, Time{0, 0}}));  

  EXPECT_FALSE(data.date_ranges[1].is_in_range(DateTime{Date{31, 1, 2021}, Time{0, 0}}));
  EXPECT_TRUE(data.date_ranges[1].is_in_range(DateTime{Date{2, 2, 2021}, Time{0, 0}}));  
  EXPECT_FALSE(data.date_ranges[1].is_in_range(DateTime{Date{5, 2, 2021}, Time{0, 0}})); 

  args = CommandLineArguments({"--file", "file1.wea", "1/2/2020", "3/2/2020", "--file", "file2.wea", "2/1/2021"});
  EXPECT_ANY_THROW(parseData(args.getArgc(), args.getArgv()));  

  args = CommandLineArguments({"--file", "file1.wea", "1/2/2020-5:00", "3/2/2020-13:00"});
  data = parseData(args.getArgc(), args.getArgv());

  EXPECT_EQ(data.filenames.size(), 1u);
  EXPECT_EQ(data.date_ranges.size(), 1u);
  EXPECT_EQ(data.filenames[0], "file1.wea");
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{2, 1, 2020}, Time{4, 0}}));
  EXPECT_TRUE(data.date_ranges[0].is_in_range(DateTime{Date{2, 1, 2020}, Time{5, 0}}));

  EXPECT_TRUE(data.date_ranges[0].is_in_range(DateTime{Date{2, 3, 2020}, Time{13, 0}}));  
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{2, 3, 2020}, Time{14, 0}}));  
 
}