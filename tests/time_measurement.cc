#include <chrono>
#include <functional>
#include <fstream>
#include <iostream>
#include <string>


template <class T>
void MeasureTime(const std::function<T(...)> & func,  
                 size_t number_of_func_cals = 10, 
                 std::string path_to_save_file = "./data/measurement_results.cvs")
{
  std::ofstream outputFile(path_to_save_file);
  
  T res = 0;

  if (outputFile.is_open())
  {
    for (size_t i = 0; i < number_of_func_cals; ++i)
    {
      auto start = std::chrono::high_resolution_clock::now();

      func();   // use the correct call of the function

      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration = end - start;
      res += duration;

    }
    outputFile << res / number_of_func_cals;
    
    outputFile.close();
  }
  else { throw std::runtime_error("Unable to open file for writing\n"); }
}
