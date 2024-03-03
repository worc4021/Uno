// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_STATISTICS_H
#define UNO_STATISTICS_H

#include <string>
#include <map>
#include "tools/Options.hpp"

class Statistics {
public:
   explicit Statistics(const Options& options);

   static int int_width;
   static int double_width;
   static int char_width;

   void add_column(std::string name, int width, int order);
   void add_statistic(std::string name, std::string value);
   void add_statistic(std::string name, int value);
   void add_statistic(std::string name, size_t value);
   void add_statistic(std::string name, double value);
   void print_header(bool first_occurrence);
   void print_current_line();
   void print_footer();
   void new_line();

private:
   size_t iteration{0};
   std::map<int, std::string> columns{};
   std::map<std::string, int> widths{};
   std::map<std::string, std::string> current_line{};

   size_t print_header_every_iterations{};
   // std::string symbol(const std::string& value) const;
   struct Symbol {
      const std::string top{""};//{"─"};
      const std::string bottom{""};//{"─"};
      const std::string mid{""};//{"-"};
      const std::string top_left{""};//{"┌"};
      const std::string top_mid{""};//{"┬"};
      const std::string top_right{""};//{"┐"};
      const std::string left_mid{""};//{"├"};
      const std::string mid_mid{""};//{"┼"};
      const std::string right_mid{""};//{"┤"};
      const std::string left{""};//{"│"};
      const std::string middle{""};//{"│"};
      const std::string right{""};//{"│"};
      const std::string bottom_left{""};//{"└"};
      const std::string bottom_mid{""};//{"┴"};
      const std::string bottom_right{""};//{"┘"};
   } symbol;
};

#endif // UNO_STATISTICS_H