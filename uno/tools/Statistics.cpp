// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <sstream>
#include <iomanip>
#include "Statistics.hpp"

// TODO move this to the option file
int Statistics::int_width = 7;
int Statistics::double_width = 17;
int Statistics::char_width = 7;

Statistics::Statistics(const Options& options): print_header_every_iterations(options.get_unsigned_int("statistics_print_header_every_iterations")) {
}

void Statistics::add_column(std::string name, int width, int order) {
   this->columns[order] = name;
   this->widths[std::move(name)] = width;
}

void Statistics::add_statistic(std::string name, std::string value) {
   this->current_line[std::move(name)] = std::move(value);
}

void Statistics::add_statistic(std::string name, int value) {
   add_statistic(std::move(name), std::to_string(value));
}

void Statistics::add_statistic(std::string name, size_t value) {
   add_statistic(std::move(name), std::to_string(value));
}

void Statistics::add_statistic(std::string name, double value) {
   std::ostringstream stream;
   stream << std::defaultfloat << std::setprecision(7) << value;
   add_statistic(std::move(name), stream.str());
}

void Statistics::print_header(bool first_occurrence) {
   /* line above */
   std::cout << (first_occurrence ? Statistics::symbol.top_left : Statistics::symbol.left_mid);
   int k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << (first_occurrence ? Statistics::symbol.top_mid : Statistics::symbol.mid_mid);
      }
      std::string header = element.second;
      for (int j = 0; j < this->widths[header]; j++) {
         std::cout << Statistics::symbol.top;
      }
      k++;
   }
   std::cout << (first_occurrence ? Statistics::symbol.top_right : Statistics::symbol.right_mid) << '\n';
   /* headers */
   std::cout << Statistics::symbol.left;
   k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbol.middle;
      }
      std::string header = element.second;
      std::cout << " " << header;
      for (int j = 0; j < this->widths[header] - static_cast<int>(header.size()) - 1; j++) {
         std::cout << " ";
      }
      k++;
   }
   std::cout << Statistics::symbol.right << '\n';
}

void Statistics::print_current_line() {
   if (this->iteration % this->print_header_every_iterations == 0) {
      this->print_header(this->iteration == 0);
   }
   std::cout << Statistics::symbol.left_mid;
   int k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbol.mid_mid;
      }
      std::string header = element.second;
      for (int j = 0; j < this->widths[header]; j++) {
         std::cout << Statistics::symbol.bottom;
      }
      k++;
   }
   std::cout << Statistics::symbol.right_mid << '\n';
   /* headers */
   std::cout << Statistics::symbol.left;
   k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbol.middle;
      }
      const std::string& header = element.second;
      int size;
      try {
         std::string value = this->current_line.at(header);
         std::cout << " " << value;
         size = 1 + static_cast<int>(value.size());
      }
      catch (const std::out_of_range&) {
         std::cout << " -";
         size = 2;
      }
      int number_spaces = (size <= this->widths[header]) ? this->widths[header] - size : 0;
      for (int j = 0; j < number_spaces; j++) {
         std::cout << " ";
      }
      k++;
   }
   std::cout << Statistics::symbol.right << '\n';
   this->iteration++;
}

void Statistics::print_footer() {
   std::cout << Statistics::symbol.bottom_left;
   int k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbol.bottom_mid;
      }
      std::string header = element.second;
      for (int j = 0; j < this->widths[header]; j++) {
         std::cout << Statistics::symbol.bottom;
      }
      k++;
   }
   std::cout << Statistics::symbol.bottom_right << '\n';
}

void Statistics::new_line() {
   this->current_line.clear();
}

// std::string Statistics::symbol(const std::string& value) const {
//     using namespace std::string_literals;
//     if (value.compare("top") || value.compare("bottom") || value.compare("mid"))
//         return "─"s;
//     if (value.compare("top-mid"))
//         return "┬"s;
//     if (value.compare("top-left"))
//         return "┌"s;
//     if (value.compare("top-right"))
//         return "┐"s;
//     if (value.compare("bottom-mid"))
//         return "┴"s;
//     if (value.compare("bottom-left"))
//         return "└"s;
//     if (value.compare("bottom-right"))
//         return "┘"s;
//     if (value.compare("left") || value.compare("right") || value.compare("middle"))
//         return "│"s;
//     if (value.compare("left-mid"))
//         return "├"s;
//     if (value.compare("mid-mid"))
//         return "┼"s;
//     if (value.compare("right-mid"))
//         return "┤"s;
//     return "";

   // using namespace std::string_literals;
   //static std::map<std::string, std::string> symbols = {
   //      {"top", u8"─"s},
   //      {"top-mid", u8"┬"s},
   //      {"top-left", u8"┌"s},
   //      {"top-right", u8"┐"s},
   //      {"bottom", u8"─"s},
   //      {"bottom-mid", u8"┴"s},
   //      {"bottom-left", u8"└"s},
   //      {"bottom-right", u8"┘"s},
   //      {"left", u8"│"s},
   //      {"left-mid", u8"├"s},
   //      {"mid", u8"─"s},
   //      {"mid-mid", u8"┼"s},
   //      {"right", u8"│"s},
   //      {"right-mid", u8"┤"s},
   //      {"middle", u8"│"s}
   //};
   //return symbols[value];
// }