//include/deal.II-translator/base/function_time.templates_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_function_time_templates_h
#define dealii_function_time_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/function_time.h>

DEAL_II_NAMESPACE_OPEN


template <typename Number>
FunctionTime<Number>::FunctionTime(const Number initial_time)
  : time(initial_time)
{}



template <typename Number>
void
FunctionTime<Number>::set_time(const Number new_time)
{
  time = new_time;
}


template <typename Number>
void
FunctionTime<Number>::advance_time(const Number delta_t)
{
  set_time(time + delta_t);
}


DEAL_II_NAMESPACE_CLOSE

#endif  /* dealii_function_time_templates_h */ 


