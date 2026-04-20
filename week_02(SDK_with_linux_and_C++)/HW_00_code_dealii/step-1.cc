// ---------------------------------------------------------------------
//
// Copyright (C) 2026 by Peter Munch
//
// This file is part of the class
// Advanced Topics in Numerics of PDEs: FEM for the Navier-Stokes Equations
// given at TU Berlin 2026.
//
// This is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// ---------------------------------------------------------------------

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

using namespace dealii;


void
first_grid()
{
  Triangulation<2> triangulation;

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4);

  std::ofstream out("grid-1.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-1.svg" << std::endl;
}


void
second_grid()
{
  Triangulation<2> triangulation;

  GridGenerator::hyper_ball_balanced(triangulation);
  triangulation.refine_global(3);

  std::ofstream out("grid-2.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-2.svg" << std::endl;
}

int
main()
{
  first_grid();
  second_grid();
}
